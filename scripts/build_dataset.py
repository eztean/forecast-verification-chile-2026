"""Build a tidy verification dataset from the raw captured JSON snapshots.

Reads data/raw/<capture-date>/ (defaults to the most recent capture) and produces
one long-format table with a row per city x valid_date x model x lead_days:

    city, valid_date, model, lead_days,
    precip_fcst_mm, gust_fcst_kmh,      # daily aggregates of the archived forecast
    precip_grid_mm, gust_grid_kmh,      # Open-Meteo archive best_match, the
                                        # gridded series every score is computed
                                        # against
    precip_era5_mm, gust_era5_kmh       # true ERA5 at 0.25 deg, for contrast

The two gridded columns are named for what they are, and the distinction is not
cosmetic. `best_match` is the archive default this project used from the start,
but the documentation says it "combines IFS HRES, ERA5 and ERA5-Land
seamlessly", and for 2026 it serves ECMWF IFS HRES at 9 km -- the operational
analysis of one of the models being verified, not an independent reanalysis.
Calling it ERA5, as earlier versions of this pipeline did, understated how
closely the arbiter is related to one of the contenders. Carrying both columns
lets the analysis separate the choice of product from the representativeness
gap between a grid cell and a rain gauge.

Hourly forecast values (precipitation_previous_dayN, wind_gusts_10m_previous_dayN)
are aggregated to daily totals (precipitation) and daily maxima (wind gusts) in
local time (America/Santiago), matching the daily observation aggregates.

Output: data/processed/<capture-date>/tidy_daily.csv (and .parquet).

Usage: uv run scripts/build_dataset.py [--capture YYYY-MM-DD]
"""

import argparse
import json
import re
from pathlib import Path

import pandas as pd

MODELS = ["best_match", "gfs_seamless", "ecmwf_ifs025", "icon_seamless"]
LEAD_DAYS = range(1, 8)
ROOT = Path(__file__).resolve().parent.parent
DATE_DIR = re.compile(r"\d{4}-\d{2}-\d{2}")


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def forecast_daily(raw: dict, model: str, city: str) -> pd.DataFrame:
    """Aggregate hourly previous-day forecasts to daily values, long format."""
    hourly = raw["hourly"]
    time = pd.to_datetime(hourly["time"])  # local time (America/Santiago)
    date = time.date
    frames = []
    for lead in LEAD_DAYS:
        df = pd.DataFrame(
            {
                "date": date,
                "precip": hourly[f"precipitation_previous_day{lead}"],
                "gust": hourly[f"wind_gusts_10m_previous_day{lead}"],
            }
        )
        # A daily total is only valid with all 24 hourly values present:
        # sum(min_count=24) yields NaN otherwise (e.g. ICON has no lead-7
        # data in the previous-runs API — all-null hours must not become 0.0).
        daily = df.groupby("date").agg(
            precip_fcst_mm=("precip", lambda s: s.sum(min_count=24)),
            gust_fcst_kmh=("gust", "max"),
            n_hours=("precip", "size"),
        )
        # Guard against partial days at the window edges, and drop lead/day
        # combinations the model never forecast (all-null precipitation).
        daily = daily[daily["n_hours"] == 24].drop(columns="n_hours")
        daily = daily.dropna(subset=["precip_fcst_mm"])
        daily = daily.reset_index().rename(columns={"date": "valid_date"})
        daily.insert(0, "city", city)
        daily.insert(2, "model", model)
        daily.insert(3, "lead_days", lead)
        frames.append(daily)
    return pd.concat(frames, ignore_index=True)


def observed_daily(raw: dict, city: str, tag: str) -> pd.DataFrame:
    daily = raw["daily"]
    return pd.DataFrame(
        {
            "city": city,
            "valid_date": pd.to_datetime(daily["time"]).date,
            f"precip_{tag}_mm": daily["precipitation_sum"],
            f"gust_{tag}_kmh": daily["wind_gusts_10m_max"],
        }
    )


def build(capture_dir: Path) -> pd.DataFrame:
    cities = sorted(
        {p.name.split("__")[0] for p in capture_dir.glob("*__observed__grid.json")}
    )
    fcst_frames, obs_frames = [], []
    for city in cities:
        for model in MODELS:
            raw = load_json(capture_dir / f"{city}__prevruns__{model}.json")
            fcst_frames.append(forecast_daily(raw, model, city))
        obs = observed_daily(load_json(capture_dir / f"{city}__observed__grid.json"),
                             city, "grid")
        # True ERA5 joined a capture that already existed, so a directory
        # without it still builds -- the era5 columns simply come out NaN.
        era5_path = capture_dir / f"{city}__observed__era5.json"
        if era5_path.exists():
            obs = obs.merge(observed_daily(load_json(era5_path), city, "era5"),
                            on=["city", "valid_date"], how="left")
        obs_frames.append(obs)
    fcst = pd.concat(fcst_frames, ignore_index=True)
    obs = pd.concat(obs_frames, ignore_index=True)
    tidy = fcst.merge(obs, on=["city", "valid_date"], how="inner")
    return tidy.sort_values(["city", "valid_date", "model", "lead_days"]).reset_index(
        drop=True
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--capture", default=None,
                        help="capture date YYYY-MM-DD (default: latest in data/raw/)")
    args = parser.parse_args()

    raw_root = ROOT / "data" / "raw"
    # Only dated directories are captures. data/raw/ also holds `stations/`,
    # which sorts after every date and would otherwise be picked as "latest",
    # failing with an unhelpful "No objects to concatenate".
    captures = sorted(p.name for p in raw_root.iterdir()
                      if p.is_dir() and DATE_DIR.fullmatch(p.name))
    if not captures:
        raise SystemExit(f"no dated capture directories under {raw_root}")
    capture = args.capture or captures[-1]
    capture_dir = raw_root / capture
    if not capture_dir.is_dir():
        raise SystemExit(f"no capture {capture}; available: {', '.join(captures)}")

    tidy = build(capture_dir)

    outdir = ROOT / "data" / "processed" / capture
    outdir.mkdir(parents=True, exist_ok=True)
    tidy.to_csv(outdir / "tidy_daily.csv", index=False)
    tidy.to_parquet(outdir / "tidy_daily.parquet", index=False)

    n_cities = tidy["city"].nunique()
    dates = tidy["valid_date"]
    print(f"capture {capture}: {len(tidy)} rows | {n_cities} cities | "
          f"{tidy['model'].nunique()} models | leads {tidy['lead_days'].min()}-"
          f"{tidy['lead_days'].max()} | valid dates {dates.min()} .. {dates.max()}")
    print(f"wrote {outdir / 'tidy_daily.csv'} and .parquet")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
