"""Capture archived weather forecasts and observations for the July 2026 frontal system in Chile.

Fetches, for a transect of Chilean cities:
  1. Previous-runs forecasts (Open-Meteo Previous Runs API): hourly precipitation and
     wind gusts as forecast 1-7 days ahead, for the recent past window.
  2. Observations (Open-Meteo Archive API): daily and hourly precipitation / wind gusts.

Raw responses are stored untouched in data/raw/<capture-date>/ as JSON, one file per
city x source x model, so the analysis is reproducible from immutable raw snapshots.
Run it repeatedly on different days; each capture goes to its own dated directory.

Usage: uv run scripts/capture_forecasts.py [--past-days N]
"""

import argparse
import datetime as dt
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

CITIES = {
    "la_serena": (-29.905, -71.249),
    # Coquimbo is deliberately kept even though every model resolves it to the
    # same grid cell as La Serena and returns a bit-identical forecast. That
    # identity is the point: the gauges inside that one cell measured 128 to
    # 214 mm during the event, so the pair makes the sub-grid representativeness
    # limit measurable instead of merely arguable. See scripts/subgrid_case.py.
    "coquimbo": (-29.953, -71.339),
    # Coquimbo region (IV) inland/coastal localities — the frontal system hit
    # this semi-arid region hardest, so it gets a denser sub-transect.
    "vicuna": (-30.032, -70.708),
    "andacollo": (-30.230, -71.085),
    "ovalle": (-30.598, -71.200),
    "monte_patria": (-30.695, -70.958),
    "combarbala": (-31.178, -71.003),
    "illapel": (-31.633, -71.166),
    "salamanca": (-31.779, -70.965),
    "los_vilos": (-31.915, -71.513),
    "valparaiso": (-33.046, -71.620),
    "santiago": (-33.447, -70.673),
    "rancagua": (-34.171, -70.741),
    "curico": (-34.983, -71.239),
    "talca": (-35.427, -71.655),
    "chillan": (-36.607, -72.104),
    "concepcion": (-36.827, -73.050),
    "temuco": (-38.736, -72.590),
    "valdivia": (-39.814, -73.246),
    "osorno": (-40.574, -73.132),
    "puerto_montt": (-41.469, -72.943),
}

MODELS = ["best_match", "gfs_seamless", "ecmwf_ifs025", "icon_seamless"]
LEAD_DAYS = range(1, 8)  # forecasts issued 1..7 days before
HOURLY_BASE_VARS = ["precipitation", "wind_gusts_10m"]
TZ = "America/Santiago"

PREVIOUS_RUNS_URL = "https://previous-runs-api.open-meteo.com/v1/forecast"
ARCHIVE_URL = "https://archive-api.open-meteo.com/v1/archive"


def fetch(url: str, params: dict) -> dict:
    full = url + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(full, timeout=60) as resp:
        data = json.loads(resp.read().decode())
    if "error" in data and data.get("error"):
        raise RuntimeError(f"API error for {full}: {data.get('reason')}")
    return data


def save(outdir: Path, name: str, payload: dict, meta: dict) -> None:
    payload["_capture"] = meta
    path = outdir / f"{name}.json"
    path.write_text(json.dumps(payload))
    print(f"  saved {path.name} ({path.stat().st_size / 1024:.0f} KiB)")


def already_captured(outdir: Path, name: str, skip: bool) -> bool:
    if skip and (outdir / f"{name}.json").exists():
        print(f"  kept {name}.json")
        return True
    return False


def observation_window(outdir: Path) -> tuple[str, str] | None:
    """The archive window an existing capture used, read back from its own JSON.

    Adding a city to a capture that already exists has to reproduce that
    capture's date window exactly, or the new city would carry different valid
    dates from every other one. The window is not guessed: it is read from the
    `_capture.params` any observation file in the directory already records.
    """
    for path in sorted(outdir.glob("*__observed__grid.json")):
        params = json.loads(path.read_text()).get("_capture", {}).get("params", {})
        if "start_date" in params and "end_date" in params:
            return params["start_date"], params["end_date"]
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--past-days", type=int, default=14,
                        help="how many past days of forecasts/observations to pull")
    parser.add_argument("--cities", default=None,
                        help="comma-separated subset of CITIES (default: all)")
    parser.add_argument("--into", default=None,
                        help="add to an existing capture directory (YYYY-MM-DD) "
                             "reusing its date window, instead of starting a new one")
    parser.add_argument("--skip-existing", action="store_true",
                        help="never overwrite a file already in the capture; use "
                             "with --into to add a city or a series without "
                             "disturbing the snapshot the published numbers rest on")
    args = parser.parse_args()

    now = dt.datetime.now(dt.timezone.utc)
    capture_tag = args.into or now.strftime("%Y-%m-%d")
    outdir = Path(__file__).resolve().parent.parent / "data" / "raw" / capture_tag
    outdir.mkdir(parents=True, exist_ok=True)
    # Each file keeps its own real capture timestamp, so a directory that gained
    # a city later is a capture *set*, not an instant, and says so.
    meta_base = {"captured_utc": now.isoformat(), "script": "capture_forecasts.py"}

    hourly_prev = ",".join(
        f"{var}_previous_day{n}" for var in HOURLY_BASE_VARS for n in LEAD_DAYS
    )
    end_date = now.astimezone(dt.timezone(dt.timedelta(hours=-4))).date()
    start_date = end_date - dt.timedelta(days=args.past_days)
    if args.into:
        window = observation_window(outdir)
        if window is None:
            raise SystemExit(f"{outdir} has no observation file to read a "
                             f"date window from")
        start_date, end_date = window
        print(f"reusing the window of capture {capture_tag}: "
              f"{start_date} .. {end_date}")

    cities = dict(CITIES)
    if args.cities:
        wanted = [c.strip() for c in args.cities.split(",")]
        unknown = [c for c in wanted if c not in CITIES]
        if unknown:
            raise SystemExit(f"unknown cities: {unknown}")
        cities = {c: CITIES[c] for c in wanted}

    failures = []
    for city, (lat, lon) in cities.items():
        print(f"{city}:")
        for model in MODELS:
            if already_captured(outdir, f"{city}__prevruns__{model}",
                                args.skip_existing):
                continue
            params = {
                "latitude": lat, "longitude": lon, "timezone": TZ,
                "hourly": hourly_prev, "models": model,
            }
            # A fresh capture walks back from today; one that joins an existing
            # capture has to ask for that capture's dates explicitly.
            if args.into:
                params |= {"start_date": str(start_date),
                           "end_date": str(end_date)}
            else:
                params |= {"past_days": args.past_days, "forecast_days": 1}
            try:
                data = fetch(PREVIOUS_RUNS_URL, params)
                save(outdir, f"{city}__prevruns__{model}",
                     data, {**meta_base, "params": params})
            except Exception as exc:  # capture what we can, report the rest
                failures.append((city, model, str(exc)))
                print(f"  FAILED prevruns {model}: {exc}", file=sys.stderr)
            time.sleep(1)

        # Two gridded observation series, and it matters which is which.
        # `best_match` is the archive default this project has always used, but
        # it is not ERA5: the docs say it "combines IFS HRES, ERA5 and ERA5-Land
        # seamlessly", and for 2026 it serves ECMWF IFS HRES at 9 km -- the
        # operational analysis of one of the very models being verified. It is
        # now requested explicitly rather than by omission, and true ERA5 at
        # 0.25 deg is captured alongside it so the two can be told apart. The
        # file suffixes say which is which: `grid` for the blend this project
        # scores against, `era5` for the reanalysis itself.
        for model, suffix in [("best_match", "grid"), ("era5", "era5")]:
            if already_captured(outdir, f"{city}__observed__{suffix}",
                                args.skip_existing):
                continue
            params = {
                "latitude": lat, "longitude": lon, "timezone": TZ,
                "start_date": str(start_date), "end_date": str(end_date),
                "models": model,
                "hourly": ",".join(HOURLY_BASE_VARS),
                "daily": "precipitation_sum,wind_gusts_10m_max",
            }
            try:
                data = fetch(ARCHIVE_URL, params)
                save(outdir, f"{city}__observed__{suffix}",
                     data, {**meta_base, "params": params})
            except Exception as exc:
                failures.append((city, f"archive/{model}", str(exc)))
                print(f"  FAILED archive {model}: {exc}", file=sys.stderr)
            time.sleep(1)

    print(f"\nCapture complete: {len(list(outdir.glob('*.json')))} files in {outdir}")
    if failures:
        print(f"{len(failures)} failures:", file=sys.stderr)
        for city, src, err in failures:
            print(f"  {city}/{src}: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
