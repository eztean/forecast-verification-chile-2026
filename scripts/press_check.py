"""Cross-check this project's gauge harvest against published DMC figures.

`data/external/media_reported_totals.csv` holds precipitation totals as
published by Chilean media and official balances during and after the July 2026
event, each row carrying its station, its accumulation window and its URL. This
script sums *our* daily series for the same station over the same window and
reports the difference.

This is a plausibility check on the harvest, never verification truth. The point
is not that the numbers agree -- it is that where they disagree, the reason is
identifiable. Two mechanisms account for essentially all of the disagreement
here, and only one of them is ours:

  * Days dropped by QC. A day whose sub-daily series is less than
    MIN_COMPLETENESS_PC complete becomes NaN rather than a partial total, so a
    month that includes such a day sums short. 31 July is incomplete for every
    station in the capture, which puts a floor under every monthly comparison.
  * Window ambiguity. "During the event" in a story filed on 19 July means
    something the story does not state precisely; those rows are matched to
    16-18 July and flagged.

Usage: uv run scripts/press_check.py [--capture YYYY-MM-DD]
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stations_common import load_gauges  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent

# Pulled out of the station name in the reference table, which carries it in
# parentheses so a human reading the CSV can see which instrument is meant.
ID_PATTERN = r"\((?:DMC|CEAZA)\s+([A-Za-z0-9]+)\)"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--capture", default="2026-07-26")
    args = parser.parse_args()

    ref = pd.read_csv(ROOT / "data" / "external" / "media_reported_totals.csv")
    if ref.empty:
        print("no reference rows; nothing to check")
        return 0
    ref["station_id"] = ref["station"].str.extract(ID_PATTERN)[0]

    daily = load_gauges(args.capture).drop_duplicates(
        subset=["station_id", "valid_date"])
    daily["sid"] = daily["station_id"].astype(str)

    rows = []
    for r in ref.itertuples():
        s = daily[(daily["sid"] == str(r.station_id))
                  & (daily["valid_date"].astype(str)
                     .between(r.window_start, r.window_end))]
        ours = s["precip_station_mm"].sum()
        n_missing = int(s["precip_station_mm"].isna().sum())
        rows.append({
            "city": r.city, "station": r.station,
            "window": f"{r.window_start}..{r.window_end}",
            "reported_mm": r.reported_mm, "ours_mm": round(float(ours), 1),
            "diff_mm": round(float(ours) - r.reported_mm, 1),
            "diff_pc": round(100 * (ours / r.reported_mm - 1), 1)
            if r.reported_mm else float("nan"),
            "days_present": int(s["precip_station_mm"].notna().sum()),
            "days_dropped_by_qc": n_missing,
            "source": r.source,
        })
    out = pd.DataFrame(rows)

    metrics = ROOT / "data" / "processed" / args.capture / "metrics"
    metrics.mkdir(parents=True, exist_ok=True)
    out.to_csv(metrics / "press_cross_check.csv", index=False)

    pd.set_option("display.width", 200)
    print(out.drop(columns=["source"]).to_string(index=False))
    complete = out[out["days_dropped_by_qc"] == 0]
    partial = out[out["days_dropped_by_qc"] > 0]
    print(f"\nrows where our series is complete over the window: {len(complete)}"
          f" -> median difference {complete['diff_pc'].median():+.1f}%, "
          f"worst {complete['diff_pc'].abs().max():.1f}%")
    print(f"rows with QC-dropped days: {len(partial)}"
          f" -> median difference {partial['diff_pc'].median():+.1f}% "
          f"(short, as dropping a day must make it)")
    print(f"\nwrote {metrics / 'press_cross_check.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
