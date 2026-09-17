"""Which conclusions survive a change of arbiter, and which do not.

The project scores every forecast against more than one truth, which raises the
obvious question of what to do when they disagree. Eyeballing two tables and
calling the difference "small" is not an answer: a gap of 0.3 mm/day is decisive
for a score whose sampling uncertainty is 0.05 and meaningless for one whose
uncertainty is 2.

So the comparison is made against each score's own sampling uncertainty. For
every metric x model x lead, the gap between two arbiters is measured in units
of the bootstrap standard error of that gap, and a verdict follows:

    robust    the arbiters agree within sampling noise
    shifted   they differ by 2-4 standard errors
    fragile   they differ by more than 4 -- the conclusion depends on the judge

The bootstrap resamples *dates*, not city-days. A frontal system hits twenty
cities on the same day, so city-days are not independent and resampling them
would understate the uncertainty by roughly the square root of the number of
cities. Both arbiters are recomputed on the same resampled dates, so the
resampling cancels out of the difference and what is left is the uncertainty of
the gap itself.

The gauge arbiter does not cover every city, so all comparisons are restricted
to the city-days both arbiters share -- otherwise the gap would mix a change of
judge with a change of sample.

Outputs robustness_by_metric.csv to data/processed/<capture>/metrics/.

Usage: uv run scripts/robustness.py [--capture YYYY-MM-DD] [--draws 2000]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify import (  # noqa: E402
    ENSEMBLE_MEMBERS,
    TRUTHS,
    TRUTH_LABELS,
    attach_gauge_truth,
)

ROOT = Path(__file__).resolve().parent.parent

THRESHOLD_MM = 10.0
SHIFTED_SIGMA = 2.0
FRAGILE_SIGMA = 4.0


def metrics(fcst: np.ndarray, obs: np.ndarray) -> dict[str, float]:
    """The scores whose arbiter-sensitivity the report discusses."""
    err = fcst - obs
    hit = ((fcst >= THRESHOLD_MM) & (obs >= THRESHOLD_MM)).sum()
    miss = ((fcst < THRESHOLD_MM) & (obs >= THRESHOLD_MM)).sum()
    false_alarm = ((fcst >= THRESHOLD_MM) & (obs < THRESHOLD_MM)).sum()
    return {
        "mae_mm": np.abs(err).mean(),
        "bias_mm": err.mean(),
        "pod": hit / (hit + miss) if hit + miss else np.nan,
        "far": false_alarm / (hit + false_alarm) if hit + false_alarm else np.nan,
    }


def paired_frame(tidy: pd.DataFrame, a: str, b: str) -> pd.DataFrame:
    """Rows where both arbiters have a value, so only the judge differs."""
    cols = [TRUTHS[a], TRUTHS[b]]
    if not set(cols).issubset(tidy.columns):
        return pd.DataFrame()
    return tidy.dropna(subset=cols)


def compare(df: pd.DataFrame, a: str, b: str, draws: int,
            rng: np.random.Generator) -> list[dict]:
    """Metric gaps between two arbiters, with a date-block bootstrap SE."""
    rows = []
    dates = np.array(sorted(df["valid_date"].unique()))
    for (model, lead), g in df.groupby(["model", "lead_days"]):
        fcst = g["precip_fcst_mm"].to_numpy()
        obs_a = g[TRUTHS[a]].to_numpy()
        obs_b = g[TRUTHS[b]].to_numpy()
        point_a, point_b = metrics(fcst, obs_a), metrics(fcst, obs_b)

        # One index array per date, so a resampled date brings all its cities.
        by_date = {d: np.flatnonzero(g["valid_date"].to_numpy() == d)
                   for d in dates}
        gaps: dict[str, list] = {k: [] for k in point_a}
        for _ in range(draws):
            pick = rng.choice(dates, size=len(dates), replace=True)
            idx = np.concatenate([by_date[d] for d in pick])
            ma = metrics(fcst[idx], obs_a[idx])
            mb = metrics(fcst[idx], obs_b[idx])
            for k in gaps:
                gaps[k].append(ma[k] - mb[k])

        for name in point_a:
            gap = point_a[name] - point_b[name]
            se = float(np.nanstd(gaps[name], ddof=1))
            sigma = abs(gap) / se if se > 0 else np.nan
            if not np.isfinite(sigma):
                verdict = "undetermined"
            elif sigma >= FRAGILE_SIGMA:
                verdict = "fragile"
            elif sigma >= SHIFTED_SIGMA:
                verdict = "shifted"
            else:
                verdict = "robust"
            # A bias that moves from -2.5 to +0.6 has not merely shifted: the
            # sentence "the models underestimate the rain" stops being true.
            # Magnitude alone cannot express that, so the sign change is its
            # own column.
            sign_flip = bool(np.isfinite(point_a[name])
                             and np.isfinite(point_b[name])
                             and point_a[name] * point_b[name] < 0)
            rows.append({
                "metric": name, "model": model, "lead_days": int(lead),
                "truth_a": a, "truth_b": b, "n": len(g),
                "n_cities": g["city"].nunique(), "n_dates": len(dates),
                f"value_{a}": point_a[name], f"value_{b}": point_b[name],
                "difference": gap, "bootstrap_se": se,
                "sigma": sigma, "verdict": verdict, "sign_flip": sign_flip,
            })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--capture", default="2026-07-26")
    parser.add_argument("--draws", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20260726)
    args = parser.parse_args()

    proc = ROOT / "data" / "processed" / args.capture
    tidy = pd.read_parquet(proc / "tidy_daily.parquet")
    tidy = attach_gauge_truth(tidy, args.capture)
    rng = np.random.default_rng(args.seed)

    rows = []
    for a, b in [("grid", "gauge"), ("grid", "era5")]:
        df = paired_frame(tidy, a, b)
        if df.empty:
            print(f"skipping {a} vs {b}: one of them is missing")
            continue
        print(f"{a} vs {b}: {df['city'].nunique()} cities, {len(df)} "
              f"forecast-days, {args.draws} bootstrap draws")
        rows += compare(df, a, b, args.draws, rng)

    out = pd.DataFrame(rows)
    path = proc / "metrics" / "robustness_by_metric.csv"
    out.to_csv(path, index=False)
    print(f"\nwrote {path}")

    for (a, b), g in out.groupby(["truth_a", "truth_b"]):
        print(f"\n=== {TRUTH_LABELS[a]} vs {TRUTH_LABELS[b]} ===")
        tally = (g.groupby(["metric", "verdict"]).size()
                 .unstack(fill_value=0))
        tally["sign_flips"] = g.groupby("metric")["sign_flip"].sum()
        print(tally.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
