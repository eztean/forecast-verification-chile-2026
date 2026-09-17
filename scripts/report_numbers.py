"""Every number quoted in the project's prose that no metrics CSV holds.

Explaining the Brier decomposition and the CRPS means working through real
cases, and those worked cases carry numbers the metrics tables never store: the
three Murphy terms, the fair-CRPS variant, the mean error of an individual
member, and two worked ensemble cases. Plus the plain descriptive facts a
write-up opens with -- how many rows the dataset has, how many cases per model
and lead, where the event window falls. Those were computed by hand once and
then had to be re-derived from memory at every refresh, which is exactly the
kind of untraceable number the project forbids.

This script recomputes all of them from the tidy dataset and writes them to
data/processed/<capture>/metrics/report_numbers.json, so the prose can be
checked against the data with a diff.

Usage: uv run scripts/report_numbers.py [--capture YYYY-MM-DD]
"""

import argparse
import datetime as dt
import json
from pathlib import Path

import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify import TRUTH_LABELS, attach_gauge_truth, by_truth  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
ENSEMBLE_MEMBERS = ["gfs_seamless", "ecmwf_ifs025", "icon_seamless"]
POOLED_LEADS = range(1, 7)  # ICON has no lead 7, so the ensemble stops at 6


def wide_members(tidy: pd.DataFrame) -> pd.DataFrame:
    members = tidy[tidy["model"].isin(ENSEMBLE_MEMBERS)]
    return members.pivot_table(
        index=["city", "valid_date", "lead_days", "precip_obs_mm"],
        columns="model", values="precip_fcst_mm").dropna().reset_index()


def murphy(wide: pd.DataFrame, thr: float) -> dict:
    """Reliability / resolution / uncertainty, pooled over leads 1-6."""
    g = wide[wide["lead_days"].isin(POOLED_LEADS)]
    p = (g[ENSEMBLE_MEMBERS] >= thr).mean(axis=1)
    o = (g["precip_obs_mm"] >= thr).astype(float)
    n, s = len(g), o.mean()
    df = pd.DataFrame({"p": p, "o": o})
    grp = df.groupby("p")["o"].agg(["size", "mean"])
    rel = (grp["size"] * (grp.index - grp["mean"]) ** 2).sum() / n
    res = (grp["size"] * (grp["mean"] - s) ** 2).sum() / n
    unc = s * (1 - s)
    bs = ((p - o) ** 2).mean()
    return {"n": int(n), "base_rate": s, "reliability": rel, "resolution": res,
            "uncertainty": unc, "brier": bs, "identity_check": rel - res + unc,
            "bss": 1 - bs / unc,
            "bins": [{"p": float(k), "n": int(r["size"]), "obs_freq": r["mean"]}
                     for k, r in grp.iterrows()]}


def crps_variants(wide: pd.DataFrame) -> dict:
    """NRG and fair CRPS by lead, plus the mean error of a single member."""
    out = {}
    for lead, g in wide.groupby("lead_days"):
        x = g[ENSEMBLE_MEMBERS].to_numpy()
        y = g["precip_obs_mm"].to_numpy()[:, None]
        m = x.shape[1]
        member_err = np.abs(x - y).mean(axis=1)
        pair = np.abs(x[:, :, None] - x[:, None, :]).sum(axis=(1, 2))
        nrg = member_err - pair / (2 * m * m)
        fair = member_err - pair / (2 * m * (m - 1))
        out[int(lead)] = {
            "n": int(len(g)), "crps_nrg_mm": float(nrg.mean()),
            "crps_fair_mm": float(fair.mean()),
            "mean_single_member_abs_error_mm": float(member_err.mean()),
        }
    return out


def worked_case(wide: pd.DataFrame, city: str, date: str, lead: int) -> dict:
    row = wide[(wide["city"] == city)
               & (wide["valid_date"].astype(str) == date)
               & (wide["lead_days"] == lead)]
    if row.empty:
        return {"city": city, "date": date, "lead": lead, "available": False}
    r = row.iloc[0]
    x = np.array([r[m] for m in ENSEMBLE_MEMBERS])
    y = float(r["precip_obs_mm"])
    m = len(x)
    member_err = float(np.abs(x - y).mean())
    pair = float(np.abs(x[:, None] - x[None, :]).sum())
    return {
        "city": city, "date": date, "lead": lead, "available": True,
        "members": {k: float(v) for k, v in zip(ENSEMBLE_MEMBERS, x)},
        "observed_mm": y, "mean_member_error_mm": member_err,
        "spread_discount_mm": pair / (2 * m * m),
        "crps_mm": member_err - pair / (2 * m * m),
        "observed_above_all_members": bool(y > x.max()),
        "observed_inside_range": bool(x.min() <= y <= x.max()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--capture", default="2026-07-26")
    args = parser.parse_args()

    proc = ROOT / "data" / "processed" / args.capture
    tidy = pd.read_parquet(proc / "tidy_daily.parquet")
    tidy = attach_gauge_truth(tidy, args.capture)

    out = {
        "capture": args.capture,
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "dataset": {
            "rows": int(len(tidy)),
            "cities": int(tidy["city"].nunique()),
            "valid_dates": [str(tidy["valid_date"].min()),
                            str(tidy["valid_date"].max())],
            "n_days": int(tidy["valid_date"].nunique()),
            "n_per_model_lead": int(
                tidy.groupby(["model", "lead_days"]).size().max()),
        },
        "truths": {},
    }

    # The one event window, derived from the gauges (stations_common) rather
    # than hand-set here, so this file and the figures cannot drift apart.
    import stations_common as sc
    est = (tidy[["city", "valid_date", "precip_gauge_mm"]].drop_duplicates()
           .dropna(subset=["precip_gauge_mm"]))
    window = sc.event_window(est)
    out["event_window"] = list(window)

    # Every prose number in the report is recomputed per arbiter, because the
    # report now quotes several of them against both and the whole point of the
    # exercise is that they differ.
    for truth, frame in by_truth(tidy):
        wide = wide_members(frame)
        cont = (frame.groupby(["model", "lead_days"])
                .apply(lambda g: pd.Series({
                    "mae": (g["precip_fcst_mm"] - g["precip_obs_mm"]).abs().mean(),
                    "bias": (g["precip_fcst_mm"] - g["precip_obs_mm"]).mean()}),
                    include_groups=False).reset_index())
        obs = frame[["city", "valid_date", "precip_obs_mm"]].drop_duplicates()
        ev = obs[obs["valid_date"].astype(str).between(*window)]
        out["truths"][truth] = {
            "label": TRUTH_LABELS[truth],
            "cities": int(frame["city"].nunique()),
            "n_forecast_days": int(len(frame)),
            "n_pooled_leads_1_6": int(
                len(wide[wide["lead_days"].isin(POOLED_LEADS)])),
            "murphy": {str(t): murphy(wide, t) for t in (1.0, 10.0, 20.0)},
            "crps": crps_variants(wide),
            "worked_cases": [
                worked_case(wide, "la_serena", "2026-07-17", 3),
                worked_case(wide, "santiago", "2026-07-17", 3),
                worked_case(wide, "la_serena", "2026-07-17", 5),
            ],
            "event_totals_mm": (
                ev.groupby("city")["precip_obs_mm"].sum()
                .round(1).sort_values(ascending=False).to_dict()),
            "mae_by_model_lead": {
                m: dict(zip(g["lead_days"].astype(int), g["mae"].round(2)))
                for m, g in cont.groupby("model")},
            "bias_by_model_lead": {
                m: dict(zip(g["lead_days"].astype(int), g["bias"].round(2)))
                for m, g in cont.groupby("model")},
        }

    path = proc / "metrics" / "report_numbers.json"
    path.write_text(json.dumps(out, indent=2, default=float))
    print(f"wrote {path}")

    d = out["dataset"]
    print(f"\n{d['cities']} cities, n per model x lead: {d['n_per_model_lead']}")
    for truth, t in out["truths"].items():
        m10 = t["murphy"]["10.0"]
        print(f"\n--- truth={truth} ({t['label']}) --- "
              f"{t['cities']} cities, {t['n_forecast_days']} forecast-days, "
              f"pooled leads 1-6: {t['n_pooled_leads_1_6']}")
        print(f"  Murphy >=10mm: REL {m10['reliability']:.3f} - RES "
              f"{m10['resolution']:.3f} + UNC {m10['uncertainty']:.3f} = "
              f"{m10['identity_check']:.3f} (BS {m10['brier']:.3f}), "
              f"BSS {m10['bss']:.3f}, base rate {m10['base_rate']:.3f}")
        print("  reliability bins: " + ", ".join(
            f"p={b['p']:.2f} n={b['n']} obs={b['obs_freq']:.3f}"
            for b in m10["bins"]))
        for lead in (1, 6):
            c = t["crps"][lead]
            print(f"  CRPS lead {lead}: NRG {c['crps_nrg_mm']:.2f}, fair "
                  f"{c['crps_fair_mm']:.2f}, single member "
                  f"{c['mean_single_member_abs_error_mm']:.2f} mm")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
