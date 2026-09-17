"""La Serena and Coquimbo: one grid cell, one forecast, many different rainfalls.

Every score in this project compares a forecast for a grid cell with an
observation, and the gap between the two is usually discussed as if it were the
model's error. Part of it is not. A cell is an areal average tens of kilometres
across, and inside a single cell the rain that actually fell can vary by a large
factor -- an error no model at this resolution can avoid, no matter how good its
physics.

Two cities in the transect make that limit measurable rather than arguable.
La Serena (-29.905, -71.249) and Coquimbo (-29.953, -71.339) are 12 km apart and
land in the same cell of all three verified models, which therefore issue
forecasts for them that are identical to the last decimal. The gauges inside
that cell do not agree with each other at all.

This script asserts the identity rather than asserting it in prose -- if a
future capture ever resolves the two cities separately, it fails loudly -- and
quantifies the spread the forecast cannot see.

Outputs subgrid_case.csv to data/processed/<capture>/metrics/ and fig9 to
figures/.

Usage: uv run scripts/subgrid_case.py [--capture YYYY-MM-DD]
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import stations_common as sc  # noqa: E402
from verify import (  # noqa: E402
    AXIS,
    ENSEMBLE_MEMBERS,
    INK,
    MODEL_COLORS,
    MODEL_LABELS,
    set_style,
)

ROOT = Path(__file__).resolve().parent.parent

PAIR = ("la_serena", "coquimbo")
# best_match is a blend that snaps to its own grid, so it is the one model that
# does tell the two cities apart. It is reported, not asserted about.
IDENTICAL_MODELS = ENSEMBLE_MEMBERS

STATION_COLOR = "#0d366b"
GRID_COLOR = "#6da7ec"


def assert_identical(tidy: pd.DataFrame) -> pd.DataFrame:
    """Check that the physical models give the two cities the same forecast."""
    pair = tidy[tidy["city"].isin(PAIR)]
    wide = pair.pivot_table(index=["model", "lead_days", "valid_date"],
                            columns="city", values="precip_fcst_mm").dropna()
    diff = (wide[PAIR[0]] - wide[PAIR[1]]).abs()
    rows = []
    for model, g in diff.groupby(level=0):
        rows.append({"model": model, "n_forecast_days": len(g),
                     "max_abs_difference_mm": float(g.max()),
                     "identical": bool(g.max() == 0.0)})
    summary = pd.DataFrame(rows)
    broken = summary[summary["model"].isin(IDENTICAL_MODELS)
                     & ~summary["identical"]]
    if not broken.empty:
        raise SystemExit(
            "the premise of this case study no longer holds -- these models now "
            f"resolve La Serena and Coquimbo separately:\n{broken.to_string()}")
    return summary


def cell_gauges(capture: str) -> pd.DataFrame:
    """Every gauge attached to either city, with its event total."""
    daily = sc.load_gauges(capture)
    quality = sc.gauge_quality(daily)
    w = daily[daily["city"].isin(PAIR)
              & daily["valid_date"].astype(str).between(*sc.EVENT_WINDOW)]
    # Several gauges sit within 30 km of both cities and therefore appear once
    # per city; the cell has one copy of each instrument, not two.
    w = w.drop_duplicates(subset=["station_id", "valid_date"])
    totals = (w.dropna(subset=["precip_station_mm"])
              .groupby(["station_id", "network"])["precip_station_mm"]
              .agg(total_mm="sum", n_days="size").reset_index())
    meta = []
    for name in ("stations_metadata.csv", "stations_metadata_dmc.csv"):
        path = ROOT / "data" / "external" / name
        if path.exists():
            meta.append(pd.read_csv(path))
    meta = (pd.concat(meta, ignore_index=True)
            .query("city in @PAIR")
            .drop_duplicates(subset="station_id")
            [["station_id", "station_name", "dist_km", "elev_m"]])
    totals["station_id"] = totals["station_id"].astype(str)
    meta["station_id"] = meta["station_id"].astype(str)
    out = totals.merge(meta, on="station_id", how="left")
    used = quality.loc[quality["used"], "station_id"].astype(str).unique()
    out["passes_qc"] = out["station_id"].isin(used)
    return out.sort_values("total_mm", ascending=False)


def fig_subgrid(gauges: pd.DataFrame, grid_total: float, fcst: dict,
                outdir: Path) -> None:
    ok = gauges[gauges["passes_qc"]]
    fig, ax = plt.subplots(figsize=(9, 4.6))
    ax.set_axisbelow(True)

    y = np.arange(len(gauges))
    colors = [STATION_COLOR if p else "#b3261e" for p in gauges["passes_qc"]]
    ax.barh(y, gauges["total_mm"], color=colors, height=0.62)
    labels = [f"{n}  ({d:.0f} km, {e:.0f} m)" if pd.notna(d) else str(n)
              for n, d, e in zip(gauges["station_name"], gauges["dist_km"],
                                 gauges["elev_m"])]
    ax.set_yticks(y, labels, fontsize=9)
    ax.invert_yaxis()

    ax.axvline(grid_total, color=GRID_COLOR, linewidth=2.5,
               label=f"gridded product: {grid_total:.0f} mm")
    for model, total in fcst.items():
        ax.axvline(total, color=MODEL_COLORS[model], linewidth=1.6,
                   linestyle="--", label=f"{MODEL_LABELS[model]}: {total:.0f} mm")

    lo, hi = ok["total_mm"].min(), ok["total_mm"].max()
    ax.axvspan(lo, hi, color=STATION_COLOR, alpha=0.07, zorder=0)
    ax.set_xlabel("Event total, 15–22 July (mm)")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title("One grid cell, one forecast, and the rain that actually fell\n"
                 "La Serena and Coquimbo share a cell: every model gives them "
                 "the same number", loc="left", fontsize=11)
    note = (f"gauges passing QC span {lo:.0f}–{hi:.0f} mm "
            f"(factor {hi / lo:.2f}) inside a single cell")
    if not gauges["passes_qc"].all():
        note += "   •   red = excluded by QC"
    ax.text(0, -0.16, note, transform=ax.transAxes, fontsize=8, color=AXIS)
    fig.tight_layout()
    fig.savefig(outdir / "fig9_subgrid_case.png", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--capture", default="2026-07-26")
    args = parser.parse_args()

    proc = ROOT / "data" / "processed" / args.capture
    tidy = pd.read_parquet(proc / "tidy_daily.parquet")
    metrics_dir = proc / "metrics"
    metrics_dir.mkdir(parents=True, exist_ok=True)

    summary = assert_identical(tidy)
    print("forecast agreement between La Serena and Coquimbo:")
    for r in summary.itertuples():
        verdict = "identical" if r.identical else "differs"
        print(f"  {r.model:14s} {verdict:9s} "
              f"max |Δ| = {r.max_abs_difference_mm:.6f} mm "
              f"over {r.n_forecast_days} forecast-days")

    gauges = cell_gauges(args.capture)
    ok = gauges[gauges["passes_qc"]]
    event = tidy[tidy["valid_date"].astype(str).between(*sc.EVENT_WINDOW)]
    grid_total = float(event[event["city"] == PAIR[0]]
                       [["valid_date", "precip_grid_mm"]].drop_duplicates()
                       ["precip_grid_mm"].sum())
    fcst = {}
    for model in IDENTICAL_MODELS:
        sel = event[(event["city"] == PAIR[0]) & (event["model"] == model)
                    & (event["lead_days"] == 1)]
        fcst[model] = float(sel["precip_fcst_mm"].sum())

    spread = {
        "n_gauges_total": len(gauges),
        "n_gauges_passing_qc": len(ok),
        "gauge_min_mm": float(ok["total_mm"].min()),
        "gauge_max_mm": float(ok["total_mm"].max()),
        "gauge_mean_mm": float(ok["total_mm"].mean()),
        "gauge_spread_factor": float(ok["total_mm"].max() / ok["total_mm"].min()),
        "grid_total_mm": grid_total,
        "model_intercity_difference_mm": 0.0,
        **{f"{m}_lead1_mm": v for m, v in fcst.items()},
    }
    pd.DataFrame([spread]).to_csv(metrics_dir / "subgrid_case.csv", index=False)
    gauges.to_csv(metrics_dir / "subgrid_gauges.csv", index=False)
    # The per-model table used to be printed and thrown away, which left the
    # one number that proves the point -- how far best_match, the interpolating
    # blend, separates two cities the physical models cannot tell apart -- with
    # nowhere to be checked against.
    summary.to_csv(metrics_dir / "subgrid_intercity.csv", index=False)

    print(f"\ngauges inside the cell: {len(ok)} passing QC of {len(gauges)}")
    print(f"  they span {spread['gauge_min_mm']:.1f}-{spread['gauge_max_mm']:.1f} mm "
          f"(factor {spread['gauge_spread_factor']:.2f}), mean "
          f"{spread['gauge_mean_mm']:.1f} mm")
    print(f"  the gridded product says {grid_total:.1f} mm")
    print(f"  the models' difference between the two cities is exactly 0.00 mm")

    set_style()
    fig_subgrid(gauges, grid_total, fcst, ROOT / "figures")
    print(f"\nmetrics -> {metrics_dir} (subgrid_case.csv, subgrid_gauges.csv)")
    print(f"figure  -> figures/fig9_subgrid_case.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
