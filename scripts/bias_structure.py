"""What structure does the gridded analysis's wet bias actually have?

The gridded product runs ~30% wetter than the rain gauges over the transect.
An earlier pass at this analysis explained the city-to-city spread of that
excess as a geographic gradient -- "coastal and foothill cells blow up above
+50%, interior valleys agree within 8-21%" -- and, two paragraphs later, said
the excess was not regional at all. Both statements cannot hold, and neither
was tested. This script tests them, and asks the question the right way round:
given that a city's gauge estimate is a mean of a handful of point
measurements inside a ~80 km^2 cell, how much of the between-city spread is
signal and how much is the sampling noise of that estimate?

Four diagnostics, all written to metrics/bias_structure_*.csv:

  C1  Intensity structure. The grid/gauge ratio in bins of the *symmetric*
      average of the two series -- binning on either one alone selects on its
      own noise and manufactures a trend. Plus the wet-day frequency vs
      wet-day intensity split, which says whether the grid rains too often or
      too hard.
  C2  Geography, tested rather than asserted: coastal vs interior, and the rank
      correlation of the excess with latitude, elevation, cell elevation,
      gauge-to-city distance and event size.
  C3  The noise floor. The within-cell scatter of gauges gives the standard
      error of each city's estimate; comparing the observed variance of
      log(grid/gauge) with the variance that scatter alone implies says how
      much real between-city structure is left to explain.
  C4  Timing. Per-city daily correlation between grid and gauges at lags of
      -2..+2 days, which is how the La Serena cell's two-day offset was found.

Usage: uv run scripts/bias_structure.py [--capture YYYY-MM-DD]
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from capture_forecasts import CITIES  # noqa: E402
from stations_common import (  # noqa: E402
    city_gauge_daily,
    event_window,
    gauge_quality,
    load_gauges,
)
from verify import (  # noqa: E402
    AXIS,
    CITY_LABELS,
    CITY_ORDER,
    INK,
    set_style,
)

ROOT = Path(__file__).resolve().parent.parent

STATION_COLOR = "#0d366b"
GRID_COLOR = "#6da7ec"

# Cities whose grid cell contains open ocean. Assigned from the coastline, not
# from the result: this is the classification the discarded "coastal cells run
# wettest" claim needs in order to be testable at all.
COASTAL = {"la_serena", "coquimbo", "los_vilos", "valparaiso", "concepcion",
           "valdivia", "puerto_montt"}

# Bin edges for the intensity structure, on the average of grid and gauge.
INTENSITY_BINS = [0, 1, 5, 15, 40, np.inf]


def spearman(a: pd.Series, b: pd.Series) -> float:
    ok = a.notna() & b.notna()
    if ok.sum() < 3:
        return np.nan
    return a[ok].rank().corr(b[ok].rank())


# ------------------------------------------------------ C1 intensity

def intensity_structure(paired: pd.DataFrame) -> pd.DataFrame:
    """Grid/gauge ratio by rainfall intensity, binned symmetrically.

    Binning on the gauge value would put every day the gauges happened to
    under-read into the low bin and make the grid look worst there; binning on
    the grid value does the mirror image. The average of the two is not immune
    but it does not tilt the answer in a known direction.
    """
    d = paired.dropna(subset=["precip_gauge_mm", "precip_grid_mm"]).copy()
    d["avg"] = (d["precip_grid_mm"] + d["precip_gauge_mm"]) / 2
    rows = []
    for lo, hi in zip(INTENSITY_BINS[:-1], INTENSITY_BINS[1:]):
        s = d[(d["avg"] >= lo) & (d["avg"] < hi)]
        if s.empty:
            continue
        gauge, grid = s["precip_gauge_mm"].sum(), s["precip_grid_mm"].sum()
        rows.append({
            "bin_lo_mm": lo, "bin_hi_mm": hi, "n": len(s),
            "gauge_mean_mm": s["precip_gauge_mm"].mean(),
            "grid_mean_mm": s["precip_grid_mm"].mean(),
            "ratio": grid / gauge if gauge > 0 else np.nan,
            "excess_mm_total": grid - gauge,
        })
    return pd.DataFrame(rows)


def frequency_intensity_split(paired: pd.DataFrame,
                              thresholds=(0.2, 1.0, 5.0)) -> pd.DataFrame:
    """Does the grid rain too often, or too hard on the days it rains?"""
    d = paired.dropna(subset=["precip_gauge_mm", "precip_grid_mm"])
    rows = []
    for thr in thresholds:
        wet_g, wet_r = d["precip_gauge_mm"] >= thr, d["precip_grid_mm"] >= thr
        rows.append({
            "wet_day_threshold_mm": thr, "n": len(d),
            "wet_freq_gauge": wet_g.mean(), "wet_freq_grid": wet_r.mean(),
            "freq_ratio": wet_r.mean() / wet_g.mean(),
            "wet_intensity_gauge_mm": d.loc[wet_g, "precip_gauge_mm"].mean(),
            "wet_intensity_grid_mm": d.loc[wet_r, "precip_grid_mm"].mean(),
            "intensity_ratio": (d.loc[wet_r, "precip_grid_mm"].mean()
                                / d.loc[wet_g, "precip_gauge_mm"].mean()),
        })
    return pd.DataFrame(rows)


# ------------------------------------------------------ C2 geography

def city_table(totals: pd.DataFrame, quality: pd.DataFrame,
               meta: pd.DataFrame) -> pd.DataFrame:
    used = quality[quality["used"]][["city", "sid"]]
    m = meta.merge(used, on=["city", "sid"])
    agg = (m.groupby("city")
           .agg(gauge_elev_m=("elev_m", "mean"),
                cell_elev_m=("cell_elev_m", "first"),
                gauge_dist_km=("dist_km", "mean"))
           .reset_index())
    t = totals.merge(agg, on="city")
    t["excess_pc"] = 100 * (t["grid_over_station"] - 1)
    t["log_ratio"] = np.log(t["grid_over_station"])
    t["lat"] = t["city"].map(lambda c: CITIES[c][0])
    t["coastal"] = t["city"].isin(COASTAL)
    return t


def geography_test(t: pd.DataFrame) -> pd.DataFrame:
    rows = [{
        "test": "coastal vs interior",
        "group_a": "coastal", "n_a": int(t["coastal"].sum()),
        "median_a_pc": t.loc[t["coastal"], "excess_pc"].median(),
        "group_b": "interior", "n_b": int((~t["coastal"]).sum()),
        "median_b_pc": t.loc[~t["coastal"], "excess_pc"].median(),
        "spearman": np.nan,
    }]
    for var in ("lat", "gauge_elev_m", "cell_elev_m", "gauge_dist_km",
                "station_mm", "n_gauges"):
        rows.append({"test": f"excess vs {var}", "group_a": "", "n_a": len(t),
                     "median_a_pc": np.nan, "group_b": "", "n_b": np.nan,
                     "median_b_pc": np.nan,
                     "spearman": spearman(t["excess_pc"], t[var])})
    return pd.DataFrame(rows)


# ----------------------------------------------------- C3 noise floor

def noise_floor(t: pd.DataFrame, quality: pd.DataFrame) -> pd.DataFrame:
    """How much of the between-city spread is gauge sampling noise?

    A city's gauge estimate is the mean of `n` point measurements drawn from a
    cell whose true rainfall varies from corner to corner. The scatter among
    gauges in the cities that have several of them measures that variation; a
    city with one gauge inherits the same uncertainty without any way to see it,
    so it is charged the pooled coefficient of variation.
    """
    used = quality[quality["used"]]
    per_city = used.groupby("city")["total_mm"].agg(["count", "mean", "std"])
    per_city["cv"] = per_city["std"] / per_city["mean"]
    cv_pooled = per_city.loc[per_city["count"] >= 3, "cv"].median()

    t = t.merge(per_city[["count"]], left_on="city", right_index=True)
    se = cv_pooled / np.sqrt(t["count"])          # SE of log(gauge mean)
    observed_var = t["log_ratio"].var(ddof=1)
    noise_var = (se ** 2).mean()
    real_var = max(observed_var - noise_var, 0.0)
    return pd.DataFrame([{
        "n_cities": len(t),
        "pooled_within_cell_gauge_cv": cv_pooled,
        "observed_sd_log_ratio": np.sqrt(observed_var),
        "noise_sd_log_ratio": np.sqrt(noise_var),
        "residual_real_sd_log_ratio": np.sqrt(real_var),
        "residual_real_sd_pc": 100 * (np.exp(np.sqrt(real_var)) - 1),
        "noise_share_of_variance": noise_var / observed_var,
    }])


# --------------------------------------------------------- C4 timing

def timing(paired: pd.DataFrame, lags=(-2, -1, 0, 1, 2)) -> pd.DataFrame:
    rows = []
    for city, g in paired.groupby("city"):
        g = g.sort_values("valid_date")
        row = {"city": city, "n_days": len(g)}
        for lag in lags:
            row[f"r_lag{lag:+d}"] = (g["precip_gauge_mm"].shift(-lag)
                                     .corr(g["precip_grid_mm"]))
        best = max(lags, key=lambda k: (row[f"r_lag{k:+d}"]
                                        if pd.notna(row[f"r_lag{k:+d}"])
                                        else -np.inf))
        row["best_lag_days"] = best
        row["gain_over_lag0"] = row[f"r_lag{best:+d}"] - row["r_lag+0"]
        rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------- figure 10

def fig_bias_structure(inten: pd.DataFrame, t: pd.DataFrame,
                       floor: pd.DataFrame, quality: pd.DataFrame,
                       outdir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6),
                             gridspec_kw={"width_ratios": [1, 1.75]})

    # Left: the excess is multiplicative and present at every intensity.
    ax = axes[0]
    labels = [f"{int(r.bin_lo_mm)}–{'' if np.isinf(r.bin_hi_mm) else int(r.bin_hi_mm)}"
              for r in inten.itertuples()]
    labels[-1] = f"≥{int(inten['bin_lo_mm'].iloc[-1])}"
    x = np.arange(len(inten))
    ax.bar(x, inten["ratio"], color=GRID_COLOR, width=0.62)
    ax.axhline(1, color=INK, linewidth=1.2, linestyle="--")
    for i, r in enumerate(inten.itertuples()):
        ax.text(i, r.ratio + 0.04, f"×{r.ratio:.2f}", ha="center", fontsize=8,
                color=INK)
        ax.text(i, 0.06, f"n={r.n}", ha="center", fontsize=8, color=SURFACE_TEXT)
    ax.set_xticks(x, labels)
    ax.set_xlabel("Daily rainfall, mean of the two series (mm)")
    ax.set_ylabel("Gridded ÷ gauges")
    ax.set_ylim(0, max(2.6, inten["ratio"].max() * 1.2))
    ax.set_title("Wet at every intensity,\nworst on light-rain days", fontsize=10)

    # Right: the per-city excess with its sampling uncertainty, ordered north
    # to south, coloured by coast vs interior. If geography explained the
    # spread, the two colours would separate. They do not.
    ax = axes[1]
    used = quality[quality["used"]]
    cv = used.groupby("city")["total_mm"].agg(["count", "mean", "std"])
    cv["cv"] = cv["std"] / cv["mean"]
    cv_pooled = float(floor["pooled_within_cell_gauge_cv"].iloc[0])
    order = [c for c in CITY_ORDER if c in set(t["city"])]
    tt = t.set_index("city").loc[order]
    y = np.arange(len(order))
    se_pc = 100 * (cv_pooled / np.sqrt(cv.loc[order, "count"].to_numpy())
                   * (1 + tt["excess_pc"].to_numpy() / 100))
    colors = [STATION_COLOR if tt.loc[c, "coastal"] else GRID_COLOR
              for c in order]
    ax.barh(y, tt["excess_pc"], color=colors, height=0.66)
    ax.errorbar(tt["excess_pc"], y, xerr=se_pc, fmt="none", ecolor=INK,
                elinewidth=1.1, capsize=3, alpha=0.75)
    ax.axvline(0, color=INK, linewidth=1.2)
    pooled = 100 * (t["grid_mm"].sum() / t["station_mm"].sum() - 1)
    ax.axvline(pooled, color=INK, linewidth=1.2, linestyle="--")
    ax.text(pooled, len(order) - 0.2, f"transect {pooled:+.0f}%", fontsize=8,
            color=INK, ha="center", va="top")
    ax.set_yticks(y, [CITY_LABELS[c] for c in order])
    ax.invert_yaxis()
    ax.set_xlabel("Gridded analysis relative to the gauges, event total (%)")
    handles = [plt.Rectangle((0, 0), 1, 1, color=STATION_COLOR),
               plt.Rectangle((0, 0), 1, 1, color=GRID_COLOR)]
    ax.legend(handles, ["cell contains ocean", "inland cell"],
              loc="lower right", fontsize=9)
    share = float(floor["noise_share_of_variance"].iloc[0])
    ax.set_title(f"No geographic ordering — and {100 * share:.0f}% of the "
                 f"city-to-city variance\nis the sampling noise of the gauge "
                 f"estimate (error bars)", fontsize=10)

    fig.tight_layout()
    fig.savefig(outdir / "fig10_bias_structure.png", bbox_inches="tight")
    plt.close(fig)


SURFACE_TEXT = "#4a4945"


# ------------------------------------------------------------- main

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--capture", default="2026-07-26")
    args = parser.parse_args()

    proc = ROOT / "data" / "processed" / args.capture
    metrics = proc / "metrics"
    paired = pd.read_csv(metrics / "station_grid_paired_daily.csv")
    totals = pd.read_csv(metrics / "event_totals_threeway.csv")
    gauges = load_gauges(args.capture)
    quality = gauge_quality(gauges)
    quality["sid"] = quality["station_id"].astype(str)
    window = event_window(city_gauge_daily(gauges, quality))

    meta = pd.concat([
        pd.read_csv(ROOT / "data" / "external" / "stations_metadata.csv"),
        pd.read_csv(ROOT / "data" / "external" / "stations_metadata_dmc.csv")])
    meta["sid"] = meta["station_id"].astype(str)

    inten = intensity_structure(paired)
    split = frequency_intensity_split(paired)
    t = city_table(totals, quality, meta)
    geo = geography_test(t)
    floor = noise_floor(t, quality)
    lags = timing(paired)

    for name, table in [("bias_structure_intensity", inten),
                        ("bias_structure_freq_intensity", split),
                        ("bias_structure_by_city", t),
                        ("bias_structure_geography", geo),
                        ("bias_structure_noise_floor", floor),
                        ("bias_structure_timing", lags)]:
        table.to_csv(metrics / f"{name}.csv", index=False)
    print(f"metrics -> {metrics} (6 CSV files)")

    set_style()
    fig_bias_structure(inten, t, floor, quality, ROOT / "figures")
    print(f"figures -> {ROOT / 'figures'} (fig10)")

    print(f"\nevent window {window[0]} .. {window[1]}")
    print("intensity structure (ratio grid/gauge by mean daily amount):")
    for r in inten.itertuples():
        hi = "+" if np.isinf(r.bin_hi_mm) else f"{r.bin_hi_mm:g}"
        print(f"  {r.bin_lo_mm:g}-{hi:>3s} mm  n={r.n:4d}  ratio {r.ratio:.2f}")
    s = split[split["wet_day_threshold_mm"] == 1.0].iloc[0]
    print(f"wet days (>=1 mm): grid rains {100 * (s.freq_ratio - 1):+.0f}% more "
          f"often and {100 * (s.intensity_ratio - 1):+.0f}% harder when it does")
    g = geo.iloc[0]
    print(f"coastal median {g.median_a_pc:+.0f}% (n={g.n_a})  vs  "
          f"interior median {g.median_b_pc:+.0f}% (n={g.n_b})")
    print("rank correlation of the excess with:")
    for r in geo.iloc[1:].itertuples():
        print(f"  {r.test:28s} {r.spearman:+.2f}")
    f = floor.iloc[0]
    print(f"between-city SD of log(grid/gauge) {f.observed_sd_log_ratio:.3f}; "
          f"gauge-sampling noise alone {f.noise_sd_log_ratio:.3f} "
          f"({100 * f.noise_share_of_variance:.0f}% of the variance); "
          f"residual real spread ±{f.residual_real_sd_pc:.0f}%")
    off = lags[lags["best_lag_days"] != 0]
    print(f"cities whose daily series correlate best at a non-zero lag: "
          f"{len(off)} of {len(lags)}")
    for _, r in off.iterrows():
        best = int(r["best_lag_days"])
        print(f"  {r['city']}: best lag {best:+d} d "
              f"(r {r['r_lag+0']:.2f} -> {r[f'r_lag{best:+d}']:.2f})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
