"""Forecast verification metrics and figures for the July 2026 Chile frontal system.

Reads the tidy dataset produced by build_dataset.py and computes:

  1. Continuous scores by model x lead: MAE, bias, RMSE for daily precipitation,
     MAE for daily max wind gusts (models that provide gusts only).
  2. Categorical scores by model x lead x threshold (precip >= 1, 10, 20 mm):
     POD, FAR, CSI, base rate, from the 2x2 contingency table.
  3. Poor man's ensemble (GFS + ECMWF + ICON as equally likely members):
     exceedance probabilities in {0, 1/3, 2/3, 1}, Brier score and skill score
     vs climatology (sample base rate), reliability table, and ensemble CRPS.

Outputs metrics as CSV to data/processed/<capture>/metrics/ and five PNG
figures to figures/.

Usage: uv run scripts/verify.py [--capture YYYY-MM-DD]
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, PowerNorm

ROOT = Path(__file__).resolve().parent.parent

# North -> south transect order (matches capture_forecasts.py).
CITY_ORDER = [
    "la_serena", "coquimbo", "vicuna", "andacollo", "ovalle", "monte_patria",
    "combarbala", "illapel", "salamanca", "los_vilos",
    "valparaiso", "santiago", "rancagua", "curico", "talca",
    "chillan", "concepcion", "temuco", "valdivia", "osorno", "puerto_montt",
]
CITY_LABELS = {
    "la_serena": "La Serena", "coquimbo": "Coquimbo",
    "vicuna": "Vicuña", "andacollo": "Andacollo",
    "ovalle": "Ovalle", "monte_patria": "Monte Patria",
    "combarbala": "Combarbalá", "illapel": "Illapel", "salamanca": "Salamanca",
    "los_vilos": "Los Vilos",
    "valparaiso": "Valparaíso", "santiago": "Santiago",
    "rancagua": "Rancagua", "curico": "Curicó", "talca": "Talca",
    "chillan": "Chillán", "concepcion": "Concepción", "temuco": "Temuco",
    "valdivia": "Valdivia", "osorno": "Osorno", "puerto_montt": "Puerto Montt",
}

# Fixed model -> color assignment (validated categorical palette, light mode).
MODEL_COLORS = {
    "best_match": "#2a78d6",
    "gfs_seamless": "#008300",
    "ecmwf_ifs025": "#e87ba4",
    "icon_seamless": "#eda100",
}
MODEL_LABELS = {
    "best_match": "best_match",
    "gfs_seamless": "GFS",
    "ecmwf_ifs025": "ECMWF IFS",
    "icon_seamless": "ICON",
}
ENSEMBLE_MEMBERS = ["gfs_seamless", "ecmwf_ifs025", "icon_seamless"]

PRECIP_THRESHOLDS_MM = [1.0, 10.0, 20.0]

# The three arbiters, and what each one actually is. No score in this project
# has a single truth any more: every table carries a `truth` column and every
# conclusion is only as strong as its agreement across the ones available.
#
#   grid   Open-Meteo archive best_match. Despite years of this pipeline
#          calling it ERA5, it is ECMWF IFS HRES at 9 km for 2026 -- the
#          operational analysis of one of the models being verified.
#   era5   True ERA5 at 0.25 deg. Still an ECMWF product, but a reanalysis
#          rather than the operational analysis, and coarser.
#   gauge  Rain gauges, averaged over every instrument within 30 km of the city
#          that passes QC. The only arbiter not made by a model, and the only
#          one that does not exist everywhere.
TRUTHS = {
    "grid": "precip_grid_mm",
    "era5": "precip_era5_mm",
    "gauge": "precip_gauge_mm",
}
TRUTH_LABELS = {
    "grid": "IFS HRES 9 km (archive best_match)",
    "era5": "ERA5 0.25°",
    "gauge": "rain gauges",
}
# Short forms, for titles and legends where the full label does not fit.
TRUTH_SHORT = {"grid": "gridded analysis", "era5": "ERA5", "gauge": "gauges"}


def attach_gauge_truth(tidy: pd.DataFrame, capture: str) -> pd.DataFrame:
    """Add the city gauge estimate to the tidy table, where one exists."""
    import stations_common as sc

    try:
        daily = sc.load_gauges(capture)
    except FileNotFoundError:
        return tidy
    est = sc.city_gauge_daily(daily, sc.gauge_quality(daily))
    return tidy.merge(est[["city", "valid_date", "precip_gauge_mm"]],
                      on=["city", "valid_date"], how="left")


def by_truth(tidy: pd.DataFrame):
    """Yield (truth_name, frame) with that truth renamed to `precip_obs_mm`.

    The scoring functions below are written against a single column name and
    stay that way: pointing them at a different arbiter is a rename here, not an
    edit in thirty places. Rows where the chosen arbiter has no value are
    dropped, so each truth is scored on exactly the city-days it covers -- which
    is why every output table also carries its own `n`.
    """
    for name, column in TRUTHS.items():
        if column not in tidy.columns or tidy[column].isna().all():
            continue
        frame = tidy.dropna(subset=[column]).copy()
        frame["precip_obs_mm"] = frame[column]
        # Gusts exist only for the gridded arbiters; a gauge network measuring
        # rainfall has nothing to say about them, so those scores come out NaN
        # rather than silently borrowing a different arbiter's wind.
        gust = f"gust_{name}_kmh"
        frame["gust_obs_kmh"] = (frame[gust] if gust in frame.columns
                                 else np.nan)
        yield name, frame


def with_truth(tables: dict, truth: str) -> dict:
    return {k: v.assign(truth=truth) for k, v in tables.items()}

SURFACE = "#fcfcfb"
GRID = "#e1e0d9"
AXIS = "#898781"
INK = "#0b0b0b"
SEQ_BLUES = LinearSegmentedColormap.from_list(
    "seq_blue",
    ["#fcfcfb", "#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
     "#256abf", "#184f95", "#0d366b"],
)


def set_style() -> None:
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE, "savefig.dpi": 150,
        "font.family": "sans-serif", "font.size": 10,
        "text.color": INK, "axes.labelcolor": INK,
        "axes.edgecolor": AXIS, "xtick.color": AXIS, "ytick.color": AXIS,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.titlesize": 11, "axes.titleweight": "bold",
        "lines.linewidth": 2, "legend.frameon": False,
    })


# ---------------------------------------------------------------- metrics

def continuous_scores(tidy: pd.DataFrame) -> pd.DataFrame:
    def scores(g: pd.DataFrame) -> pd.Series:
        perr = g["precip_fcst_mm"] - g["precip_obs_mm"]
        gerr = (g["gust_fcst_kmh"] - g["gust_obs_kmh"]).dropna()
        return pd.Series({
            "n": len(g),
            "precip_mae_mm": perr.abs().mean(),
            "precip_bias_mm": perr.mean(),
            "precip_rmse_mm": np.sqrt((perr ** 2).mean()),
            "gust_mae_kmh": gerr.abs().mean() if len(gerr) else np.nan,
        })

    return (tidy.groupby(["model", "lead_days"])
            .apply(scores, include_groups=False).reset_index())


def categorical_scores(tidy: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for thr in PRECIP_THRESHOLDS_MM:
        fcst = tidy["precip_fcst_mm"] >= thr
        obs = tidy["precip_obs_mm"] >= thr
        df = tidy[["model", "lead_days"]].assign(
            hit=fcst & obs, miss=~fcst & obs, false_alarm=fcst & ~obs, obs=obs)
        agg = df.groupby(["model", "lead_days"]).agg(
            n=("obs", "size"), hits=("hit", "sum"), misses=("miss", "sum"),
            false_alarms=("false_alarm", "sum"), base_rate=("obs", "mean"))
        agg["pod"] = agg["hits"] / (agg["hits"] + agg["misses"])
        agg["far"] = (agg["false_alarms"]
                      / (agg["hits"] + agg["false_alarms"]).replace(0, np.nan))
        agg["csi"] = agg["hits"] / (agg["hits"] + agg["misses"]
                                    + agg["false_alarms"])
        rows.append(agg.reset_index().assign(threshold_mm=thr))
    return pd.concat(rows, ignore_index=True)


def ensemble_table(tidy: pd.DataFrame) -> pd.DataFrame:
    """One row per city x valid_date x lead with the member forecasts as columns."""
    members = tidy[tidy["model"].isin(ENSEMBLE_MEMBERS)]
    wide = members.pivot_table(
        index=["city", "valid_date", "lead_days", "precip_obs_mm"],
        columns="model", values="precip_fcst_mm").dropna().reset_index()
    return wide


def ensemble_brier(wide: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Brier score / skill score by lead x threshold, plus reliability table."""
    brier_rows, rel_rows = [], []
    for thr in PRECIP_THRESHOLDS_MM:
        prob = (wide[ENSEMBLE_MEMBERS] >= thr).mean(axis=1)
        event = (wide["precip_obs_mm"] >= thr).astype(float)
        df = pd.DataFrame({"lead_days": wide["lead_days"],
                           "prob": prob, "event": event})
        for lead, g in df.groupby("lead_days"):
            bs = ((g["prob"] - g["event"]) ** 2).mean()
            base = g["event"].mean()
            bs_ref = (base * (1 - base))  # climatology Brier score
            brier_rows.append({
                "threshold_mm": thr, "lead_days": lead, "n": len(g),
                "base_rate": base, "brier": bs, "brier_climatology": bs_ref,
                "bss": 1 - bs / bs_ref if bs_ref > 0 else np.nan,
            })
        # Reliability pooled over leads: 4 possible forecast probabilities.
        rel = df.groupby("prob").agg(n=("event", "size"),
                                     obs_freq=("event", "mean")).reset_index()
        rel_rows.append(rel.assign(threshold_mm=thr))
    return pd.DataFrame(brier_rows), pd.concat(rel_rows, ignore_index=True)


def ensemble_crps(wide: pd.DataFrame) -> pd.DataFrame:
    """Empirical CRPS of the 3-member ensemble, averaged by lead."""
    x = wide[ENSEMBLE_MEMBERS].to_numpy()          # (n, m)
    y = wide["precip_obs_mm"].to_numpy()[:, None]  # (n, 1)
    term1 = np.abs(x - y).mean(axis=1)
    term2 = np.abs(x[:, :, None] - x[:, None, :]).mean(axis=(1, 2))
    crps = term1 - 0.5 * term2
    df = pd.DataFrame({"lead_days": wide["lead_days"], "crps": crps,
                       "abs_err_best_single": np.abs(x - y).min(axis=1)})
    return (df.groupby("lead_days")
            .agg(n=("crps", "size"), crps_mm=("crps", "mean")).reset_index())


# ---------------------------------------------------------------- figures

def fig_event_overview(tidy: pd.DataFrame, outdir: Path,
                       window: tuple[str, str] | None = None) -> None:
    """City x day heatmap of what the rain gauges measured.

    Earlier versions drew the gridded analysis here, which made the overview
    figure disagree with every gauge total quoted later in the report by up to
    70%. What fell is what the instruments caught; the gridded field is a model
    product and is shown as such, in the figures that compare the two.
    """
    obs = (tidy[["city", "valid_date", "precip_obs_mm"]]
           .drop_duplicates()
           .pivot(index="city", columns="valid_date", values="precip_obs_mm"))
    obs = obs.reindex([c for c in CITY_ORDER if c in obs.index])
    fig, ax = plt.subplots(figsize=(9, 1.2 + 0.28 * len(obs.index)))
    # Square-root colour scale: daily rainfall is strongly right-skewed, and on
    # a linear scale one 174 mm cell flattens the whole frontal passage into
    # pale blue. The tick labels stay in mm.
    mesh = ax.pcolormesh(np.arange(obs.shape[1] + 1), np.arange(obs.shape[0] + 1),
                         obs.to_numpy(), cmap=SEQ_BLUES,
                         norm=PowerNorm(0.5, vmin=0,
                                        vmax=float(np.nanmax(obs.to_numpy()))),
                         edgecolors=SURFACE, linewidth=2)
    ax.set_yticks(np.arange(len(obs.index)) + 0.5,
                  [CITY_LABELS[c] for c in obs.index])
    ax.set_xticks(np.arange(obs.shape[1]) + 0.5,
                  [d.strftime("%d %b") for d in obs.columns], rotation=45,
                  ha="right")
    ax.invert_yaxis()  # north at the top
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(mesh, ax=ax, label="Measured precipitation (mm/day)")
    cbar.outline.set_visible(False)
    ax.set_title("July 2026 frontal system: daily precipitation measured by "
                 "rain gauges\n(city mean of the CEAZA-Met and DMC gauges "
                 "within 30 km)")
    if window is not None:
        # Mark the event window on the date axis, since every event total in
        # the report is summed over exactly these columns.
        cols = [str(pd.Timestamp(d).date()) for d in obs.columns]
        lo = min(i for i, d in enumerate(cols) if d >= window[0])
        hi = max(i for i, d in enumerate(cols) if d <= window[1])
        ax.plot([lo, hi + 1], [len(obs.index) + 0.35] * 2, color=INK,
                linewidth=2.5, clip_on=False)
        ax.text((lo + hi + 1) / 2, len(obs.index) + 0.95, "event window",
                ha="center", va="center", fontsize=8, color=INK)
    fig.tight_layout()
    fig.savefig(outdir / "fig1_event_overview.png")
    plt.close(fig)


def _lead_axis(ax) -> None:
    ax.set_xticks(list(range(1, 8)))
    ax.set_xlabel("Forecast lead (days)")


def fig_error_vs_lead(cont: pd.DataFrame, outdir: Path,
                      truths: tuple[str, ...] = ("grid", "gauge")) -> None:
    """MAE and bias against lead time, one row per arbiter.

    Precipitation MAE and bias are the two scores that moved most when the
    arbiter changed, so they are drawn against both -- side by side, at a shared
    y scale, because the comparison is the result.
    """
    truths = tuple(t for t in truths if t in set(cont["truth"]))
    panels = [("precip_mae_mm", "MAE, precipitation (mm/day)"),
              ("precip_bias_mm", "Bias, precipitation (mm/day)"),
              ("gust_mae_kmh", "MAE, max wind gusts (km/h)")]
    fig, axes = plt.subplots(len(truths), 3, figsize=(11, 3.5 * len(truths)),
                             sharex=True, squeeze=False)
    for row, truth in enumerate(truths):
        sub = cont[cont["truth"] == truth]
        for ax, (col, title) in zip(axes[row], panels):
            for model, g in sub.groupby("model"):
                g = g.dropna(subset=[col])
                if g.empty:
                    continue
                ax.plot(g["lead_days"], g[col], marker="o", markersize=5,
                        color=MODEL_COLORS[model], label=MODEL_LABELS[model])
            _lead_axis(ax)
            if row == 0:
                ax.set_title(title)
        axes[row][1].axhline(0, color=AXIS, linewidth=1, linestyle="--")
        n = int(sub["n"].max()) if "n" in sub.columns else None
        axes[row][0].set_ylabel(f"vs {TRUTH_LABELS[truth]}", fontweight="bold")
        if n:
            axes[row][0].text(0.03, 0.94, f"n = {n} forecast-days per lead",
                              transform=axes[row][0].transAxes, fontsize=8,
                              color=AXIS, va="top")
        if sub["gust_mae_kmh"].isna().all():
            ax = axes[row][2]
            ax.text(0.5, 0.5, "rain gauges do not\nmeasure wind gusts",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=9, color=AXIS, style="italic")
            ax.set_yticks([])
    # Same y range per column, or the two rows cannot be compared by eye.
    for col in range(3):
        drawn = [r for r in range(len(truths)) if axes[r][col].has_data()]
        if len(drawn) < 2:
            continue
        lims = [axes[r][col].get_ylim() for r in drawn]
        lo, hi = min(l for l, _ in lims), max(h for _, h in lims)
        for r in drawn:
            axes[r][col].set_ylim(lo, hi)
    axes[0][0].legend(loc="lower right", fontsize=9)
    axes[0][2].set_xlabel("")
    axes[0][2].text(0.03, 0.03, "ECMWF previous runs do not include gusts",
                    transform=axes[0][2].transAxes, fontsize=8, color=AXIS)
    fig.suptitle("Forecast error grows with lead time — and how much of it is "
                 "the model depends on the arbiter", fontweight="bold")
    fig.tight_layout()
    fig.savefig(outdir / "fig2_error_vs_lead.png")
    plt.close(fig)


def fig_event_total_vs_lead(tidy: pd.DataFrame, outdir: Path,
                            event_dates: list, cities: list[str]) -> None:
    """Event-total forecasts converging on what the gauges measured.

    Both reference lines are drawn. The gauges are what fell and are the line a
    forecast should be judged against; the gridded analysis is drawn beside them
    because the distance between the two is the size of the arbiter problem, and
    a reader who only saw the gridded line would conclude the models were drier
    than they were.
    """
    ev = tidy[tidy["valid_date"].isin(event_dates)]
    fig, axes = plt.subplots(1, len(cities), figsize=(11.5, 4.0), sharex=True)
    for ax, city in zip(axes, cities):
        g = ev[ev["city"] == city]
        days = g[["valid_date", "precip_gauge_mm", "precip_grid_mm"]
                 ].drop_duplicates()
        gauge_total = days["precip_gauge_mm"].sum()
        grid_total = days["precip_grid_mm"].sum()
        for model, gm in g.groupby("model"):
            tot = gm.groupby("lead_days")["precip_fcst_mm"].sum()
            ax.plot(tot.index, tot.to_numpy(), marker="o", markersize=5,
                    color=MODEL_COLORS[model], label=MODEL_LABELS[model])
        ax.axhline(gauge_total, color=INK, linewidth=1.8, linestyle="--")
        ax.text(7.1, gauge_total, f"gauges\n{gauge_total:.0f}", fontsize=8,
                color=INK, va="center")
        ax.axhline(grid_total, color=AXIS, linewidth=1.4, linestyle=":")
        ax.text(7.1, grid_total, f"grid\n{grid_total:.0f}", fontsize=8,
                color=AXIS, va="center")
        ax.set_xlim(0.5, 8.6)
        _lead_axis(ax)
        ax.set_title(CITY_LABELS[city])
    axes[0].set_ylabel("Event-total precipitation (mm)")
    axes[0].legend(loc="lower left", fontsize=8)
    dates_str = (f"{min(event_dates).strftime('%d')}–"
                 f"{max(event_dates).strftime('%d %b %Y')}")
    fig.suptitle(f"How the event total ({dates_str}) looked N days ahead",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(outdir / "fig3_event_total_vs_lead.png")
    plt.close(fig)


def fig_categorical(cat: pd.DataFrame, outdir: Path, thr: float = 10.0,
                    truths: tuple[str, ...] = ("gauge", "grid")) -> None:
    truths = tuple(t for t in truths if t in set(cat["truth"]))
    sel = cat[cat["threshold_mm"] == thr]
    panels = [("pod", "POD (hit rate)"), ("far", "FAR (false-alarm ratio)"),
              ("csi", "CSI (critical success index)")]
    fig, axes = plt.subplots(len(truths), 3, figsize=(11, 3.5 * len(truths)),
                             sharex=True, sharey=True, squeeze=False)
    for row, truth in enumerate(truths):
        sub = sel[sel["truth"] == truth]
        for ax, (col, title) in zip(axes[row], panels):
            for model, g in sub.groupby("model"):
                ax.plot(g["lead_days"], g[col], marker="o", markersize=5,
                        color=MODEL_COLORS[model], label=MODEL_LABELS[model])
            _lead_axis(ax)
            ax.set_ylim(0, 1)
            if row == 0:
                ax.set_title(title)
        base = sub["base_rate"].iloc[0] if len(sub) else np.nan
        axes[row][0].set_ylabel(f"vs {TRUTH_LABELS[truth]}", fontweight="bold")
        axes[row][1].text(0.03, 0.92, f"base rate {base:.3f}",
                          transform=axes[row][1].transAxes, fontsize=8,
                          color=AXIS, va="top")
    axes[0][0].legend(loc="lower left", fontsize=8)
    fig.suptitle(f"Detection of heavy-rain days (≥ {thr:.0f} mm) vs lead time",
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig(outdir / "fig4_categorical_skill.png")
    plt.close(fig)


TRUTH_MARKERS = {"gauge": ("o", "#0d366b"), "grid": ("s", "#6da7ec"),
                 "era5": ("^", "#eda100")}


def fig_reliability(rel: pd.DataFrame, brier: pd.DataFrame, outdir: Path,
                    truths: tuple[str, ...] = ("gauge", "grid")) -> None:
    """Reliability of the same forecasts read against two different arbiters.

    Both curves come from one set of ensemble probabilities; only the definition
    of "it rained" changes between them. Drawing them on one panel is the point:
    the ordering survives the change of judge and the calibration does not.
    """
    truths = tuple(t for t in truths if t in set(rel["truth"]))
    thresholds = [1.0, 10.0]
    fig, axes = plt.subplots(1, len(thresholds), figsize=(9.5, 4.6), sharey=True)
    for ax, thr in zip(axes, thresholds):
        ax.plot([0, 1], [0, 1], color=AXIS, linewidth=1, linestyle="--")
        bits = []
        for truth in truths:
            g = rel[(rel["threshold_mm"] == thr) & (rel["truth"] == truth)]
            marker, color = TRUTH_MARKERS[truth]
            sizes = 40 + 260 * g["n"] / g["n"].max()
            ax.scatter(g["prob"], g["obs_freq"], s=sizes, color=color,
                       marker=marker, zorder=3, edgecolors=SURFACE,
                       linewidth=2, label=TRUTH_LABELS[truth])

            ax.plot(g["prob"], g["obs_freq"], color=color, linewidth=1.2,
                    zorder=2)
            if truth == truths[0]:
                for _, r in g.iterrows():
                    ax.annotate(f"n={r['n']:.0f}", (r["prob"], r["obs_freq"]),
                                textcoords="offset points", xytext=(9, -14),
                                fontsize=8, color=AXIS)
            b = brier[(brier["threshold_mm"] == thr) & (brier["truth"] == truth)]
            if len(b):
                bs = (b["brier"] * b["n"]).sum() / b["n"].sum()
                bits.append(f"{TRUTH_SHORT[truth]} {bs:.3f}")
        ax.set_title(f"≥ {thr:.0f} mm/day\nBrier: " + ", ".join(bits),
                     fontsize=10)
        ax.set_xlabel("Forecast probability (3-model ensemble)")
        ax.set_xticks([0, 1 / 3, 2 / 3, 1], ["0", "1/3", "2/3", "1"])
        ax.set_xlim(-0.08, 1.08)
        ax.set_ylim(-0.05, 1.05)
    axes[0].set_ylabel("Observed frequency")
    axes[0].legend(loc="upper left", fontsize=9)
    fig.suptitle("The poor man's ensemble is reliable — and how reliable "
                 "depends on the arbiter", fontweight="bold")
    fig.tight_layout()
    fig.savefig(outdir / "fig5_reliability.png")
    plt.close(fig)


# ---------------------------------------------------------------- main

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--capture", default=None,
                        help="capture date YYYY-MM-DD (default: latest)")
    args = parser.parse_args()

    proc_root = ROOT / "data" / "processed"
    capture = args.capture or sorted(
        p.name for p in proc_root.iterdir() if p.is_dir())[-1]
    tidy = pd.read_parquet(proc_root / capture / "tidy_daily.parquet")
    tidy = attach_gauge_truth(tidy, capture)

    metrics_dir = proc_root / capture / "metrics"
    metrics_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = ROOT / "figures"
    fig_dir.mkdir(exist_ok=True)

    # Every score, once per arbiter. The tables are stacked rather than written
    # side by side so that any of them can be read with a single `truth` filter.
    parts: dict[str, list] = {}
    coverage = []
    for truth, frame in by_truth(tidy):
        wide = ensemble_table(frame)
        brier, rel = ensemble_brier(wide)
        tables = {
            "continuous_by_lead": continuous_scores(frame),
            "categorical_by_lead": categorical_scores(frame),
            "ensemble_brier": brier,
            "ensemble_reliability": rel,
            "ensemble_crps": ensemble_crps(wide),
        }
        for name, table in with_truth(tables, truth).items():
            parts.setdefault(name, []).append(table)
        coverage.append((truth, frame["city"].nunique(), len(frame)))

    for name, tabs in parts.items():
        pd.concat(tabs, ignore_index=True).to_csv(
            metrics_dir / f"{name}.csv", index=False)
    print(f"metrics -> {metrics_dir} ({len(parts)} CSV files)")
    for truth, ncity, nrow in coverage:
        print(f"  truth={truth:6s} {TRUTH_LABELS[truth]:38s} "
              f"{ncity:2d} cities, {nrow} forecast-days")

    # One event window for the whole project, derived from the gauges rather
    # than hand-set, and covering the frontal passage end to end. Earlier
    # versions cut the case-study figures to a 3-day core (16-18 July) while
    # the gauge tables used a wider window, so the same city carried two
    # different "event totals" in one report.
    import stations_common as sc
    est = tidy[["city", "valid_date", "precip_gauge_mm"]
               ].drop_duplicates().dropna(subset=["precip_gauge_mm"])
    window = sc.event_window(est)
    event_dates = sorted({d for d in tidy["valid_date"].unique()
                          if window[0] <= str(pd.Timestamp(d).date())
                          <= window[1]})
    print(f"event window (majority-wet run around the peak): "
          f"{window[0]} .. {window[1]} ({len(event_dates)} days)")

    cont = pd.concat(parts["continuous_by_lead"], ignore_index=True)
    cat = pd.concat(parts["categorical_by_lead"], ignore_index=True)
    brier = pd.concat(parts["ensemble_brier"], ignore_index=True)
    rel = pd.concat(parts["ensemble_reliability"], ignore_index=True)

    # Figures show what the instruments measured. The gridded field is a model
    # product; where it still appears it is labelled as one and drawn against
    # the gauges, never in place of them.
    gauge = tidy.assign(precip_obs_mm=tidy["precip_gauge_mm"])

    set_style()
    fig_event_overview(gauge, fig_dir, window)
    fig_error_vs_lead(cont, fig_dir)
    fig_event_total_vs_lead(tidy, fig_dir, event_dates,
                            ["la_serena", "ovalle", "santiago", "concepcion"])
    fig_categorical(cat, fig_dir)
    fig_reliability(rel, brier, fig_dir)
    print(f"figures -> {fig_dir} (5 PNG files)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
