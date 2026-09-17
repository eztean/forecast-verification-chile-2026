"""Empirical error bands and coverage for the July 2026 Chile frontal system.

Answers the question the deterministic model output never states: a forecast of
80 mm is 80 mm plus or minus *what*, and did the observation land inside?

Everything here is IN-SAMPLE: the bands are fitted and evaluated on the same
event. Pooled coverage of an empirical quantile band is nominal by construction
and is not a finding; what is informative is the *width* of the band, its
asymmetry, and where it fails when evaluated on strata that did not define it.
Out-of-sample calibration is a separate project.

Sign convention: the error is e = observed - forecast, so that a band is applied
around a forecast as [f + q_lo, f + q_hi]. NOTE this is the opposite sign of
`precip_bias_mm` in verify.py, which reports forecast - observed.

Computes:

  B1  Coverage of +/- k*MAE bands (k = 1, 2) per model x lead, with Wilson 95%
      intervals, against the Gaussian reference (for a normal error, MAE =
      sigma*sqrt(2/pi), so |e| <= 1*MAE covers 57.5% and 2*MAE covers 88.9% --
      not the 68/95 that belong to +/- sigma).
  B2  Empirical quantile bands (q05, q16, q50, q84, q95) per model x lead, plus
      coverage of those bands evaluated on strata that did not define them
      (forecast-amount bin, city group).
  B3  Bands conditioned on the forecast amount -- the object closest to
      "80 +/- 20": the band around a large forecast is far wider and more
      asymmetric than the marginal one.
  B4  Rank histogram of the 3-member poor man's ensemble (ties randomised,
      fixed seed) and the spread-error relationship.
  B5  figures/fig8_uncertainty_cone.png -- the forecast for the peak day as it
      converged from lead 7 to lead 1, with its band.

Outputs CSVs to data/processed/<capture>/metrics/ and one PNG to figures/.

Usage: uv run scripts/coverage_analysis.py [--capture YYYY-MM-DD]
"""

import argparse
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify import (  # noqa: E402  (path set above)
    AXIS,
    CITY_LABELS,
    ENSEMBLE_MEMBERS,
    INK,
    MODEL_COLORS,
    MODEL_LABELS,
    attach_gauge_truth,
    by_truth,
    set_style,
)

ROOT = Path(__file__).resolve().parent.parent

# Forecast-amount strata. The band a user actually needs depends on how much
# rain is being forecast, so this is the stratification that matters most.
FCST_BINS = [-0.01, 0.05, 10.0, np.inf]
FCST_BIN_LABELS = ["f = 0", "0 < f < 10", "f >= 10"]

# City strata: the event peaked in the interior valleys of Coquimbo, where
# sub-grid orography is strongest, so coast and interior are kept apart.
CITY_GROUPS = {
    "la_serena": "Coquimbo coast", "los_vilos": "Coquimbo coast",
    "vicuna": "Coquimbo interior", "andacollo": "Coquimbo interior",
    "ovalle": "Coquimbo interior", "monte_patria": "Coquimbo interior",
    "combarbala": "Coquimbo interior", "illapel": "Coquimbo interior",
    "salamanca": "Coquimbo interior",
    "valparaiso": "Centre-south", "santiago": "Centre-south",
    "rancagua": "Centre-south", "curico": "Centre-south",
    "talca": "Centre-south", "chillan": "Centre-south",
    "concepcion": "Centre-south", "temuco": "Centre-south",
    "valdivia": "Centre-south", "osorno": "Centre-south",
    "puerto_montt": "Centre-south",
}

LEAD_GROUPS = {"1-3": [1, 2, 3], "4-7": [4, 5, 6, 7]}

RANK_SEED = 20260726  # reproducible tie-breaking in the rank histogram


# ------------------------------------------------------------------ helpers

def wilson(k: int, n: int, z: float = 1.959964) -> tuple[float, float]:
    """Wilson 95% interval for a binomial proportion (k successes of n)."""
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z / denom * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, centre - half), min(1.0, centre + half))


def norm_cdf(x: float) -> float:
    return 0.5 * (1 + math.erf(x / math.sqrt(2)))


def gaussian_kmae_reference(k: float) -> float:
    """P(|e| <= k*MAE) for a zero-mean normal error, where MAE = sigma*sqrt(2/pi)."""
    return 2 * norm_cdf(k * math.sqrt(2 / math.pi)) - 1


def prepare(tidy: pd.DataFrame) -> pd.DataFrame:
    df = tidy.copy()
    df["error_mm"] = df["precip_obs_mm"] - df["precip_fcst_mm"]
    df["fcst_bin"] = pd.cut(df["precip_fcst_mm"], FCST_BINS,
                            labels=FCST_BIN_LABELS)
    df["city_group"] = df["city"].map(CITY_GROUPS)
    return df


# ------------------------------------------------------- B1 coverage of k*MAE

def coverage_kmae(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (model, lead), g in df.groupby(["model", "lead_days"]):
        mae = g["error_mm"].abs().mean()
        for k in (1, 2):
            inside = int((g["error_mm"].abs() <= k * mae).sum())
            lo, hi = wilson(inside, len(g))
            rows.append({
                "model": model, "lead_days": lead, "n": len(g),
                "mae_mm": mae, "k": k, "band_half_width_mm": k * mae,
                "coverage": inside / len(g), "wilson_lo": lo, "wilson_hi": hi,
                "gaussian_reference": gaussian_kmae_reference(k),
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------- B2 quantile bands

QUANTILES = {"q05": 0.05, "q16": 0.16, "q50": 0.50, "q84": 0.84, "q95": 0.95}


def quantile_bands(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (model, lead), g in df.groupby(["model", "lead_days"]):
        e = g["error_mm"]
        row = {"model": model, "lead_days": lead, "n": len(g),
               "mean_error_mm": e.mean()}
        row.update({name: e.quantile(q) for name, q in QUANTILES.items()})
        row["width68_mm"] = row["q84"] - row["q16"]
        row["width90_mm"] = row["q95"] - row["q05"]
        rows.append(row)
    return pd.DataFrame(rows)


def _attach_band(df: pd.DataFrame, bands: pd.DataFrame) -> pd.DataFrame:
    """Tag every case with the band of its own (model, lead) and hit/miss flags."""
    cols = ["model", "lead_days", "q05", "q16", "q84", "q95"]
    out = df.merge(bands[cols], on=["model", "lead_days"], how="left")
    out["in68"] = (out["error_mm"] >= out["q16"]) & (out["error_mm"] <= out["q84"])
    out["in90"] = (out["error_mm"] >= out["q05"]) & (out["error_mm"] <= out["q95"])
    return out


def _coverage_rows(g: pd.DataFrame, keys: dict) -> dict:
    n = len(g)
    k68, k90 = int(g["in68"].sum()), int(g["in90"].sum())
    lo68, hi68 = wilson(k68, n)
    lo90, hi90 = wilson(k90, n)
    return {**keys, "n": n,
            "coverage68": k68 / n, "wilson68_lo": lo68, "wilson68_hi": hi68,
            "coverage90": k90 / n, "wilson90_lo": lo90, "wilson90_hi": hi90}


def coverage_stratified(tagged: pd.DataFrame) -> pd.DataFrame:
    """Coverage of the pooled per-lead bands within strata that did not define them.

    Every case is scored against the band of its own lead, so leads can then be
    pooled freely. The 'all cases' rows are nominal by construction (the bands
    were fitted on exactly that sample) and are printed only as a control.
    """
    rows = []
    for model, gm in tagged.groupby("model"):
        for lead_group, leads in {"all": list(range(1, 8)), **LEAD_GROUPS}.items():
            gl = gm[gm["lead_days"].isin(leads)]
            if gl.empty:
                continue
            rows.append(_coverage_rows(gl, {
                "model": model, "lead_group": lead_group,
                "stratum_type": "all", "stratum": "all cases"}))
            for stratum_type, col in [("forecast_bin", "fcst_bin"),
                                      ("city_group", "city_group")]:
                for stratum, gs in gl.groupby(col, observed=True):
                    rows.append(_coverage_rows(gs, {
                        "model": model, "lead_group": lead_group,
                        "stratum_type": stratum_type, "stratum": str(stratum)}))
    return pd.DataFrame(rows)


# ------------------------------------------- B3 bands by forecast amount

def conditional_bands(df: pd.DataFrame) -> pd.DataFrame:
    """Error quantiles by forecast-amount bin, per lead and per lead group.

    Per-lead cells fall to n ~ 66 in the heaviest bin, which is thin for a q05 /
    q95 (the 4th and 63rd order statistics), so the lead groups 1-3 and 4-7 are
    the rows meant for reporting; the per-lead rows are kept with their n so the
    aggregation can be checked rather than trusted.
    """
    rows = []
    groups = {str(lead): [lead] for lead in range(1, 8)} | LEAD_GROUPS
    for model, gm in df.groupby("model"):
        for lead_group, leads in groups.items():
            gl = gm[gm["lead_days"].isin(leads)]
            for fbin, g in gl.groupby("fcst_bin", observed=True):
                if g.empty:
                    continue
                e = g["error_mm"]
                row = {"model": model, "lead_group": lead_group,
                       "forecast_bin": str(fbin), "n": len(g),
                       "mean_forecast_mm": g["precip_fcst_mm"].mean(),
                       "mean_observed_mm": g["precip_obs_mm"].mean(),
                       "mean_error_mm": e.mean(), "mae_mm": e.abs().mean()}
                row.update({name: e.quantile(q) for name, q in QUANTILES.items()})
                row["width68_mm"] = row["q84"] - row["q16"]
                row["width90_mm"] = row["q95"] - row["q05"]
                rows.append(row)
    return pd.DataFrame(rows)


# ------------------------------------- B4 rank histogram and spread-error

def _wide_members(df: pd.DataFrame) -> pd.DataFrame:
    members = df[df["model"].isin(ENSEMBLE_MEMBERS)]
    return members.pivot_table(
        index=["city", "valid_date", "lead_days", "precip_obs_mm"],
        columns="model", values="precip_fcst_mm").dropna().reset_index()


def rank_histogram(wide: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rank of the observation among the 3 members, ties randomised.

    Daily precipitation has a large atom at zero: whenever the observation and
    one or more members are all 0.0 the rank is not defined, and the standard
    fix is to draw it uniformly from the tied block (Hamill 2001). The fraction
    of cases affected is reported alongside, and a 'wet' subset (any of the four
    values above 0.1 mm) is tabulated separately because the zero atom otherwise
    dominates the histogram.
    """
    rng = np.random.default_rng(RANK_SEED)
    x = wide[ENSEMBLE_MEMBERS].to_numpy()
    y = wide["precip_obs_mm"].to_numpy()[:, None]
    n_below = (x < y).sum(axis=1)
    n_tied = (x == y).sum(axis=1)
    rank = n_below + 1 + rng.integers(0, n_tied + 1)  # 1..m+1

    w = wide.assign(rank=rank, tied=n_tied > 0,
                    ens_min=x.min(axis=1), ens_max=x.max(axis=1),
                    ens_mean=x.mean(axis=1))
    w["wet"] = np.maximum(w["ens_max"], w["precip_obs_mm"]) > 0.1
    w["below"] = w["precip_obs_mm"] < w["ens_min"]
    w["above"] = w["precip_obs_mm"] > w["ens_max"]
    w["inside"] = ~w["below"] & ~w["above"]

    rows = []
    for subset, gs in [("all", w), ("wet", w[w["wet"]])]:
        for lead, g in gs.groupby("lead_days"):
            n = len(g)
            counts = g["rank"].value_counts().reindex(range(1, 5), fill_value=0)
            row = {"subset": subset, "lead_days": lead, "n": n,
                   "tied_frac": g["tied"].mean()}
            row.update({f"rank{r}_frac": counts[r] / n for r in range(1, 5)})
            row["below_min_frac"] = g["below"].mean()
            row["inside_frac"] = g["inside"].mean()
            row["above_max_frac"] = g["above"].mean()
            lo, hi = wilson(int(g["inside"].sum()), n)
            row["inside_wilson_lo"], row["inside_wilson_hi"] = lo, hi
            rows.append(row)
    return pd.DataFrame(rows), w


def spearman(a: pd.Series, b: pd.Series) -> float:
    """Spearman rho as Pearson on mid-ranks (pandas defers this one to scipy)."""
    return a.rank().corr(b.rank())


def spread_error(w: pd.DataFrame) -> pd.DataFrame:
    """Spearman correlation between ensemble range and |error of the ensemble mean|.

    Over all cases this correlation is inflated to ~0.9 by the dry/wet contrast
    alone: on a dry day both the spread and the error are 0, so the statistic
    mostly measures whether it rained. The 'wet' subset (any member or the
    observation above 0.1 mm) is the one that says whether spread carries
    information about the error *given* that something was going to happen.
    """
    rows = []
    for subset, gs in [("all", w), ("wet", w[w["wet"]])]:
        for lead, g in gs.groupby("lead_days"):
            spread = g["ens_max"] - g["ens_min"]
            abs_err = (g["ens_mean"] - g["precip_obs_mm"]).abs()
            rows.append({
                "subset": subset, "lead_days": lead, "n": len(g),
                "mean_spread_mm": spread.mean(),
                "mean_abs_error_mm": abs_err.mean(),
                "spearman": spearman(spread, abs_err),
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------- B5 the cone

def fig_uncertainty_cone(df: pd.DataFrame, bands: pd.DataFrame,
                         cond: pd.DataFrame, outdir: Path,
                         cities: list[str], targets: dict, model: str) -> None:
    """Forecast for one day as it converged from lead 7 to lead 1, with its band.

    Two bands are drawn: the marginal one of B2 (all cases of that lead) and the
    conditional one of B3 for the forecast-amount bin each forecast actually
    falls in. The gap between them is the point of the figure -- a band fitted
    on a sample that is 40% dry days is far too narrow for a forecast of 80 mm.

    Each panel uses that city's own wettest day rather than one date shared by
    all of them. The transect's peak day is not every city's peak day, and
    pinning the figure to it once put La Serena on a date when the gauges in
    that cell were still two days short of their maximum.
    """
    fig, axes = plt.subplots(1, len(cities), figsize=(10, 4.4), sharex=True)
    colour = MODEL_COLORS[model]
    for ax, city in zip(axes, cities):
        target_date = targets[city]
        g = (df[(df["city"] == city) & (df["valid_date"] == target_date)
                & (df["model"] == model)]
             .sort_values("lead_days"))
        b = bands[bands["model"] == model].set_index("lead_days")
        leads = g["lead_days"].to_numpy()
        f = g["precip_fcst_mm"].to_numpy()
        obs = g["precip_obs_mm"].iloc[0]

        lo90 = np.clip(f + b.loc[leads, "q05"].to_numpy(), 0, None)
        hi90 = f + b.loc[leads, "q95"].to_numpy()
        lo68 = np.clip(f + b.loc[leads, "q16"].to_numpy(), 0, None)
        hi68 = f + b.loc[leads, "q84"].to_numpy()

        # Conditional band: pick, per lead, the row for the bin this forecast is in.
        cm = cond[(cond["model"] == model)].set_index(
            ["lead_group", "forecast_bin"])
        clo, chi = [], []
        for lead, fv in zip(leads, f):
            group = next(k for k, v in LEAD_GROUPS.items() if lead in v)
            fbin = str(pd.cut([fv], FCST_BINS, labels=FCST_BIN_LABELS)[0])
            r = cm.loc[(group, fbin)]
            clo.append(max(0.0, fv + r["q05"]))
            chi.append(fv + r["q95"])

        ax.fill_between(leads, lo90, hi90, color=colour, alpha=0.13, linewidth=0,
                        label="marginal band, 90%")
        ax.fill_between(leads, lo68, hi68, color=colour, alpha=0.25, linewidth=0,
                        label="marginal band, 68%")
        ax.plot(leads, clo, color=colour, linewidth=1.2, linestyle=":",
                label="band conditioned on\nforecast amount, 90%")
        ax.plot(leads, chi, color=colour, linewidth=1.2, linestyle=":")
        ax.plot(leads, f, marker="o", markersize=6, color=colour,
                label=f"{MODEL_LABELS[model]} forecast")
        ax.axhline(obs, color=INK, linewidth=1.5, linestyle="--")
        ax.text(0.98, obs, f"observed {obs:.0f} mm", fontsize=8, color=INK,
                va="bottom", ha="right")
        ax.set_xlim(7.4, 0.6)  # lead 7 on the left: time runs towards the event
        ax.set_ylim(bottom=0)
        ax.set_xticks(list(range(1, 8)))
        ax.set_xlabel("Forecast lead (days before the event)")
        ax.set_title(f"{CITY_LABELS[city]} — {target_date:%d %b}")
    axes[0].set_ylabel("Daily precipitation (mm)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=8,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(
        "The cone of uncertainty: each city's wettest day as it was forecast, "
        "with its empirical error bands (in-sample)", fontweight="bold")
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(outdir / "fig8_uncertainty_cone.png", bbox_inches="tight")
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

    # The error bands are the project's headline answer to "80 mm plus or minus
    # what", so they are computed against every arbiter: a band calibrated on a
    # gridded product that runs wet is a band around the wrong number.
    parts: dict[str, list] = {}
    # The cone figure is drawn against the gauges: it answers "did what fell
    # land inside the band?", and what fell is what the instruments caught.
    # The bands themselves are published for every arbiter, as everything else
    # in this project is.
    FIG_TRUTH = "gauge"
    fig_df = fig_bands = fig_cond = None
    for truth, frame in by_truth(tidy):
        df = prepare(frame)
        bands = quantile_bands(df)
        cond = conditional_bands(df)
        ranks, w = rank_histogram(_wide_members(df))
        tables = {
            "coverage_kmae": coverage_kmae(df),
            "quantile_bands": bands,
            "coverage_stratified": coverage_stratified(_attach_band(df, bands)),
            "conditional_bands": cond,
            "rank_histogram": ranks,
            "spread_error": spread_error(w),
        }
        for name, table in tables.items():
            parts.setdefault(name, []).append(table.assign(truth=truth))
        if truth == FIG_TRUTH:
            fig_df, fig_bands, fig_cond = df, bands, cond

    for name, tabs in parts.items():
        pd.concat(tabs, ignore_index=True).to_csv(
            metrics_dir / f"{name}.csv", index=False)
    print(f"metrics -> {metrics_dir} ({len(parts)} CSV files)")

    # Each city's own wettest day, chosen from the data rather than assumed.
    obs = fig_df[["city", "valid_date", "precip_obs_mm"]].drop_duplicates()
    cone_cities = ["la_serena", "combarbala"]
    targets = {c: obs[obs["city"] == c].set_index("valid_date")
               ["precip_obs_mm"].idxmax() for c in cone_cities}
    print("wettest day at the gauges: "
          + ", ".join(f"{c} {d}" for c, d in targets.items()))

    set_style()
    fig_uncertainty_cone(fig_df, fig_bands, fig_cond, fig_dir,
                         cone_cities, targets, "ecmwf_ifs025")
    print(f"figures -> {fig_dir / 'fig8_uncertainty_cone.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
