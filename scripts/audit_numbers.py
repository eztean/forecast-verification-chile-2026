"""Check the numbers quoted in the deliverables against the metrics CSVs.

Every number this project states in prose is supposed to be traceable to a CSV
produced by this pipeline. That rule is easy to state and easy to break: during
one rewrite, a table of twenty-one cities was partly reconstructed from an
earlier intermediate output and eight of its rows were wrong. Nothing in the
pipeline noticed.

So the check is a script. Each entry below names a claim, the number the prose
states, and the value recomputed from the CSVs -- optionally followed by the
number of decimals the prose rounded to, since "+30%" is an honest way to write
30.4 and the check should know that. The script prints a line per claim and
exits non-zero if any of them disagree.

The section tags on the claim names (§2.1, §5.1, ...) group the checks by the
part of the analysis they belong to, in the order the results are presented.
This covers the load-bearing numbers, the ones a reader would quote, not every
digit written down.

Usage: uv run scripts/audit_numbers.py [--capture YYYY-MM-DD]
"""

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
TOL = 0.006  # relative; prose rounds to 2-3 significant figures


def build_checks(m: Path, numbers: dict) -> list[tuple]:
    cont = pd.read_csv(m / "continuous_by_lead.csv")
    cat = pd.read_csv(m / "categorical_by_lead.csv")
    rel = pd.read_csv(m / "ensemble_reliability.csv")
    brier = pd.read_csv(m / "ensemble_brier.csv")
    crps = pd.read_csv(m / "ensemble_crps.csv")
    daily = pd.read_csv(m / "station_vs_grid_daily.csv")
    totals = pd.read_csv(m / "event_totals_threeway.csv").set_index("city")
    rob = pd.read_csv(m / "robustness_by_metric.csv")
    sub = pd.read_csv(m / "subgrid_case.csv").iloc[0]
    inter = pd.read_csv(m / "subgrid_intercity.csv")
    qual = pd.read_csv(m / "gauge_quality_check.csv")
    inten = pd.read_csv(m / "bias_structure_intensity.csv")
    split = pd.read_csv(m / "bias_structure_freq_intensity.csv")
    geo = pd.read_csv(m / "bias_structure_geography.csv")
    floor = pd.read_csv(m / "bias_structure_noise_floor.csv").iloc[0]
    lags = pd.read_csv(m / "bias_structure_timing.csv")
    strat = pd.read_csv(m / "coverage_stratified.csv")
    condb = pd.read_csv(m / "conditional_bands.csv")
    press = pd.read_csv(m / "press_cross_check.csv")

    def mae(model, lead, truth):
        r = cont[(cont.model == model) & (cont.lead_days == lead)
                 & (cont.truth == truth)]
        return float(r["precip_mae_mm"].iloc[0])

    def bias(model, lead, truth):
        r = cont[(cont.model == model) & (cont.lead_days == lead)
                 & (cont.truth == truth)]
        return float(r["precip_bias_mm"].iloc[0])

    def reliab(prob, truth):
        r = rel[(rel.threshold_mm == 10.0) & (rel.truth == truth)
                & (abs(rel.prob - prob) < 0.01)]
        return 100 * float(r["obs_freq"].iloc[0])

    def pooled(col, gridded):
        r = daily[(daily.city == "(pooled)") & (daily.gridded == gridded)]
        return float(r[col].iloc[0])

    def verdicts(metric, verdict):
        r = rob[(rob.truth_a == "grid") & (rob.truth_b == "gauge")
                & (rob.metric == metric)]
        return int((r.verdict == verdict).sum())

    def bss(thr, lead, truth):
        return float(brier[(brier.threshold_mm == thr) & (brier.truth == truth)
                           & (brier.lead_days == lead)].bss.iloc[0])

    def cond(group, lo_or_hi):
        r = condb[(condb.truth == "gauge") & (condb.model == "ecmwf_ifs025")
                  & (condb.lead_group == group)
                  & (condb.forecast_bin == "f >= 10")]
        return float(r[lo_or_hi].iloc[0])

    def cover(fbin, level):
        r = strat[(strat.truth == "gauge") & (strat.model == "ecmwf_ifs025")
                  & (strat.lead_group == "all") & (strat.stratum == fbin)]
        return 100 * float(r[f"coverage{level}"].iloc[0])

    COQ = ["la_serena", "coquimbo", "vicuna", "andacollo", "ovalle",
           "monte_patria", "combarbala", "illapel", "salamanca", "los_vilos"]
    gauge_tot = totals["station_mm"]

    def transect(col):
        return float(totals[col].sum())

    def signed(col, subset=None):
        d = totals if subset is None else totals.loc[subset]
        return float((d[col] - d["station_mm"]).mean())

    ga, gr = numbers["truths"]["gauge"], numbers["truths"]["grid"]
    return [
        # 2.1 / 2.4 -- dataset and gauge network
        ("§2.1 filas del dataset", 11907, numbers["dataset"]["rows"]),
        ("§2.1 n por modelo-plazo", 441, numbers["dataset"]["n_per_model_lead"]),
        ("§2.1 ciudades", 21, numbers["dataset"]["cities"]),
        ("§2.3 ventana inicio", "2026-07-16", numbers["event_window"][0]),
        ("§2.3 ventana fin", "2026-07-21", numbers["event_window"][1]),
        ("§2.4 pares ciudad-estacion", 80, len(qual)),
        ("§2.4 pares que pasan QC", 69, int(qual["used"].sum())),
        ("§2.4 pares ciudad-dia", 438, int(
            daily[(daily.city == "(pooled)")
                  & (daily.gridded == "grid")].n.iloc[0])),
        ("§2.4 ciudades con un pluviometro", 5, int(
            (qual[qual["used"]].groupby("city").size() == 1).sum())),
        ("§2.4 excluidos por desacuerdo", 2, int(
            (qual["exclusion"] == "disagrees_with_peers").sum())),
        ("§2.4 excluidos por cobertura", 9, int(
            (qual["exclusion"] == "insufficient_coverage").sum())),
        # 4 -- event totals at the gauges
        ("§4 Combarbala pluviometros", 304.0, float(gauge_tot["combarbala"])),
        ("§4 La Serena pluviometros", 182.2, float(gauge_tot["la_serena"])),
        ("§4 La Serena grillado", 310.2, float(totals.loc["la_serena", "grid_mm"])),
        ("§4 Rancagua pluviometros", 93.0, float(gauge_tot["rancagua"])),
        ("§4 Temuco pluviometros", 27.5, float(gauge_tot["temuco"])),
        # 5.1 -- error growth, both arbiters
        ("§5.1 MAE best_match p1 gauge", 4.70, mae("best_match", 1, "gauge")),
        ("§5.1 MAE ECMWF p1 gauge", 5.07, mae("ecmwf_ifs025", 1, "gauge")),
        ("§5.1 MAE GFS p1 gauge", 5.10, mae("gfs_seamless", 1, "gauge")),
        ("§5.1 MAE ICON p1 gauge", 6.02, mae("icon_seamless", 1, "gauge")),
        ("§5.1 MAE GFS p3 gauge", 5.97, mae("gfs_seamless", 3, "gauge")),
        ("§5.1 MAE GFS p7 gauge", 9.10, mae("gfs_seamless", 7, "gauge")),
        ("§5.1 MAE ECMWF p1 grid", 3.63, mae("ecmwf_ifs025", 1, "grid")),
        ("§5.1 MAE GFS p1 grid", 6.08, mae("gfs_seamless", 1, "grid")),
        ("§5.1 sesgo ECMWF p1 gauge", 3.12, bias("ecmwf_ifs025", 1, "gauge")),
        ("§5.1 sesgo best_match p1 gauge", 2.08, bias("best_match", 1, "gauge")),
        ("§5.1 sesgo GFS p1 gauge", 0.72, bias("gfs_seamless", 1, "gauge")),
        ("§5.1 sesgo ICON p1 gauge", 0.56, bias("icon_seamless", 1, "gauge")),
        ("§5.1 sesgo ECMWF p1 grid", -0.15, bias("ecmwf_ifs025", 1, "grid")),
        ("§5.1 sesgo GFS p1 grid", -2.55, bias("gfs_seamless", 1, "grid")),
        ("§5.1 sesgo GFS p5 gauge", -5.42, bias("gfs_seamless", 5, "gauge")),
        ("§5.1 sesgo GFS p7 gauge", -6.52, bias("gfs_seamless", 7, "gauge")),
        # 5.2 -- THE headline: transect event totals
        ("§5.2 total transecto pluv", 3392, transect("station_mm"), 0),
        ("§5.2 total transecto grid", 4445, transect("grid_mm"), 0),
        ("§5.2 total transecto ECMWF p1", 4457,
         transect("ecmwf_ifs025_lead1_mm"), 0),
        ("§5.2 total transecto GFS p1", 3638,
         transect("gfs_seamless_lead1_mm"), 0),
        ("§5.2 total transecto ICON p1", 3575,
         transect("icon_seamless_lead1_mm"), 0),
        ("§5.2 error medio ECMWF todas", 50.7, signed("ecmwf_ifs025_lead1_mm")),
        ("§5.2 error medio ECMWF Coquimbo", 74.9,
         signed("ecmwf_ifs025_lead1_mm", COQ)),
        ("§5.2 error medio ICON Coquimbo", -1.0,
         signed("icon_seamless_lead1_mm", COQ), 1),
        ("§5.2 error medio GFS Coquimbo", -15.1,
         signed("gfs_seamless_lead1_mm", COQ)),
        ("§5.2 error medio grid Coquimbo", 80.6, signed("grid_mm", COQ)),
        # 5.3 -- categorical
        ("§5.3 tasa base gauge", 0.292, float(
            cat[(cat.threshold_mm == 10.0) & (cat.truth == "gauge")]
            .base_rate.iloc[0])),
        ("§5.3 tasa base grid", 0.342, float(
            cat[(cat.threshold_mm == 10.0) & (cat.truth == "grid")]
            .base_rate.iloc[0])),
        ("§5.3 POD GFS p7 gauge", 0.33, float(
            cat[(cat.threshold_mm == 10.0) & (cat.truth == "gauge")
                & (cat.model == "gfs_seamless") & (cat.lead_days == 7)]
            .pod.iloc[0])),
        ("§5.3 FAR ECMWF p1 gauge", 0.192, float(
            cat[(cat.threshold_mm == 10.0) & (cat.truth == "gauge")
                & (cat.model == "ecmwf_ifs025") & (cat.lead_days == 1)]
            .far.iloc[0])),
        ("§5.3 FAR ECMWF p1 grid", 0.066, float(
            cat[(cat.threshold_mm == 10.0) & (cat.truth == "grid")
                & (cat.model == "ecmwf_ifs025") & (cat.lead_days == 1)]
            .far.iloc[0])),
        # 5.4 -- the ensemble
        ("§5.4 reliability 0/3 gauge", 2.4, reliab(0.0, "gauge"), 1),
        ("§5.4 reliability 1/3 gauge", 31.8, reliab(1 / 3, "gauge")),
        ("§5.4 reliability 2/3 gauge", 64.0, reliab(2 / 3, "gauge")),
        ("§5.4 reliability 3/3 gauge", 86.4, reliab(1.0, "gauge")),
        ("§5.4 reliability 2/3 grid", 78.1, reliab(2 / 3, "grid")),
        ("§5.4 reliability 3/3 grid", 94.6, reliab(1.0, "grid")),
        ("§5.4 n plazos 1-6 gauge", 2628, ga["n_pooled_leads_1_6"]),
        ("§5.4 Murphy REL gauge", 0.004, ga["murphy"]["10.0"]["reliability"]),
        ("§5.4 Murphy RES gauge", 0.123, ga["murphy"]["10.0"]["resolution"]),
        ("§5.4 Murphy UNC gauge", 0.207, ga["murphy"]["10.0"]["uncertainty"]),
        ("§5.4 Murphy BS gauge", 0.088, ga["murphy"]["10.0"]["brier"]),
        ("§5.4 Murphy BSS gauge", 0.574, ga["murphy"]["10.0"]["bss"]),
        ("§5.4 Murphy RES grid", 0.154, gr["murphy"]["10.0"]["resolution"]),
        ("§5.4 Murphy BSS grid", 0.664, gr["murphy"]["10.0"]["bss"]),
        ("§5.4 BSS >=1mm p1 gauge", 0.797, bss(1.0, 1, "gauge")),
        ("§5.4 BSS >=10mm p1 gauge", 0.735, bss(10.0, 1, "gauge")),
        ("§5.4 BSS >=10mm p6 gauge", 0.328, bss(10.0, 6, "gauge")),
        ("§5.4 BSS >=20mm p6 gauge", 0.206, bss(20.0, 6, "gauge")),
        ("§5.4 CRPS p1 gauge", 3.45, float(
            crps[(crps.truth == "gauge") & (crps.lead_days == 1)]
            .crps_mm.iloc[0])),
        ("§5.4 CRPS p6 gauge", 6.02, float(
            crps[(crps.truth == "gauge") & (crps.lead_days == 6)]
            .crps_mm.iloc[0])),
        ("§5.4 miembro individual p1", 5.40,
         ga["crps"]["1"]["mean_single_member_abs_error_mm"]),
        ("§5.4 miembro individual p6", 9.05,
         ga["crps"]["6"]["mean_single_member_abs_error_mm"]),
        ("§3.4 CRPS fair p1 gauge", 2.48, ga["crps"]["1"]["crps_fair_mm"]),
        # 5.5 -- error bands
        ("§5.5 cobertura 68% f=0", 96.7, cover("f = 0", 68)),
        ("§5.5 cobertura 68% 0<f<10", 77.0, cover("0 < f < 10", 68)),
        ("§5.5 cobertura 68% f>=10", 27.8, cover("f >= 10", 68)),
        ("§5.5 cobertura 90% f>=10", 75.2, cover("f >= 10", 90)),
        ("§5.5 banda cond. p1-3 q05", -35.9, cond("1-3", "q05")),
        ("§5.5 banda cond. p1-3 q95", 18.2, cond("1-3", "q95")),
        ("§5.5 banda cond. p4-7 q05", -44.7, cond("4-7", "q05")),
        ("§5.5 banda cond. p4-7 q95", 32.7, cond("4-7", "q95")),
        # 6.2 -- the wet bias and its structure
        ("§6.2 sesgo grid vs pluv (%)", 30, pooled("bias_pc", "grid"), 0),
        ("§6.2 sesgo grid vs pluv (mm)", 3.27, pooled("bias_mm", "grid")),
        ("§6.2 media diaria grid", 14.01, pooled("gridded_mean_mm", "grid")),
        ("§6.2 media diaria pluv", 10.75, pooled("station_mean_mm", "grid")),
        ("§6.2 r diario grid", 0.91, pooled("pearson_r", "grid")),
        ("§6.2 r sqrt grid", 0.95, pooled("pearson_r_sqrt", "grid")),
        ("§6.2 sesgo era5 vs pluv (%)", 33, pooled("bias_pc", "era5"), 0),
        ("§6.2 sesgo era5 vs pluv (mm)", 3.52, pooled("bias_mm", "era5")),
        ("§6.2 r diario era5", 0.90, pooled("pearson_r", "era5")),
        ("§6.2 ratio 0-1 mm", 2.13, float(inten.ratio.iloc[0])),
        ("§6.2 ratio 5-15 mm", 1.66, float(inten.ratio.iloc[2])),
        ("§6.2 ratio 15-40 mm", 1.29, float(inten.ratio.iloc[3])),
        ("§6.2 ratio >=40 mm", 1.25, float(inten.ratio.iloc[4])),
        ("§6.2 exceso de frecuencia (%)", 10, 100 * (float(
            split[split.wet_day_threshold_mm == 1.0].freq_ratio.iloc[0]) - 1), 0),
        ("§6.2 exceso de intensidad (%)", 19, 100 * (float(
            split[split.wet_day_threshold_mm == 1.0]
            .intensity_ratio.iloc[0]) - 1), 0),
        ("§6.2 mediana costera (%)", 8, float(geo.median_a_pc.iloc[0]), 0),
        ("§6.2 mediana interior (%)", 33, float(geo.median_b_pc.iloc[0]), 0),
        ("§6.2 CV dentro de celda", 0.19,
         float(floor["pooled_within_cell_gauge_cv"])),
        ("§6.2 SD observada log-ratio", 0.235,
         float(floor["observed_sd_log_ratio"])),
        ("§6.2 SD de ruido log-ratio", 0.130,
         float(floor["noise_sd_log_ratio"])),
        ("§6.2 ruido / varianza (%)", 30,
         100 * float(floor["noise_share_of_variance"]), 0),
        ("§6.2 dispersion real residual (%)", 22,
         float(floor["residual_real_sd_pc"]), 0),
        ("§6.2 Ovalle exceso (%)", 44, 100 * (float(
            totals.loc["ovalle", "grid_over_station"]) - 1), 0),
        ("§6.2 Illapel exceso (%)", 7, 100 * (float(
            totals.loc["illapel", "grid_over_station"]) - 1), 0),
        # 6.3 -- robustness verdicts
        ("§6.3 sesgo desplazado", 27, verdicts("bias_mm", "shifted")),
        ("§6.3 sesgo cambia de signo", 21, int(
            rob[(rob.truth_a == "grid") & (rob.truth_b == "gauge")
                & (rob.metric == "bias_mm")].sign_flip.sum())),
        ("§6.3 MAE robusto", 21, verdicts("mae_mm", "robust")),
        ("§6.3 POD robusto", 19, verdicts("pod", "robust")),
        ("§6.3 grid vs era5 robustos", 106, int(
            (rob[(rob.truth_a == "grid")
                 & (rob.truth_b == "era5")].verdict == "robust").sum())),
        # 6.4 -- the sub-grid case
        ("§6.4 pluviometros en la celda", 7, int(sub["n_gauges_passing_qc"])),
        ("§6.4 minimo de la celda", 128.4, float(sub["gauge_min_mm"])),
        ("§6.4 maximo de la celda", 214.4, float(sub["gauge_max_mm"])),
        ("§6.4 factor de dispersion", 1.67, float(sub["gauge_spread_factor"])),
        ("§6.4 media de los pluviometros", 179.5, float(sub["gauge_mean_mm"])),
        ("§6.4 grillado en la celda", 310.2, float(sub["grid_total_mm"])),
        ("§6.4 ECMWF en la celda", 322.5, float(sub["ecmwf_ifs025_lead1_mm"])),
        ("§6.4 GFS en la celda", 131.5, float(sub["gfs_seamless_lead1_mm"])),
        ("§6.4 ICON en la celda", 216.8, float(sub["icon_seamless_lead1_mm"])),
        ("§6.4 diferencia entre ciudades", 0.0,
         float(sub["model_intercity_difference_mm"])),
        ("§6.4 combinaciones fecha-plazo", 147, int(
            inter.loc[inter.model == "ecmwf_ifs025", "n_forecast_days"].iloc[0])),
        ("§6.4 best_match entre ciudades", 29.6, float(
            inter.loc[inter.model == "best_match",
                      "max_abs_difference_mm"].iloc[0])),
        # 6.5 -- the two-day offset
        ("§6.5 ciudades sin desfase", 19, int((lags.best_lag_days == 0).sum())),
        ("§6.5 ciudades con desfase", 2, int((lags.best_lag_days != 0).sum())),
        # 7 -- press cross-check
        ("§7 filas de referencia", 12, len(press)),
        ("§7 filas con ventana completa", 5, int(
            (press.days_dropped_by_qc == 0).sum())),
        ("§7 mediana ventana completa (%)", -2.1, float(
            press[press.days_dropped_by_qc == 0].diff_pc.median())),
        ("§7 peor caso ventana completa (%)", 5.5, float(
            press[press.days_dropped_by_qc == 0].diff_pc.abs().max())),
        ("§7 mediana mensual (%)", -9.7, float(
            press[press.days_dropped_by_qc > 0].diff_pc.median())),
    ]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--capture", default="2026-07-26")
    args = parser.parse_args()

    m = ROOT / "data" / "processed" / args.capture / "metrics"
    numbers = json.loads((m / "report_numbers.json").read_text())
    checks = build_checks(m, numbers)

    bad = []
    for check in checks:
        label, quoted, actual = check[:3]
        digits = check[3] if len(check) > 3 else None
        if isinstance(quoted, str):
            ok = quoted == actual
            shown = str(actual)
        elif digits is not None:
            ok = round(actual, digits) == quoted
            shown = f"{actual:.4g}"
        else:
            scale = max(abs(quoted), 1e-9)
            ok = (abs(quoted - actual) / scale <= TOL
                  or abs(quoted - actual) < 0.005)
            shown = f"{actual:.4g}"
        mark = "ok  " if ok else "FAIL"
        print(f"  {mark} {label:38s} prosa {quoted!s:>10}  datos {shown:>10}")
        if not ok:
            bad.append((label, quoted, actual))

    print(f"\n{len(checks) - len(bad)}/{len(checks)} cifras verificadas")
    if bad:
        print("\nDESAJUSTES:")
        for label, quoted, actual in bad:
            print(f"  {label}: la prosa dice {quoted}, los datos dicen {actual}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
