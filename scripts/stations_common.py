"""Shared station-to-city matching and QC, used by every gauge network.

Two networks feed this project's independent arbiter -- CEAZA-Met in the
Coquimbo region and the DMC national network for the rest of the transect --
and they must be matched to cities and quality-controlled by exactly the same
rules, or the "gauge" truth would mean something different north and south of
Los Vilos. This module holds those rules; the network modules hold only what is
specific to their own web service.

The station table each network hands in is network-neutral: one row per station
with `station_id`, `station_name`, `lat`, `lon`, `elev_m`.
"""

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent

# A station is accepted as a proxy for a city below MAX_DIST_KM, and flagged as
# a loose match (reported, used, caveated) between PREFERRED and MAX.
PREFERRED_DIST_KM = 15.0
MAX_DIST_KM = 30.0
MAX_DELTA_ELEV_M = 300.0
MIN_COMPLETENESS_PC = 90.0  # sub-daily completeness below this makes a day NaN

# The frontal passage. Earlier versions of this pipeline carried two hand-set
# windows -- a 3-day "core" (16-18 July) and a wider 8-day episode (15-22) --
# and used them interchangeably in prose, so the same city had two different
# "event totals". There is now one window, and it is derived from the data by
# `event_window()` rather than typed in. These constants are the value that
# derivation returns for the final capture; `event_window()` asserts it.
EVENT_WINDOW = ("2026-07-16", "2026-07-21")

# Rule for `event_window()`: a day belongs to the event if a majority of the
# transect's cities recorded measurable rain at their gauges, and the window is
# the contiguous run of such days containing the wettest day. The rule is stated
# on the gauges rather than on the gridded product so that the window is not set
# by the arbiter under test, and the majority requirement is what separates a
# transect-wide frontal passage from local rain at one end of it.
EVENT_WET_DAY_MM = 1.0
EVENT_MAJORITY_FRAC = 0.5

# A gauge differing from its peers by more than this factor is treated as an
# instrument problem rather than as weather, and left out of the city estimate.
# The ratio is reported per gauge so the choice can be second-guessed.
OUTLIER_FACTOR = 2.0

# Convicting a gauge takes at least this many peers to convict it with. With a
# single peer, two gauges that disagree accuse each other symmetrically and
# nothing in the data says which one is wrong -- excluding both would throw away
# the city's only measurement to avoid choosing. Those cities keep their
# estimate and are flagged `gauges_disagree` instead, with the spread between
# the gauges left visible as the honest uncertainty. Adjudicating with the
# gridded product is not an option: it is the thing under test.
MIN_PEERS_TO_EXCLUDE = 2

# Below this many gauges reporting, a city-day has no gauge estimate at all.
MIN_GAUGES_PER_DAY = 1

# A gauge reporting fewer than this share of the event's days is too thin to
# judge or to average. It is deliberately looser than MIN_COMPLETENESS_PC: that
# one asks whether a single day is trustworthy, this one whether a series is,
# and over an eight-day window a 90% bar would discard a gauge for missing one
# day. The comparison itself uses each gauge's mean over the days it did
# report, so moderate gaps do not bias the verdict.
MIN_GAUGE_COVERAGE_PC = 60.0


def haversine_km(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    r = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlam = math.radians(lon2 - lon1)
    a = (math.sin(dphi / 2) ** 2
         + math.cos(p1) * math.cos(p2) * math.sin(dlam / 2) ** 2)
    return 2 * r * math.asin(math.sqrt(a))


def model_cell_elevations(capture: str, cities) -> dict[str, float]:
    """Elevation of the Open-Meteo grid cell each city was sampled from.

    The station-to-model elevation gap is the physically meaningful one here:
    it is a direct measure of how much sub-grid orography the ~9-31 km cell has
    smoothed away, which is the representativeness question this task is about.
    """
    out = {}
    for city in cities:
        path = ROOT / "data" / "raw" / capture / f"{city}__observed__grid.json"
        out[city] = json.loads(path.read_text()).get("elevation", float("nan"))
    return out


def match_stations(stations: pd.DataFrame, cell_elev: dict, cities: dict,
                   network: str) -> pd.DataFrame:
    """Every station within MAX_DIST_KM of each city, nearest first.

    The nearest one is the city's `primary` gauge and the rest are `neighbour`s.
    Both roles are kept and both are used: the city's rainfall is estimated from
    all of its gauges (see station_validation.py), and the spread among them is
    what tells a single under-catching instrument apart from a whole
    neighbourhood reading low. Only the second is evidence about the gridded
    product.
    """
    rows = []
    for city, (lat, lon) in cities.items():
        cand = stations.assign(
            dist_km=[haversine_km(lat, lon, r.lat, r.lon)
                     for r in stations.itertuples()])
        cand = cand[cand["dist_km"] <= MAX_DIST_KM].sort_values("dist_km")
        for rank, st in enumerate(cand.itertuples()):
            delta = st.elev_m - cell_elev[city]
            notes = []
            if st.dist_km > PREFERRED_DIST_KM:
                notes.append(f"loose match: {st.dist_km:.0f} km from the city")
            if abs(delta) > MAX_DELTA_ELEV_M:
                notes.append(f"elevation gap {delta:+.0f} m vs the model cell")
            rows.append({
                "city": city, "role": "primary" if rank == 0 else "neighbour",
                "station_id": st.station_id, "station_name": st.station_name,
                "network": network, "lat": st.lat, "lon": st.lon,
                "elev_m": st.elev_m, "cell_elev_m": cell_elev[city],
                "dist_km": round(st.dist_km, 1),
                "delta_elev_m": round(delta, 0),
                "notes": "; ".join(notes),
            })
    return pd.DataFrame(rows)


def load_gauges(capture: str) -> pd.DataFrame:
    """Every gauge-day from both networks, in one frame."""
    proc = ROOT / "data" / "processed" / capture
    frames = []
    for name in ("stations_daily.csv", "stations_daily_dmc.csv"):
        path = proc / name
        if path.exists():
            frames.append(pd.read_csv(path))
    if not frames:
        raise FileNotFoundError(f"no station data in {proc}")
    daily = pd.concat(frames, ignore_index=True)
    daily["valid_date"] = pd.to_datetime(daily["valid_date"]).dt.date
    # A station can serve two cities; the pair is the unit, not the station.
    return daily.drop_duplicates(subset=["city", "station_id", "valid_date"])


def gauge_quality(daily: pd.DataFrame,
                  window: tuple[str, str] = EVENT_WINDOW) -> pd.DataFrame:
    """Judge every gauge against its peers in the same city.

    The earlier version of this check compared the single nearest gauge with the
    median of its neighbours, which asked whether a city's chosen instrument was
    odd. Now that a city's rainfall is estimated from all of its gauges, the
    question changes: every gauge has to earn its place in that average, so each
    is compared with the median of the *other* gauges in its city and the ones
    that disagree by more than OUTLIER_FACTOR are left out -- provided there are
    enough peers to convict (see MIN_PEERS_TO_EXCLUDE). Nothing is deleted: the
    verdict is a column, and the ratio behind it is reported so the threshold
    can be argued with.

    A gauge that is merely missing days would look dry for a trivial reason, so
    coverage is checked first and thin series are excluded before the ratios are
    computed rather than being allowed to drag the peer median around.
    """
    w = daily[daily["valid_date"].astype(str).between(*window)]
    n_days = w["valid_date"].nunique()
    per_gauge = (w.groupby(["city", "station_id", "network"], dropna=False)
                 ["precip_station_mm"]
                 .agg(total_mm="sum", n_reported="count").reset_index())
    per_gauge["coverage_pc"] = 100.0 * per_gauge["n_reported"] / n_days
    # The ratio is built on the daily mean, not the total, so a gauge that
    # missed two days is not mistaken for a gauge that stayed dry.
    per_gauge["mean_mm"] = (per_gauge["total_mm"]
                            / per_gauge["n_reported"].replace(0, pd.NA))
    thin = per_gauge["coverage_pc"] < MIN_GAUGE_COVERAGE_PC

    # Each gauge is compared with the median of the *others*, not of the whole
    # city including itself. With only two gauges an all-inclusive median sits
    # exactly between them, so their ratio to it can never exceed 2 however far
    # apart they are -- which would let Los Vilos keep a gauge reading 85.7 mm
    # next to one reading 238.4. Leaving the gauge out fixes that and changes
    # almost nothing where there are several.
    good = per_gauge[~thin]
    loo, n_peers = [], []
    for row in per_gauge.itertuples():
        others = good[(good["city"] == row.city)
                      & (good["station_id"] != row.station_id)]["mean_mm"]
        loo.append(others.median() if len(others) else np.nan)
        n_peers.append(len(others))
    per_gauge["peer_median_mean_mm"] = loo
    per_gauge["n_peers"] = n_peers
    ratio = per_gauge["mean_mm"] / per_gauge["peer_median_mean_mm"]
    per_gauge["ratio_to_peers"] = ratio
    outlier = (ratio > OUTLIER_FACTOR) | (ratio < 1 / OUTLIER_FACTOR)
    convictable = per_gauge["n_peers"] >= MIN_PEERS_TO_EXCLUDE
    per_gauge["exclusion"] = ""
    per_gauge.loc[thin, "exclusion"] = "insufficient_coverage"
    per_gauge.loc[~thin & convictable & outlier.fillna(False),
                  "exclusion"] = "disagrees_with_peers"
    per_gauge["used"] = per_gauge["exclusion"] == ""
    # A two-gauge city whose gauges disagree keeps both and says so.
    per_gauge["gauges_disagree"] = (~thin & ~convictable
                                    & outlier.fillna(False))
    return per_gauge


def city_gauge_daily(daily: pd.DataFrame, quality: pd.DataFrame) -> pd.DataFrame:
    """A city's daily rainfall as the mean of the gauges that passed QC.

    A single nearest gauge is a point measurement standing in for a city, and
    the choice of point moves the answer: at La Serena the nearest gauge is also
    the lowest of the seven inside that grid cell. Averaging the gauges is both
    steadier and closer in kind to the cell average it is being compared with,
    and the spread among them is carried alongside as the honest uncertainty of
    the estimate rather than being thrown away.
    """
    keep = quality.loc[quality["used"], ["city", "station_id"]]
    used = daily.merge(keep, on=["city", "station_id"], how="inner")
    used = used.dropna(subset=["precip_station_mm"])
    est = (used.groupby(["city", "valid_date"])["precip_station_mm"]
           .agg(precip_gauge_mm="mean", n_gauges="size",
                gauge_sd_mm=lambda s: s.std(ddof=0))
           .reset_index())
    return est[est["n_gauges"] >= MIN_GAUGES_PER_DAY]


def event_window(city_est: pd.DataFrame,
                 expect: tuple[str, str] | None = EVENT_WINDOW
                 ) -> tuple[str, str]:
    """The event's first and last day, derived from the gauges (see above).

    Returns the contiguous run of majority-wet days containing the wettest day
    of the transect. Passing `expect` asserts the result, so a future capture
    that moved the window could not silently rescale every event total in the
    report; pass `expect=None` to derive without asserting.
    """
    day = (city_est.assign(wet=city_est["precip_gauge_mm"] >= EVENT_WET_DAY_MM)
           .groupby("valid_date")
           .agg(frac_wet=("wet", "mean"), total=("precip_gauge_mm", "sum"))
           .sort_index())
    majority = day["frac_wet"] >= EVENT_MAJORITY_FRAC
    days = list(day.index)
    peak = day["total"].idxmax()
    if not majority[peak]:
        raise ValueError(f"wettest day {peak} is not a majority-wet day")
    lo = hi = days.index(peak)
    while lo > 0 and majority.iloc[lo - 1]:
        lo -= 1
    while hi < len(days) - 1 and majority.iloc[hi + 1]:
        hi += 1
    got = (str(days[lo]), str(days[hi]))
    if expect is not None and got != tuple(expect):
        raise AssertionError(
            f"event window derived from the gauges is {got}, but the pipeline "
            f"is pinned to {tuple(expect)}. If the capture changed, update "
            f"EVENT_WINDOW and re-run every downstream script -- every event "
            f"total in the report is summed over this window.")
    return got


def qc_daily(series: pd.DataFrame, precip_col: str,
             completeness_col: str) -> pd.DataFrame:
    """Flag unusable daily totals without ever dropping or editing a value.

    A day whose sub-daily completeness is too low, or whose total is negative,
    keeps its row and its flag; only the precipitation becomes NaN. Deleting the
    row instead would silently shorten the series and let a gauge with a bad
    week look like a gauge that agrees with its neighbours.
    """
    out = series.copy()
    out["qc_flags"] = ""
    low = out[completeness_col] < MIN_COMPLETENESS_PC
    out.loc[low, "qc_flags"] = "incomplete_day"
    negative = out[precip_col] < 0
    out.loc[negative, "qc_flags"] += "negative_value;"
    out["precip_station_mm"] = out[precip_col].where(~(low | negative))
    out["_usable"] = ~(low | negative)
    return out
