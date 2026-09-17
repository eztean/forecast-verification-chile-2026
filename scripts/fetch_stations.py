"""Harvest rain-gauge observations from the CEAZA-Met network (Coquimbo region).

The gridded series used everywhere else in this project does not assimilate
Chilean rain gauges, and it comes from the same ECMWF system as one of the
models being verified. Gauges are the only arbiter here that a model did not
make, and this script fetches the northern half of them; fetch_stations_dmc.py
covers the rest of the transect and stations_common.py holds the matching and
QC rules both of them share.

CEAZA-Met serves a public CSV web service (http://www.ceazamet.cl/ws/) with no
registration; the caller identifies itself by e-mail in the `user` parameter.
Three calls are used:

    fn=GetListaEstaciones   station catalogue (code, name, lat, lon, elevation)
    fn=GetListaSensores     sensors of one station (the gauge is "Precipitación")
    fn=GetSerieSensor       the series itself, with interv=dia for daily totals

Daily aggregation: the service returns timestamps in GMT-4, which is Chilean
standard time and therefore local time in July (no DST between April and
September), and `interv=dia` sums a midnight-to-midnight local calendar day.
That is exactly the convention build_dataset.py uses for the forecasts and for
the gridded observations, so the series are directly comparable with no
realignment -- unlike the DMC feed, which arrives in UTC. Verified against the
hourly series: La Serena [CEAZA] on 16 July sums to 6.4 mm by hour and reports
6.4 mm as its daily total.

Raw responses are written verbatim to data/raw/stations/ceazamet/ alongside a
fetch_manifest.json recording every URL, its parameters and the UTC timestamp of
the request, so the harvest is reproducible and auditable like the API captures.

Usage: uv run scripts/fetch_stations.py [--start 2026-07-01] [--end 2026-07-27]
"""

import argparse
import datetime as dt
import io
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from capture_forecasts import CITIES  # noqa: E402  (path set above)
from stations_common import (  # noqa: E402
    MAX_DIST_KM,
    ROOT,
    match_stations,
    model_cell_elevations,
    qc_daily,
)

WS = "http://www.ceazamet.cl/ws/pop_ws.php"
NETWORK = "ceazamet"
# The service asks callers to identify themselves; no registration or key.
USER = "eztean.morales@gmail.com"
USER_AGENT = ("forecast-verification-chile-2026/0.1 "
              "(research script; contact: eztean.morales@gmail.com)")

CHUNK_DAYS = 20             # GetSerieSensor truncates any response to 24 records
RETRIES = 4                 # the service returns sporadic 502s

OUTDIR = ROOT / "data" / "raw" / "stations" / NETWORK


def ws_get(fn: str, params: dict, manifest: list, label: str,
           retries: int = RETRIES) -> str:
    """Call the web service, record the call in the manifest, return the body."""
    full = {"fn": fn, "user": USER, "p_cod": NETWORK, **params}
    url = WS + "?" + urllib.parse.urlencode(full)
    requested = dt.datetime.now(dt.timezone.utc).isoformat()
    # The service rejects urllib's default user-agent with a 403, and returns a
    # sporadic 502 under load, so every call is retried with a growing backoff.
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=90) as resp:
                raw = resp.read()
            break
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
            if attempt == retries - 1:
                manifest.append({"label": label, "fn": fn, "url": url,
                                 "params": full, "requested_utc": requested,
                                 "attempts": retries, "error": str(exc)})
                raise
            wait = 4 * (attempt + 1)
            print(f"    {label}: {exc} -- retrying in {wait}s", file=sys.stderr)
            time.sleep(wait)
    try:
        body = raw.decode("utf-8")
    except UnicodeDecodeError:  # older records in the catalogue are latin-1
        body = raw.decode("latin-1")
    manifest.append({"label": label, "fn": fn, "url": url, "params": full,
                     "requested_utc": requested, "attempts": attempt + 1,
                     "bytes": len(body)})
    time.sleep(1)  # be polite to a public service
    return body


def parse_csv(body: str, expect: list[str]) -> pd.DataFrame:
    """Parse a CEAZA CSV body.

    Every line is prefixed with '#' except the data: the last comment line is
    the header. Columns are read from that header rather than assumed, because
    the service appends fields of its own beyond the ones requested via c0..cN
    and every row carries a trailing comma.
    """
    lines = body.splitlines()
    header = next((ln for ln in lines
                   if ln.startswith("#") and ln.lstrip("#").startswith(expect[0])),
                  None)
    rows = [ln for ln in lines if ln and not ln.startswith("#") and "," in ln]
    if header is None or not rows:
        return pd.DataFrame(columns=expect)
    # The catalogue header repeats e_nombre, so duplicates get a suffix and the
    # first occurrence -- the one the c0..cN request asked for -- keeps the name.
    names, seen = [], {}
    for col in header.lstrip("#").rstrip(",").split(","):
        col = col.strip()
        seen[col] = seen.get(col, 0) + 1
        names.append(col if seen[col] == 1 else f"{col}.{seen[col] - 1}")
    df = pd.read_csv(io.StringIO("\n".join(rows)), header=None,
                     names=names, usecols=range(len(names)), index_col=False)
    missing = [c for c in expect if c not in df.columns]
    if missing:
        raise RuntimeError(f"CEAZA response is missing columns {missing}; "
                           f"got {list(df.columns)}")
    return df


def find_gauge(station_id: str, manifest: list) -> tuple[str, str] | None:
    body = ws_get("GetListaSensores", {"e_cod": station_id}, manifest,
                  f"sensors/{station_id}")
    (OUTDIR / f"sensors_{station_id}.csv").write_text(body)
    sensors = parse_csv(body, ["e_cod", "s_cod", "tf_nombre", "um_notacion",
                               "s_altura", "s_ultima_lectura"])
    gauge = sensors[sensors["tf_nombre"].astype(str)
                    .str.strip().str.lower().str.startswith("precipitaci")]
    if gauge.empty:
        return None
    row = gauge.iloc[0]
    return str(row["s_cod"]), str(row["um_notacion"]).strip()


def _fetch_range(s_cod: str, a: dt.date, b: dt.date,
                 manifest: list, out: list) -> None:
    """Fetch one date range, bisecting around ranges the service cannot serve.

    Besides the 24-row cap, GetSerieSensor answers 502 to some particular
    (sensor, range) combinations however often they are retried -- Vicuña's
    gauge over 21-27 July is one -- while both halves of the same range succeed.
    So a failing range is split rather than abandoned, and only a single day
    that still fails is recorded as missing, leaving the rest of the series
    intact.
    """
    try:
        body = ws_get("GetSerieSensor",
                      {"s_cod": s_cod, "fecha_inicio": str(a),
                       "fecha_fin": str(b), "interv": "dia"},
                      manifest, f"series/{s_cod}/{a}_{b}",
                      retries=2 if a < b else RETRIES)
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
        if a == b:
            print(f"    series/{s_cod}: {a} unavailable ({exc}) -- left as a gap",
                  file=sys.stderr)
            return
        mid = a + (b - a) // 2
        _fetch_range(s_cod, a, mid, manifest, out)
        _fetch_range(s_cod, mid + dt.timedelta(days=1), b, manifest, out)
        return
    (OUTDIR / f"series_{s_cod}_daily_{a}_{b}.csv").write_text(body)
    chunk = parse_csv(body, ["s_cod", "ultima_lectura", "prom", "data_pc"])
    if not chunk.empty:
        out.append(chunk)


def fetch_series(s_cod: str, start: str, end: str,
                 manifest: list) -> pd.DataFrame:
    """Daily totals over [start, end], paged around the service's 24-row cap.

    GetSerieSensor silently truncates every response to 24 records whatever the
    requested range, so a longer window has to be walked in chunks and the
    result de-duplicated on the date.
    """
    first, last = dt.date.fromisoformat(start), dt.date.fromisoformat(end)
    chunks: list = []
    cursor = first
    while cursor <= last:
        stop = min(cursor + dt.timedelta(days=CHUNK_DAYS - 1), last)
        _fetch_range(s_cod, cursor, stop, manifest, chunks)
        cursor = stop + dt.timedelta(days=1)
    if not chunks:
        return pd.DataFrame()
    df = pd.concat(chunks, ignore_index=True)
    df["valid_date"] = pd.to_datetime(df["ultima_lectura"]).dt.date
    return (df.drop_duplicates(subset="valid_date")
            .sort_values("valid_date").reset_index(drop=True))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", default="2026-07-01")
    parser.add_argument("--end", default="2026-07-27")
    parser.add_argument("--capture", default="2026-07-26",
                        help="capture whose grid-cell elevations are compared against")
    args = parser.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    manifest: list = []

    body = ws_get("GetListaEstaciones",
                  {"c0": "e_cod", "c1": "e_nombre", "c2": "e_lat", "c3": "e_lon",
                   "c4": "e_altitud"}, manifest, "catalogue")
    (OUTDIR / "estaciones.csv").write_text(body)
    stations = parse_csv(body, ["e_cod", "e_nombre", "e_lat", "e_lon",
                                "e_altitud"]).dropna(subset=["e_lat"])
    # The matcher is shared with the DMC network, so CEAZA's column names are
    # translated to the network-neutral schema it expects.
    stations = stations.rename(columns={
        "e_cod": "station_id", "e_nombre": "station_name", "e_lat": "lat",
        "e_lon": "lon", "e_altitud": "elev_m"})
    print(f"catalogue: {len(stations)} stations")

    matches = match_stations(stations,
                             model_cell_elevations(args.capture, CITIES),
                             CITIES, "CEAZA-Met")
    print(f"{matches['city'].nunique()} cities have a station within "
          f"{MAX_DIST_KM:.0f} km ({len(matches)} station-city pairs)")

    # A station can be the neighbour of two cities (Huintil serves Illapel and
    # Salamanca), so both the gauge lookup and the series are cached by code.
    gauge_cache: dict = {}
    series_cache: dict = {}
    frames, kept = [], []
    for row in matches.itertuples():
        if row.station_id not in gauge_cache:
            gauge_cache[row.station_id] = find_gauge(row.station_id, manifest)
        gauge = gauge_cache[row.station_id]
        if gauge is None:
            print(f"  {row.city}/{row.role}: {row.station_name} has no rain gauge")
            continue
        s_cod, unit = gauge
        if s_cod not in series_cache:
            series_cache[s_cod] = fetch_series(s_cod, args.start, args.end,
                                               manifest)
        series = series_cache[s_cod].copy()
        if series.empty:
            print(f"  {row.city}/{row.role}: {s_cod} returned no data")
            continue
        series = qc_daily(series, "prom", "data_pc")
        frames.append(pd.DataFrame({
            "city": row.city, "role": row.role, "station_id": row.station_id,
            "network": NETWORK, "valid_date": series["valid_date"],
            "precip_station_mm": series["precip_station_mm"],
            "completeness_pc": series["data_pc"], "qc_flags": series["qc_flags"],
        }))
        rec = row._asdict()
        rec.pop("Index", None)
        kept.append({**rec, "sensor_id": s_cod, "unit": unit,
                     "n_days": len(series),
                     "n_usable": int(series["_usable"].sum())})
        print(f"  {row.city}/{row.role}: {row.station_name} ({s_cod}, {unit}) "
              f"-> {len(series)} days, "
              f"{int((~series['_usable']).sum())} unusable")

    daily = pd.concat(frames, ignore_index=True)
    meta = pd.DataFrame(kept)

    ext = ROOT / "data" / "external"
    ext.mkdir(parents=True, exist_ok=True)
    meta.to_csv(ext / "stations_metadata.csv", index=False)

    procdir = ROOT / "data" / "processed" / args.capture
    procdir.mkdir(parents=True, exist_ok=True)
    daily.to_csv(procdir / "stations_daily.csv", index=False)

    (OUTDIR / "fetch_manifest.json").write_text(json.dumps({
        "network": NETWORK, "service": WS,
        "fetched_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "window": {"start": args.start, "end": args.end},
        "timezone": "GMT-4 (= America/Santiago in July, no DST)",
        "daily_convention": "interv=dia, midnight-to-midnight local calendar day",
        "calls": manifest,
    }, indent=2))

    usable = meta[meta["n_usable"] > 0]
    print(f"\n{usable[usable['role'] == 'primary']['city'].nunique()} cities with a "
          f"usable primary gauge, {usable['station_id'].nunique()} distinct gauges, "
          f"{len(daily)} station-days")
    print(f"raw + manifest -> {OUTDIR}")
    print(f"metadata -> {ext / 'stations_metadata.csv'}")
    print(f"daily -> {procdir / 'stations_daily.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
