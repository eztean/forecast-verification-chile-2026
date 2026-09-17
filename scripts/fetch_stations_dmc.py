"""Harvest rain-gauge observations from the DMC national network.

CEAZA-Met covers only the Coquimbo region, so the eleven cities from Valparaíso
to Puerto Montt had no independent arbiter at all -- they were verified against
a gridded product alone. The Dirección Meteorológica de Chile is the only source
that closes that gap: every other candidate was checked and ruled out
(GHCN-Daily carries a single Chilean station for 2026 and no precipitation at
all; Chilean METARs report hourly precipitation as 0.00 straight through the
event, i.e. they do not measure it).

The service needs a personal token (free self-service registration at
climatologia.meteochile.gob.cl); `DMC_API_USER` and `DMC_API_TOKEN` are read
from a git-ignored .env at the project root.

Three findings from probing the API shape everything below:

  1. Only `Automática` and `Mixta` stations answer `getDatosRecientesEma`;
     `Convencional` ones return an HTML error page. That still leaves at least
     one serving station within 30 km of all twenty cities except Los Vilos,
     which CEAZA already covers.

  2. `aguaCaidaDelMinuto` is one minute of rain sampled every fifteen, exactly
     as documented. Summing it over July at Quinta Normal gives 10.6 mm for a
     month that actually recorded 155. It is unusable as a daily total and is
     not used here.

  3. `aguaCaida6Horas` is a running accumulation that resets on fixed six-hour
     UTC blocks, and the sample at HH:00 closes the previous block rather than
     opening the next. Differencing it within each block recovers a true
     15-minute incremental series, which is what this script does. The result
     was checked three ways at Quinta Normal for July 2026: 155.4 mm by
     increments, 155.1 mm by summing the per-block maxima, and an event total
     of 132.6 mm against 139.1 mm from the gridded product.

`aguaCaida24Horas` is deliberately unused: it accumulates 12:01-12:00 UTC, a
synoptic window that cannot be reconciled with the local calendar day that
build_dataset.py and CEAZA both use.

Timestamps arrive in UTC (the payload says so) and are shifted by a fixed -4 h.
Chile observes no DST between April and September, so that is local time for
the whole window; the script refuses windows that leave that period.

Raw responses are written verbatim (gzipped -- one station-month is 2.6 MB of
JSON that compresses 23x) alongside a fetch_manifest.json recording every URL
and the UTC timestamp of the request, matching the CEAZA harvest.

Usage: uv run scripts/fetch_stations_dmc.py [--year 2026] [--month 7]
"""

import argparse
import datetime as dt
import gzip
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
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

BASE = "https://climatologia.meteochile.gob.cl/application"
NETWORK = "dmc"
USER_AGENT = ("forecast-verification-chile-2026/0.1 "
              "(research script; contact: eztean.morales@gmail.com)")

# Only these answer getDatosRecientesEma; Convencional stations return HTML.
SERVING_TYPES = {"Automática", "Mixta"}
UTC_OFFSET_H = -4          # America/Santiago, April-September (no DST)
SAMPLES_PER_DAY = 96       # 15-minute records
RETRIES = 3

OUTDIR = ROOT / "data" / "raw" / "stations" / NETWORK


def load_credentials() -> tuple[str, str]:
    """Read DMC_API_USER / DMC_API_TOKEN from the environment or .env."""
    env = dict(os.environ)
    dotenv = ROOT / ".env"
    if dotenv.exists():
        for line in dotenv.read_text().splitlines():
            line = line.strip().removeprefix("export ").strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, val = line.partition("=")
            env.setdefault(key.strip(), val.strip().strip("'\""))
    try:
        return env["DMC_API_USER"], env["DMC_API_TOKEN"]
    except KeyError:
        raise SystemExit(
            "missing DMC_API_USER / DMC_API_TOKEN.\n"
            "Register at https://climatologia.meteochile.gob.cl/application/"
            "usuario/registroUsuario and put them in a git-ignored .env "
            "at the project root.")


def get(path: str, user: str, token: str, manifest: list,
        label: str) -> bytes | None:
    """Call the service, record the call, return the raw body (None on error)."""
    url = f"{BASE}/{path}?" + urllib.parse.urlencode(
        {"usuario": user, "token": token})
    requested = dt.datetime.now(dt.timezone.utc).isoformat()
    # The URL carries the token, so the manifest records the path only -- it is
    # committed to a public repository.
    entry = {"label": label, "path": path, "requested_utc": requested}
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    for attempt in range(RETRIES):
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                raw = resp.read()
            manifest.append({**entry, "attempts": attempt + 1,
                             "bytes": len(raw)})
            time.sleep(1)  # be polite to a public service
            return raw
        except (urllib.error.HTTPError, urllib.error.URLError,
                TimeoutError) as exc:
            if attempt == RETRIES - 1:
                manifest.append({**entry, "attempts": RETRIES,
                                 "error": str(exc)})
                print(f"    {label}: {exc}", file=sys.stderr)
                return None
            time.sleep(4 * (attempt + 1))
    return None


def parse_mm(value) -> float | None:
    if value is None:
        return None
    m = re.match(r"\s*(-?[\d.]+)", str(value))
    return float(m.group(1)) if m else None


def catalogue(user: str, token: str, manifest: list,
              reprocess: bool = False) -> pd.DataFrame:
    cached = OUTDIR / "catastro.json"
    if reprocess and cached.exists():
        raw = cached.read_bytes()
    else:
        raw = get("geoservicios/getCatastroEstacionesGeo", user, token,
                  manifest, "catalogue")
        if raw is None:
            raise SystemExit("could not fetch the DMC station catalogue")
        cached.write_bytes(raw)
    feats = json.loads(raw.decode("utf-8"))["features"]
    rows = []
    for f in feats:
        p = f["features"]["properties"]
        if p.get("tipoEstacion") not in SERVING_TYPES:
            continue
        try:
            rows.append({
                "station_id": str(p["CodigoNacional"]),
                "station_name": str(p["nombreEstacion"]).strip(),
                "lat": float(p["latitud"]), "lon": float(p["longitud"]),
                "elev_m": float(p["altitud"]),
                "tipo": p.get("tipoEstacion"),
                "clasificacion": p.get("tipoClasificacion"),
            })
        except (KeyError, TypeError, ValueError):
            continue
    return pd.DataFrame(rows)


def closing_block(t: dt.datetime) -> tuple[dt.date, int]:
    """The six-hour UTC block a sample reports on, under the closing convention.

    Stepping back one minute puts the HH:00 reading in the block it closes
    rather than the one it would otherwise open.
    """
    tt = t - dt.timedelta(minutes=1)
    return tt.date(), tt.hour // 6


def is_block_boundary(t: dt.datetime) -> bool:
    return t.minute == 0 and t.hour % 6 == 0


def daily_from_month(payload: dict) -> tuple[pd.DataFrame, dict]:
    """Local-day rainfall totals rebuilt from the six-hour running counter.

    Stations do not agree on what the reading at a block boundary means. Quinta
    Normal reports the closing total of the block that just ended (06:00 reads
    14.4, then 06:15 reads 0.3); Carriel Sur has already reset (11:45 reads 0.2,
    12:00 reads 0.0). Assuming either convention for all of them silently
    corrupts the other half of the network, so the convention is detected from
    the data one boundary at a time: a boundary reading that has dropped is a
    reset and opens the new block, and one that has not is a close.

    Under the reset convention the rain that fell in the last quarter hour of a
    block is not recoverable from this field -- the counter was already back to
    zero when it was sampled. That is a bounded undercount of up to four
    fifteen-minute slivers a day, and it is counted and reported rather than
    hidden.

    Returns the daily frame and a diagnostics dict.
    """
    datos = payload.get("datosEstaciones", {}).get("datos") or []
    samples = []
    for r in datos:
        try:
            t = dt.datetime.strptime(r["momento"], "%Y-%m-%d %H:%M:%S")
        except (KeyError, ValueError):
            continue
        samples.append((t, parse_mm(r.get("aguaCaida6Horas"))))
    samples.sort()

    totals: dict[dt.date, float] = defaultdict(float)
    counted: dict[dt.date, int] = defaultdict(int)
    prev_block = prev_value = None
    diag = {"boundary_resets": 0, "boundary_closes": 0, "glitches": 0}
    for t, value in samples:
        local_day = (t + dt.timedelta(hours=UTC_OFFSET_H)).date()
        if value is None:
            continue
        dropped = prev_value is not None and value < prev_value - 0.05
        if is_block_boundary(t) and dropped:
            # Reset convention: this sample opens the next block.
            diag["boundary_resets"] += 1
            block = (t.date(), t.hour // 6)
            step = value
        else:
            if is_block_boundary(t):
                diag["boundary_closes"] += 1
            block = closing_block(t)
            step = value if block != prev_block else value - prev_value
            if step < -0.05:
                # A drop inside a block is a genuine glitch, not a convention.
                diag["glitches"] += 1
                step = 0.0
        prev_block, prev_value = block, value
        totals[local_day] += max(step, 0.0)
        counted[local_day] += 1

    daily = pd.DataFrame({
        "valid_date": sorted(totals),
        "precip_mm": [round(totals[d], 2) for d in sorted(totals)],
        "completeness_pc": [100.0 * counted[d] / SAMPLES_PER_DAY
                            for d in sorted(totals)],
    })
    return daily, diag


def fetch_station_month(station_id: str, year: int, month: int, user: str,
                        token: str, manifest: list,
                        reprocess: bool) -> pd.DataFrame | None:
    path = OUTDIR / f"ema_{station_id}_{year}{month:02d}.json.gz"
    if reprocess:
        # Rebuild the daily series from the archived response, so a change to
        # the aggregation can be re-run without touching the service again.
        if not path.exists():
            return None
        raw = gzip.open(path, "rb").read()
    else:
        label = f"series/{station_id}/{year}-{month:02d}"
        raw = get(
            f"servicios/getDatosRecientesEma/{station_id}/{year}/{month:02d}",
            user, token, manifest, label)
        if raw is None:
            return None
    try:
        payload = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError):
        # Convencional stations answer with an HTML error page.
        print(f"    {station_id}: not an EMA station (non-JSON response)")
        return None
    if "datosEstaciones" not in payload:
        print(f"    {station_id}: {payload.get('mensaje') or 'no data'}")
        return None
    if not reprocess:
        with gzip.open(path, "wb") as fh:
            fh.write(raw)
    daily, diag = daily_from_month(payload)
    resets, closes = diag["boundary_resets"], diag["boundary_closes"]
    convention = ("resets at HH:00" if resets > closes else "closes at HH:00")
    note = f"    {station_id}: {convention} ({resets} resets / {closes} closes)"
    if diag["glitches"]:
        note += f", {diag['glitches']} in-block drops clamped"
    print(note)
    return daily


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--year", type=int, default=2026)
    parser.add_argument("--month", type=int, default=7)
    parser.add_argument("--capture", default="2026-07-26",
                        help="capture whose grid-cell elevations are compared against")
    parser.add_argument("--reprocess", action="store_true",
                        help="rebuild the daily series from the archived raw "
                             "responses without calling the service again")
    args = parser.parse_args()

    if not 4 <= args.month <= 8:
        raise SystemExit(
            f"month {args.month} may fall in Chilean DST; this script assumes a "
            f"fixed UTC{UTC_OFFSET_H:+d} offset and only covers April-August.")

    user, token = ("", "") if args.reprocess else load_credentials()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    manifest: list = []

    stations = catalogue(user, token, manifest, args.reprocess)
    print(f"catalogue: {len(stations)} serving stations "
          f"({'/'.join(sorted(SERVING_TYPES))})")

    matches = match_stations(stations,
                             model_cell_elevations(args.capture, CITIES),
                             CITIES, "DMC")
    print(f"{matches['city'].nunique()} cities have a station within "
          f"{MAX_DIST_KM:.0f} km ({len(matches)} station-city pairs, "
          f"{matches['station_id'].nunique()} distinct stations)")

    # A station can serve two cities, so the month is downloaded once per code.
    series_cache: dict = {}
    frames, kept = [], []
    for row in matches.itertuples():
        if row.station_id not in series_cache:
            print(f"  {row.station_id} {row.station_name[:40]}")
            series_cache[row.station_id] = fetch_station_month(
                row.station_id, args.year, args.month, user, token, manifest,
                args.reprocess)
        series = series_cache[row.station_id]
        if series is None or series.empty:
            continue
        if series["precip_mm"].isna().all():
            print(f"    {row.station_id}: no rain gauge on this station")
            continue
        series = qc_daily(series, "precip_mm", "completeness_pc")
        frames.append(pd.DataFrame({
            "city": row.city, "role": row.role, "station_id": row.station_id,
            "network": NETWORK, "valid_date": series["valid_date"],
            "precip_station_mm": series["precip_station_mm"],
            "completeness_pc": series["completeness_pc"].round(1),
            "qc_flags": series["qc_flags"],
        }))
        rec = row._asdict()
        rec.pop("Index", None)
        kept.append({**rec, "sensor_id": row.station_id, "unit": "mm",
                     "n_days": len(series),
                     "n_usable": int(series["_usable"].sum())})

    if not frames:
        raise SystemExit("no DMC station returned usable data")
    daily = pd.concat(frames, ignore_index=True)
    meta = pd.DataFrame(kept)

    ext = ROOT / "data" / "external"
    ext.mkdir(parents=True, exist_ok=True)
    meta.to_csv(ext / "stations_metadata_dmc.csv", index=False)

    procdir = ROOT / "data" / "processed" / args.capture
    procdir.mkdir(parents=True, exist_ok=True)
    daily.to_csv(procdir / "stations_daily_dmc.csv", index=False)

    (OUTDIR / "fetch_manifest.json").write_text(json.dumps({
        "network": NETWORK, "service": BASE,
        "fetched_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "window": {"year": args.year, "month": args.month},
        "timezone": f"payload is UTC; shifted {UTC_OFFSET_H:+d} h to "
                    f"America/Santiago (no DST April-September)",
        "daily_convention": "sum of 15-min increments differenced from "
                            "aguaCaida6Horas, over the local calendar day",
        "unused_fields": {
            "aguaCaidaDelMinuto": "1 minute sampled every 15; not an accumulation",
            "aguaCaida24Horas": "accumulates 12:01-12:00 UTC, not a local day",
        },
        "calls": manifest,
    }, indent=2))

    usable = meta[meta["n_usable"] > 0]
    print(f"\n{usable['city'].nunique()} cities with a usable gauge, "
          f"{usable['station_id'].nunique()} distinct gauges, "
          f"{len(daily)} station-days")
    print(f"raw + manifest -> {OUTDIR}")
    print(f"metadata -> {ext / 'stations_metadata_dmc.csv'}")
    print(f"daily -> {procdir / 'stations_daily_dmc.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
