# External reference data (manual)

This directory holds *manually curated* reference values used to cross-check the ERA5 observations — it is the only part of the data tree that is not produced by a script, and every row must carry its source.

- `media_reported_totals.csv` — precipitation totals as reported by Chilean media and official communications (DMC, ONEMI/SENAPRED, municipal statements) during and after the July 2026 event. Media numbers are heterogeneous (different accumulation windows, station vs. city-wide values, sometimes unsourced), so each row records the exact accumulation window and the source URL. These values are used for *plausibility checks only*, never as verification truth.

Columns of `media_reported_totals.csv`:

| column | meaning |
|---|---|
| `city` | snake_case city id matching `scripts/capture_forecasts.py` |
| `window_start`, `window_end` | accumulation window (local dates, inclusive) |
| `reported_mm` | precipitation total as reported |
| `station` | station name if the report cites one (e.g. DMC La Florida), else empty |
| `source` | outlet or institution |
| `url` | link to the report |
| `accessed` | date the value was transcribed |
| `notes` | anything odd (rounding, conflicting values, etc.) |
