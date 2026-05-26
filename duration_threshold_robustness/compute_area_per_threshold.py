"""
Robustness check: green-area trend under different duration-day thresholds.

For each year (2002-2023), in two duration datasets:
  1. First_last/Land  -> first-to-last NDVI>0.2 day method
  2. all_whole_map/Land -> daily-count method (Fig. S19 alt metric)
count pixels whose duration exceeds each of {1, 2, 3} days, convert to km^2
(pixel = 500 m), and run Mann-Kendall + Sen's slope + linear regression.

Output:
  duration_threshold_trend_results.csv  -- one row per (dataset, threshold, year)
  duration_threshold_trend_summary.csv  -- one row per (dataset, threshold)
"""

from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
import pymannkendall as mk
from scipy import stats

ROOT = Path(r"G:/Hangkai/Anttarctic Vegetation Dynamic/Version_2_data/Antarctica_Vegetation_Duration")

DATASETS = {
    "first_last": {
        "dir":  ROOT / "First_last" / "Land",
        "pattern": "Duration_{year}.tif",
    },
    "daily_count": {
        "dir":  ROOT / "all_whole_map" / "Land",
        "pattern": "Land_{year}.tif",
    },
}

YEARS = list(range(2002, 2024))
THRESHOLDS = list(range(1, 15))          # pixel kept if duration > T days (1..14)
PIXEL_AREA_KM2 = (500.0 * 500.0) / 1e6   # 0.25 km^2 per 500 m pixel

def area_for_year(tif_path: Path, t_days: int) -> float:
    with rasterio.open(tif_path) as ds:
        arr = ds.read(1)
    return float(np.sum(arr > t_days)) * PIXEL_AREA_KM2

rows = []
for ds_name, cfg in DATASETS.items():
    for year in YEARS:
        tif = cfg["dir"] / cfg["pattern"].format(year=year)
        if not tif.exists():
            print(f"missing: {tif}")
            continue
        for t in THRESHOLDS:
            area = area_for_year(tif, t)
            rows.append({
                "dataset": ds_name,
                "threshold_days": t,
                "year": year,
                "area_km2": area,
            })
            print(f"{ds_name} t>{t}d {year}: {area:,.1f} km^2")

df = pd.DataFrame(rows)
df.to_csv(ROOT / "duration_threshold_trend_results.csv", index=False)

summary_rows = []
for (ds_name, t), grp in df.groupby(["dataset", "threshold_days"]):
    grp = grp.sort_values("year")
    yrs = grp["year"].to_numpy(dtype=float)
    vals = grp["area_km2"].to_numpy(dtype=float)

    mk_res = mk.original_test(vals)
    sen = mk.sens_slope(vals)
    lr  = stats.linregress(yrs, vals)

    summary_rows.append({
        "dataset": ds_name,
        "threshold_days": t,
        "n_years": len(vals),
        "first_year_area_km2": vals[0],
        "last_year_area_km2":  vals[-1],
        "delta_km2":           vals[-1] - vals[0],
        "pct_change":          (vals[-1] - vals[0]) / vals[0] * 100.0 if vals[0] else np.nan,
        "mean_km2":            vals.mean(),
        "MK_trend":            mk_res.trend,
        "MK_p":                mk_res.p,
        "MK_z":                mk_res.z,
        "Sens_slope_km2_per_yr": sen.slope,
        "OLS_slope_km2_per_yr":  lr.slope,
        "OLS_intercept":         lr.intercept,
        "OLS_r2":                lr.rvalue ** 2,
        "OLS_p":                 lr.pvalue,
    })

summary = pd.DataFrame(summary_rows)
summary.to_csv(ROOT / "duration_threshold_trend_summary.csv", index=False)

with pd.option_context("display.max_columns", None, "display.width", 200):
    print("\n=== SUMMARY ===")
    print(summary.to_string(index=False))
