"""Stage 1: per-pixel OLS-detrend each climate variable, run all 256 detrending
coalitions through the trained RF, and save per-fold yearly continental and
regional area time series.

Math (per variable v, per pixel p):
    full(p, y)        = spatial(p) + temporal(y) + residual(p, y)
    slope_p           = OLS slope of full(p, ·) vs y
    mean_slope        = mean over pixels of slope_p
    y_mean            = midpoint of the period

    new_temporal(y)   = temporal(y)   − mean_slope             × (y − y_mean)
    new_residual(p,y) = residual(p,y) − (slope_p − mean_slope) × (y − y_mean)
    new_spatial       = spatial   (unchanged; centered detrending preserves pixel mean)

For each subset S ⊆ {8 climate variables} we feed the RF the perturbed
features for variables in S (originals for the rest), aggregate point-wise
predictions to a yearly continental area Σ(pred × land_area), and save it.

Outputs (under OUTPUT_DIR):
    <Region>_per_fold.csv       columns: region, scenario, fold, year, area
    <Region>_detrend_curves.csv columns: region, scenario, year, mean_area, std_area
"""

import gc
import time
import warnings
from itertools import combinations  # not used directly, kept for parity
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.exceptions import InconsistentVersionWarning

from paths import DATA, MODELS_DIR, OUTPUT_DIR

warnings.filterwarnings("ignore", category=InconsistentVersionWarning)
warnings.filterwarnings("ignore", message="X does not have valid feature names")


VARS = [
    "temperature_2m",
    "uv_radiation",
    "solar_radiation",
    "volumetric_soil_water",
    "runoff",
    "X10m_wind_speed",
    "precipitation",
    "vapor_pressure_deficit",
]
N_VARS = len(VARS)
N_SCENARIOS = 1 << N_VARS  # 256

REGIONS = {
    "Continental": (MODELS_DIR / "Continental_results_permutation", None),
    "East":        (MODELS_DIR / "Regional_results_permutation" / "East",      "East"),
    "West":        (MODELS_DIR / "Regional_results_permutation" / "West",      "West"),
    "Peninsula":   (MODELS_DIR / "Regional_results_permutation" / "Peninsula", "Peninsula"),
}


def patch_for_modern_sklearn(model):
    """RFs saved with sklearn ≤1.3 lack monotonic_cst on each tree, which
    sklearn 1.7 expects during predict. Patch to None."""
    if hasattr(model, "estimators_"):
        for tree in model.estimators_:
            if not hasattr(tree, "monotonic_cst"):
                tree.monotonic_cst = None
    if not hasattr(model, "monotonic_cst"):
        model.monotonic_cst = None
    return model


def scenario_name(subset_vars):
    if not subset_vars:
        return "baseline"
    return ",".join(sorted(subset_vars))


def load_baseline(region_filter):
    """Reconstruct the un-detrended baseline DataFrame from the detrended_<v>.csv
    files. Each detrended_<v>.csv has all-baseline columns *except* its own
    <v>_temporal and <v>_residual; we recover full baseline by taking
    detrended_temperature_2m.csv as the starting point and overlaying the
    temperature_2m_temporal/_residual from detrended_uv_radiation.csv (where
    those columns are at their original baseline values)."""
    print("[load] Constructing baseline DataFrame...")
    df = pd.read_csv(DATA / "detrended_temperature_2m.csv")
    df_uv = pd.read_csv(
        DATA / "detrended_uv_radiation.csv",
        usecols=["temperature_2m_temporal", "temperature_2m_residual"],
    )
    df["temperature_2m_temporal"] = df_uv["temperature_2m_temporal"].values
    df["temperature_2m_residual"] = df_uv["temperature_2m_residual"].values
    del df_uv
    gc.collect()

    if region_filter is not None:
        mask = (df["Regions"].values == region_filter)
        df = df.loc[mask].reset_index(drop=True)

    print(f"[load] Baseline rows: {len(df):,}")
    return df


def compute_perpixel_detrend(baseline_df):
    """Compute per-variable per-pixel centered detrending caches.

    Returns:
        cache              -- dict v -> (new_temporal, new_residual) as float32 arrays
                              aligned to baseline_df rows.
        mean_slope_per_var -- dict v -> mean across-pixels slope used for
                              the temporal redistribution (for reporting).
    """
    pixel  = baseline_df["pixel"].to_numpy()
    year   = baseline_df["year"].to_numpy()
    y_mean = (year.min() + year.max()) / 2.0
    y_centered = (year - y_mean).astype(np.float64)

    unique_pixels, pixel_idx = np.unique(pixel, return_inverse=True)
    n_pixels = len(unique_pixels)

    sum_yc2_per_pixel  = np.zeros(n_pixels, dtype=np.float64)
    np.add.at(sum_yc2_per_pixel, pixel_idx, y_centered ** 2)

    cache = {}
    mean_slope_per_var = {}

    for v in VARS:
        spatial  = baseline_df[f"{v}_spatial"].to_numpy(dtype=np.float64)
        temporal = baseline_df[f"{v}_temporal"].to_numpy(dtype=np.float64)
        residual = baseline_df[f"{v}_residual"].to_numpy(dtype=np.float64)
        full = spatial + temporal + residual

        # Per-pixel mean of full
        mean_full_pixel = np.zeros(n_pixels, dtype=np.float64)
        count_pixel     = np.zeros(n_pixels, dtype=np.float64)
        np.add.at(mean_full_pixel, pixel_idx, full)
        np.add.at(count_pixel, pixel_idx, 1.0)
        mean_full_pixel = mean_full_pixel / count_pixel

        # Per-pixel OLS slope: numer/denom
        full_centered_per_pixel = full - mean_full_pixel[pixel_idx]
        numer_per_pixel = np.zeros(n_pixels, dtype=np.float64)
        np.add.at(numer_per_pixel, pixel_idx, y_centered * full_centered_per_pixel)
        slope_p_unique = numer_per_pixel / sum_yc2_per_pixel

        slope_p = slope_p_unique[pixel_idx]
        mean_slope = float(slope_p_unique.mean())
        mean_slope_per_var[v] = mean_slope

        new_temporal = temporal - mean_slope * y_centered
        new_residual = residual - (slope_p - mean_slope) * y_centered

        cache[v] = (new_temporal.astype(np.float32),
                    new_residual.astype(np.float32))

        # Sanity: per-pixel slope of new_full should be ~0
        new_full = spatial + new_temporal + new_residual
        nf_centered = new_full - mean_full_pixel[pixel_idx]
        new_numer = np.zeros(n_pixels, dtype=np.float64)
        np.add.at(new_numer, pixel_idx, y_centered * nf_centered)
        new_slope_p_max = float(np.abs(new_numer / sum_yc2_per_pixel).max())
        print(f"   {v:<26} mean_slope_p={mean_slope:>12.4g}  "
              f"max|new_slope_p|={new_slope_p_max:.3e}")

    return cache, mean_slope_per_var


def run_region(region_name, model_dir, region_filter):
    t0 = time.time()
    print(f"\n=== {region_name} ===")

    m1 = patch_for_modern_sklearn(joblib.load(model_dir / "RF_model_fold1.joblib"))
    feat_cols = list(m1.feature_names_in_)
    del m1
    gc.collect()

    baseline_df = load_baseline(region_filter)
    n_rows = len(baseline_df)

    print("[detrend] Computing per-pixel centered detrending per variable...")
    cache, mean_slope_per_var = compute_perpixel_detrend(baseline_df)

    X_baseline = baseline_df[feat_cols].to_numpy(dtype=np.float32)
    year_arr = baseline_df["year"].to_numpy()
    land_arr = baseline_df["land_area"].to_numpy(dtype=np.float64)
    del baseline_df
    gc.collect()

    idx_for_var = {v: (feat_cols.index(f"{v}_temporal"),
                       feat_cols.index(f"{v}_residual")) for v in VARS}

    rows_per_fold = []
    for fold in range(1, 6):
        tf = time.time()
        print(f"[{region_name}] Loading fold {fold}...")
        model = patch_for_modern_sklearn(joblib.load(model_dir / f"RF_model_fold{fold}.joblib"))
        try:
            model.set_params(n_jobs=-1)
        except Exception:
            pass

        # Working copy of X — modify in place per scenario, restore after each.
        X = X_baseline.copy()

        for mask_int in range(N_SCENARIOS):
            subset = [VARS[i] for i in range(N_VARS) if mask_int & (1 << i)]
            sname = scenario_name(subset)

            saved = []
            for v in subset:
                t_idx, r_idx = idx_for_var[v]
                t_new, r_new = cache[v]
                saved.append((t_idx, X[:, t_idx].copy()))
                saved.append((r_idx, X[:, r_idx].copy()))
                X[:, t_idx] = t_new
                X[:, r_idx] = r_new

            pred = model.predict(X)
            ya = pd.Series(pred * land_arr).groupby(year_arr).sum()
            for yr, a in ya.items():
                rows_per_fold.append({
                    "region": region_name,
                    "scenario": sname,
                    "fold": fold,
                    "year": int(yr),
                    "area": float(a),
                })

            for col_idx, orig in saved:
                X[:, col_idx] = orig

            if (mask_int + 1) % 32 == 0 or mask_int == 0:
                print(f"   fold={fold}  {mask_int + 1:3d}/{N_SCENARIOS}  "
                      f"elapsed {time.time() - tf:5.0f}s")

        del model, X
        gc.collect()
        print(f"[{region_name}] fold {fold} done in {time.time() - tf:.1f}s")

    df_pf = pd.DataFrame(rows_per_fold)
    df_pf.to_csv(OUTPUT_DIR / f"{region_name}_per_fold.csv", index=False)

    df_summ = (
        df_pf.groupby(["region", "scenario", "year"])["area"]
        .agg(["mean", "std"])
        .reset_index()
        .rename(columns={"mean": "mean_area", "std": "std_area"})
    )
    df_summ.to_csv(OUTPUT_DIR / f"{region_name}_detrend_curves.csv", index=False)
    print(f"[{region_name}] DONE in {time.time() - t0:.1f}s "
          f"({len(df_pf):,} per-fold rows, {len(df_summ):,} summary rows)")


def main():
    t0 = time.time()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for name, (md, rf) in REGIONS.items():
        run_region(name, md, rf)
    print(f"\nALL DONE in {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
