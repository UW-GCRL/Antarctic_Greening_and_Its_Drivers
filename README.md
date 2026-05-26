# Antarctic Greening and Its Drivers

This repository contains the code, notebooks, and figure-generation scripts used in the study of Antarctic greenness (MODIS NDVI) and its environmental drivers. The workflow covers (i) MODIS-based greenness mapping, (ii) ERA5 environmental-driver preprocessing, (iii) Random Forest modeling at continental and biogeographic-regional scales, (iv) Shapley-based driver attribution using per-pixel mean-centered detrending counterfactuals, and (v) every robustness, sensitivity, and validation analysis reported in the manuscript and supplementary materials.

> The workflow was originally developed and run with **Python 3.8**.

---

## Repository layout

### Core data preparation
- `Antarctic_NDVI_DIY_NBAR.ipynb` — NDVI from MCD43A4 NBAR reflectance.
- `Antarctic_greening_mapping.ipynb` — Continental greenness maps (Fig. 1A, Fig. 2).
- `First_Last_NDVI.js` — Google Earth Engine script generating First/Last DOY of NDVI > 0.2.
- `MOD10A1_Antarctica_MinSnowCover_clean.js` — GEE script for MODIS snow-cover masking.
- `MODIS_Area_Rock&DEM&Lat_Map.ipynb` — Auxiliary maps (rock / DEM / latitude masks).
- `Download ERA5.ipynb`, `Process ERA5 data to tif.ipynb` — ERA5 acquisition and conversion to GeoTIFFs.
- `1-preprocess_rf_data.R` — Environmental-driver preprocessing and spatial / temporal / residual decomposition (RF inputs).

### Random Forest modeling
- `RF_continent.py` — Random Forest at the continental scale (5-fold spatiotemporal CV).
- `RF_region.py` — Random Forest at the biogeographic-regional scale (East / West / Peninsula).
- `Visualization of AG RF results.ipynb` — RF performance hexbin plots (Fig. S6, Fig. S7) and permutation-feature importance (Data S1).

### Duration and greenness-area analyses
- `Duration using last and first.ipynb` — First-and-last-day greenness duration (Fig. 1C).
- `SM_Duration_Day_Count.ipynb` — Daily-count duration metric (Fig. S24).

### Shapley-based driver attribution (`shapley_attribution/`)
The data-driven attribution pipeline used in Fig. 1D, Fig. 3, Fig. S22, and Fig. S23. Per-pixel mean-centered detrending of the eight climate variables, all 256 detrending coalitions evaluated through the trained RF models, exact Shapley values, and Mann–Kendall significance per scenario. See `shapley_attribution/README.md` for usage.
- `shapley_attribution/build_detrend_curves_perpixel.py` — Build the per-fold, per-scenario predicted green-area curves.
- `shapley_attribution/shapley_continental.py` — Exact Shapley values, continental scale.
- `shapley_attribution/analyze_regional.py` — Per-region Shapley values + sanity checks.
- `shapley_attribution/compute_mk_continental.py` — Per-scenario MK *p*-values (continental).
- `shapley_attribution/compute_mk_all_regions.py` — Per-scenario MK *p*-values for **all 4 regions** (data behind Fig. S23).
- `shapley_attribution/build_var_trend_table.py` — Area-weighted long-term variable trend table.
- `shapley_attribution/build_shapley_tables_with_std.py` — Region × variable Shapley tables.
- `shapley_attribution/plot_shapley_tables.py`, `shapley_attribution/plot_combined_attribution_figure.py` — Publication figures.
- `shapley_attribution/paths.py`, `shapley_attribution/requirements.txt` — Env-var path config and Python deps.

### Manuscript figure-generation scripts (`figures/`)
One script per main-text or SI figure that depends on the new Shapley pipeline.
- `figures/plot_main_fig3_shapley.py` — **Fig. 3** (4-panel regional Shapley bars with Σφ).
- `figures/build_continental_climate_anomalies.py` — **Fig. 1** climate sparkline CSVs (annual continental area-weighted means, anomalies, and spatial spread for the eight climate variables).
- `figures/plot_continental_detrend_counterfactuals.py` — **Fig. S22** (modeled continental green-area trajectories under representative detrending counterfactuals).
- `figures/plot_uav_three_site_subpixel.py` — **Fig. S8** (UAV three-site sub-pixel moss extent in the Windmill Islands; expects the cached `*_uav_ndvi.npz` files from the AAS 4046 dataset).
- `figures/plot_duration_threshold_robustness.py` — **Fig. S11 and Fig. S12** (green-area and mean-duration trends across the ≥1..≥14 d duration-threshold sweep).
- `figures/plot_significance_matrix.py` — **Fig. S23** (per-region MK significance pattern across the 256 detrending counterfactuals).

### Duration-threshold robustness data pipeline (`duration_threshold_robustness/`)
- `duration_threshold_robustness/compute_area_per_threshold.py` — Per-year green-area under thresholds ≥1..≥14 d (data behind Fig. S11).
- `duration_threshold_robustness/compute_mean_duration_per_threshold.py` — Per-year mean greenness duration under thresholds ≥1..≥14 d (data behind Fig. S12).

### Supplementary analyses (one notebook per SI figure)
- `Supplementary _materials_mapping_area.ipynb` — **Fig. S1** (per-year hexagon area maps).
- `SM_Trend robustness without certain years.ipynb` — **Fig. S2 and Fig. S5** (±2σ / ±3σ outlier robustness for area and duration).
- `SM_NDVI_Sensitivity_New_Stable_Disappeared_Greenness_of_Hot_Spot_Region.ipynb` — **Fig. S3** (NDVI threshold sensitivity of stable / new / disappeared green area).
- `Supplementary _materials_mapping_duration.ipynb` — **Fig. S4** (per-year hexagon duration maps).
- `Supplementary _materials_Sensitivity_test.ipynb` — **Fig. S9 and Fig. S10** (NDVI threshold sensitivity for area and mean duration).
- `SM_Trend_of_Key_Environmental_Variables.ipynb` — **Fig. S13 and Fig. S14** (environmental-variable spatial distributions and pixel-wise trends over greenness pixels).
- `Validate_MODIS_Green_Fraction_Via_Sentinel2.ipynb`, `Visualization of Validate_MODIS_With_Without_Green_Via_Sentinel2.ipynb` — **Fig. S15** (MODIS–Sentinel-2 cross-sensor agreement).
- `SM_MODIS_AREA_ROCK_83_DEM.ipynb` — **Fig. S16** (MODIS green area under rock / latitude / DEM masks).
- `Visualization of MODIS LC and CC.ipynb`, `SM_Land_Cover.ipynb` — **Fig. S17** (MOD08/MYD08 land-observation and cloud-fraction trends).
- `Simulation of LC Map.ipynb` — **Fig. S18** (synthetic DART land-cover mosaics).
- `DEM_Generation.ipynb`, `SM_DEM.ipynb` — **Fig. S19** (synthetic DART DEMs).
- `SM_DART_Visualization_Ice_vs_Minimum_Vegetation.ipynb` — **Fig. S20** (DART minimum vegetation fraction under ice cover; Data S2).
- `SM_DART_Visualization_View_Az_Angle.ipynb` — **Fig. S21** (DART NDVI vs azimuth and view geometry; Data S3).

---

## Data sources

### MODIS
MODIS reflectance can be accessed via Google Earth Engine (GEE) following the Methods in the manuscript:
- **`MCD43A4.061 MODIS Nadir BRDF-Adjusted Reflectance Daily 500m`**

### ERA5 environmental drivers
Acquisition and preprocessing are implemented in:
- `Download ERA5.ipynb`, `Process ERA5 data to tif.ipynb`
- `1-preprocess_rf_data.R` (decomposition into spatial / temporal / residual components)

### UAV (AAS 4046, Windmill Islands)
The `figures/plot_uav_three_site_subpixel.py` script consumes cached `*_uav_ndvi.npz` files derived from the Australian Antarctic Data Centre AAS 4046 orthomosaics. Those orthomosaics are published separately by the AAD and are not redistributed here.

---

## Citation

If you use this code or any derived product, please cite the associated manuscript.

---

## Contact

- Min Chen — min.chen@wisc.edu
- Hangkai You — hyou34@wisc.edu

---

*This README was summarized by Claude Code and reviewed by Hangkai You.*
