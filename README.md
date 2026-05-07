# Antarctic Greening and Its Drivers

This repository contains code and notebooks for mapping Antarctic greenness (MODIS NDVI) and quantifying its environmental drivers using Random Forest modeling and detrending-based attribution analyses.

> The workflow was originally developed and run with **Python 3.8**.

---

## What’s in this repo

- `RF_continent.py` — Random Forest modeling at the continental scale  
- `RF_region.py` — Random Forest modeling at the regional scale  
- `Detrend_continent.py` — Detrending-based attribution (continent, legacy anchor-at-2002)  
- `Detrend_region.py` — Detrending-based attribution (regions, legacy anchor-at-2002)  
- `shapley_attribution/` — **Shapley-based driver attribution pipeline** (data-driven; supersedes the legacy detrending). Per-pixel-centered detrending of the eight climate variables, all 256 detrending coalitions evaluated through the trained RF, and exact Shapley values fairly partition each variable's contribution to the modeled continental and regional greening trend. See `shapley_attribution/README.md` for usage.
- `1-preprocess_rf_data.R` — Environmental driver preprocessing + decomposition (used as RF inputs)
- `Antarctic_NDVI_DIY_NBAR.ipynb` — NDVI workflow based on MODIS NBAR reflectance  
- `Antarctic_greening_mapping.ipynb` — Greening mapping / product generation  
- `MODIS_Area_Rock&DEM&Lat_Map.ipynb` — Supporting maps (area/rock/DEM/lat)  
- `Duration using last and first.ipynb` — Duration calculation
- `Download ERA5.ipynb` — ERA5 download helper notebook  
- `Process ERA5 data to tif.ipynb` — Convert/process ERA5 data into GeoTIFFs
- `MOD10A1_Antarctica_MinSnowCover_clean.js` — GEE script for MODIS snow cover product handling

**These notebooks support sensitivity analyses, robustness checks, and visualizations used in the manuscript and supplementary materials:**
- `SM_DART_Visualization_Ice_vs_Minimum_Vegetation.ipynb`
- `SM_DART_Visualization_View_Az_Angle.ipynb`
- `SM_Land_Cover.ipynb`
- `DEM_Generation.ipynb`
- `SM_DEM.ipynb`
- `Visualization of MODIS LC and CC.ipynb`
- `Simulation of LC Map.ipynb`
- `SM_Duration_Day_Count.ipynb`
- `Supplementary _materials_mapping_duration.ipynb`
- `SM_Trend robustness without certain years.ipynb`
- `SM_Trend_of_Key_Environmental_Variables.ipynb`
- `SM_MODIS_AREA_ROCK_83_DEM.ipynb`
- `Supplementary _materials_mapping_area.ipynb`
- `SM_NDVI_Sensitivity_New_Stable_Disappeared_Greenness_of_Hot_Spot_Region.ipynb`
- `Supplementary _materials_Sensitivity_test.ipynb`
- `Validate_MODIS_Green_Fraction_Via_Sentinel2.ipynb`
- `Visualization of Validate_MODIS_With_Without_Green_Via_Sentinel2.ipynb`
- `Visualization of AG RF results.ipynb`

---

## Data sources

### MODIS
MODIS reflectance can be accessed via Google Earth Engine (GEE) following the Methods in the manuscript:
- **`MCD43A4.061 MODIS Nadir BRDF-Adjusted Reflectance Daily 500m`**

### Environmental drivers
Environmental driver preprocessing & decomposition are implemented in:
- `1-preprocess_rf_data.R`

ERA5 download / processing notebooks are provided in this repo (see above).

---

## Citation

If you use this code or derived products, please cite the associated manuscript.
---

## Contact

- Min Chen — min.chen@wisc.edu  
- Hangkai You — hyou34@wisc.edu  

---
