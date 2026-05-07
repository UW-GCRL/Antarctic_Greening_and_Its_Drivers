# Shapley-based driver attribution of the Antarctic greening trend

This folder contains the data-driven Shapley-attribution pipeline that
quantifies how the 2002–2023 trend in each climate variable contributes to the
modeled continental and regional greening trend in Antarctica.

The approach treats the trained Random Forest plus the area-weighted regional
aggregation as a black-box predictor mapping per-pixel climate state to a
yearly greening-area time series. We perturb the inputs by per-pixel
linear-detrending each variable (centered on the period mean year, preserving
each pixel's mean and interannual variability), re-feed through the model,
and aggregate. Across all 2⁸ = 256 variable subsets, exact Shapley values
fairly partition the model's trend-reduction across the eight climate
variables.

This pipeline supersedes the earlier anchor-at-2002 detrending (see
`Detrend_continent.py` and `Detrend_region.py` at the repo root) which pushed
inputs outside the model's training distribution.

---

## Pipeline overview

| Stage | Script | Purpose |
|---|---|---|
| 1 | `build_detrend_curves_perpixel.py` | Per-pixel detrend + 256 model evaluations × 5 folds × 4 regions; produces `<Region>_per_fold.csv` and `<Region>_detrend_curves.csv`. |
| 2a | `shapley_continental.py` | Exact Shapley over all 256 subsets, Continental scale. |
| 2b | `analyze_regional.py` | Same for each of the four regions, plus baseline/all-8 Mann–Kendall p-values. |
| 3 | `compute_mk_continental.py` | Per-scenario Mann–Kendall p-value and Sen's slope; optional comparison to a legacy MK file via `AGID_OLD_MK_FILE`. |
| 4 | `build_var_trend_table.py` | Region × variable area-weighted per-pixel trends with MK significance stars. |
| 5 | `build_shapley_tables_with_std.py` | Region × variable Shapley tables (absolute km²/yr and % of climate-explained), mean ± std across folds. |
| 6 | `plot_shapley_tables.py` | 2-panel figure of the Shapley tables. |
| 7 | `plot_combined_attribution_figure.py` | 3-panel publication figure (variable trends + Shapley absolute + Shapley %). |

Running stages 1–2 is enough for the headline numbers; stages 3–7 produce
supplementary tables and figures.

---

## Configuration

All scripts read paths from `paths.py`, which prefers environment variables.
Set the following before running:

```bash
# Folder with the per-variable detrended_<v>.csv decomposition files
# (output of 1-preprocess_rf_data.R at the repo root).
export AGID_DATA="/path/to/cleared_data"

# Folder with the trained RF .joblib files in the layout:
#   $AGID_MODELS_DIR/Continental_results_permutation/RF_model_fold{1..5}.joblib
#   $AGID_MODELS_DIR/Regional_results_permutation/{East,West,Peninsula}/RF_model_fold{1..5}.joblib
export AGID_MODELS_DIR="/path/to/models"

# Output folder. Defaults to ./DetrendCurves_PerPixelCentered next to this code.
export AGID_OUTPUT_DIR="/path/to/output"

# OPTIONAL: legacy MK file to enable NEW-vs-OLD significance comparison
# in compute_mk_continental.py. Skip if you don't have it.
export AGID_OLD_MK_FILE="/path/to/detrend_mk_pvalues_avg_std.csv"
```

On Windows PowerShell:
```powershell
$env:AGID_DATA = "G:\path\to\cleared_data"
$env:AGID_MODELS_DIR = "G:\path\to\models"
$env:AGID_OUTPUT_DIR = "G:\path\to\output"
```

---

## Running

After setting environment variables and `pip install -r requirements.txt`:

```bash
# Stage 1: ~1.5–2 hr full pipeline (Continental fold ~30 min × 5 folds; regional faster)
python build_detrend_curves_perpixel.py

# Stage 2: seconds
python shapley_continental.py
python analyze_regional.py

# Stages 3–7: seconds–minutes, in any order after stages 1–2
python compute_mk_continental.py
python build_var_trend_table.py
python build_shapley_tables_with_std.py
python plot_shapley_tables.py
python plot_combined_attribution_figure.py
```

---

## Outputs

In `$AGID_OUTPUT_DIR` (default: `./DetrendCurves_PerPixelCentered`):

- `<Region>_per_fold.csv` — long-format yearly area for each (region, scenario, fold, year)
- `<Region>_detrend_curves.csv` — mean ± std across folds per scenario per year
- `Continental_shapley.csv` and `<Region>_shapley.csv` — per-region Shapley attribution
- `Continental_mk_perpixel.csv` — Mann–Kendall per scenario (Continental)
- `variable_trends_per_region.csv` — area-weighted variable trends per region
- `shapley_table_abs_with_std.csv`, `shapley_table_pct_with_std.csv` — region × variable wide tables
- `Variable_trends_per_region.{png,pdf}`
- `Shapley_attribution_table.{png,pdf}`
- `Combined_attribution_figure.{png,pdf}`

---

## Reading the Shapley sign

Each Shapley value `φ_v` is in km²/yr.
- **Positive φ** = the variable's 2002–2023 trend has, on net, *helped* the
  greening trend; if its trend were absent, predicted greening would slow by
  that amount.
- **Negative φ** = the variable's trend has, on net, *slowed* the greening
  trend.

The sign is **not** the variable's own direction of change. It is the product
of (variable's trend direction) × (variable's marginal effect on vegetation).
Example: UV radiation declined over 2002–2023 and harms vegetation
(negative effect); the two negatives multiply to a positive Shapley value,
meaning UV decline has helped greening.

By the Shapley efficiency property, `Σ_v φ_v = baseline_slope − all-detrended_slope`,
i.e. the climate-explained portion of each region's greening rate.

---

## Method note

This is a model-based / data-driven driver attribution, not formal
Detection & Attribution in the Hasselmann/IPCC sense. The Shapley values
reflect what the random forest's learned response surface attributes to
counterfactual per-pixel detrending, not direct causal effects from a
structural climate model. The RF underestimates the observed continental
trend (yearly R² = 0.98 but slope is dampened by ~23%, a well-known
property of tree ensembles), so the Shapley sum represents a *lower bound*
on the climate-driven fraction of observed greening.
