"""Path configuration for the Shapley-attribution pipeline.

All scripts in this folder import from this module so that input/output
locations are configurable via environment variables. Set them once in your
shell (or .env), or edit the defaults below.

Environment variables (recommended):
    AGID_DATA           — folder containing the per-variable detrended_<v>.csv
                          decomposition files (output of 1-preprocess_rf_data.R).
                          Required.
    AGID_MODELS_DIR     — folder containing trained RF .joblib model files
                          organized as:
                              <AGID_MODELS_DIR>/Continental_results_permutation/RF_model_fold{1..5}.joblib
                              <AGID_MODELS_DIR>/Regional_results_permutation/East/RF_model_fold{1..5}.joblib
                              <AGID_MODELS_DIR>/Regional_results_permutation/West/RF_model_fold{1..5}.joblib
                              <AGID_MODELS_DIR>/Regional_results_permutation/Peninsula/RF_model_fold{1..5}.joblib
                          Required for build_detrend_curves_perpixel.py.
    AGID_OUTPUT_DIR     — folder where outputs are written. Optional; defaults
                          to ./DetrendCurves_PerPixelCentered alongside this code.
    AGID_OLD_MK_FILE    — (optional) path to a legacy detrend_mk_pvalues_avg_std.csv
                          if you want compute_mk_continental.py to print a
                          "new vs old" significance comparison. Skip this var
                          to skip the comparison.
"""

import os
from pathlib import Path

DATA = Path(os.environ.get(
    "AGID_DATA",
    "data/cleared_data",  # placeholder default; override via env var
))

MODELS_DIR = Path(os.environ.get(
    "AGID_MODELS_DIR",
    "data/models",
))

OUTPUT_DIR = Path(os.environ.get(
    "AGID_OUTPUT_DIR",
    Path(__file__).resolve().parent / "DetrendCurves_PerPixelCentered",
))

OLD_MK_FILE = os.environ.get("AGID_OLD_MK_FILE")  # optional; None if not set
OLD_MK_FILE = Path(OLD_MK_FILE) if OLD_MK_FILE else None

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
