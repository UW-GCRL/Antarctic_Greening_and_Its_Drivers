"""Compute Mann-Kendall p-value + Sen's slope per scenario per fold for
each region (Continental, East, West, Peninsula) under the per-pixel
mean-centered detrending. Saves one CSV per region in the same folder as
the per-fold area files.

Each row of the output CSV is one (scenario) aggregate across the 5
cross-validation folds:
    scenario, mean_p_value, std_p_value, mean_sen_slope, std_sen_slope

These feed the regional significance-matrix figure (formerly Fig. 4 in
main, now Fig. S23).
"""

from pathlib import Path
import numpy as np
import pandas as pd
import pymannkendall as mk

NEW = Path(r"G:\Hangkai\Antarctica_Mapping_Data\New\DetrendCurves_PerPixelCentered")

REGIONS = ["Continental", "East", "West", "Peninsula"]


def per_scenario_mk(per_fold_csv):
    pf = pd.read_csv(per_fold_csv)
    rows = []
    for sname, grp in pf.groupby("scenario"):
        pvals, slopes = [], []
        for fold, g in grp.groupby("fold"):
            g = g.sort_values("year")
            result = mk.original_test(g["area"].values)
            pvals.append(result.p)
            slopes.append(result.slope)
        rows.append({
            "scenario":        sname,
            "mean_p_value":    float(np.mean(pvals)),
            "std_p_value":     float(np.std(pvals)),
            "mean_sen_slope":  float(np.mean(slopes)),
            "std_sen_slope":   float(np.std(slopes)),
        })
    return pd.DataFrame(rows)


for region in REGIONS:
    pf_path = NEW / f"{region}_per_fold.csv"
    out_path = NEW / f"{region}_mk_perpixel.csv"
    if not pf_path.exists():
        print(f"[skip] {pf_path} not found")
        continue
    df = per_scenario_mk(pf_path)
    df.to_csv(out_path, index=False)
    n_sig = int((df.mean_p_value < 0.05).sum())
    print(f"[{region:<11}] {len(df):3d} scenarios  |  {n_sig} significant (p<0.05)  |  saved {out_path.name}")
