"""Compute Mann-Kendall p-values per (scenario, fold) for the Continental
per-fold output, then average across the 5 folds. Optionally compares with a
legacy detrend_mk_pvalues_avg_std.csv if AGID_OLD_MK_FILE is set, to highlight
how significance status changes between methods."""

import numpy as np
import pandas as pd
import pymannkendall as mk

from paths import OUTPUT_DIR, OLD_MK_FILE

VARS = ["temperature_2m", "uv_radiation", "solar_radiation",
        "volumetric_soil_water", "runoff", "X10m_wind_speed",
        "precipitation", "vapor_pressure_deficit"]
ALL_NAME = ",".join(sorted(VARS))


def main():
    pf = pd.read_csv(OUTPUT_DIR / "Continental_per_fold.csv")

    rows = []
    for sname, grp in pf.groupby("scenario"):
        pvals, slopes = [], []
        for fold, g in grp.groupby("fold"):
            g = g.sort_values("year")
            result = mk.original_test(g["area"].values)
            pvals.append(result.p)
            slopes.append(result.slope)
        rows.append({
            "scenario": sname,
            "mean_p_value": np.mean(pvals),
            "std_p_value": np.std(pvals),
            "mean_sen_slope": np.mean(slopes),
            "std_sen_slope": np.std(slopes),
        })
    df_new = pd.DataFrame(rows)
    out_path = OUTPUT_DIR / "Continental_mk_perpixel.csv"
    df_new.to_csv(out_path, index=False)
    print(f"Saved: {out_path}  ({len(df_new)} scenarios)")

    if OLD_MK_FILE is None or not OLD_MK_FILE.exists():
        print("\n(Skipping NEW-vs-OLD comparison: set AGID_OLD_MK_FILE to enable.)")
        return

    old = pd.read_csv(OLD_MK_FILE)
    old_clean = old.copy()

    # Old uses ", " separator and original variable order; new uses ","
    # alphabetical. Normalize both to alphabetical-comma-joined.
    def _norm(s):
        parts = [p.strip() for p in str(s).split(",")]
        return ",".join(sorted(parts))

    old_clean["scenario"] = old_clean["variables"].apply(_norm)

    merged = df_new.merge(
        old_clean[["scenario", "mean_p_value", "std_p_value", "mean_slope"]],
        on="scenario", how="left", suffixes=("_new", "_old"),
    )

    print()
    print(f"{'scenario':<55} {'p_new':>8} {'p_old':>8} {'sig_new':>8} {'sig_old':>8} {'changed':>11}")
    print("-" * 110)
    SIG = 0.05
    n_changed = n_no_longer_sig = n_newly_sig = 0
    for _, r in merged.iterrows():
        sname = r["scenario"]
        pn, po = r["mean_p_value_new"], r["mean_p_value_old"]
        sig_n = "YES" if pn < SIG else "no"
        sig_o = ("YES" if po < SIG else "no") if pd.notna(po) else "?"
        changed = ""
        if pd.notna(po):
            if (pn < SIG) != (po < SIG):
                changed = "**CHANGED**"
                n_changed += 1
                if pn >= SIG and po < SIG:
                    n_no_longer_sig += 1
                elif pn < SIG and po >= SIG:
                    n_newly_sig += 1
        label = "all-8" if sname == ALL_NAME else (sname if len(sname) < 50 else sname[:47] + "...")
        po_str = f"{po:.4f}" if pd.notna(po) else "  N/A "
        print(f"{label:<55} {pn:>8.4f} {po_str:>8} {sig_n:>8} {sig_o:>8} {changed:>11}")

    print()
    print(f"Total scenarios where significance status changed: {n_changed}")
    print(f"  No-longer-significant (was sig, now not): {n_no_longer_sig}")
    print(f"  Newly-significant (was not sig, now sig): {n_newly_sig}")
    print(f"  Total compared: {sum(merged['mean_p_value_old'].notna())}")


if __name__ == "__main__":
    main()
