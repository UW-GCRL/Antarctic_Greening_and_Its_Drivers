import pandas as pd
import numpy as np
import joblib
import itertools
import pymannkendall as mk

# === Paths ===
data_path = '/home/hyou34/RF/decomp_df.csv'
detrend_path_template = '/home/hyou34/RF/decomposed/detrended_{}.csv'
regions = ['West', 'East', 'Peninsula']
model_base_dir = '/home/hyou34/RF/Regional_results_permutation'

# === Load full base dataset ===
df_all = pd.read_csv(data_path)
base_variables = ['temperature_2m', 'solar_radiation', 'uv_radiation', 'volumetric_soil_water',
                  'runoff', 'X10m_wind_speed', 'precipitation', 'vapor_pressure_deficit']

# === Get all variable combinations ===
all_combinations = list(itertools.chain(*[itertools.combinations(base_variables, r) for r in range(1, len(base_variables) + 1)]))

# === Loop over regions ===
for region in regions:
    print(f"\nProcessing region: {region}")
    trend_results = []

    # Build region-specific model paths
    model_paths = [
        f'{model_base_dir}/{region}/RF_model_fold{i+1}.joblib'
        for i in range(5)
    ]

    # Filter region data
    df_original = df_all[df_all['Regions'] == region].copy()

    for combination in all_combinations:
        yearly_areas = []

        for model_path in model_paths:
            rf_model = joblib.load(model_path)
            df_full = df_original.copy()

            # Load detrended variables
            for var in combination:
                detrended = pd.read_csv(detrend_path_template.format(var))
                for component in ['temporal', 'residual']:
                    col_name = f"{var}_{component}"
                    if col_name in detrended.columns:
                        df_full[col_name] = detrended[col_name]

            # Predict
            predictors = [col for col in rf_model.feature_names_in_ if col in df_full.columns]
            X = df_full[predictors]
            df_full['predicted_ratio'] = rf_model.predict(X)
            df_full['predicted_area'] = df_full['land_area'] * df_full['predicted_ratio']

            # Group by year and store for averaging
            yearly_area = df_full.groupby('year')['predicted_area'].sum().reset_index()
            yearly_areas.append(yearly_area.set_index('year')['predicted_area'])

        # Combine folds into a DataFrame: each column is a fold
        combined_df = pd.concat(yearly_areas, axis=1)
        combined_df.columns = [f"fold_{i+1}" for i in range(len(yearly_areas))]

        # Compute mean across folds
        combined_df['mean_predicted_area'] = combined_df.mean(axis=1)

        # Apply Mann-Kendall test on mean
        result = mk.original_test(combined_df['mean_predicted_area'])

        # Store single result
        trend_results.append({
            'Region': region,
            'Variables': ', '.join(combination),
            'Mean_p_value': result.p,
            'Mean_slope': result.slope
        })

    # === Save per-region results ===
    trend_df = pd.DataFrame(trend_results)
    out_path = f'/home/hyou34/RF/detrend_mk_pvalues_avg_std_{region}.csv'
    trend_df.to_csv(out_path, index=False)
    print(f"Saved {region} results to {out_path}")