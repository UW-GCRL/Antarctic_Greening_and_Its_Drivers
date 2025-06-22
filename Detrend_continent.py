import pandas as pd
import numpy as np
import joblib
import itertools
import pymannkendall as mk

# === Paths ===
data_path = '/home/hyou34/RF/decomp_df.csv'
detrend_path_template = '/home/hyou34/RF/decomposed/detrended_{}.csv'
model_paths = [
    fr'/home/hyou34/RF/Continental_results_permutation/RF_model_fold{i+1}.joblib'
    for i in range(5)
]

# === Load full base dataset ===
df_original = pd.read_csv(data_path)
base_variables = ['temperature_2m', 'solar_radiation', 'uv_radiation', 'volumetric_soil_water',
                  'runoff', 'X10m_wind_speed', 'precipitation', 'vapor_pressure_deficit']

# === Get all variable combinations ===
all_combinations = list(itertools.chain(*[itertools.combinations(base_variables, r) for r in range(1, len(base_variables) + 1)]))

# === Store results ===
trend_results = []

for combination in all_combinations:
    yearly_areas = []

    for model_path in model_paths:
        rf_model = joblib.load(model_path)
        df_full = df_original.copy()

        # Load detrended variables for current combination
        for var in combination:
            detrended = pd.read_csv(detrend_path_template.format(var))
            for component in ['temporal', 'residual']:
                col_name = f"{var}_{component}"
                if col_name in detrended.columns:
                    df_full[col_name] = detrended[col_name]

        # Predict using model (only valid columns)
        predictors = [col for col in rf_model.feature_names_in_ if col in df_full.columns]
        X = df_full[predictors]
        df_full['predicted_ratio'] = rf_model.predict(X)
        df_full['predicted_area'] = df_full['land_area'] * df_full['predicted_ratio']

        # Aggregate to yearly predicted area
        yearly_area = df_full.groupby('year')['predicted_area'].sum().reset_index()
        yearly_areas.append(yearly_area.set_index('year')['predicted_area'])

    # Combine all folds into one DataFrame (columns: fold1, fold2, ...)
    combined_df = pd.concat(yearly_areas, axis=1)
    combined_df.columns = [f"fold_{i+1}" for i in range(len(yearly_areas))]

    # Compute mean predicted area across folds for each year
    combined_df['mean_predicted_area'] = combined_df.mean(axis=1)

    # Run Mann-Kendall test on the mean series
    result = mk.original_test(combined_df['mean_predicted_area'])

    # Store results
    trend_results.append({
        'variables': ', '.join(combination),
        'mean_p_value': result.p,
        'mean_slope': result.slope
    })

# Convert to DataFrame and save
trend_df = pd.DataFrame(trend_results)
trend_df.to_csv('/home/hyou34/RF/detrend_mk_pvalues_avg_std.csv', index=False)

print("Saved detrend analysis with averaged p-values.")