import pandas as pd
import numpy as np
import os
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import r2_score, mean_squared_error
import joblib

# === Paths ===
data_path = '/home/hyou34/RF/decomp_df.csv'
base_output_path = '/home/hyou34/RF/Regional_results_permutation'
os.makedirs(base_output_path, exist_ok=True)

# === Load full dataset ===
df = pd.read_csv(data_path)

# === Define regions to loop over ===
regions = ['East', 'West', 'Peninsula']

# === Loop through each region ===
for region in regions:
    print(f"\n===== Processing Region: {region} =====")

    # === Create output folder for this region ===
    region_output = os.path.join(base_output_path, region)
    os.makedirs(region_output, exist_ok=True)

    # === Filter and clean data ===
    region_df = df[df['Regions'] == region].copy()
    region_df.drop(['Num', 'Regions', 'Subregions', 'pixel', 'year', 'latitude', 'longitude', 'land_area'], axis=1, inplace=True)

    # === Separate features and target ===
    X = region_df.drop('vegetation_area_ratio', axis=1).reset_index(drop=True)
    y = region_df['vegetation_area_ratio'].reset_index(drop=True)

    # === Setup 5-fold cross-validation ===
    cv = KFold(n_splits=5, shuffle=True, random_state=42)
    r2_scores = []
    rmse_scores = []
    perm_scores_list = []
    all_preds = []  # to store predictions and true values from all folds

    for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X)):
        print(f"\n--- Fold {fold_idx + 1} ---")
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        model = RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=8)
        model.fit(X_train, y_train)

        y_pred = model.predict(X_test)
        r2 = r2_score(y_test, y_pred)
        rmse = mean_squared_error(y_test, y_pred, squared=False)
        r2_scores.append(r2)
        rmse_scores.append(rmse)

        # === Save test predictions for this fold ===
        fold_pred_df = pd.DataFrame({
            'Index': test_idx,
            'Fold': fold_idx + 1,
            'Observed': y_test.values,
            'Predicted': y_pred
        })
        all_preds.append(fold_pred_df)

        # === Permutation importance ===
        perm_result = permutation_importance(model, X_test, y_test, n_repeats=10, random_state=42, n_jobs=8)
        perm_scores_list.append(perm_result.importances_mean)

        # === Save model ===
        joblib.dump(model, os.path.join(region_output, f"RF_model_fold{fold_idx + 1}.joblib"))

    # === Save all fold test predictions ===
    all_preds_df = pd.concat(all_preds).sort_values(by='Index')
    all_preds_df.to_csv(os.path.join(region_output, "all_test_predictions.csv"), index=False)

    # === Compute average and std of evaluation metrics ===
    mean_r2 = np.mean(r2_scores)
    std_r2 = np.std(r2_scores)
    mean_rmse = np.mean(rmse_scores)
    std_rmse = np.std(rmse_scores)

    print(f"\n[Performance] Mean R2: {mean_r2:.4f} +/- {std_r2:.4f}")
    print(f"[Performance] Mean RMSE: {mean_rmse:.4f} +/- {std_rmse:.4f}")

    # === Save metrics to file ===
    with open(os.path.join(region_output, "cv_scores.txt"), "w") as f:
        f.write(f"Mean R2: {mean_r2:.4f} +/- {std_r2:.4f}\n")
        f.write(f"Mean RMSE: {mean_rmse:.4f} +/- {std_rmse:.4f}\n")

    # === Compute average permutation importance ===
    perm_scores_array = np.vstack(perm_scores_list)
    perm_mean = perm_scores_array.mean(axis=0)
    perm_std = perm_scores_array.std(axis=0)

    perm_df = pd.DataFrame({
        'Feature': X.columns,
        'Importance_Mean': perm_mean,
        'Importance_STD': perm_std
    }).sort_values(by='Importance_Mean', ascending=False)

    # === Save permutation importance to CSV ===
    perm_df.to_csv(os.path.join(region_output, "permutation_importance.csv"), index=False)

    print(f"Saved cross-validation results for {region} to {region_output}")