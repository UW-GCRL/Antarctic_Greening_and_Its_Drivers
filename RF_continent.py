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
output_path = '/home/hyou34/RF/Continental_results_permutation'
os.makedirs(output_path, exist_ok=True)

# === Load full dataset ===
df = pd.read_csv(data_path)

# === Drop irrelevant columns ===
df_clean = df.drop(['Num', 'Regions', 'Subregions', 'pixel', 'year', 'latitude', 'longitude', 'land_area'], axis=1)

# === Separate features and target ===
X = df_clean.drop('vegetation_area_ratio', axis=1).reset_index(drop=True)
y = df_clean['vegetation_area_ratio'].reset_index(drop=True)

# === Setup 5-fold cross-validation ===
cv = KFold(n_splits=5, shuffle=True, random_state=42)

# === Storage for evaluation results ===
r2_scores = []
rmse_scores = []
perm_scores_list = []
all_preds = []  # to collect predicted and observed values for all test data

# === Cross-validation loop ===
for fold_idx, (train_idx, test_idx) in enumerate(cv.split(X)):
    print(f"\n--- Fold {fold_idx + 1} ---")

    # Split training and test data
    X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
    y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

    # Train Random Forest model
    model = RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=8)
    model.fit(X_train, y_train)

    # Predict and evaluate
    y_pred = model.predict(X_test)
    r2 = r2_score(y_test, y_pred)
    rmse = mean_squared_error(y_test, y_pred, squared=False)
    r2_scores.append(r2)
    rmse_scores.append(rmse)

    # Save test predictions for later plotting
    fold_pred_df = pd.DataFrame({
        'Index': test_idx,
        'Fold': fold_idx + 1,
        'Observed': y_test.values,
        'Predicted': y_pred
    })
    all_preds.append(fold_pred_df)

    # Compute permutation importance
    perm_result = permutation_importance(model, X_test, y_test, n_repeats=10, random_state=42, n_jobs=8)
    perm_scores_list.append(perm_result.importances_mean)

    # Save the model for this fold
    joblib.dump(model, os.path.join(output_path, f"RF_model_fold{fold_idx + 1}.joblib"))

# === Combine all predictions and save to CSV ===
all_preds_df = pd.concat(all_preds).sort_values(by='Index')
all_preds_df.to_csv(os.path.join(output_path, "all_test_predictions.csv"), index=False)

# === Compute mean and std of performance metrics ===
mean_r2 = np.mean(r2_scores)
std_r2 = np.std(r2_scores)
mean_rmse = np.mean(rmse_scores)
std_rmse = np.std(rmse_scores)

print(f"\n[Performance] Mean R2: {mean_r2:.4f} +/- {std_r2:.4f}")
print(f"[Performance] Mean RMSE: {mean_rmse:.4f} +/- {std_rmse:.4f}")

# Save metrics to file
with open(os.path.join(output_path, "cv_scores.txt"), "w") as f:
    f.write(f"Mean R2: {mean_r2:.4f} +/- {std_r2:.4f}\n")
    f.write(f"Mean RMSE: {mean_rmse:.4f} +/- {std_rmse:.4f}\n")

# === Average permutation importance across folds ===
perm_scores_array = np.vstack(perm_scores_list)
perm_mean = perm_scores_array.mean(axis=0)
perm_std = perm_scores_array.std(axis=0)

perm_df = pd.DataFrame({
    'Feature': X.columns,
    'Importance_Mean': perm_mean,
    'Importance_STD': perm_std
}).sort_values(by='Importance_Mean', ascending=False)

# Save permutation importance to CSV
perm_df.to_csv(os.path.join(output_path, "permutation_importance.csv"), index=False)

print(f"Saved full continental CV results and predictions to {output_path}")