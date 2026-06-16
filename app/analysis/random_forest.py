"""
Random Forest classification for identifying discriminant taxa between groups.
"""
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_predict, StratifiedKFold
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.preprocessing import LabelEncoder


def run_random_forest(count_df: pd.DataFrame, meta_df: pd.DataFrame,
                      sid_col: str, group_col: str,
                      n_estimators: int = 500, top_n: int = 30,
                      cv_folds: int = 5, random_state: int = 42) -> dict:
    """Run Random Forest classification and return feature importances.

    Parameters
    ----------
    count_df : DataFrame
        Features (rows) x samples (columns).
    meta_df : DataFrame
        Metadata.
    sid_col : str
        Sample ID column.
    group_col : str
        Group/class column.
    n_estimators : int
        Number of trees.
    top_n : int
        Number of top important features to report.
    cv_folds : int
        Number of cross-validation folds.
    random_state : int
        Random seed.

    Returns
    -------
    dict with keys: importance_df, accuracy, cv_accuracy, classification_report,
                    confusion_matrix, class_labels, predictions
    """
    # Align samples
    meta_indexed = meta_df.set_index(meta_df[sid_col].astype(str))
    common = sorted(set(count_df.columns) & set(meta_indexed.index))

    if len(common) < 6:
        raise ValueError(f"Need at least 6 matching samples, found {len(common)}")

    X = count_df[common].T  # samples x features
    y_raw = meta_indexed.loc[common, group_col].values

    # Encode labels
    le = LabelEncoder()
    y = le.fit_transform(y_raw)
    class_labels = le.classes_.tolist()

    # Check minimum samples per class for CV
    unique, counts = np.unique(y, return_counts=True)
    min_class_size = counts.min()
    actual_folds = min(cv_folds, min_class_size)
    if actual_folds < 2:
        raise ValueError(
            f"Smallest group has {min_class_size} sample(s). "
            f"Need at least 2 samples per group for cross-validation."
        )

    # Train full model
    rf = RandomForestClassifier(
        n_estimators=n_estimators,
        random_state=random_state,
        n_jobs=-1,
    )
    rf.fit(X, y)

    # Feature importances
    importances = rf.feature_importances_
    feature_names = X.columns.tolist()
    imp_df = pd.DataFrame({
        "Feature": feature_names,
        "Importance": importances,
    }).sort_values("Importance", ascending=False)

    # Cross-validation
    cv = StratifiedKFold(n_splits=actual_folds, shuffle=True, random_state=random_state)
    y_pred_cv = cross_val_predict(rf, X, y, cv=cv)
    cv_accuracy = accuracy_score(y, y_pred_cv)

    # Full-model predictions and metrics
    y_pred_full = rf.predict(X)
    full_accuracy = accuracy_score(y, y_pred_full)

    report = classification_report(y, y_pred_cv, target_names=class_labels, output_dict=True)
    cm = confusion_matrix(y, y_pred_cv)

    return {
        "importance_df": imp_df.head(top_n),
        "full_importance_df": imp_df,
        "accuracy": float(full_accuracy),
        "cv_accuracy": float(cv_accuracy),
        "cv_folds": actual_folds,
        "classification_report": report,
        "confusion_matrix": cm.tolist(),
        "class_labels": class_labels,
        "predictions": le.inverse_transform(y_pred_cv).tolist(),
        "actual": y_raw.tolist(),
    }
