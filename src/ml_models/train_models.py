from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.linear_model import LogisticRegression
from xgboost import XGBClassifier
from sklearn.linear_model import LogisticRegression
from lightgbm import LGBMClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.metrics import classification_report, accuracy_score, f1_score
from sklearn.exceptions import ConvergenceWarning
import os
import json
os.environ["LIGHTGBM_VERBOSE"] = "0"
import pandas as pd
import warnings

def get_models_and_grids():
    models = {
        "random_forest": (RandomForestClassifier(), {

    'n_estimators': [50, 100, 200],         # Number of trees in the forest
    'max_depth': [None, 10, 20],            # Maximum depth of the tree
    'min_samples_split': [2, 5, 10],        # Minimum samples to split an internal node
    'min_samples_leaf': [1, 2, 4], 
    'oob_score': [True]         # Minimum samples at a leaf node
}
        ),
        "lightgbm": (LGBMClassifier(verbosity=-1), {
           
                # tree complexity
                "num_leaves":    [31, 63],       # size of individual trees
                "max_depth":     [-1, 10],       # -1 = no limit, or shallow trees

                # learning
                "learning_rate": [0.01, 0.1],    # slow vs. fast
                "n_estimators":  [100, 200],     # number of boosting rounds
            }),
             "svm": 
            (SVC(random_state=42, probability=False),  # drop probability=True for speed
            {
                "C":      [0.1, 1],     # just three scales
                "kernel": ["linear", "rbf"],# two kernels
                "gamma":  ["scale"]         # only the default
            }),
            "mlp": (
            MLPClassifier(
                random_state=42,
                max_iter=500,      # keep training time reasonable
                early_stopping=True
            ),
            {
                # only two sizes instead of four
                "hidden_layer_sizes": [(50,), (100,)],
                # the two most common activations
                "activation": ["relu", "tanh"],
                # two L2 penalties
                "alpha": [1e-4, 1e-3],
                # adaptive lets the solver adjust the learning rate
                "learning_rate": ["adaptive"]
            }),
            "xgboost": (XGBClassifier(tree_method="hist", eval_metric="mlogloss"), {
           "max_depth": [3,6],
           "learning_rate": [0.05,0.1],
           "n_estimators": [100,200]
       }),
         "logistic": (
            LogisticRegression(
                random_state=42,
                solver="saga",      # supports l1 & l2 with multinomial
                max_iter=2000
            ),
            {
                "C":       [0.01, 0.1, 1, 10],      # inverse regularization strength
                "penalty": ["l2", "l1"],           # L2 ridge vs. L1 lasso
                # note: with solver='saga' you can mix l1 & l2
            }
        ),   
    }
        # TabTransformer would be added later via a custom wrapper

    return models

def evaluate_model(model, X_test, y_test, cancer_type, feature_type, model_name):
    # predict (suppress UserWarning about feature names etc)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        warnings.simplefilter("ignore", FutureWarning)
        warnings.simplefilter("ignore", ConvergenceWarning)
        preds = model.predict(X_test)

    acc   = accuracy_score(y_test, preds)
    f1_w  = f1_score(y_test, preds, average="weighted", zero_division=0)
    # per-class metrics
    report = classification_report(y_test, preds, output_dict=True, zero_division=0)

    # create directory
    out_dir = os.path.join("results", cancer_type, feature_type)
    os.makedirs(out_dir, exist_ok=True)

    # save CSV
    df = pd.DataFrame(report).T
    df["support"] = df["support"].astype(int)
    df.to_csv(os.path.join(out_dir, f"{model_name}_class_report.csv"))

    return {
        "accuracy":    acc,
        "f1_weighted": f1_w,
        "report":      report
    }


def train_all_models(X_train, y_train, X_test, y_test, cancer_type, feature_type):
    # Assign column names if not a DataFrame
    if not isinstance(X_train, pd.DataFrame):
        feature_names = [f"f{i}" for i in range(X_train.shape[1])]
        X_train = pd.DataFrame(X_train, columns=feature_names)
        X_test = pd.DataFrame(X_test, columns=feature_names)

    models = get_models_and_grids()


    results = {}

    for name, (model, param_grid) in get_models_and_grids().items():
        print(f"🔍 Training model: {name} on {feature_type}")
        grid = GridSearchCV(model, param_grid, cv=2, n_jobs=-1, scoring='accuracy')
        grid.fit(X_train, y_train)
        best = grid.best_estimator_

        prefix = f"results/{cancer_type}/{feature_type}/{name}"
        results[name] = evaluate_model(
            best,
            X_test, y_test,
            cancer_type=cancer_type,
            feature_type=feature_type,
            model_name=name
        )

    return results
