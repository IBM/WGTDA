from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, mean_squared_error, log_loss
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from xgboost import XGBClassifier
from sklearn.neural_network import MLPClassifier
import numpy as np

def train_and_evaluate_svm(train_features, train_labels, test_features, test_labels, param_grid=None, cv=5, random_state=42):
    """
    Train and evaluate an SVM model with Grid Search for hyperparameter tuning.

    Parameters:
    - train_features: np.ndarray, training data features.
    - train_labels: np.ndarray, labels for training data.
    - test_features: np.ndarray, testing data features.
    - test_labels: np.ndarray, labels for testing data.
    - param_grid: dict, parameter grid for hyperparameter tuning.
    - cv: int, number of cross-validation folds.
    - random_state: int, random seed.

    Returns:
    - best_svm_model: Best trained SVM model.
    - metrics: Dictionary containing accuracy, MSE, log loss, and other metrics.
    """
    if param_grid is None:
        param_grid = {
            'C': [0.1, 1, 10, 100],
            'kernel': ['linear', 'rbf', 'poly'],
            'gamma': ['scale', 'auto']
        }

    # Grid Search for hyperparameter tuning
    grid_search = GridSearchCV(
        SVC(probability=True, random_state=random_state, verbose=False),
        param_grid,
        cv=cv,
        scoring='accuracy'
    )
    grid_search.fit(train_features, train_labels)

    best_svm_model = grid_search.best_estimator_
    print("Best parameters:", grid_search.best_params_)

    # Evaluate the model
    train_predictions = best_svm_model.predict(train_features)
    test_predictions = best_svm_model.predict(test_features)
    test_probabilities = best_svm_model.predict_proba(test_features)

    train_accuracy = accuracy_score(train_labels, train_predictions)
    accuracy = accuracy_score(test_labels, test_predictions)
    mse = mean_squared_error(test_labels, np.argmax(test_probabilities, axis=1))
    logloss = log_loss(test_labels, test_probabilities)
    classification_rep = classification_report(test_labels, test_predictions)
    conf_matrix = confusion_matrix(test_labels, test_predictions)
    print(f"Training Accuracy: {train_accuracy:.4f}")
    print(f"Test Accuracy: {accuracy:.4f}")
    print(f"Mean Squared Error: {mse:.4f}")
    print(f"Log Loss: {logloss:.4f}")
    print("\nClassification Report:\n", classification_rep)
    print("\nConfusion Matrix:\n", conf_matrix)

    metrics = {'accuracy': accuracy, 'mse': mse, 'logloss': logloss}

    return best_svm_model, metrics

def train_and_evaluate_xgboost(train_features, train_labels, test_features, test_labels, param_grid=None, cv=5, random_state=42):
    """
    Train and evaluate an XGBoost model with Grid Search for hyperparameter tuning.

    Parameters:
    - train_features: np.ndarray, training data features.
    - train_labels: np.ndarray, labels for training data.
    - test_features: np.ndarray, testing data features.
    - test_labels: np.ndarray, labels for testing data.
    - param_grid: dict, parameter grid for hyperparameter tuning.
    - cv: int, number of cross-validation folds.
    - random_state: int, random seed.

    Returns:
    - best_xgb_model: Best trained XGBoost model.
    - metrics: Dictionary containing accuracy, MSE, log loss, and other metrics.
    """
    if param_grid is None:
        param_grid = {
            'n_estimators': [50, 100, 200],
            'max_depth': [3, 5, 7],
            'learning_rate': [0.01, 0.1, 0.2],
            'subsample': [0.8, 1.0]
        }

    # Grid Search for hyperparameter tuning
    grid_search = GridSearchCV(
        XGBClassifier(use_label_encoder=False, eval_metric='mlogloss', random_state=random_state),
        param_grid,
        cv=cv,
        scoring='accuracy'
    )
    grid_search.fit(train_features, train_labels)

    best_xgb_model = grid_search.best_estimator_
    print("Best parameters:", grid_search.best_params_)

    # Evaluate the model
    train_predictions = best_xgb_model.predict(train_features)
    test_predictions = best_xgb_model.predict(test_features)
    test_probabilities = best_xgb_model.predict_proba(test_features)

    accuracy = accuracy_score(test_labels, test_predictions)
    mse = mean_squared_error(test_labels, np.argmax(test_probabilities, axis=1))
    logloss = log_loss(test_labels, test_probabilities)
    classification_rep = classification_report(test_labels, test_predictions)
    conf_matrix = confusion_matrix(test_labels, test_predictions)
    train_accuracy = accuracy_score(train_labels, train_predictions)

    print(f"Training Accuracy: {train_accuracy:.4f}")
    print(f"Test Accuracy: {accuracy:.4f}")
    print(f"Mean Squared Error: {mse:.4f}")
    print(f"Log Loss: {logloss:.4f}")
    print("\nClassification Report:\n", classification_rep)
    print("\nConfusion Matrix:\n", conf_matrix)

    metrics = {'accuracy': accuracy, 'mse': mse, 'logloss': logloss}

    return best_xgb_model, metrics

def train_and_evaluate_nn(train_features, train_labels, test_features, test_labels, param_grid=None, cv=5, random_state=42):
    """
    Train and evaluate a Neural Network (MLP) with Grid Search for hyperparameter tuning.

    Parameters:
    - train_features: np.ndarray, training data features.
    - train_labels: np.ndarray, labels for training data.
    - test_features: np.ndarray, testing data features.
    - test_labels: np.ndarray, labels for testing data.
    - param_grid: dict, parameter grid for hyperparameter tuning.
    - cv: int, number of cross-validation folds.
    - random_state: int, random seed.

    Returns:
    - best_nn_model: Best trained Neural Network model.
    - metrics: Dictionary containing accuracy, MSE, log loss, and other metrics.
    """
    if param_grid is None:
        param_grid = {
            'hidden_layer_sizes': [(50,), (100,), (50, 50), (100, 50)],
            'activation': ['relu', 'tanh', 'logistic'],
            'alpha': [0.0001, 0.001, 0.01],
            'learning_rate': ['constant', 'adaptive']
        }

    # Grid Search for hyperparameter tuning
    grid_search = GridSearchCV(
        MLPClassifier(max_iter=500, random_state=random_state),
        param_grid,
        cv=cv,
        scoring='accuracy'
    )
    grid_search.fit(train_features, train_labels)

    best_nn_model = grid_search.best_estimator_
    print("Best parameters:", grid_search.best_params_)

    # Evaluate the model
    train_predictions = best_nn_model.predict(train_features)
    test_predictions = best_nn_model.predict(test_features)
    test_probabilities = best_nn_model.predict_proba(test_features)

    accuracy = accuracy_score(test_labels, test_predictions)
    mse = mean_squared_error(test_labels, np.argmax(test_probabilities, axis=1))
    logloss = log_loss(test_labels, test_probabilities)
    classification_rep = classification_report(test_labels, test_predictions)
    conf_matrix = confusion_matrix(test_labels, test_predictions)
    train_accuracy = accuracy_score(train_labels, train_predictions)

    print(f"Training Accuracy: {train_accuracy:.4f}")
    print(f"Test Accuracy: {accuracy:.4f}")
    print(f"Mean Squared Error: {mse:.4f}")
    print(f"Log Loss: {logloss:.4f}")
    print("\nClassification Report:\n", classification_rep)
    print("\nConfusion Matrix:\n", conf_matrix)

    metrics = {'accuracy': accuracy, 'mse': mse, 'logloss': logloss}

    return best_nn_model, metrics


from lightgbm import LGBMClassifier

def train_and_evaluate_lightgbm(train_features, train_labels, test_features, test_labels, param_grid=None, cv=5, random_state=42):
    if param_grid is None:
        param_grid = {
            'n_estimators': [50, 100, 200],
            'max_depth': [-1, 10, 20],
            'learning_rate': [0.01, 0.1, 0.2],
            'num_leaves': [31, 50, 100],
        }

    # Grid Search for hyperparameter tuning
    grid_search = GridSearchCV(
        LGBMClassifier(random_state=random_state, verbosity=-1),
        param_grid,
        cv=cv,
        scoring='accuracy'
    )
    grid_search.fit(train_features, train_labels)

    best_model = grid_search.best_estimator_
    print("Best parameters:", grid_search.best_params_)

    # Evaluate the model
    train_predictions = best_model.predict(train_features)
    test_predictions = best_model.predict(test_features)
    test_probabilities = best_model.predict_proba(test_features)
    accuracy = accuracy_score(test_labels, test_predictions)
    report = classification_report(test_labels, test_predictions)
    conf_matrix = confusion_matrix(test_labels, test_predictions)
    train_accuracy = accuracy_score(train_labels, train_predictions)
    mse = mean_squared_error(test_labels, np.argmax(test_probabilities, axis=1))
    logloss = log_loss(test_labels, test_probabilities)

    print(f"Training Accuracy: {train_accuracy:.4f}")
    print(f"Test Accuracy: {accuracy:.4f}")
    print(f"Mean Squared Error: {mse:.4f}")
    print(f"Log Loss: {logloss:.4f}")
    print("\nClassification Report:\n", report)
    print("\nConfusion Matrix:\n", conf_matrix)

    return best_model, {'accuracy': accuracy, 'classification_report': report, 'confusion_matrix': conf_matrix}
