#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 12:30:05 2025

@author: qaraa001
"""

import numpy as np
import pandas as pd
import os
from boruta import BorutaPy
from xgboost import XGBClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import make_scorer, accuracy_score
from sklearn.preprocessing import LabelEncoder
import numpy 
import random
from copy import deepcopy
from scipy.special import gamma
from sklearn.model_selection import train_test_split
from numpy import abs, zeros, log10, where, arctanh, tanh
from numpy.random import uniform, standard_cauchy
import numpy as np
import random
import math
import time
import os



def iPSOeL(objf, dim, lb, ub, Max_iter):
    # Algorithm parameters
    PopSize = 20
    Vmax = 6
    wMax = 0.9
    wMin = 0.2
    c1 = 2
    c2 = 2
    
    # Initialize population
    pos = np.zeros((PopSize, dim))
    vel = np.zeros((PopSize, dim))
    Cost = np.full(PopSize, float("inf"))
    pBestScore = np.full(PopSize, float("inf"))
    pBest = np.zeros((PopSize, dim))
    gBest = np.zeros(dim)
    gBestScore = float("inf")
    
    # Initialize position and velocity
    for i in range(dim):
        pos[:, i] = np.random.uniform(0, 1, PopSize) * (ub[i] - lb[i]) + lb[i]
        vel[:, i] = np.random.uniform(-Vmax, Vmax, PopSize)
    start_time = time.time()
    # Evaluate initial population
    for i in range(PopSize):
        Cost[i] = objf(*pos[i, :])
        pBestScore[i] = Cost[i]
        pBest[i, :] = pos[i, :].copy()
        
        if Cost[i] < gBestScore:
            gBestScore = Cost[i]
            gBest = pos[i, :].copy()
    
    # Initialize additional variables for iPSOeL
    uF = np.zeros(PopSize)
    uCR = np.zeros(PopSize)
    convergence_curve = np.zeros(Max_iter)
    
    # Main optimization loop
    for t in range(Max_iter):
        # Rank population
        SmellOrder = np.sort(Cost)
        SmellIndex = np.argsort(Cost)
        
        # Update uF and uCR based on ranking
        for i in range(PopSize):
            uF[SmellIndex[i]] = i / PopSize
            uCR[SmellIndex[i]] = i / PopSize
        
        # Get best solutions
        Best_X = pos[SmellIndex[0], :]
        Best_X2 = pos[SmellIndex[1], :]
        Best_X3 = pos[SmellIndex[2], :]
        
        # Inertia weight update
        w = wMax - t * ((wMax - wMin) / Max_iter)
        
        for i in range(PopSize):
            # Select two random individuals
            rand = random.sample([x for x in range(PopSize) if x != i], 2)
            r1, r2 = rand[0], rand[1]
            
            # Generate ro1 and ro2
            if random.random() < 0.5:
                ro1 = uF[r1] + 0.1 * random.random()
                ro2 = uF[r2] + 0.1 * random.random()
            else:
                ro1 = uF[i] + 0.1 * random.random()
                ro2 = uF[i] + 0.1 * random.random()
            
            # Update velocity and position
            for j in range(dim):
                cof1 = random.random()
                cof2 = random.random()
                
                # Standard PSO velocity update
                vel[i, j] = (w * vel[i, j] + 
                             c1 * cof1 * (pBest[i, j] - pos[i, j]) + 
                             c2 * cof2 * (gBest[j] - pos[i, j]))
                
                # Velocity clamping
                if vel[i, j] > Vmax:
                    vel[i, j] = Vmax
                if vel[i, j] < -Vmax:
                    vel[i, j] = -Vmax
                
                # Position update
                pos[i, j] = pos[i, j] + vel[i, j]
            
            # iPSOeL specific operations
            Xnew = np.zeros(dim)
            jrand = np.random.randint(0, dim)
            CR = uCR[i] + 0.1 * random.random()
            Xm = (pos[r1, :] + pos[r2, :]) / 2  # Midpoint
            
            for j in range(dim):
                if random.random() < CR or j == jrand:
                    Xnew[j] = Xm[j]
                else:
                    Xnew[j] = pos[i, j]
            
            # Special position update (similar to eMPA)
            r = np.random.rand()
            if r < 0.2:
                ids_except_current = [x for x in range(PopSize) if x != i]
                id_1, id_2 = random.sample(ids_except_current, 2)
                
                Xnew1 = gBest + 0.618 * (pos[id_1, :] - pos[id_2, :])
                Xnew2 = pos[i, :] - 0.618 * (pos[id_1, :] - pos[id_2, :])
                alpha = 0.7
                Xnew = alpha * Xnew1 + (1 - alpha) * Xnew2
            
            # Boundary check
            Xnew = np.clip(Xnew, lb, ub)
            Xnew_Cost = objf(*Xnew)
            
            # Update if better solution found
            if Cost[i] > Xnew_Cost:
                Cost[i] = Xnew_Cost
                pos[i, :] = Xnew
                
                if Cost[i] < pBestScore[i]:
                    pBestScore[i] = Cost[i]
                    pBest[i, :] = Xnew.copy()
                
                if Cost[i] < gBestScore:
                    gBestScore = Cost[i]
                    gBest = Xnew.copy()
        
    #     # Store convergence
    #     convergence_curve[t] = gBestScore
    #     print(f"Current Iteration: {t+1}, Best Cost: {gBestScore}")
    
    # return gBest, gBestScore


        convergence_curve[t] = gBestScore
        end_time = time.time()  # End timing
        run_time = end_time - start_time     
        # print("Current Iteration + :  ", l+1)
        print(f"Current Iteration: {t+1}, Best Fitness: {gBestScore:.4f}")
    return gBest, gBestScore, convergence_curve, run_time 

from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.metrics import accuracy_score
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from xgboost import XGBClassifier
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import (
    accuracy_score, f1_score, recall_score, precision_score,
    roc_auc_score, confusion_matrix, classification_report
)

def train_xgboost_model_kfold(X, y,class_names, k=5):
    """Train XGBoost model for binary classification with K-Fold cross validation"""
    
    print(f"Training XGBoost model for binary classification with {k}-Fold Cross Validation...")
    
    # Convert to numpy arrays and ensure proper shape
    X = np.array(X)
    y = np.array(y).flatten()  # Ensure y is 1D
    
    print(f"Data shapes - X: {X.shape}, y: {y.shape}")
    print(f"Unique labels in y: {np.unique(y)}")
    print(f"Class distribution: {np.bincount(y)}")
    
    def objective_function_kfold(learning_rate, max_depth, colsample_bytree, reg_alpha, reg_lambda):
        params = {
            'objective': 'binary:logistic',  # Changed for binary classification
            'max_depth': int(max_depth),
            'learning_rate': learning_rate,
            'gamma': 0.0,
            'colsample_bytree': colsample_bytree,
            'reg_alpha': reg_alpha,
            'reg_lambda': reg_lambda,
            'random_state': random_seed,
            'eval_metric': 'logloss'  # Changed to logloss for binary
        }
        
        # Stratified K-Fold Cross Validation
        skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=random_seed)
        fold_accuracies = []
        
        for train_idx, val_idx in skf.split(X, y):
            X_train_fold, X_val_fold = X[train_idx], X[val_idx]
            y_train_fold, y_val_fold = y[train_idx], y[val_idx]
            
            dtrain = xgb.DMatrix(X_train_fold, label=y_train_fold)
            dval = xgb.DMatrix(X_val_fold, label=y_val_fold)
            evallist = [(dtrain, 'train'), (dval, 'val')]
            
            model = xgb.train(
                params,
                dtrain,
                num_boost_round=200,
                evals=evallist,
                early_stopping_rounds=20,
                verbose_eval=False
            )
            
            # Make predictions on validation fold
            y_pred_proba = model.predict(dval)
            y_pred = (y_pred_proba > 0.5).astype(int)  # Threshold for binary classification
            fold_accuracy = accuracy_score(y_val_fold, y_pred)
            fold_accuracies.append(fold_accuracy)
        
        # Return negative mean accuracy for minimization
        mean_accuracy = np.mean(fold_accuracies)
        return -mean_accuracy

    param_bounds = {
        'learning_rate': np.arange(0.001, 1, 0.1),
        'max_depth': np.arange(5, 16, 1),
        'colsample_bytree': np.arange(0.1, 1, 0.1),     
        'reg_alpha': np.arange(0.0000001, 100, 0.01), 
        'reg_lambda': np.arange(0.0000001, 100, 0.01)
    }

    # Run optimization (assuming PSOL is your optimization algorithm)
    dim = len(param_bounds)
    lb = [param_bounds[k][0] for k in param_bounds]
    ub = [param_bounds[k][-1] for k in param_bounds]
    best_params_vector, best_score, convergence_curve, run_time = iPSOeL(objective_function_kfold, 
                                                                        dim, lb, ub, Max_iter=20)
    
    optimization_results = {
        'convergence_curve': convergence_curve,
        'run_time': run_time,
        'best_score': best_score
    }

    best_params = {k: best_params_vector[i] for i, k in enumerate(param_bounds)}
    best_params['max_depth'] = int(best_params['max_depth'])
    best_params['gamma'] = 0.0
    
    print("\nBest Hyperparameters Found:")
    print(best_params)
    print(f"Best Cross-Validation Accuracy (negative): {best_score:.4f}")
    print(f"Best Cross-Validation Accuracy: {-best_score:.4f}")
    
    # Train final model on full dataset with best parameters
    final_params = {
        'objective': 'binary:logistic',  # Binary objective
        'max_depth': best_params['max_depth'],
        'learning_rate': best_params['learning_rate'],
        'gamma': 0.0,
        'colsample_bytree': best_params['colsample_bytree'],
        'reg_alpha': best_params['reg_alpha'],
        'reg_lambda': best_params['reg_lambda'],
        'random_state': random_seed,
        'eval_metric': 'logloss'  # Binary evaluation metric
    }
    
    print(f"Training final model on full dataset...")
    print(f"Final training data shapes - X: {X.shape}, y: {y.shape}")
    
    dtrain_full = xgb.DMatrix(X, label=y)
    
    final_model = xgb.train(
        final_params,
        dtrain_full,
        num_boost_round=200,
        verbose_eval=10
    )
    
    return final_model, best_params, -best_score, optimization_results

def evaluate_xgboost_model_kfold(X, y, best_params, k=5, dataset_name=""):
    """Evaluate XGBoost binary classification model with k-fold cross validation"""
    
    skf = StratifiedKFold(n_splits=k, shuffle=True, random_state=random_seed)
    fold_results = []
    
    for fold, (train_idx, val_idx) in enumerate(skf.split(X, y), 1):
        X_train_fold, X_val_fold = X[train_idx], X[val_idx]
        y_train_fold, y_val_fold = y[train_idx], y[val_idx]
        
        # Train model with best parameters for binary classification
        params = {
            'objective': 'binary:logistic',  # Changed to binary objective
            'max_depth': best_params['max_depth'],
            'learning_rate': best_params['learning_rate'],
            'gamma': best_params['gamma'],
            'colsample_bytree': best_params['colsample_bytree'],
            'reg_alpha': best_params['reg_alpha'],
            'reg_lambda': best_params['reg_lambda'],
            'random_state': random_seed,
            'eval_metric': 'logloss'  # Changed to logloss for binary
        }
        
        dtrain = xgb.DMatrix(X_train_fold, label=y_train_fold)
        dval = xgb.DMatrix(X_val_fold, label=y_val_fold)
        
        model = xgb.train(
            params,
            dtrain,
            num_boost_round=200,
            verbose_eval=False
        )
        
        # Make predictions for binary classification
        y_pred_proba = model.predict(dval)
        y_pred = (y_pred_proba > 0.5).astype(int)  # Apply threshold for binary
        
        # Calculate metrics for binary classification
        accuracy = accuracy_score(y_val_fold, y_pred)
        f1 = f1_score(y_val_fold, y_pred, average='binary')  # Changed to binary average
        precision = precision_score(y_val_fold, y_pred, average='binary')  # Changed to binary average
        recall = recall_score(y_val_fold, y_pred, average='binary')  # Changed to binary average
        roc_auc = roc_auc_score(y_val_fold, y_pred_proba)  # Added ROC AUC
        
        # Additional binary-specific metrics
        tn, fp, fn, tp = confusion_matrix(y_val_fold, y_pred).ravel()
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        
        fold_results.append({
            'fold': fold,
            'accuracy': accuracy,
            'f1_score': f1,
            'precision': precision,
            'recall': recall,
            'roc_auc': roc_auc,
            'specificity': specificity,
            'y_true': y_val_fold,
            'y_pred': y_pred,
            'y_pred_proba': y_pred_proba
        })
    
    # Calculate mean metrics
    mean_accuracy = np.mean([r['accuracy'] for r in fold_results])
    std_accuracy = np.std([r['accuracy'] for r in fold_results])
    mean_f1 = np.mean([r['f1_score'] for r in fold_results])
    mean_precision = np.mean([r['precision'] for r in fold_results])
    mean_recall = np.mean([r['recall'] for r in fold_results])
    mean_roc_auc = np.mean([r['roc_auc'] for r in fold_results])
    mean_specificity = np.mean([r['specificity'] for r in fold_results])
    
    # Print comprehensive results
    print(f"\n{'='*50}")
    print(f"Binary Classification Evaluation Results")
    print(f"{'='*50}")
    print(f"Dataset: {dataset_name}")
    print(f"Number of folds: {k}")
    print(f"Mean Accuracy: {mean_accuracy:.4f} ± {std_accuracy:.4f}")
    print(f"Mean F1 Score: {mean_f1:.4f}")
    print(f"Mean Precision: {mean_precision:.4f}")
    print(f"Mean Recall (Sensitivity): {mean_recall:.4f}")
    print(f"Mean Specificity: {mean_specificity:.4f}")
    print(f"Mean ROC AUC: {mean_roc_auc:.4f}")
    
    # Print per-fold results
    print(f"\nPer-fold results:")
    for result in fold_results:
        print(f"Fold {result['fold']}: "
              f"Acc={result['accuracy']:.4f}, "
              f"F1={result['f1_score']:.4f}, "
              f"AUC={result['roc_auc']:.4f}")
    
    return {
        'fold_results': fold_results,
        'mean_accuracy': mean_accuracy,
        'std_accuracy': std_accuracy,
        'mean_f1': mean_f1,
        'mean_precision': mean_precision,
        'mean_recall': mean_recall,
        'mean_roc_auc': mean_roc_auc,
        'mean_specificity': mean_specificity
    }

# Additional function to plot ROC curve for binary classification
def plot_roc_curve_binary(evaluation_results, dataset_name=""):
    """Plot ROC curve for binary classification results"""
    from sklearn.metrics import roc_curve, auc
    
    # Collect all predictions across folds
    all_y_true = np.concatenate([r['y_true'] for r in evaluation_results['fold_results']])
    all_y_pred_proba = np.concatenate([r['y_pred_proba'] for r in evaluation_results['fold_results']])
    
    # Calculate overall ROC curve
    fpr, tpr, _ = roc_curve(all_y_true, all_y_pred_proba)
    roc_auc = auc(fpr, tpr)
    
    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (AUC = {roc_auc:.3f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Random classifier')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'ROC Curve - {dataset_name}\nBinary Classification')
    plt.legend(loc="lower right")
    plt.grid(True, alpha=0.3)
    plt.show()
    
    return roc_auc

# Function to plot confusion matrix for binary classification
def plot_confusion_matrix_binary(evaluation_results, dataset_name=""):
    """Plot confusion matrix for binary classification results"""
    # Collect all predictions across folds
    all_y_true = np.concatenate([r['y_true'] for r in evaluation_results['fold_results']])
    all_y_pred = np.concatenate([r['y_pred'] for r in evaluation_results['fold_results']])
    
    cm = confusion_matrix(all_y_true, all_y_pred)
    
    plt.figure(figsize=(6, 5))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                xticklabels=['Class 0', 'Class 1'],
                yticklabels=['Class 0', 'Class 1'])
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    plt.title(f'Confusion Matrix - {dataset_name}\nBinary Classification')
    plt.show()
    
    return cm


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score
import xgboost as xgb
from sklearn.model_selection import train_test_split, StratifiedKFold
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import (
    accuracy_score, f1_score, recall_score, precision_score,
    roc_auc_score, confusion_matrix, classification_report
)
import random
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
import warnings
warnings.filterwarnings("ignore")

# Set random seed
random_seed = 100
np.random.seed(random_seed)
random.seed(random_seed)






# Load the dataset
try:
    X = pd.read_csv("DEseq_matched_Genes_expression_data_trian_validation.csv", index_col=0)
    X = X.T
    y = pd.read_csv("TCGA-LUAD_LUSC_labels_combined.csv")
    unique_count = y['label'].nunique()
    print(f"Number of unique labels: {unique_count}")
    # Label encoding
    label_encoder = LabelEncoder()
    y = label_encoder.fit_transform(y.iloc[:, 0])
    class_names = label_encoder.classes_
    print(" X shape : ", X.shape)
    print(" y shape : ", y.shape)
    import joblib 
    joblib.dump(label_encoder, 'label_encoder.pkl')

except FileNotFoundError:
    print("X_features.csv or y_labels.csv not found.")
    exit(1)

# ------------------------ Standardization ------------------------ #
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

dataset_name = "NSCLC" 

# Number of repetitions for k-fold
kfold_number = 5
N_RUNS = 20
all_models = []
all_optimization_results = []
all_best_params = []
all_kfold_results = []

print(f"\n{'='*60}")
print(f"RUNNING {N_RUNS} REPETITIONS OF K-FOLD CROSS VALIDATION")
print(f"{'='*60}")

for run in range(1, N_RUNS + 1):
    print(f"\nRun {run}/{N_RUNS}")
    print("-" * 40)
    
    # Train the model for this run with k-fold CV
    model, best_params, best_cv_accuracy, optimization_results = train_xgboost_model_kfold(
        X_scaled, y, class_names, k=kfold_number
    )
    
    # Store results
    all_models.append(model)
    all_best_params.append(best_params)
    all_optimization_results.append(optimization_results)
    
    # Evaluate with k-fold CV
    kfold_results = evaluate_xgboost_model_kfold(
        X_scaled, 
        y, 
        best_params, 
        k=kfold_number, 
        dataset_name=f"{dataset_name}_run_{run}"
    )
    # Plot additional visualizations
    # plot_roc_curve_binary(evaluation_results, "GSE99701 Binary")
    # plot_confusion_matrix_binary(evaluation_results, "GSE99701 Binary")
    all_kfold_results.append(kfold_results)
    
    print(f"Run {run} completed - Best CV Accuracy: {best_cv_accuracy:.4f}")

# Save convergence data for all runs
all_convergence_data = []
all_run_stats = []
all_kfold_stats = []

for run_idx, (results, kfold_results) in enumerate(zip(all_optimization_results, all_kfold_results), 1):
    # Save individual run convergence curve
    convergence_df = pd.DataFrame({
        'Iteration': range(len(results['convergence_curve'])),
        'Fitness_Value': results['convergence_curve'],
        'Run': [run_idx] * len(results['convergence_curve'])
    })
    all_convergence_data.append(convergence_df)
    
    # Save run statistics
    run_stats_df = pd.DataFrame({
        'Run_Number': [run_idx],
        'Run_Time_Seconds': [results['run_time']],
        'Best_Score': [results['best_score']],
        'Best_CV_Accuracy': [-results['best_score']],  # Convert back to accuracy
        'Dataset': [dataset_name]
    })
    all_run_stats.append(run_stats_df)
    
    # Save k-fold results
    kfold_stats_df = pd.DataFrame({
        'Run_Number': [run_idx],
        'Mean_CV_Accuracy': [kfold_results['mean_accuracy']],
        'Std_CV_Accuracy': [kfold_results['std_accuracy']],
        'Mean_F1_Score': [kfold_results['mean_f1']],
        'Mean_Precision': [kfold_results['mean_precision']],
        'Mean_Recall': [kfold_results['mean_recall']],
        'Dataset': [dataset_name]
    })
    all_kfold_stats.append(kfold_stats_df)
    
    # Save iteration-level data for this run
    # iteration_df = pd.DataFrame(results['iteration_data'])
    # iteration_df.to_csv(f'kfold_iteration_data_run_{run_idx}_{dataset_name}.csv', index=False)
    
    # Save k-fold detailed results
    kfold_details_df = pd.DataFrame(kfold_results['fold_results'])
    kfold_details_df.to_csv(f'kfold_detailed_results_run_{run_idx}_{dataset_name}.csv', index=False)

# Combine all convergence data
combined_convergence_df = pd.concat(all_convergence_data, ignore_index=True)
combined_convergence_df.to_csv(f'kfold_convergence_curve_all_runs_{dataset_name}.csv', index=False)

# Combine all run statistics
combined_stats_df = pd.concat(all_run_stats, ignore_index=True)
combined_stats_df.to_csv(f'kfold_optimization_results_all_runs_{dataset_name}.csv', index=False)

# Combine all k-fold statistics
combined_kfold_stats_df = pd.concat(all_kfold_stats, ignore_index=True)
combined_kfold_stats_df.to_csv(f'kfold_cv_results_all_runs_{dataset_name}.csv', index=False)

# Calculate and print summary statistics
print(f"\n{'='*60}")
print("SUMMARY STATISTICS ACROSS ALL K-FOLD RUNS")
print(f"{'='*60}")
print(f"Mean Best CV Accuracy: {combined_stats_df['Best_CV_Accuracy'].mean():.4f} ± {combined_stats_df['Best_CV_Accuracy'].std():.4f}")
print(f"Mean K-Fold CV Accuracy: {combined_kfold_stats_df['Mean_CV_Accuracy'].mean():.4f} ± {combined_kfold_stats_df['Mean_CV_Accuracy'].std():.4f}")
print(f"Mean Run Time: {combined_stats_df['Run_Time_Seconds'].mean():.2f} ± {combined_stats_df['Run_Time_Seconds'].std():.2f} seconds")
print(f"Best CV Accuracy Overall: {combined_stats_df['Best_CV_Accuracy'].max():.4f}")

# Plot convergence curves for all runs
plt.figure(figsize=(12, 8))
for run_idx, results in enumerate(all_optimization_results, 1):
    plt.plot(results['convergence_curve'], 
             label=f'Run {run_idx} (Final: {-results["best_score"]:.4f})', 
             alpha=0.7, linewidth=2)

plt.xlabel('Iteration')
plt.ylabel('Fitness Value (Negative Accuracy)')
plt.title(f'K-Fold Convergence Curves for {N_RUNS} Runs - {dataset_name}')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig(f'kfold_convergence_curves_all_runs_{dataset_name}.png', dpi=300, bbox_inches='tight')
plt.show()

# Plot k-fold CV accuracy across runs
plt.figure(figsize=(10, 6))
plt.bar(range(1, N_RUNS + 1), combined_kfold_stats_df['Mean_CV_Accuracy'], 
        yerr=combined_kfold_stats_df['Std_CV_Accuracy'], 
        capsize=5, alpha=0.7)
plt.axhline(y=combined_kfold_stats_df['Mean_CV_Accuracy'].mean(), 
           color='red', linestyle='--', label=f'Mean: {combined_kfold_stats_df["Mean_CV_Accuracy"].mean():.4f}')
plt.xlabel('Run Number')
plt.ylabel('K-Fold CV Accuracy')
plt.title(f'K-Fold CV Accuracy Across {N_RUNS} Runs - {dataset_name}')
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig(f'kfold_accuracy_across_runs_{dataset_name}.png', dpi=300, bbox_inches='tight')
plt.show()

# Select the best model based on k-fold CV performance
best_run_idx = combined_kfold_stats_df['Mean_CV_Accuracy'].idxmax()
best_model = all_models[best_run_idx]
best_params = all_best_params[best_run_idx]

print(f"\nBest model from Run {best_run_idx + 1}:")
print(f"K-Fold CV Accuracy: {combined_kfold_stats_df.loc[best_run_idx, 'Mean_CV_Accuracy']:.4f}")
print(f"Best Parameters: {best_params}")

# Save final best model (optional)
import joblib
joblib.dump(best_model, f'best_xgboost_model_kfold_{dataset_name}.pkl')
joblib.dump(best_params, f'best_params_kfold_{dataset_name}.pkl')

print(f"\nAll {N_RUNS} k-fold runs completed successfully!")
print(f"Results saved to CSV files with prefix 'kfold_*_{dataset_name}.csv'")
