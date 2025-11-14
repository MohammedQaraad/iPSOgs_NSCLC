import pandas as pd
import numpy as np
import joblib
from sklearn.preprocessing import LabelEncoder, MinMaxScaler
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score, f1_score, recall_score, precision_score, roc_auc_score, roc_curve, auc, matthews_corrcoef
import matplotlib.pyplot as plt
import seaborn as sns

# Load the saved model and label encoder
optimized_xgb = joblib.load('optimized_xgboost_model.pkl')
label_encoder = joblib.load('label_encoder.pkl')

# Load test data
X_test = pd.read_csv("GSE81089_matched_Genes_expression_data_test.csv", index_col=0) # with batch effect removal
X_test = X_test.T
y_test = pd.read_csv("GSE81089_LUAD_LUSC_Labels_filtered.csv")

# Preprocess test data
X_test_col = X_test.columns

# Check for missing values
if X_test.isna().sum().sum() > 0:
    X_test = X_test.fillna(X_test.mean())  # Simple imputation
    print("Filled missing values with column means")

# Encode test labels
y_test_encoded = label_encoder.transform(y_test['Subtype'])

# Verify shapes
print("X_test shape:", X_test.shape)
print("y_test shape:", y_test_encoded.shape)

# Check if binary classification
if len(label_encoder.classes_) != 2:
    raise ValueError("This code is designed for binary classification. Found {} classes.".format(len(label_encoder.classes_)))

scaler = MinMaxScaler()
X_test_scaled = scaler.fit_transform(X_test)

# Make predictions
y_pred = optimized_xgb.predict(X_test_scaled)
y_pred_proba = optimized_xgb.predict_proba(X_test_scaled)[:, 1]  # Probability of positive class

# Evaluate predictions
print("\n" + "="*50)
print("TEST SET EVALUATION RESULTS")
print("="*50)

# Calculate and print metrics
accuracy = accuracy_score(y_test_encoded, y_pred)
precision = precision_score(y_test_encoded, y_pred)
recall = recall_score(y_test_encoded, y_pred)
f1 = f1_score(y_test_encoded, y_pred)
mcc = matthews_corrcoef(y_test_encoded, y_pred)
roc_auc = roc_auc_score(y_test_encoded, y_pred_proba)

print("\nTest Accuracy: {:.4f}".format(accuracy))
print("Precision: {:.4f}".format(precision))
print("Recall: {:.4f}".format(recall))
print("F1 Score: {:.4f}".format(f1))
print("MCC: {:.4f}".format(mcc))
print("ROC AUC: {:.4f}".format(roc_auc))

print("\nClassification Report (Test Set):")
print(classification_report(y_test_encoded, y_pred, target_names=label_encoder.classes_))

print("\nConfusion Matrix (Test Set):")
cm = confusion_matrix(y_test_encoded, y_pred)
print(cm)

# Plotting ROC curve (binary classification)
plt.figure(figsize=(6, 5))
fpr, tpr, _ = roc_curve(y_test_encoded, y_pred_proba)
roc_auc = auc(fpr, tpr)

plt.plot(fpr, tpr, color='blue', lw=2,
         label=f'ROC curve (AUC = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], 'k--', lw=2, label='Random guessing')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate', fontsize=10)
plt.ylabel('True Positive Rate', fontsize=10)
plt.title('Receiver Operating Characteristic (ROC) Curve', fontsize=12)
plt.legend(loc="lower right", fontsize=8)
plt.grid(True)
plt.tight_layout()
plt.savefig('GSE81089_test_ROC_curve.png', dpi=300)
plt.show()

# Plotting Normalized Confusion Matrix
plt.figure(figsize=(6, 5))
cm_normalized = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
sns.heatmap(cm_normalized, annot=True, fmt=".2f", cmap="Blues",
            xticklabels=label_encoder.classes_, yticklabels=label_encoder.classes_,
            annot_kws={"size": 10})
plt.xlabel('Predicted Label', fontsize=10)
plt.ylabel('True Label', fontsize=10)
plt.title('Normalized Confusion Matrix', fontsize=12)
plt.tight_layout()
plt.savefig('GSE81089_test_Normalized_Confusion_Matrix.png', dpi=300)
plt.show()

# Plotting Feature Importance (Top 10)
if hasattr(optimized_xgb, 'feature_importances_'):
    importances = optimized_xgb.feature_importances_
    feature_names = X_test_col

    feature_importances = pd.Series(importances, index=feature_names)

    top_n = 10
    top_features = feature_importances.nlargest(top_n)

    plt.figure(figsize=(8, 6))
    sns.barplot(x=top_features.values, y=top_features.index, palette='viridis')
    plt.xlabel('Feature Importance', fontsize=10)
    plt.ylabel('Feature Name', fontsize=10)
    plt.title(f'Top {top_n} Feature Importances', fontsize=12)
    plt.tight_layout()
    plt.savefig(f'Top_{top_n}_Feature_Importances.png', dpi=300)
    plt.show()
else:
    print("\nFeature importances not available for this model type.")