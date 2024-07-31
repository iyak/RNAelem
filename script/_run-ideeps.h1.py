from sklearn.metrics import roc_auc_score
import gzip
import sys


seq  = sys.argv[1]
pred = sys.argv[2]


# Initialize lists to hold labels and predictions
y_true = []
y_pred = []


# Read gzipped FASTA file and populate y_true based on sequence headers
with gzip.open(seq, "rt") as fasta_file:
    for line in fasta_file:
        line = line.strip()
        if line.startswith(">"):
            if "mark:neg" in line:
                y_true.append(0)
            else:
                y_true.append(1)


# Read prediction scores and populate y_pred
with open(pred, "r") as pred_file:
    for line in pred_file:
        y_pred.append(float(line.strip()))


# Calculate and print AUC
if len(y_true) != len(y_pred):
    print("Number of labels and predictions do not match.")
else:
    auc = roc_auc_score(y_true, y_pred)
    print(auc)
