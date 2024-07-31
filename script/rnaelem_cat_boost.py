from catboost import CatBoostRegressor
from lightgbm import LGBMRegressor
from pathlib import Path
from sklearn.ensemble import RandomForestRegressor
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
import numpy as np
import sys


def train_and_eval(X_train, y_train, X_test, y_test, save_dir):
    imputer = SimpleImputer(strategy="constant")
    X_train = imputer.fit_transform(X_train)
    X_test = imputer.transform(X_test)
    model = CatBoostRegressor(iterations=1000, learning_rate=0.05, depth=10, verbose=False)

    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    auc = roc_auc_score(y_test, y_pred)
    with open(save_dir / "auc.txt", "w") as fo:
        fo.write(str(auc) + "\n")

    model.save_model(save_dir / "model.cbm", format="cbm")

    # Save feature importance
    feature_importance = model.get_feature_importance(type="FeatureImportance")
    import pandas as pd
    feature_names = [str(i+1) for i in range(X_train.shape[1])]
    importance_df = pd.DataFrame({"Feature": feature_names, "Importance": feature_importance})
    importance_df.sort_values(by="Importance", ascending=False, inplace=True)
    importance_df.to_csv(save_dir / "features.csv", sep="\t", index=False)


def get_train_data(root_dir: Path, sub_dir="elem_out_w50"):
    n = 200
    x = dict()
    y = dict()
    for k in ["0", "1"]:
        for i in range(n):
            for ox in [0, 1]:
                fname = ["negative.raw", "positive.raw"]
                raw = root_dir / sub_dir / f"cv-{k}" / "test" / f"pattern-{i+1}" / fname[ox]
                try:
                    with open(raw) as f:
                        lines = f.readlines()
                        for s, p in zip(lines[0::10], lines[6::10]):
                            s += str(ox)
                            if s not in x:
                                x[s] = [np.nan for _ in range(n)]
                            x[s][i] = float(p.strip().split(": ")[1])
                            if s not in y:
                                y[s] = ox
                except:
                    pass
    x_data = []
    y_data = []
    for key in set(x.keys()) | set(y.keys()):
        x_data.append(x[key])
        y_data.append(y[key])
    return np.array(x_data), np.array(y_data)


def get_test_data(root_dir: Path, sub_dir="scan_out_w50", dname=("negative", "positive")):
    n = 200
    x = dict()
    y = dict()
    for i in range(n):
        for ox in [0, 1]:
            raw = root_dir / sub_dir / dname[ox] / f"pattern-{i+1}" / "scan.raw"
            try:
                with open(raw) as f:
                    lines = f.readlines()
                    for s, p in zip(lines[0::10], lines[6::10]):
                        s += str(ox)
                        if s not in x:
                            x[s] = [np.nan for _ in range(n)]
                        x[s][i] = float(p.strip().split(": ")[1])
                        if s not in y:
                            y[s] = ox
            except:
                pass
    x_data = []
    y_data = []
    for key in set(x.keys()) | set(y.keys()):
        x_data.append(x[key])
        y_data.append(y[key])
    return np.array(x_data), np.array(y_data)


if "__main__" == __name__:
    data_dir = Path(sys.argv[1])
    save_dir = Path(sys.argv[2])
    save_dir.mkdir(parents=True, exist_ok=True)

    X_train, y_train = get_train_data(data_dir)
    X_test, y_test = get_test_data(data_dir)

    train_and_eval(X_train, y_train, X_test, y_test, save_dir)
