"""
Fusion model: logistic regression combining structural + conservation + physicochemical features.

Features per variant:
  - rsa (0..1)
  - conservation (0..1)
  - delta_hydrophobicity
  - delta_volume
  - grantham_distance
  - blosum62_score
  - secondary_structure one-hot (H, E, C)
  - abs(delta_hydrophobicity), abs(delta_volume) for magnitude

Model: sklearn LogisticRegression (L2) or GradientBoosting.
Simple, defensible; not an overcomplicated deep model.
"""

from __future__ import annotations

import math
import pickle
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

from .physicochemical import delta_hydrophobicity, delta_volume, grantham_distance, blosum62_score
from .structure import ResidueFeatures

FEATURE_NAMES = [
    "rsa",
    "conservation",
    "delta_hydro",
    "delta_vol",
    "grantham",
    "blosum62",
    "abs_delta_hydro",
    "abs_delta_vol",
    "buried",  # 1 if rsa<0.2 else 0
    "ss_H",
    "ss_E",
    "ss_C",
]

SS_ORDER = ["H", "E", "C"]


def featurize_variant(
    wt: str,
    mut: str,
    residue_feature: Optional[ResidueFeatures],
    conservation: float,
) -> np.ndarray:
    """Build feature vector for a single variant."""
    wt = wt.upper()
    mut = mut.upper()
    rsa = residue_feature.rsa if residue_feature else 0.5
    ss = residue_feature.ss if residue_feature else "C"
    d_hydro = delta_hydrophobicity(wt, mut)
    d_vol = delta_volume(wt, mut)
    grantham = grantham_distance(wt, mut)
    blosum = float(blosum62_score(wt, mut))
    buried = 1.0 if rsa < 0.2 else 0.0
    ss_vec = [1.0 if ss == s else 0.0 for s in SS_ORDER]
    feats = np.array(
        [rsa, conservation, d_hydro, d_vol, grantham, blosum, abs(d_hydro), abs(d_vol), buried] + ss_vec,
        dtype=float,
    )
    return feats


def featurize_batch(variants: List[Dict[str, Any]]) -> Tuple[np.ndarray, List[str]]:
    """
    Variants: list of dicts with keys wt, mut, residue_feature, conservation.
    Returns (X, names)
    """
    X = []
    for v in variants:
        feats = featurize_variant(v["wt"], v["mut"], v.get("residue_feature"), v.get("conservation", 0.5))
        X.append(feats)
    return np.array(X), FEATURE_NAMES


def train_logistic_regression(
    X: np.ndarray, y: np.ndarray, C: float = 1.0, max_iter: int = 1000
) -> Pipeline:
    pipe = Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(C=C, max_iter=max_iter, solver="lbfgs")),
    ])
    pipe.fit(X, y)
    return pipe


def train_gradient_boosting(X: np.ndarray, y: np.ndarray) -> GradientBoostingClassifier:
    clf = GradientBoostingClassifier(random_state=42)
    clf.fit(X, y)
    return clf


def evaluate(y_true: np.ndarray, y_score: np.ndarray) -> Dict[str, float]:
    """AUROC / AUPRC."""
    auroc = roc_auc_score(y_true, y_score) if len(set(y_true)) > 1 else 0.5
    auprc = average_precision_score(y_true, y_score) if len(set(y_true)) > 1 else float(np.mean(y_true))
    return {"auroc": float(auroc), "auprc": float(auprc)}


def evaluate_with_baselines(
    X: np.ndarray, y: np.ndarray, feature_names: List[str], model
) -> Dict[str, Dict[str, float]]:
    """Compare fusion model vs single-feature baselines (conservation-only, hydro-only)."""
    # Fusion
    y_score = model.predict_proba(X)[:, 1] if hasattr(model, "predict_proba") else model.decision_function(X)
    results = {"fusion": evaluate(y, y_score)}

    # Baselines: train logistic on single feature
    for feat_name in ["conservation", "abs_delta_hydro"]:
        if feat_name in feature_names:
            idx = feature_names.index(feat_name)
            X_single = X[:, idx].reshape(-1, 1)
            pipe = Pipeline([("scaler", StandardScaler()), ("clf", LogisticRegression(max_iter=1000))])
            pipe.fit(X_single, y)
            y_s = pipe.predict_proba(X_single)[:, 1]
            results[f"baseline_{feat_name}"] = evaluate(y, y_s)
    return results


@dataclass
class ModelArtifacts:
    model: Any
    feature_names: List[str]
    scaler: Optional[Any] = None


def save_artifacts(model, feature_names: List[str], path: str | Path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump({"model": model, "feature_names": feature_names}, f)


def load_artifacts(path: str | Path):
    with open(path, "rb") as f:
        data = pickle.load(f)
    return data["model"], data["feature_names"]
