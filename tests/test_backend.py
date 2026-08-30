import os
import tempfile
import pickle
from pathlib import Path

import numpy as np
from fastapi.testclient import TestClient

from backend.app import create_app
from backend.model_store import ModelStore
from data_pipeline.model import train_logistic_regression, featurize_batch, FEATURE_NAMES
from data_pipeline.structure import ResidueFeatures


def _make_dummy_model(tmp_path):
    # Train a tiny model to use as artifact
    variants = []
    labels = []
    for i in range(30):
        rsa = 0.1 if i % 2 == 0 else 0.9
        cons = 0.9 if i % 2 == 0 else 0.1
        wt, mut = ("R", "W") if i % 2 == 0 else ("S", "T")
        label = 1 if i % 2 == 0 else 0
        rf = ResidueFeatures(chain_id="A", resseq=i+1, resname="ALA", rsa=rsa, ss="H", burial_count=10, ca_coord=(0,0,0))
        variants.append({"wt": wt, "mut": mut, "residue_feature": rf, "conservation": cons})
        labels.append(label)
    X, _ = featurize_batch(variants)
    y = np.array(labels)
    model = train_logistic_regression(X, y)
    p = tmp_path / "model.pkl"
    with open(p, "wb") as f:
        pickle.dump({"model": model, "feature_names": FEATURE_NAMES}, f)
    return p


def test_health_not_loaded_by_default():
    os.environ.pop("MODEL_RELEASE_APPROVED", None)
    os.environ.pop("APPROVED_ARTIFACT_REVISION", None)
    from backend.model_store import store
    store.reset()
    store.try_load()
    app = create_app()
    client = TestClient(app)
    resp = client.get("/health")
    assert resp.status_code == 200
    data = resp.json()
    assert data["model_loaded"] is False


def test_readiness_503_when_not_loaded():
    os.environ.pop("MODEL_RELEASE_APPROVED", None)
    os.environ.pop("APPROVED_ARTIFACT_REVISION", None)
    from backend.model_store import store
    store.reset()
    store.try_load()
    app = create_app()
    client = TestClient(app)
    resp = client.get("/readiness")
    assert resp.status_code == 503
    assert resp.json()["model_loaded"] is False


def test_predict_503_when_gate_closed():
    os.environ.pop("MODEL_RELEASE_APPROVED", None)
    os.environ.pop("APPROVED_ARTIFACT_REVISION", None)
    from backend.model_store import store
    store.reset()
    store.try_load()
    app = create_app()
    client = TestClient(app)
    resp = client.post("/predict", json={"protein": "P53", "position": 175, "wt": "R", "mut": "H"})
    assert resp.status_code == 503
    assert "not yet released" in resp.json()["detail"]


def test_release_gate_loads_when_approved(tmp_path):
    p = _make_dummy_model(tmp_path)
    os.environ["MODEL_RELEASE_APPROVED"] = "true"
    os.environ["APPROVED_ARTIFACT_REVISION"] = "rev-123"
    os.environ["MODEL_ARTIFACT_PATH"] = str(p)
    from backend.model_store import store
    store.reset()
    ok = store.try_load()
    assert ok is True
    assert store.is_loaded is True
    assert store.revision == "rev-123"

    app = create_app()
    client = TestClient(app)
    # health should reflect loaded
    resp = client.get("/health")
    assert resp.json()["model_loaded"] is True
    # readiness 200
    resp = client.get("/readiness")
    assert resp.status_code == 200
    # predict should work
    resp = client.post("/predict", json={"protein": "P53", "position": 175, "wt": "R", "mut": "H", "chain": "A"})
    assert resp.status_code == 200
    data = resp.json()
    assert 0 <= data["score"] <= 1
    assert data["label"] in ("damaging", "benign")
    assert "features" in data
    # cleanup
    os.environ.pop("MODEL_RELEASE_APPROVED", None)
    os.environ.pop("APPROVED_ARTIFACT_REVISION", None)
    os.environ.pop("MODEL_ARTIFACT_PATH", None)
    store.reset()
    store.try_load()


def test_auth_stub_permissive_when_no_service_account():
    os.environ.pop("FIREBASE_SERVICE_ACCOUNT_JSON", None)
    from backend.auth import verify_bearer_token
    # Should not raise, returns None
    result = verify_bearer_token(authorization=None)
    assert result is None


def test_auth_stub_strict_with_valid_token(tmp_path):
    sa_path = tmp_path / "sa.json"
    sa_path.write_text('{"project_id": "test"}')
    os.environ["FIREBASE_SERVICE_ACCOUNT_JSON"] = str(sa_path)
    from backend.auth import verify_bearer_token
    # Valid test token
    result = verify_bearer_token(authorization="Bearer valid-test-token")
    assert result is not None
    assert result["uid"] == "test-user"
    # Invalid token should raise 401
    import pytest
    from fastapi import HTTPException
    with pytest.raises(HTTPException) as exc:
        verify_bearer_token(authorization="Bearer bad-token")
    assert exc.value.status_code == 401
    # Missing header should raise
    with pytest.raises(HTTPException):
        verify_bearer_token(authorization=None)
    os.environ.pop("FIREBASE_SERVICE_ACCOUNT_JSON", None)
