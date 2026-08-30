"""
Fail-closed model artifact store.

Model artifacts are NOT loaded/served unless explicit env vars are set:
  MODEL_RELEASE_APPROVED=true
  APPROVED_ARTIFACT_REVISION=<non-empty>

This is the developer's standard release-gate pattern.
"""

from __future__ import annotations

import os
import pickle
from pathlib import Path
from typing import Optional


class ModelStore:
    def __init__(self):
        self._model = None
        self._feature_names = None
        self._revision: Optional[str] = None
        self._loaded: bool = False

    def try_load(self) -> bool:
        """Attempt to load model if release gate passes. Returns True if loaded."""
        approved = os.getenv("MODEL_RELEASE_APPROVED", "").lower() == "true"
        revision = os.getenv("APPROVED_ARTIFACT_REVISION", "").strip()
        if not approved or not revision:
            self._loaded = False
            self._model = None
            self._revision = None
            return False
        # Gate passed: try to load artifact file
        artifact_path = os.getenv("MODEL_ARTIFACT_PATH", "artifacts/model.pkl")
        p = Path(artifact_path)
        if not p.exists():
            # Gate passed but no artifact file — still not loaded (honest)
            self._loaded = False
            self._revision = revision
            return False
        try:
            with open(p, "rb") as f:
                data = pickle.load(f)
            self._model = data.get("model") if isinstance(data, dict) else data
            self._feature_names = data.get("feature_names") if isinstance(data, dict) else None
            self._revision = revision
            self._loaded = True
            return True
        except Exception:
            self._loaded = False
            return False

    @property
    def is_loaded(self) -> bool:
        return self._loaded

    @property
    def revision(self) -> Optional[str]:
        return self._revision

    @property
    def model(self):
        return self._model

    @property
    def feature_names(self):
        return self._feature_names

    def reset(self):
        self._model = None
        self._feature_names = None
        self._revision = None
        self._loaded = False


# Global singleton
store = ModelStore()
