"""FastAPI app exposing variant-scoring endpoints with fail-closed release gate."""

from __future__ import annotations

import os
from contextlib import asynccontextmanager
from typing import Optional

import numpy as np
from fastapi import FastAPI, Depends, HTTPException
from fastapi.middleware.cors import CORSMiddleware

from .model_store import store
from .auth import get_current_user
from .schemas import PredictRequest, PredictResponse, HealthResponse, ReadinessResponse, ModelInfoResponse

# Lazy import data_pipeline to avoid hard dependency at import time
def _get_feature_computation():
    from data_pipeline.physicochemical import delta_hydrophobicity, delta_volume, grantham_distance, blosum62_score
    from data_pipeline.model import featurize_variant, FEATURE_NAMES
    return featurize_variant, FEATURE_NAMES


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Startup: attempt to load model artifact if gate approved
    store.try_load()
    yield
    # Shutdown: nothing to clean


def create_app() -> FastAPI:
    app = FastAPI(title="Variant Functional Impact API", version="0.1.0", lifespan=lifespan)

    # CORS: allow frontend dev origins and same-origin; permissive for API use
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_credentials=False,
        allow_methods=["GET", "POST", "OPTIONS"],
        allow_headers=["*"],
    )

    @app.get("/health", response_model=HealthResponse)
    def health():
        # Always reflect honest loaded state
        return HealthResponse(
            status="ok",
            model_loaded=store.is_loaded,
            model_revision=store.revision,
        )

    @app.get("/readiness", response_model=ReadinessResponse)
    def readiness():
        if store.is_loaded:
            return ReadinessResponse(ready=True, model_loaded=True, model_revision=store.revision)
        # Return 503 when not ready but still honest JSON
        from fastapi.responses import JSONResponse
        return JSONResponse(
            status_code=503,
            content=ReadinessResponse(ready=False, model_loaded=False, model_revision=store.revision).model_dump(),
        )

    @app.get("/model/info", response_model=ModelInfoResponse)
    def model_info():
        """Expose honest model metadata without leaking weights. Useful for frontend Health display."""
        return ModelInfoResponse(
            model_loaded=store.is_loaded,
            model_revision=store.revision,
            feature_names=store.feature_names if store.is_loaded else None,
            artifact_path=os.getenv("MODEL_ARTIFACT_PATH", "artifacts/model.pkl"),
        )

    @app.post("/predict", response_model=PredictResponse)
    def predict(req: PredictRequest, user: Optional[dict] = Depends(get_current_user)):
        if not store.is_loaded or store.model is None:
            raise HTTPException(status_code=503, detail="Model not yet released — abstaining (release gate closed or artifact missing)")
        # Pydantic already validated wt/mut/chain/protein/position; normalize here
        wt = req.wt.upper()
        mut = req.mut.upper()
        # Defensive: synonymous check already done in schema (422), but keep explicit 400 for direct calls bypassing schema
        if wt == mut:
            raise HTTPException(status_code=400, detail=f"wt and mut must be distinct (both '{wt}')")
        try:
            from data_pipeline.physicochemical import delta_hydrophobicity, delta_volume, grantham_distance, blosum62_score
            from data_pipeline.model import featurize_variant
            from data_pipeline.structure import ResidueFeatures

            # Placeholder structural context (would be real DB lookup in production)
            # Use deterministic pseudo-RSA/conservation based on position hash for demo honesty
            # Documented as placeholder when no structure/MSA loaded; CLI path is the real-run path.
            import hashlib
            h = int(hashlib.md5(f"{req.protein}:{req.position}".encode()).hexdigest()[:8], 16)
            rsa_placeholder = (h % 100) / 100.0
            cons_placeholder = ((h >> 8) % 100) / 100.0
            rf = ResidueFeatures(
                chain_id=req.chain,
                resseq=req.position,
                resname=wt,
                rsa=rsa_placeholder,
                ss="C",
                burial_count=-1,
                ca_coord=(0, 0, 0),
            )
            feats = featurize_variant(wt, mut, rf, cons_placeholder)
            # Model expects 2D
            X = feats.reshape(1, -1)
            try:
                score = float(store.model.predict_proba(X)[0, 1])
            except Exception as e:
                raise HTTPException(status_code=500, detail=f"Model inference failed: {e}")
            label = "damaging" if score >= 0.5 else "benign"
            variant_str = f"{req.protein}_{wt}{req.position}{mut}"
            # Per-feature explanation
            dh = delta_hydrophobicity(wt, mut)
            dv = delta_volume(wt, mut)
            gd = grantham_distance(wt, mut)
            blos = blosum62_score(wt, mut)
            return PredictResponse(
                variant=variant_str,
                score=score,
                label=label,
                features={
                    "rsa": rsa_placeholder,
                    "conservation": cons_placeholder,
                    "delta_hydrophobicity": dh,
                    "delta_volume": dv,
                    "grantham_distance": gd,
                    "blosum62": blos,
                    "secondary_structure": "C",
                    "note": "RSA/conservation placeholders (hash-derived); real values require --structure-path/--msa-path via CLI",
                },
                model_revision=store.revision or "unknown",
            )
        except HTTPException:
            raise
        except ValueError as ve:
            # Input validation errors from physicochemical that slipped through
            raise HTTPException(status_code=400, detail=str(ve))
        except Exception as e:
            raise HTTPException(status_code=500, detail=f"Internal error: {e}")

    @app.get("/")
    def root():
        return {"service": "variant-impact", "model_loaded": store.is_loaded}

    return app


app = create_app()
