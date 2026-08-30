from __future__ import annotations
from typing import Optional
from pydantic import BaseModel, ConfigDict, Field


class PredictRequest(BaseModel):
    protein: str = Field(default="UNKNOWN", description="Protein / UniProt ID / gene")
    position: int = Field(..., ge=1, description="1-indexed residue position")
    wt: str = Field(..., min_length=1, max_length=1, description="Wild-type AA single letter")
    mut: str = Field(..., min_length=1, max_length=1, description="Mutant AA single letter")
    chain: str = Field(default="A", description="Chain ID")


class PredictResponse(BaseModel):
    model_config = ConfigDict(protected_namespaces=())
    variant: str
    score: float = Field(..., ge=0, le=1)
    label: str  # "damaging" or "benign"
    features: dict
    model_revision: str


class HealthResponse(BaseModel):
    model_config = ConfigDict(protected_namespaces=())
    status: str
    model_loaded: bool
    model_revision: Optional[str] = None


class ReadinessResponse(BaseModel):
    model_config = ConfigDict(protected_namespaces=())
    ready: bool
    model_loaded: bool
    model_revision: Optional[str] = None
