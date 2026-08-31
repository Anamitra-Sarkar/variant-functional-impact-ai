from __future__ import annotations

import re
from typing import Optional

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

_VALID_AA = set("ARNDCQEGHILKMFPSTWYV")
_CHAIN_RE = re.compile(r"^[A-Za-z0-9]$")


class PredictRequest(BaseModel):
    protein: str = Field(default="UNKNOWN", description="Protein / UniProt ID / gene")
    position: int = Field(..., ge=1, le=100000, description="1-indexed residue position")
    wt: str = Field(..., min_length=1, max_length=1, description="Wild-type AA single letter")
    mut: str = Field(..., min_length=1, max_length=1, description="Mutant AA single letter")
    chain: str = Field(default="A", description="Chain ID single character A-Z/0-9")

    @field_validator("protein")
    @classmethod
    def validate_protein(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("protein must be non-empty")
        if len(v) > 64:
            raise ValueError("protein identifier too long (max 64)")
        # Allow alphanumeric, dash, underscore only for safety (UniProt/gene)
        if not re.match(r"^[A-Za-z0-9_\-]+$", v):
            raise ValueError("protein must contain only letters, numbers, underscore or dash")
        return v

    @field_validator("wt", "mut")
    @classmethod
    def validate_aa(cls, v: str) -> str:
        v = v.strip().upper()
        if v not in _VALID_AA:
            raise ValueError(f"Invalid amino acid '{v}': must be one of {''.join(sorted(_VALID_AA))}")
        return v

    @field_validator("chain")
    @classmethod
    def validate_chain(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("chain must be non-empty")
        if len(v) != 1:
            raise ValueError("chain must be a single character")
        if not _CHAIN_RE.match(v):
            raise ValueError("chain must be alphanumeric A-Z or 0-9")
        return v.upper()

    @model_validator(mode="after")
    def validate_distinct(self):
        if self.wt and self.mut and self.wt == self.mut:
            raise ValueError(f"wt and mut must be distinct (both '{self.wt}') — synonymous variant has no impact to score")
        return self


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


class ModelInfoResponse(BaseModel):
    model_config = ConfigDict(protected_namespaces=())
    model_loaded: bool
    model_revision: Optional[str] = None
    feature_names: Optional[list[str]] = None
    artifact_path: Optional[str] = None
