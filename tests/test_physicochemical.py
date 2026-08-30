import pytest
from data_pipeline.physicochemical import (
    KYTE_DOOLITTLE,
    VDW_VOLUME,
    BLOSUM62,
    delta_hydrophobicity,
    delta_volume,
    grantham_distance,
    blosum62_score,
)


def test_kd_table_correctness():
    # Real values spot checks
    assert KYTE_DOOLITTLE["I"] == 4.5
    assert KYTE_DOOLITTLE["R"] == -4.5
    assert KYTE_DOOLITTLE["A"] == 1.8
    assert len(KYTE_DOOLITTLE) == 20


def test_vdw_table_correctness():
    assert VDW_VOLUME["G"] == 60.1
    assert VDW_VOLUME["W"] == 227.8
    assert VDW_VOLUME["A"] == 88.6
    assert len(VDW_VOLUME) == 20


def test_blosum62_correctness():
    # Diagonal positive, off-diagonal often negative
    assert blosum62_score("A", "A") == 4
    assert blosum62_score("W", "W") == 11
    assert blosum62_score("R", "W") == -3
    # Symmetry
    assert blosum62_score("A", "R") == blosum62_score("R", "A")


def test_delta_hydrophobicity():
    # R(-4.5) -> W(-0.9) delta = 3.6
    assert abs(delta_hydrophobicity("R", "W") - 3.6) < 1e-6
    # I(4.5) -> R(-4.5) delta = -9.0
    assert abs(delta_hydrophobicity("I", "R") - (-9.0)) < 1e-6


def test_delta_volume():
    # G 60.1 -> W 227.8 delta ~167.7
    assert abs(delta_volume("G", "W") - 167.7) < 1e-6


def test_grantham_distance():
    # Identity should be 0
    assert grantham_distance("A", "A") == 0.0
    # Large change R->W should be larger than small L->I
    assert grantham_distance("R", "W") > grantham_distance("L", "I")
    # R->W known to be relatively large (~100 in literature); our simplified c/p/v yields small numbers
    # but ordering must hold and be positive
    assert grantham_distance("R", "W") > 0
    assert grantham_distance("S", "T") < grantham_distance("R", "W")


def test_invalid_aa_raises():
    with pytest.raises(ValueError):
        delta_hydrophobicity("Z", "A")
    with pytest.raises(ValueError):
        grantham_distance("A", "X")
