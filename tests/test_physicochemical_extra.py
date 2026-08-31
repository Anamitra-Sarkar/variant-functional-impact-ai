"""Extra physicochemical tests: lower-case handling, BLOSUM symmetry broader, grantham robustness."""
import pytest
from data_pipeline.physicochemical import delta_hydrophobicity, delta_volume, grantham_distance, blosum62_score, KYTE_DOOLITTLE, VDW_VOLUME, BLOSUM62

def test_lower_case_inputs():
    assert delta_hydrophobicity("r", "w") == delta_hydrophobicity("R", "W")
    assert delta_volume("g", "w") == delta_volume("G", "W")
    assert grantham_distance("r", "w") == grantham_distance("R", "W")
    assert blosum62_score("a", "r") == blosum62_score("A", "R")

def test_blosum_symmetry_all():
    # Sample symmetry across a few pairs
    pairs = [("A","R"), ("W","Y"), ("L","I"), ("D","E"), ("K","R")]
    for a,b in pairs:
        assert blosum62_score(a,b) == blosum62_score(b,a)

def test_blosum_diagonal_range():
    for aa in "ARNDCQEGHILKMFPSTWYV":
        score = blosum62_score(aa, aa)
        assert 4 <= score <= 11, f"diagonal {aa} {score} out of 4..11"

def test_grantham_zero_identity_all():
    for aa in "ARNDCQEGHILKMFPSTWYV":
        assert grantham_distance(aa, aa) == 0.0

def test_grantham_ordering():
    # R->W large, S->T small, L->I very small
    assert grantham_distance("R","W") > grantham_distance("S","T")
    assert grantham_distance("S","T") > grantham_distance("L","I") or grantham_distance("L","I") < grantham_distance("S","T")
    # Hydrophobic -> charged should be large
    assert grantham_distance("I","R") > grantham_distance("I","L")

def test_delta_range():
    # Hydrophobicity delta range should be within -9 to 9
    for wt in KYTE_DOOLITTLE:
        for mut in KYTE_DOOLITTLE:
            d = delta_hydrophobicity(wt, mut)
            assert -9.0 <= d <= 9.0
    # Volume delta range within -167..167
    for wt in VDW_VOLUME:
        for mut in VDW_VOLUME:
            d = delta_volume(wt, mut)
            assert -170 <= d <= 170

def test_invalid_raises_case_insensitive():
    with pytest.raises(ValueError):
        delta_hydrophobicity("z", "a")
    with pytest.raises(ValueError):
        blosum62_score("X", "A")
