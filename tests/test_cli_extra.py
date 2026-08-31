"""CLI parsing edge cases."""
import pytest
from data_pipeline.cli import parse_variant_string

def test_simple():
    assert parse_variant_string("R175H") == ("UNKNOWN",175,"R","H")
    assert parse_variant_string("r175h") == ("UNKNOWN",175,"R","H")

def test_gene_prefix():
    assert parse_variant_string("P53_R175H") == ("P53",175,"R","H")
    assert parse_variant_string("BRCA1_A100V") == ("BRCA1",100,"A","V")

def test_four_part():
    # TEST_2_A_V -> protein TEST, pos2, wt A, mut V
    assert parse_variant_string("TEST_2_A_V") == ("TEST",2,"A","V")

def test_with_dash_gene():
    # gene with dash? parse keeps protein as UNKNOWN? Actually "GENE-1_R100H" -> split "_" gives ["GENE-1","R100H"] -> protein GENE-1
    assert parse_variant_string("GENE-1_R100H")[0] == "GENE-1"

def test_invalid_raises():
    with pytest.raises(ValueError):
        parse_variant_string("bad")
    with pytest.raises(ValueError):
        parse_variant_string("P53_R175")
    with pytest.raises(ValueError):
        parse_variant_string("P53_175H")
    with pytest.raises(ValueError):
        parse_variant_string("")

def test_extra_whitespace():
    assert parse_variant_string("  P53_R175H  ") == ("P53",175,"R","H")
