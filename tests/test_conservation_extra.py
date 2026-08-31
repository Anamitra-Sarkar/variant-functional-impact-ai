"""Hardened conservation parser tests: real-file quirks (comments, header variants, whitespace, ambiguous AA, gaps)."""
import pytest
from pathlib import Path
from data_pipeline.conservation import parse_fasta, conservation_scores, column_entropy, get_conservation_for_position, LOG2_20

def write(path, content):
    Path(path).write_text(content)

def test_comment_lines_ignored(tmp_path):
    p = tmp_path / "msa.fasta"
    p.write_text("; this is a comment\n# another comment\n>seq1\nAAAA\n>seq2\nAAAA\n")
    recs = parse_fasta(p)
    assert len(recs) == 2
    assert recs[0][1] == "AAAA"

def test_leading_whitespace_header(tmp_path):
    p = tmp_path / "msa.fasta"
    p.write_text("  >seq1\nAAAA\n  >seq2\nAAAA\n")
    recs = parse_fasta(p)
    assert len(recs) == 2
    assert recs[0][0] == "seq1"

def test_empty_header_autonamed(tmp_path):
    p = tmp_path / "msa.fasta"
    p.write_text(">\nAAAA\n>seq2\nAAAA\n")
    recs = parse_fasta(p)
    assert recs[0][0].startswith("seq")

def test_spaced_and_lowercase_sequences(tmp_path):
    p = tmp_path / "msa.fasta"
    p.write_text(">seq1\na a a a\n>seq2\naaaa\n>seq3\nAaAa\n")
    recs = parse_fasta(p)
    assert all(s == "AAAA" for _, s in recs)

def test_multiline_sequence(tmp_path):
    p = tmp_path / "msa.fasta"
    p.write_text(">seq1\nAA\nAA\n>seq2\nAAAA\n")
    recs = parse_fasta(p)
    assert recs[0][1] == "AAAA"
    assert len(recs[0][1]) == 4

def test_differing_lengths_raises_with_detail(tmp_path):
    p = tmp_path / "msa.fasta"
    p.write_text(">seq1\nAAAA\n>seq2\nAAA\n")
    with pytest.raises(ValueError, match="differing lengths"):
        parse_fasta(p)

def test_empty_file_raises(tmp_path):
    p = tmp_path / "empty.fasta"
    p.write_text("")
    with pytest.raises(ValueError, match="No records"):
        parse_fasta(p)

def test_only_comments_raises(tmp_path):
    p = tmp_path / "only.fasta"
    p.write_text("; comment\n# comment\n\n")
    with pytest.raises(ValueError, match="No records"):
        parse_fasta(p)

def test_sequence_before_header_raises(tmp_path):
    p = tmp_path / "bad.fasta"
    p.write_text("AAAA\n>seq1\nAAAA\n")
    with pytest.raises(ValueError, match="sequence data before first header"):
        parse_fasta(p)

def test_empty_sequence_raises(tmp_path):
    p = tmp_path / "bad2.fasta"
    p.write_text(">seq1\n\n>seq2\nAAAA\n")
    # seq1 empty, should raise about empty sequence
    with pytest.raises(ValueError):
        parse_fasta(p)

def test_ambiguous_AA_ignored_in_entropy():
    # Column with B/Z/X/* should be treated as gaps -> entropy 0 if all ambiguous
    assert column_entropy("BBBB") == 0.0
    assert column_entropy("XXXX") == 0.0
    assert column_entropy("****") == 0.0
    # Mixed: A + B (ambiguous) => only A counted => entropy 0
    assert column_entropy("ABAB") == 0.0
    # A,R,N are valid, B ignored => column "ARBN" -> residues ARN => variable
    h = column_entropy("ARBN")
    assert h > 0

def test_gap_chars_both_dash_and_dot():
    # '-' and '.' both are gaps
    assert column_entropy("----") == 0.0
    assert column_entropy("....") == 0.0
    assert column_entropy("-.-.") == 0.0
    # A with gaps -> conserved
    assert column_entropy("A-A.") == 0.0
    # Mixed variable with dot gaps
    h = column_entropy("A.R.N")
    assert h > 0

def test_case_insensitive_entropy():
    assert column_entropy("aaaa") == column_entropy("AAAA")
    assert column_entropy("arnd") == column_entropy("ARND")

def test_conservation_scores_with_gaps_and_ambiguous(tmp_path):
    p = tmp_path / "msa.fasta"
    p.write_text(">seq1\nA-A.\n>seq2\nA-A.\n>seq3\nA-A.\n")
    scores = conservation_scores(p)
    # All columns either A or gaps -> conserved 1.0 for non-gap columns, gap col entropy 0 -> cons 1.0 as well (all gaps -> 0 entropy)
    assert scores[1]["conservation"] == 1.0
    assert scores[2]["conservation"] == 1.0  # all gaps -> entropy 0 -> cons 1
    assert scores[3]["conservation"] == 1.0

def test_conservation_scores_ambiguous(tmp_path):
    p = tmp_path / "mixed.fasta"
    p.write_text(">s1\nAXAA\n>s2\nARAA\n>s3\nANAA\n>s4\nADAA\n")
    scores = conservation_scores(p)
    # col2 has X,R,N,D -> X ignored => R,N,D => variable (3 distinct, entropy ~1.585) -> cons ~0.633 <0.7
    assert scores[2]["conservation"] < 0.7

def test_get_conservation_position_validation(tmp_path):
    p = tmp_path / "msa.fasta"
    p.write_text(">s1\nAAAA\n>s2\nAAAA\n")
    with pytest.raises(ValueError, match="position must be integer"):
        get_conservation_for_position(p, 0)
    with pytest.raises(ValueError, match="out of range"):
        get_conservation_for_position(p, 5)
    # offset test
    assert get_conservation_for_position(p, 1, chain_offset=1) == get_conservation_for_position(p, 2)
    with pytest.raises(ValueError, match="out of range"):
        get_conservation_for_position(p, 4, chain_offset=2)

def test_blank_lines_and_trailing_whitespace(tmp_path):
    p = tmp_path / "msa.fasta"
    p.write_text(">seq1   \n AAAA  \n\n>seq2 \nAAAA\n\n")
    recs = parse_fasta(p)
    assert len(recs) == 2
    assert recs[0][1] == "AAAA"
