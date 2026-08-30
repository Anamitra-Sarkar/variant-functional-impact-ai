import tempfile
from pathlib import Path
from data_pipeline.conservation import conservation_scores, column_entropy, LOG2_20


def write_msa(path, records):
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n{seq}\n")


def test_entropy_conserved_vs_variable():
    # Fully conserved column -> low entropy, high conservation
    assert column_entropy("AAAAA") == 0.0
    # Fully variable (5 distinct AAs equally) -> higher entropy
    h_variable = column_entropy("ARNDC")
    h_conserved = column_entropy("AAAAA")
    assert h_variable > h_conserved
    # Max entropy for 20 equally represented would be log2(20)
    h_two = column_entropy("AB" * 10)  # but B not in alphabet, filter? use valid
    # With 2 equally represented AAs, entropy =1
    assert abs(column_entropy("ARARARAR") - 1.0) < 1e-6


def test_conservation_scoring(tmp_path):
    # Synthetic MSA: col1 fully conserved, col2 variable, col3 gaps
    msa_path = tmp_path / "test.fasta"
    # 4 sequences, 5 columns
    # Col1: AAAA (conserved) -> cons ~1
    # Col2: ARND (variable) -> cons low
    # Col3: AAAA (conserved)
    # Col4: AARR (semi-conserved)
    # Col5: AAAA (conserved)
    records = [
        ("seq1", "AAAAA"),
        ("seq2", "ARAAA"),
        ("seq3", "ANAAA"),
        ("seq4", "ADAAA"),
    ]
    # Actually above: col2 = ARND? Let's design:
    records = [
        ("seq1", "AARAA"),
        ("seq2", "ARA AA".replace(" ", "")),
        ("seq3", "ANA AA".replace(" ", "")),
        ("seq4", "ADA AA".replace(" ", "")),
    ]
    # Simpler: construct directly
    records = [
        ("seq1", "AAAAA"),
        ("seq2", "ARAAA"),
        ("seq3", "ANAAA"),
        ("seq4", "ADAAA"),
    ]
    # This gives: col1 AAAA cons 1, col2 A,R,N,D variable, col3 AAAA cons1, etc.
    write_msa(msa_path, records)
    scores = conservation_scores(msa_path)
    # col1 fully conserved
    assert scores[1]["conservation"] == 1.0
    assert scores[1]["entropy"] == 0.0
    # col2 variable (4 distinct AAs: A,R,N,D) -> entropy=2.0, IC~2.32, cons~0.537
    # So threshold <0.6 is appropriate (not <0.5)
    assert scores[2]["conservation"] < 0.6
    assert scores[2]["entropy"] > 1.5
    # col3 conserved again
    assert scores[3]["conservation"] == 1.0


def test_gap_handling(tmp_path):
    msa_path = tmp_path / "gaps.fasta"
    records = [
        ("seq1", "A-A"),
        ("seq2", "A-A"),
        ("seq3", "A-A"),
    ]
    write_msa(msa_path, records)
    scores = conservation_scores(msa_path)
    # Column 2 all gaps -> entropy 0, conservation 1? But gaps ignored -> no residues -> entropy 0
    assert scores[2]["entropy"] == 0.0
