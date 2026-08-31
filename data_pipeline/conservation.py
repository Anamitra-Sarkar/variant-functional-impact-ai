"""
Real MSA-entropy conservation scoring.

Computes per-column Shannon entropy and information content from a FASTA MSA.
Real, correctly-implemented information theory; no fabricated math.

References:
- Shannon 1948, Schneider & Stephens 1990, ConSurf (Ashkenazy et al. 2016)
- phyloP/phastCons (Pollard et al. 2010, Siepel et al. 2005) documented as alternative source
"""

from __future__ import annotations

import math
from collections import Counter
from pathlib import Path

AMINO_ACIDS = set("ARNDCQEGHILKMFPSTWYV")
ALPHABET_SIZE = 20
LOG2_20 = math.log2(20)  # ~4.3219
GAP_CHARS = set("-.")


def parse_fasta(fasta_path: str | Path) -> list[tuple[str, str]]:
    """Parse FASTA file, returns list of (header, sequence).

    Robust to real-file quirks:
    - Comment lines starting with ';' or '#' are ignored (common in some FASTA/MSA exports)
    - Header variants: '>' may have leading whitespace, header may be empty (auto-named)
    - Sequence lines: whitespace inside sequence is stripped, lower-case is upper-cased
    - Multiline sequences, blank lines, and trailing whitespace are handled
    - Validates all sequences have same aligned length (MSA invariant)
    """
    path = Path(fasta_path)
    if not path.exists():
        raise FileNotFoundError(f"MSA file not found: {path}")
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []
    # Count headers for auto-naming empty headers
    auto_idx = 0
    with open(path, encoding="utf-8", errors="strict") as f:
        for raw_line in f:
            # Handle CRLF and preserve detection of header with leading spaces
            stripped = raw_line.strip()
            if not stripped:
                continue
            # Comment lines (some tools emit ';' comments before FASTA)
            if stripped.startswith(";") or stripped.startswith("#"):
                continue
            # Header detection: allow leading whitespace before '>'
            lstrip = raw_line.lstrip()
            if lstrip.startswith(">"):
                # Save previous record
                if header is not None:
                    seq = "".join(seq_parts).upper().replace(" ", "").replace("\t", "")
                    # Remove any remaining internal whitespace
                    seq = "".join(seq.split())
                    if not seq:
                        raise ValueError(f"Empty sequence for header '{header}' in {fasta_path}")
                    records.append((header, seq))
                # Parse new header
                h = lstrip[1:].strip()
                if not h:
                    auto_idx += 1
                    h = f"seq{auto_idx}"
                header = h
                seq_parts = []
            else:
                if header is None:
                    raise ValueError(f"FASTA parse error in {fasta_path}: sequence data before first header line: '{stripped[:60]}'")
                # Strip all whitespace inside sequence line (some exports have spaced alignments)
                cleaned = "".join(stripped.split())
                if cleaned:
                    seq_parts.append(cleaned)
        if header is not None:
            seq = "".join(seq_parts).upper().replace(" ", "").replace("\t", "")
            seq = "".join(seq.split())
            if not seq:
                raise ValueError(f"Empty sequence for header '{header}' in {fasta_path}")
            records.append((header, seq))
    if not records:
        raise ValueError(f"No records found in {fasta_path} (empty or only comments/blank lines)")
    # Validate all same length (MSA) – report lengths and offending headers
    lengths = {len(s) for _, s in records}
    if len(lengths) > 1:
        detail = ", ".join(f"{h}={len(s)}" for h, s in records)
        raise ValueError(f"MSA sequences have differing lengths: {lengths} ({detail})")
    if len(records[0][1]) == 0:
        raise ValueError(f"MSA has zero-length alignment in {fasta_path}")
    return records


def column_entropy(column: str) -> float:
    """
    Shannon entropy H = -Σ p_a * log2(p_a) for a column string.
    Gaps ('-', '.') and ambiguous/non-canonical codes (B,Z,X,*,J,O,U) are ignored.
    If all gaps/ambiguous, returns 0 (fully non-informative column).
    Input is case-insensitive; column chars are upper-cased before evaluation.
    """
    # Filter gaps and non-canonical; upper-case first
    residues = [c.upper() for c in column if c not in GAP_CHARS]
    # Keep only canonical 20 AAs
    residues = [c for c in residues if c in AMINO_ACIDS]
    if not residues:
        return 0.0
    n = len(residues)
    counts = Counter(residues)
    h = 0.0
    for cnt in counts.values():
        p = cnt / n
        h -= p * math.log2(p)
    return h


def conservation_scores(msa_path: str | Path) -> dict[int, dict]:
    """
    Compute per-column conservation metrics.
    Returns dict position (1-indexed) -> {entropy, information_content, normalized_conservation}

    - entropy: H in [0, log2(20)] bits
    - information_content: IC = log2(20) - H (0..log2(20)), higher = more conserved
    - normalized_conservation: IC / log2(20) in [0,1], 1 = fully conserved, 0 = fully variable
    """
    records = parse_fasta(msa_path)
    seqs = [s for _, s in records]
    aln_len = len(seqs[0])
    result: dict[int, dict] = {}
    for col_idx in range(aln_len):
        column = "".join(seq[col_idx] for seq in seqs)
        h = column_entropy(column)
        ic = LOG2_20 - h
        norm = ic / LOG2_20 if LOG2_20 else 0.0
        # Clamp to [0,1]
        norm = max(0.0, min(1.0, norm))
        result[col_idx + 1] = {
            "entropy": h,
            "information_content": ic,
            "conservation": norm,
            "column": column,
        }
    return result


def get_conservation_for_position(msa_path: str | Path, position: int, chain_offset: int = 0) -> float:
    """
    Convenience: get normalized conservation for a 1-indexed protein position.
    chain_offset allows mapping if MSA columns != protein positions (e.g. offset).
    Validates position >=1; raises ValueError with available range on miss.
    """
    if not isinstance(position, int) or position < 1:
        raise ValueError(f"position must be integer >=1, got {position}")
    scores = conservation_scores(msa_path)
    # Simple mapping: MSA column = position + offset
    col = position + chain_offset
    if col not in scores:
        raise ValueError(f"Position {position} (col {col}) out of range (MSA length {len(scores)}, valid 1..{len(scores)})")
    return float(scores[col]["conservation"])
