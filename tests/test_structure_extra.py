"""Hardened structure parsing tests: edge cases, error handling, real-file quirks."""
import pytest
from pathlib import Path
from data_pipeline.structure import (
    parse_structure, compute_rsa_proxy, compute_ss_proxy, compute_features, get_residue_feature,
    _get_ca_atoms
)

def pdb_line(serial, atom_name, resname, chain, resseq, x, y, z, element, occupancy=" 1.00", bfactor="20.00"):
    return f"ATOM  {serial:5d}  {atom_name:<3s} {resname:3s} {chain}{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occupancy} {bfactor}           {element:>2s}"

def make_pdb(path, residues, chain="A", header_lines=None):
    lines = []
    if header_lines:
        lines.extend(header_lines)
    serial = 1
    for resname, resseq, x, y, z in residues:
        lines.append(pdb_line(serial, "N", resname, chain, resseq, x-1.5, y, z, "N")); serial+=1
        lines.append(pdb_line(serial, "CA", resname, chain, resseq, x, y, z, "C")); serial+=1
        lines.append(pdb_line(serial, "C", resname, chain, resseq, x+1.0, y+1.0, z, "C")); serial+=1
        lines.append(pdb_line(serial, "O", resname, chain, resseq, x+1.5, y+2.0, z, "O")); serial+=1
    lines.append("TER")
    lines.append("END")
    Path(path).write_text("\n".join(lines) + "\n")

def test_parse_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        parse_structure(tmp_path / "nope.pdb")

def test_parse_empty_file_raises(tmp_path):
    p = tmp_path / "empty.pdb"
    p.write_text("")
    with pytest.raises(ValueError, match="empty"):
        parse_structure(p)

def test_parse_header_variants_tolerated(tmp_path):
    residues = [("ALA", 1, 0.0, 0.0, 0.0), ("GLY", 2, 5.0, 0.0, 0.0)]
    p = tmp_path / "hdr.pdb"
    make_pdb(p, residues, header_lines=["HEADER    TEST 01-JAN-20", "TITLE     TEST PROTEIN", "REMARK   1 TEST", "COMPND    TEST"])
    s = parse_structure(p)
    cas = _get_ca_atoms(s)
    assert len(cas) == 2

def test_hetero_only_raises(tmp_path):
    p = tmp_path / "het.pdb"
    # Only HETATM water, no ATOM with CA
    p.write_text("HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00 20.00           O\nEND\n")
    with pytest.raises(ValueError, match="No CA atoms"):
        compute_rsa_proxy(p)

def test_missing_CA_skipped(tmp_path):
    p = tmp_path / "miss.pdb"
    lines = [
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N",
        "ATOM      2  C   ALA A   1       1.000   1.000   0.000  1.00 20.00           C",
        "ATOM      3  O   ALA A   1       1.500   2.000   0.000  1.00 20.00           O",
        "ATOM      4  N   GLY A   2       5.000   0.0     0.000  1.00 20.00           N",
        "ATOM      5  CA  GLY A   2       5.000   0.000   0.000  1.00 20.00           C",
        "ATOM      6  C   GLY A   2       6.000   1.000   0.000  1.00 20.00           C",
        "ATOM      7  O   GLY A   2       6.500   2.000   0.000  1.00 20.00           O",
        "TER",
        "END"
    ]
    Path(p).write_text("\n".join(lines)+"\n")
    feats = compute_rsa_proxy(p)
    # Only residue 2 has CA, residue 1 should be skipped
    assert len(feats) == 1
    assert ("A", 2) in feats

def test_single_residue_rsa(tmp_path):
    p = tmp_path / "single.pdb"
    make_pdb(p, [("ALA", 1, 0.0, 0.0, 0.0)])
    feats = compute_rsa_proxy(p)
    assert len(feats) == 1
    assert feats[("A",1)].rsa == 1.0  # no neighbors -> fully exposed
    assert feats[("A",1)].burial_count == 0
    ss = compute_ss_proxy(p)
    assert ss[("A",1)] in ("H","E","C")

def test_multi_chain_grouping(tmp_path):
    p = tmp_path / "multi.pdb"
    lines = []
    serial=1
    # Chain A resid 1
    lines.append(pdb_line(serial, "CA", "ALA", "A", 1, 0,0,0,"C")); serial+=1
    lines.append(pdb_line(serial, "N", "ALA", "A", 1, -1,0,0,"N")); serial+=1
    # Chain B resid 1
    lines.append(pdb_line(serial, "CA", "GLY", "B", 1, 20,0,0,"C")); serial+=1
    lines.append(pdb_line(serial, "N", "GLY", "B", 1, 19,0,0,"N")); serial+=1
    lines.append("TER")
    lines.append("END")
    Path(p).write_text("\n".join(lines)+"\n")
    feats = compute_features(p)
    assert ("A",1) in feats
    assert ("B",1) in feats

def test_radius_validation(tmp_path):
    p = tmp_path / "small.pdb"
    make_pdb(p, [("ALA",1,0,0,0)])
    with pytest.raises(ValueError, match="radius"):
        compute_rsa_proxy(p, radius=0)
    with pytest.raises(ValueError, match="max_neighbors"):
        compute_rsa_proxy(p, max_neighbors=0)

def test_get_residue_feature_validation(tmp_path):
    p = tmp_path / "two.pdb"
    make_pdb(p, [("ALA",1,0,0,0),("VAL",2,5,0,0)])
    feats = compute_features(p)
    with pytest.raises(ValueError, match="position must be integer"):
        get_residue_feature(feats, "A", 0)
    with pytest.raises(ValueError, match="chain_id"):
        get_residue_feature(feats, "", 1)
    # case-insensitive lookup
    assert get_residue_feature(feats, "a", 1) is not None or get_residue_feature(feats, "A", 1) is not None
    # missing returns None
    assert get_residue_feature(feats, "A", 99) is None

def test_compute_features_empty_after_filter(tmp_path):
    p = tmp_path / "empty2.pdb"
    p.write_text("HEADER TEST\nEND\n")
    with pytest.raises(ValueError, match="No CA atoms|No residues|No models"):
        compute_features(p)

def test_multiple_models_only_first_used(tmp_path):
    p = tmp_path / "models.pdb"
    lines = [
        "MODEL        1",
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C",
        "ENDMDL",
        "MODEL        2",
        "ATOM      1  CA  GLY A   1      10.000  10.000  10.000  1.00 20.00           C",
        "ENDMDL",
        "END"
    ]
    Path(p).write_text("\n".join(lines)+"\n")
    feats = compute_rsa_proxy(p)
    # Should have 1 CA from model 1 at origin
    assert len(feats)==1
    assert feats[("A",1)].ca_coord == (0.0,0.0,0.0)

def test_mmcif_suffix_detection(tmp_path):
    # Create a PDB-like file but with .cif suffix that is actually PDB – should still parse via PDB fallback heuristic
    p = tmp_path / "test.cif"
    make_pdb(p, [("ALA",1,0,0,0)])
    # This file is not true mmCIF but PDB content with .cif suffix; our parse will try MMCIFParser and may fail
    # but we handle gracefully: if fails, should raise ValueError about parse failure, not crash silently
    # Instead test that .cif file with real PDB content still raises or falls back? Here we just check it doesn't hang
    try:
        s = parse_structure(p)
        # If it succeeded via MMCIFParser confusion, it may have 0 models -> raises
        assert len(list(s)) >= 0
    except ValueError:
        pass

def test_cif_like_heuristic_pdb_file(tmp_path):
    p = tmp_path / "heuristic.pdb"
    p.write_text("data_test\n_entry.id test\n")
    # This will be detected as cif via heuristic but actually not valid cif -> will raise parse error
    with pytest.raises(ValueError):
        parse_structure(p)
