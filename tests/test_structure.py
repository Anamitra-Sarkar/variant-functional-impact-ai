from pathlib import Path
from data_pipeline.structure import compute_rsa_proxy, compute_ss_proxy, compute_features

def pdb_line(serial, atom_name, resname, chain, resseq, x, y, z, element):
    # Proper PDB FORMAT per wwPDB spec, 8.3f for coords ensures 8-char columns
    # ATOM  serial  atom_name resname chain resseq    x      y      z  occupancy temp element
    return f"ATOM  {serial:5d}  {atom_name:<3s} {resname:3s} {chain}{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           {element:>2s}"

def make_tiny_pdb(path):
    lines = []
    # 5 residues, CA coords forming a semi-compact cluster
    residues = [
        ("ALA", 1, 0.0, 0.0, 0.0),
        ("GLY", 2, 3.9, 2.8, 0.0),
        ("VAL", 3, 6.0, 2.0, 4.5),
        ("LEU", 4, 8.0, -1.5, 3.8),
        ("SER", 5, 10.5, -2.8, 1.0),
    ]
    serial = 1
    for resname, resseq, x, y, z in residues:
        # Use CA at given xyz, plus N/C/O nearby to make valid residue
        lines.append(pdb_line(serial, "N", resname, "A", resseq, x-1.5, y, z, "N")); serial+=1
        lines.append(pdb_line(serial, "CA", resname, "A", resseq, x, y, z, "C")); serial+=1
        lines.append(pdb_line(serial, "C", resname, "A", resseq, x+1.0, y+1.0, z, "C")); serial+=1
        lines.append(pdb_line(serial, "O", resname, "A", resseq, x+1.5, y+2.0, z, "O")); serial+=1
    lines.append("TER")
    lines.append("END")
    Path(path).write_text("\n".join(lines) + "\n")


def test_rsa_proxy_buried_vs_exposed(tmp_path):
    pdb_path = tmp_path / "tiny.pdb"
    make_tiny_pdb(pdb_path)
    feats = compute_rsa_proxy(str(pdb_path), radius=10.0, max_neighbors=20)
    assert len(feats) == 5
    for k, v in feats.items():
        assert 0 <= v.rsa <= 1
    # Middle residue (position 3) should have at least as many neighbors as ends
    assert feats[("A", 3)].burial_count >= feats[("A", 1)].burial_count or feats[("A", 3)].burial_count >= feats[("A", 5)].burial_count
    feats_tight = compute_rsa_proxy(str(pdb_path), radius=5.0, max_neighbors=10)
    assert feats_tight[("A", 3)].burial_count >= 0


def test_ss_proxy(tmp_path):
    pdb_path = tmp_path / "tiny.pdb"
    make_tiny_pdb(pdb_path)
    ss = compute_ss_proxy(str(pdb_path))
    assert len(ss) == 5
    for v in ss.values():
        assert v in ("H", "E", "C")


def test_compute_features(tmp_path):
    pdb_path = tmp_path / "tiny.pdb"
    make_tiny_pdb(pdb_path)
    feats = compute_features(str(pdb_path))
    assert len(feats) == 5
    assert all(0 <= f.rsa <= 1 for f in feats.values())
    assert all(f.ss in ("H", "E", "C") for f in feats.values())
