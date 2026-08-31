"""
Real protein structure parsing with Bio.PDB, RSA/secondary-structure computation.

Production path: DSSP via Bio.PDB.DSSP (requires mkdssp binary).
Sandbox proxy: neighbor-count-based burial estimate (real geometry, documented).

Usage:
  from data_pipeline.structure import parse_structure, compute_rsa_proxy, compute_ss_proxy
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

try:
    from Bio.PDB import PDBParser, MMCIFParser
    from Bio.PDB.DSSP import DSSP
    BIOPDB_AVAILABLE = True
except ImportError:
    BIOPDB_AVAILABLE = False

# Max ASA per residue (Tien et al. 2013 / Miller et al.) for RSA normalization
# Used when DSSP ASA is available. For proxy we don't need this, but store for DSSP path.
MAX_ASA_TIEN: Dict[str, float] = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "Q": 225.0, "E": 223.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}


@dataclass
class ResidueFeatures:
    chain_id: str
    resseq: int
    resname: str  # 3-letter
    rsa: float  # 0..1 (0 buried, 1 exposed)
    ss: str  # H / E / C (helix/sheet/coil)
    burial_count: int  # raw neighbor count (for debugging)
    ca_coord: Tuple[float, float, float]


def parse_structure(structure_path: str | Path):
    """Parse PDB or mmCIF file using Bio.PDB. Returns Bio.PDB Structure object.

    Handles real-file quirks:
    - Detects mmCIF by suffix .cif/.mmcif (also .CIF case-insensitive)
    - PDBParser is tolerant to missing optional columns (occupancy/B-factor/element)
      and multiple header variants (HEADER, TITLE, REMARK, etc. are ignored by Biopython)
    - Supports .ent (PDB flatfile alternative suffix)
    - Validates file is non-empty; raises informative error if no atoms parsed
    """
    if not BIOPDB_AVAILABLE:
        raise ImportError("biopython not installed; pip install biopython")
    path = Path(structure_path)
    if not path.exists():
        raise FileNotFoundError(f"Structure file not found: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"Structure file is empty: {path}")
    suffix = path.suffix.lower()
    # Also handle .ent (common PDB mirror suffix) and compressed-like names
    if suffix in (".cif", ".mmcif", ".ent"):
        # .ent files are PDB format, not CIF
        is_cif = suffix in (".cif", ".mmcif")
    else:
        is_cif = suffix in (".cif", ".mmcif")
        # Fallback heuristic: inspect first non-comment line for mmCIF keyword
        if not is_cif:
            try:
                with open(path, encoding="utf-8", errors="ignore") as fh:
                    for line in fh:
                        s = line.strip()
                        if not s or s.startswith("#"):
                            continue
                        if s.startswith("data_") or s.startswith("_entry."):
                            is_cif = True
                        break
            except Exception:
                pass
    if is_cif:
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("protein", str(path))
    except Exception as e:
        raise ValueError(f"Failed to parse structure file {path}: {e}") from e
    # Validate at least one model with content
    if len(list(structure)) == 0:
        raise ValueError(f"No models found in structure file: {path}")
    return structure


def _get_ca_atoms(structure) -> list:
    """Collect all CA atoms with their residue info (first model only).

    Filters:
    - Only standard residues (residue.id[0] == ' ') – excludes hetero/water (HETATM)
    - Requires CA atom present; handles altloc implicitly via Biopython's default
    - Preserves insertion code info: residue.id = (hetfield, resseq, icode); we key by
      (chain.id, resseq) so insertion variants at same resseq are treated as one
      (first encountered). This matches DSSP's resseq-centric lookup.
    """
    cas = []
    for model in structure:
        for chain in model:
            for residue in chain:
                # Skip hetero residues / water (hetfield != ' ')
                if residue.id[0] != " ":
                    continue
                if "CA" not in residue:
                    continue
                try:
                    coord = residue["CA"].get_coord()
                except Exception:
                    continue
                cas.append((chain.id, residue.id[1], residue.resname, coord, residue))
        break  # only first model
    return cas


def compute_rsa_proxy(structure_path: str | Path, radius: float = 10.0, max_neighbors: int = 20) -> Dict[Tuple[str, int], ResidueFeatures]:
    """
    Neighbor-count-based RSA proxy (real geometry).
    For each CA, count neighboring CAs within `radius` Å.
    burial_fraction = min(count / max_neighbors, 1.0)
    RSA = 1 - burial_fraction
    This is a documented simplified geometric proxy when DSSP/mkdssp is not available.
    Production DSSP path would use actual SASA via DSSP.

    Returns dict keyed by (chain_id, resseq).
    Raises ValueError if structure has no CA atoms (e.g. hetero-only file).
    """
    if radius <= 0:
        raise ValueError(f"radius must be >0, got {radius}")
    if max_neighbors <= 0:
        raise ValueError(f"max_neighbors must be >0, got {max_neighbors}")
    structure = parse_structure(structure_path)
    cas = _get_ca_atoms(structure)
    if not cas:
        raise ValueError(f"No CA atoms found in structure file: {structure_path} (empty or hetero-only)")
    result: Dict[Tuple[str, int], ResidueFeatures] = {}
    for idx, (chain_id, resseq, resname, coord, _res) in enumerate(cas):
        count = 0
        for jdx, (_, _, _, coord2, _) in enumerate(cas):
            if idx == jdx:
                continue
            dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(coord, coord2)))
            if dist < radius:
                count += 1
        burial = min(count / max_neighbors, 1.0) if max_neighbors > 0 else 0.0
        rsa = 1.0 - burial
        # Clamp to [0,1] defensively
        rsa = max(0.0, min(1.0, rsa))
        result[(chain_id, resseq)] = ResidueFeatures(
            chain_id=chain_id,
            resseq=resseq,
            resname=resname,
            rsa=rsa,
            ss="C",  # placeholder; updated by compute_ss_proxy
            burial_count=count,
            ca_coord=(float(coord[0]), float(coord[1]), float(coord[2])),
        )
    return result


def compute_ss_proxy(structure_path: str | Path) -> Dict[Tuple[str, int], str]:
    """
    Simplified secondary-structure proxy based on local CA geometry.
    NOT a replacement for DSSP; documented as proxy for sandbox.

    Heuristic:
    - Helix: CA(i) - CA(i+3) ~5.5 Å and CA(i)-CA(i+4) ~6.5 Å (alpha helix pitch)
    - Sheet: extended distances >9 Å alternating
    - Otherwise coil.

    Returns dict (chain_id, resseq) -> 'H'|'E'|'C'.
    Production path uses DSSP 8-state -> 3-state.
    Raises ValueError if no CA atoms found.
    """
    structure = parse_structure(structure_path)
    cas = _get_ca_atoms(structure)
    if not cas:
        raise ValueError(f"No CA atoms found in structure file: {structure_path}")
    # Group by chain
    from collections import defaultdict
    by_chain: Dict[str, list] = defaultdict(list)
    for entry in cas:
        by_chain[entry[0]].append(entry)

    result: Dict[Tuple[str, int], str] = {}
    for chain_id, entries in by_chain.items():
        entries.sort(key=lambda x: x[1])
        n = len(entries)
        for i in range(n):
            # Default to coil
            ss = "C"
            if i + 3 < n:
                d1 = math.sqrt(sum((a - b) ** 2 for a, b in zip(entries[i][3], entries[i+3][3])))
                # Helix characteristic: i to i+3 ~5.5, i to i+4 ~6.3
                if 4.5 < d1 < 6.5:
                    if i + 4 < n:
                        d2 = math.sqrt(sum((a - b) ** 2 for a, b in zip(entries[i][3], entries[i+4][3])))
                        if 5.0 < d2 < 7.5:
                            ss = "H"
                    else:
                        ss = "H"
                # Extended sheet: larger separation
                elif d1 > 8.5:
                    # Mark as strand if consistently extended
                    ss = "E"
            result[(chain_id, entries[i][1])] = ss
    return result


def compute_features(structure_path: str | Path, use_dssp: bool = False, dssp_executable: str = "mkdssp") -> Dict[Tuple[str, int], ResidueFeatures]:
    """
    Main entry: compute RSA + SS per residue.
    If use_dssp=True and mkdssp available, uses DSSP; otherwise uses proxies.

    DSSP path (production, requires mkdssp):
      DSSP returns ASA; RSA = ASA / MAX_ASA[aa]. SS 8-state -> 3-state.
    Proxy path (sandbox):
      neighbor-count RSA + geometric SS proxy.
    Raises ValueError if no CA atoms found.
    """
    if use_dssp and BIOPDB_AVAILABLE:
        try:
            structure = parse_structure(structure_path)
            # Ensure at least one model with content
            try:
                model = next(iter(structure))
            except StopIteration:
                raise ValueError(f"No models in structure: {structure_path}")
            dssp = DSSP(model, str(structure_path), dssp=dssp_executable)
            # dssp keys: (chain_id, (resseq, icode))
            result: Dict[Tuple[str, int], ResidueFeatures] = {}
            cas = { (c, r): coord for c, r, _, coord, _ in _get_ca_atoms(structure) }
            for key, val in dssp.property_dict.items():
                # val layout: (aa, ss, acc, phi, psi, ...) but check Bio.PDB version
                # Use dssp dict: dssp[key] = (dssp_index, aa, ss, rsa, phi, psi, ...)
                # Simpler: access via dssp[key]
                dssp_row = dssp[key]
                aa = dssp_row[1]
                ss8 = dssp_row[2]
                asa = dssp_row[3]
                chain_id = key[0]
                resseq = key[1][1]
                # Map 8-state to 3-state
                if ss8 in ("H", "G", "I"):
                    ss3 = "H"
                elif ss8 in ("E", "B"):
                    ss3 = "E"
                else:
                    ss3 = "C"
                max_asa = MAX_ASA_TIEN.get(aa, 150.0)
                rsa = min(asa / max_asa, 1.0) if max_asa else 0.0
                rsa = max(0.0, min(1.0, rsa))
                coord = cas.get((chain_id, resseq), (0, 0, 0))
                # 3-letter resname not directly available; use aa
                result[(chain_id, resseq)] = ResidueFeatures(
                    chain_id=chain_id,
                    resseq=resseq,
                    resname=aa,
                    rsa=rsa,
                    ss=ss3,
                    burial_count=-1,
                    ca_coord=(float(coord[0]), float(coord[1]), float(coord[2])) if hasattr(coord, '__iter__') else (0, 0, 0),
                )
            if result:
                return result
        except Exception as e:
            # Fall through to proxy – log but do not hide programming errors for empty files
            msg = str(e)
            if "No CA atoms" in msg or "No models" in msg or "empty" in msg.lower():
                raise
            print(f"DSSP failed ({e}), falling back to geometric proxy")
            pass

    # Proxy path
    rsa_dict = compute_rsa_proxy(structure_path)
    ss_dict = compute_ss_proxy(structure_path)
    for key, ss in ss_dict.items():
        if key in rsa_dict:
            rsa_dict[key].ss = ss
    if not rsa_dict:
        raise ValueError(f"No residues with features computed for {structure_path}")
    return rsa_dict


def get_residue_feature(features: Dict[Tuple[str, int], ResidueFeatures], chain_id: str, position: int) -> Optional[ResidueFeatures]:
    """Lookup by chain and 1-indexed position. Chain lookup is case-sensitive (PDB convention); also tries upper."""
    if not isinstance(position, int) or position < 1:
        raise ValueError(f"position must be integer >=1, got {position}")
    if not chain_id or not isinstance(chain_id, str):
        raise ValueError(f"chain_id must be non-empty string, got {chain_id!r}")
    feat = features.get((chain_id, position))
    if feat is None and chain_id.upper() != chain_id:
        feat = features.get((chain_id.upper(), position))
    if feat is None and chain_id.lower() != chain_id:
        feat = features.get((chain_id.lower(), position))
    return feat
