"""
CLI entry points for real-run procedure (Kaggle/Modal).

Examples:
  python -m data_pipeline.cli --structure-path AF-P04637.pdb --msa-path orthologs.fasta --train --labels clinvar_labels.csv
  python -m data_pipeline.cli --structure-path AF-P04637.pdb --msa-path orthologs.fasta --predict --variant P53_R175H
  python -m data_pipeline.cli --structure-path tests/fixtures/tiny.pdb --msa-path tests/fixtures/tiny.fasta --predict --variant TEST_2_A_V

For synthetic verification without real downloads, labels can be generated via --synthetic-labels.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

from .structure import compute_features
from .conservation import conservation_scores
from .model import featurize_variant, featurize_batch, train_logistic_regression, evaluate_with_baselines, save_artifacts, load_artifacts, FEATURE_NAMES
from .physicochemical import KYTE_DOOLITTLE, VDW_VOLUME


def parse_variant_string(s: str) -> tuple[str, int, str, str]:
    """Parse 'P53_R175H' or 'TEST_2_A_V' or 'R175H' -> (protein, pos, wt, mut)"""
    # Try formats: GENE_WTposMUT or WTposMUT
    s = s.strip()
    protein = "UNKNOWN"
    var_part = s
    if "_" in s:
        # Could be GENE_R175H or GENE_2_A_V
        parts = s.split("_")
        if len(parts) == 2 and len(parts[1]) >= 2:
            protein, var_part = parts[0], parts[1]
        elif len(parts) == 4:
            # TEST_2_A_V -> protein TEST, pos 2, wt A, mut V
            protein = parts[0]
            var_part = f"{parts[2]}{parts[1]}{parts[3]}"
        else:
            # Fallback: last part is variant
            var_part = parts[-1]
            protein = "_".join(parts[:-1])
    # Now var_part like R175H
    import re
    m = re.match(r"^([A-Z])(\d+)([A-Z])$", var_part.strip().upper())
    if not m:
        raise ValueError(f"Cannot parse variant string: {s} (expected WTposMUT like R175H)")
    wt, pos, mut = m.group(1), int(m.group(2)), m.group(3)
    return protein, pos, wt, mut


def main():
    parser = argparse.ArgumentParser(description="Variant impact pipeline CLI")
    parser.add_argument("--structure-path", type=str, help="Path to PDB/mmCIF structure file")
    parser.add_argument("--msa-path", type=str, help="Path to FASTA MSA file")
    parser.add_argument("--chain", type=str, default="A", help="Chain ID (default A)")
    parser.add_argument("--labels", type=str, help="CSV with columns: variant, label (0/1)")
    parser.add_argument("--train", action="store_true", help="Train fusion model")
    parser.add_argument("--predict", action="store_true", help="Predict single variant")
    parser.add_argument("--variant", type=str, help="Variant string e.g. P53_R175H or R175H")
    parser.add_argument("--model-out", type=str, default="artifacts/model.pkl", help="Model output path")
    parser.add_argument("--model-in", type=str, help="Model input path for predict")
    parser.add_argument("--use-dssp", action="store_true", help="Use DSSP if mkdssp available")
    parser.add_argument("--synthetic-labels", action="store_true", help="Generate synthetic injected-signal labels for verification")
    args = parser.parse_args()

    if not args.structure_path and not args.msa_path:
        parser.print_help()
        sys.exit(0)

    # Load features if paths provided
    features = None
    cons_scores = None
    if args.structure_path:
        print(f"[info] Loading structure: {args.structure_path}")
        features = compute_features(args.structure_path, use_dssp=args.use_dssp)
        print(f"[info] Computed features for {len(features)} residues")
    if args.msa_path:
        print(f"[info] Loading MSA: {args.msa_path}")
        cons_scores = conservation_scores(args.msa_path)
        print(f"[info] MSA length {len(cons_scores)} columns")

    if args.train:
        # Need labels
        if args.synthetic_labels:
            # Generate synthetic labels: buried+conserved+large physchem delta -> damaging
            # For demo, create grid of synthetic variants
            print("[info] Generating synthetic injected-signal labels")
            variants = []
            labels = []
            # Use residues 1..len(features) or MSA length
            n = len(cons_scores) if cons_scores else (len(features) if features else 20)
            import random
            random.seed(42)
            aas = list(KYTE_DOOLITTLE.keys())
            for i in range(100):
                pos = (i % n) + 1
                rsa_val = 0.1 if i % 3 == 0 else 0.7  # alternate buried/exposed
                cons_val = 0.9 if i % 3 == 0 else 0.2
                # large vs small change
                if i % 2 == 0:
                    wt, mut = "R", "W"  # large Grantham, large hydro change
                    label = 1 if (rsa_val < 0.3 and cons_val > 0.7) else 0
                else:
                    wt, mut = "S", "T"  # small change
                    label = 0
                # But for injected signal we want correlation: damaging = buried+conserved+large
                # Simplify: label = 1 if i%3==0 and i%2==0 else 0
                # Actually generate deterministically
                is_damaging_context = (rsa_val < 0.3 and cons_val > 0.7)
                is_large_change = (wt == "R" and mut == "W")
                label = 1 if (is_damaging_context and is_large_change) else 0
                # Build mock residue feature
                from data_pipeline.structure import ResidueFeatures
                rf = ResidueFeatures(chain_id=args.chain, resseq=pos, resname="ALA", rsa=rsa_val, ss="H" if pos % 2 == 0 else "C", burial_count=15 if rsa_val < 0.3 else 2, ca_coord=(0, 0, 0))
                var = {"wt": wt, "mut": mut, "residue_feature": rf, "conservation": cons_val}
                variants.append(var)
                labels.append(label)
            X, _ = featurize_batch(variants)
            y = np.array(labels)
            print(f"[info] Synthetic: {sum(labels)} damaging / {len(labels)-sum(labels)} benign")
            model = train_logistic_regression(X, y)
            results = evaluate_with_baselines(X, y, FEATURE_NAMES, model)
            print(json.dumps(results, indent=2))
            save_artifacts(model, FEATURE_NAMES, args.model_out)
            print(f"[info] Saved model to {args.model_out}")
        elif args.labels:
            # Load real labels CSV
            rows = []
            with open(args.labels) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    rows.append(row)
            variants = []
            labels = []
            for row in rows:
                var_str = row.get("variant") or row.get("Variant") or row.get("id")
                label = int(row.get("label") or row.get("Label") or row.get("pathogenic"))
                protein, pos, wt, mut = parse_variant_string(var_str)
                chain = row.get("chain", args.chain)
                rf = features.get((chain, pos)) if features else None
                cons = cons_scores.get(pos, {}).get("conservation", 0.5) if cons_scores else 0.5
                variants.append({"wt": wt, "mut": mut, "residue_feature": rf, "conservation": cons})
                labels.append(label)
            X, _ = featurize_batch(variants)
            y = np.array(labels)
            model = train_logistic_regression(X, y)
            results = evaluate_with_baselines(X, y, FEATURE_NAMES, model)
            print(json.dumps(results, indent=2))
            save_artifacts(model, FEATURE_NAMES, args.model_out)
            print(f"[info] Saved model to {args.model_out}")
        else:
            print("[error] --train requires --labels or --synthetic-labels")
            sys.exit(1)

    if args.predict:
        if not args.variant:
            print("[error] --predict requires --variant")
            sys.exit(1)
        protein, pos, wt, mut = parse_variant_string(args.variant)
        print(f"[info] Variant: {protein} {wt}{pos}{mut} chain {args.chain}")
        rf = features.get((args.chain, pos)) if features else None
        if rf:
            print(f"[info] RSA={rf.rsa:.3f} SS={rf.ss} resname={rf.resname}")
        else:
            print(f"[warn] No structure feature for chain {args.chain} pos {pos}; using defaults")
        cons = cons_scores.get(pos, {}).get("conservation", 0.5) if cons_scores else 0.5
        print(f"[info] Conservation={cons:.3f}")

        model_path = args.model_in or args.model_out
        try:
            model, _ = load_artifacts(model_path)
            X, _ = featurize_batch([{"wt": wt, "mut": mut, "residue_feature": rf, "conservation": cons}])
            score = float(model.predict_proba(X)[0, 1])
            print(json.dumps({"variant": args.variant, "score": score, "rsa": rf.rsa if rf else None, "conservation": cons, "ss": rf.ss if rf else None}, indent=2))
        except FileNotFoundError:
            # No model artifact; show per-feature explanation without fusion score
            from data_pipeline.physicochemical import delta_hydrophobicity, delta_volume, grantham_distance
            dh = delta_hydrophobicity(wt, mut)
            dv = delta_volume(wt, mut)
            gd = grantham_distance(wt, mut)
            print(json.dumps({"variant": args.variant, "score": None, "note": "no model artifact; showing per-feature deltas", "rsa": rf.rsa if rf else None, "conservation": cons, "delta_hydro": dh, "delta_vol": dv, "grantham": gd}, indent=2))


if __name__ == "__main__":
    main()
