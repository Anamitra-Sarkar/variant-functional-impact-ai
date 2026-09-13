"""
Real training run for Variant Functional Impact AI, on Modal CPU.

Builds a real, modest-sized labeled dataset from public sources documented in
docs/data_sources.md:
  - AlphaFold DB structures (per-UniProt PDB) for RSA/secondary-structure proxy
    features (data_pipeline.structure)
  - UniProt REST orthologs, aligned with MAFFT, for conservation scores
    (data_pipeline.conservation)
  - ClinVar variant_summary.txt.gz for real pathogenic/benign missense labels

Scope: a modest, well-characterized gene panel (not genome-wide) so the whole
pipeline runs on CPU in one Modal call. Real data throughout; no fabricated
labels or features. Trains the same fusion logistic regression already
implemented in data_pipeline/model.py, then pushes model.pkl to HF Hub.

Run: modal run --detach data_pipeline/train_real.py
"""
import modal

image = (
    modal.Image.debian_slim(python_version="3.11")
    .apt_install("mafft", "git")
    .pip_install(
        "biopython>=1.83", "numpy", "scikit-learn", "pandas", "requests",
        "huggingface_hub>=0.25.0",
    )
)

app = modal.App("vfi-real-train", image=image)
vol = modal.Volume.from_name("vfi-artifacts", create_if_missing=True)

REPO = "https://github.com/Anamitra-Sarkar/variant-functional-impact-ai.git"
HF_REPO_ID = "bhumika-tewari-282006/variant-functional-impact-ai"

# UniProt accession -> gene symbol, for a modest well-studied disease-gene panel
GENE_PANEL = {
    "TP53": "P04637", "BRCA1": "P38398", "BRCA2": "P51587", "PTEN": "P60484",
    "MLH1": "P40692", "MSH2": "P43246", "APC": "P25054", "VHL": "P40337",
    "RB1": "P06400", "KRAS": "P01116", "EGFR": "P00533", "LDLR": "P01130",
    "CFTR": "P13569", "G6PD": "P11413", "PAH": "P00439", "MYH7": "P12883",
    "FBN1": "P35555", "ATM": "Q13315", "RET": "P07949", "SMAD4": "Q13485",
}

ORTHOLOG_ORGANISMS = [9606, 10090, 10116, 9615, 9913, 9823, 7955, 8364, 9031]  # human,mouse,rat,dog,cow,pig,zebrafish,frog,chicken


@app.function(timeout=14400, cpu=4.0, memory=8192, volumes={"/art": vol}, secrets=[modal.Secret.from_name("bhumika-hf-token")])
def train_and_export() -> dict:
    import gzip
    import json
    import re
    import subprocess
    import sys
    from pathlib import Path

    import numpy as np
    import pandas as pd
    import requests

    subprocess.run(["git", "clone", "--depth", "1", REPO, "/repo"], check=True)
    sys.path.insert(0, "/repo")
    from data_pipeline.model import featurize_variant, train_logistic_regression, evaluate_with_baselines, save_artifacts, FEATURE_NAMES
    from data_pipeline.structure import compute_rsa_proxy, compute_ss_proxy
    from data_pipeline.conservation import conservation_scores

    RAW = Path("/art/raw")
    RAW.mkdir(parents=True, exist_ok=True)

    # ---- Step 1: ClinVar variant_summary (real labels) ----
    cv_path = RAW / "variant_summary.txt.gz"
    if not cv_path.exists():
        print("downloading ClinVar variant_summary.txt.gz", flush=True)
        url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
        with requests.get(url, stream=True, timeout=1800) as r:
            r.raise_for_status()
            with open(cv_path, "wb") as f:
                for chunk in r.iter_content(1 << 20):
                    f.write(chunk)
        print(f"  {cv_path} {cv_path.stat().st_size} bytes", flush=True)
    vol.commit()

    print("parsing ClinVar for our gene panel, missense, clear P/B calls", flush=True)
    protein_change_re = re.compile(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})")
    aa3to1 = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
        "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
        "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
        "Tyr": "Y", "Val": "V",
    }
    pathogenic_terms = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
    benign_terms = {"Benign", "Likely benign", "Benign/Likely benign"}

    rows = []
    genes_upper = set(GENE_PANEL.keys())
    with gzip.open(cv_path, "rt", encoding="utf-8", errors="ignore") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        need = ["GeneSymbol", "ClinicalSignificance", "Name", "Assembly", "Type"]
        for c in need:
            if c not in col:
                raise RuntimeError(f"ClinVar column missing: {c}; header={header[:20]}")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(col.values()):
                continue
            gene = parts[col["GeneSymbol"]]
            if gene not in genes_upper:
                continue
            if parts[col["Assembly"]] != "GRCh38":
                continue
            if parts[col["Type"]] != "single nucleotide variant":
                continue
            sig = parts[col["ClinicalSignificance"]]
            if sig in pathogenic_terms:
                label = 1
            elif sig in benign_terms:
                label = 0
            else:
                continue
            m = protein_change_re.search(parts[col["Name"]])
            if not m:
                continue
            wt3, pos, mut3 = m.group(1), int(m.group(2)), m.group(3)
            if wt3 not in aa3to1 or mut3 not in aa3to1:
                continue
            rows.append({"gene": gene, "pos": pos, "wt": aa3to1[wt3], "mut": aa3to1[mut3], "label": label})

    df = pd.DataFrame(rows).drop_duplicates(subset=["gene", "pos", "wt", "mut"])
    print(f"real ClinVar missense variants with clear P/B labels: {len(df)}", flush=True)
    print(df["gene"].value_counts().to_dict(), flush=True)
    print(f"label balance: {df['label'].value_counts().to_dict()}", flush=True)
    if len(df) < 30:
        raise RuntimeError(f"too few real labeled variants ({len(df)}) -- expand gene panel")
    df.to_csv("/art/clinvar_panel_labels.csv", index=False)
    vol.commit()

    # ---- Step 2: per-gene structure (AlphaFold) + conservation (orthologs + MAFFT) ----
    struct_dir = RAW / "structures"
    struct_dir.mkdir(exist_ok=True)
    msa_dir = RAW / "msa"
    msa_dir.mkdir(exist_ok=True)

    residue_features_by_gene = {}
    conservation_by_gene = {}

    for gene, uniprot in GENE_PANEL.items():
        if gene not in set(df["gene"]):
            continue
        pdb_path = struct_dir / f"{uniprot}.pdb"
        if not pdb_path.exists():
            api_r = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}", timeout=60)
            pdb_url = None
            if api_r.status_code == 200 and api_r.json():
                pdb_url = api_r.json()[0].get("pdbUrl")
            if not pdb_url:
                print(f"  [warn] no AlphaFold structure for {gene} ({uniprot}), skipping structure features", flush=True)
                continue
            r = requests.get(pdb_url, timeout=120)
            if r.status_code != 200:
                print(f"  [warn] AlphaFold fetch failed for {gene} ({uniprot}): {r.status_code}", flush=True)
                continue
            pdb_path.write_bytes(r.content)
        try:
            rsa_map = compute_rsa_proxy(str(pdb_path))
            ss_map = compute_ss_proxy(str(pdb_path))
            for k, v in ss_map.items():
                if k in rsa_map:
                    rsa_map[k].ss = v
            residue_features_by_gene[gene] = rsa_map
            print(f"  {gene}: {len(rsa_map)} residues with structure features", flush=True)
        except Exception as e:
            print(f"  [warn] structure featurization failed for {gene}: {e}", flush=True)

        msa_path = msa_dir / f"{gene}.fasta"
        if not msa_path.exists():
            try:
                seqs = []
                for org in ORTHOLOG_ORGANISMS:
                    q = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene}+AND+organism_id:{org}+AND+reviewed:true&format=fasta&size=1"
                    r = requests.get(q, timeout=60)
                    if r.status_code == 200 and r.text.strip():
                        seqs.append(r.text.strip())
                raw_fasta = msa_dir / f"{gene}_raw.fasta"
                raw_fasta.write_text("\n".join(seqs) + "\n")
                if len(seqs) >= 2:
                    aligned = subprocess.run(
                        ["mafft", "--quiet", "--auto", str(raw_fasta)],
                        capture_output=True, text=True, timeout=600,
                    )
                    msa_path.write_text(aligned.stdout)
                    print(f"  {gene}: aligned {len(seqs)} ortholog sequences", flush=True)
                else:
                    print(f"  [warn] too few orthologs for {gene} ({len(seqs)}), skipping conservation", flush=True)
            except Exception as e:
                print(f"  [warn] MSA build failed for {gene}: {e}", flush=True)

        if msa_path.exists() and msa_path.stat().st_size > 0:
            try:
                conservation_by_gene[gene] = conservation_scores(str(msa_path))
            except Exception as e:
                print(f"  [warn] conservation scoring failed for {gene}: {e}", flush=True)

    vol.commit()

    # ---- Step 3: featurize every labeled variant with whatever real features exist ----
    X_list, y_list, meta = [], [], []
    for _, row in df.iterrows():
        gene, pos, wt, mut, label = row["gene"], int(row["pos"]), row["wt"], row["mut"], int(row["label"])
        rsa_map = residue_features_by_gene.get(gene, {})
        residue_feature = None
        for (chain, resseq), rf in rsa_map.items():
            if resseq == pos:
                residue_feature = rf
                break
        cons_map = conservation_by_gene.get(gene, {})
        conservation = cons_map.get(pos, {}).get("conservation", 0.5) if isinstance(cons_map.get(pos), dict) else 0.5
        feats = featurize_variant(wt, mut, residue_feature, conservation)
        X_list.append(feats)
        y_list.append(label)
        meta.append({"gene": gene, "pos": pos, "wt": wt, "mut": mut})

    X = np.array(X_list)
    y = np.array(y_list)
    print(f"final real training set: X={X.shape}, positives={int(y.sum())}, negatives={int((1 - y).sum())}", flush=True)

    # ---- Step 4: train + honest evaluation with a held-out split ----
    from sklearn.model_selection import train_test_split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=42, stratify=y)
    model = train_logistic_regression(X_train, y_train)
    metrics = evaluate_with_baselines(X_test, y_test, FEATURE_NAMES, model)
    print(f"REAL_METRICS: {json.dumps(metrics, indent=2)}", flush=True)

    final_model = train_logistic_regression(X, y)  # refit on all real data for serving
    art_path = Path("/art/model.pkl")
    save_artifacts(final_model, FEATURE_NAMES, art_path)
    (Path("/art/eval.json")).write_text(json.dumps({
        "n_variants": len(df), "n_genes": int(df["gene"].nunique()),
        "held_out_metrics": metrics, "gene_panel": list(GENE_PANEL.keys()),
    }, indent=2))
    vol.commit()

    print("pushing model.pkl to HF Hub", flush=True)
    import os
    from huggingface_hub import HfApi
    hf_token = os.environ.get("HF_TOKEN") or os.environ.get("HUGGINGFACE_TOKEN") or os.environ.get("token")
    api = HfApi(token=hf_token)
    api.create_repo(repo_id=HF_REPO_ID, repo_type="model", exist_ok=True, private=True)
    sha = api.upload_file(path_or_fileobj=str(art_path), path_in_repo="model.pkl", repo_id=HF_REPO_ID, repo_type="model")
    print(f"PUSHED: {sha}", flush=True)

    return {"n_variants": len(df), "metrics": metrics}


@app.local_entrypoint()
def main():
    import json
    print(json.dumps(train_and_export.remote(), indent=2, default=str))
