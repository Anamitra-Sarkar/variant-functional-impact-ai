"""
Real training run for Variant Functional Impact AI, on Kaggle (CPU-only,
internet-enabled kernel). See train_real.py (Modal version) for full
rationale/citations -- identical logic, adapted to run standalone in
/kaggle/working with no GPU needed.

Expects: a Kaggle Dataset providing hf_token.txt (HF write token), attached
as kernel input, OR HF push is skipped and model.pkl is just saved to
/kaggle/working/model.pkl for later pickup.
"""
import glob
import gzip
import json
import os
import re
import subprocess
import sys
from pathlib import Path

WORK = Path("/kaggle/working")
os.chdir(WORK)

print("=== Step 0: system deps (mafft) ===", flush=True)
subprocess.run(["bash", "-lc", "apt-get update -qq && apt-get install -y -qq mafft"], check=True)
subprocess.run([sys.executable, "-m", "pip", "install", "-q", "huggingface_hub>=0.25.0"], check=True)

print("=== Step 1: clone repo ===", flush=True)
subprocess.run(["git", "clone", "--depth", "1",
                 "https://github.com/Anamitra-Sarkar/variant-functional-impact-ai.git", "/kaggle/working/repo"],
                check=True)
sys.path.insert(0, "/kaggle/working/repo")

from data_pipeline.model import featurize_variant, train_logistic_regression, evaluate_with_baselines, save_artifacts, FEATURE_NAMES  # noqa: E402
from data_pipeline.structure import compute_rsa_proxy, compute_ss_proxy  # noqa: E402
from data_pipeline.conservation import conservation_scores  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import requests  # noqa: E402

GENE_PANEL = {
    "TP53": "P04637", "BRCA1": "P38398", "BRCA2": "P51587", "PTEN": "P60484",
    "MLH1": "P40692", "MSH2": "P43246", "APC": "P25054", "VHL": "P40337",
    "RB1": "P06400", "KRAS": "P01116", "EGFR": "P00533", "LDLR": "P01130",
    "CFTR": "P13569", "G6PD": "P11413", "PAH": "P00439", "MYH7": "P12883",
    "FBN1": "P35555", "ATM": "Q13315", "RET": "P07949", "SMAD4": "Q13485",
}
ORTHOLOG_ORGANISMS = [9606, 10090, 10116, 9615, 9913, 9823, 7955, 8364, 9031]

RAW = WORK / "raw"
RAW.mkdir(exist_ok=True)

print("=== Step 2: ClinVar variant_summary (real labels) ===", flush=True)
cv_path = RAW / "variant_summary.txt.gz"
if not cv_path.exists():
    url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
    with requests.get(url, stream=True, timeout=1800) as r:
        r.raise_for_status()
        with open(cv_path, "wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
    print(f"  {cv_path} {cv_path.stat().st_size} bytes", flush=True)

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
    for c in ["GeneSymbol", "ClinicalSignificance", "Name", "Assembly", "Type"]:
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
df.to_csv(WORK / "clinvar_panel_labels.csv", index=False)

print("=== Step 3: per-gene structure + conservation ===", flush=True)
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
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb"
        r = requests.get(url, timeout=120)
        if r.status_code != 200:
            print(f"  [warn] no AlphaFold structure for {gene} ({uniprot})", flush=True)
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
                print(f"  [warn] too few orthologs for {gene} ({len(seqs)})", flush=True)
        except Exception as e:
            print(f"  [warn] MSA build failed for {gene}: {e}", flush=True)

    if msa_path.exists() and msa_path.stat().st_size > 0:
        try:
            conservation_by_gene[gene] = conservation_scores(str(msa_path))
        except Exception as e:
            print(f"  [warn] conservation scoring failed for {gene}: {e}", flush=True)

print("=== Step 4: featurize + train + evaluate ===", flush=True)
X_list, y_list = [], []
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

X = np.array(X_list)
y = np.array(y_list)
print(f"final real training set: X={X.shape}, positives={int(y.sum())}, negatives={int((1 - y).sum())}", flush=True)

from sklearn.model_selection import train_test_split  # noqa: E402
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=42, stratify=y)
model = train_logistic_regression(X_train, y_train)
metrics = evaluate_with_baselines(X_test, y_test, FEATURE_NAMES, model)
print(f"REAL_METRICS: {json.dumps(metrics, indent=2)}", flush=True)

final_model = train_logistic_regression(X, y)
art_path = WORK / "model.pkl"
save_artifacts(final_model, FEATURE_NAMES, art_path)
(WORK / "eval.json").write_text(json.dumps({
    "n_variants": len(df), "n_genes": int(df["gene"].nunique()),
    "held_out_metrics": metrics, "gene_panel": list(GENE_PANEL.keys()),
}, indent=2))
print("=== done -- model.pkl + eval.json in /kaggle/working ===", flush=True)

token_candidates = glob.glob("/kaggle/input/**/token.txt", recursive=True) + glob.glob("/kaggle/input/**/hf_token.txt", recursive=True)
if token_candidates:
    print("=== Step 5: push model.pkl to HF Hub ===", flush=True)
    from huggingface_hub import HfApi
    hf_token = open(token_candidates[0]).read().strip()
    api = HfApi(token=hf_token)
    REPO_ID = "bhumika-tewari-282006/variant-functional-impact-ai"
    api.create_repo(repo_id=REPO_ID, repo_type="model", exist_ok=True, private=True)
    sha = api.upload_file(path_or_fileobj=str(art_path), path_in_repo="model.pkl", repo_id=REPO_ID, repo_type="model")
    print(f"PUSHED: {sha}", flush=True)
else:
    print("no HF token dataset attached; skipping HF push (pick up model.pkl from kernel Output instead)", flush=True)
