# Architecture

## Overview
Structure/evolution-aware AI for per-variant missense functional impact prediction.
Single-variant scoring (in the spirit of SIFT, PolyPhen-2, AlphaMissense), not
whole-exome multi-gene prioritization.

## Components

### 1. `data_pipeline/` — Feature extraction + fusion model

```
data_pipeline/
  __init__.py
  physicochemical.py   # Kyte-Doolittle hydrophobicity + vdW volumes (real values)
  structure.py         # Bio.PDB parsing, RSA proxy, secondary-structure proxy, stability delta
  conservation.py      # MSA Shannon-entropy / information-content conservation scoring
  model.py             # LogisticRegression fusion model (structural + conservation + physicochemical)
  evaluate.py          # AUROC/AUPRC vs baselines
  cli.py               # CLI: --structure-path, --msa-path, --train, --predict
```

**Structural features (per mutated residue):**
- **RSA** — Relative Solvent Accessibility. Production path: DSSP via `Bio.PDB.DSSP`
  (`mkdssp` binary) on real PDB/AlphaFold structures. Sandbox proxy: neighbor-count
  burial estimate (counts CA neighbors within 10 Å, maps to 0..1 burial fraction,
  then RSA = 1 - burial). Real geometric computation, documented DSSP alternative.
  Parser hardened for real-file quirks: empty file, hetero-only, missing CA, insertion codes, multiple models, header/comment variants, optional PDB columns (B-factor/occupancy) tolerance.
- **Secondary structure** — DSSP 8-state → 3-state (H helix, E sheet, C loop) when
  DSSP available; proxy uses backbone geometry heuristics (CA-CA-CA angle / hydrogen-bond
  proxy is not attempted; we use a documented simplified rule: local CA distance pattern,
  clearly labeled as proxy). Proxy hardened for single-residue chains and multi-chain grouping.
- **Stability proxy** — Δvolume, Δhydrophobicity, Grantham-like distance between WT and mutant
  residues using real per-AA tables. No learned energy function; honest physicochemical delta.

**Conservation:**
- Real Shannon-entropy and information-content scoring from MSA.
- Input: FASTA MSA file (orthologs/paralogs). Per-column entropy:
  `H = -Σ p_a * log2(p_a)`, conservation = `1 - H / log2(20)` or `IC = log2(20) - H`.
  Implemented for real given an input alignment file. Parser hardened for real MSA quirks: comment lines (`;`, `#`), leading-whitespace headers, empty headers (auto-named), spaced alignments, lower-case sequences, blank lines, empty-sequence detection, differing-length validation with header detail, ambiguous AA (`B/Z/X/*`) ignored, both gap chars (`-`, `.`).

**Fusion model:**
- LogisticRegression (L2) or GradientBoosting (sklearn). Features: RSA, SS one-hot,
  conservation, Δvolume, Δhydrophobicity, Grantham distance, BLOSUM62 penalty.
  Outputs continuous damaging score in [0,1].

### 2. `backend/` — FastAPI

```
backend/
  __init__.py
  app.py        # FastAPI app factory (lifespan, CORS, fail-closed gate)
  model_store.py # Fail-closed release gate
  auth.py       # Firebase-auth-shaped bearer token stub
  schemas.py    # Pydantic request/response models with strict validators
```

**Release gate (fail-closed):**
- Model artifacts are NOT loaded/served unless BOTH env vars are set:
  `MODEL_RELEASE_APPROVED=true` AND `APPROVED_ARTIFACT_REVISION=<non-empty>`.
- `GET /health` always returns 200 with `model_loaded: bool`.
- `GET /readiness` returns 200 if loaded, 503 if not (JSON body still honest).
- `GET /model/info` returns `model_loaded`, `model_revision`, `feature_names` (when loaded).
- `GET /` returns service identity + `model_loaded` flag.
- `POST /predict` returns 503 with honest "model not released" if gate closed.
- `POST /predict` validates input strictly via Pydantic: `wt`/`mut` must be canonical 20 AAs (422), `wt != mut` (422), `chain` single alphanumeric (422), `protein` alphanumeric/underscore/dash (422), `position` 1..100000 (422). Malformed JSON returns 422, synonymous `wt==mut` returns 422, invalid AA returns 422 (never 500).
- CORS middleware enabled for frontend dev origins; startup uses `lifespan` (not deprecated `on_event`).
- Frontend shows abstention banner when backend reports `model_loaded=false`.

**Auth stub:**
- Reads `FIREBASE_SERVICE_ACCOUNT_JSON` env var (path to JSON). If absent, auth is
  permissive in dev (documented). If present, verifies Bearer token via
  `google.oauth2.id_token` or mocked verifier. Unit-tested with mocked verifier.
- Dependency `get_current_user` raises 401 on invalid token when strict mode.

### 3. `frontend/` — React + Vite + TypeScript

```
frontend/
  src/
    App.tsx        # Main dashboard: form, health banner, score display (accessible)
    api.ts         # fetchHealth / predictVariant (handles 422 detail arrays + 503)
    main.tsx
  index.html
  vite.config.ts
  package.json
  tsconfig.json
```

- Variant input: protein (UniProt ID / gene), position (1-indexed), WT, mutant, chain.
- All form controls have associated `<label htmlFor>` + `id`, `aria-label`, `aria-required`, `aria-live` regions for errors/results, `role="alert"` on banners, `role="progressbar"` for score bar, skip-link for keyboard navigation.
- Responsive: `flex-wrap` + `@media (max-width: 640px)` stacks rows vertically; works at 320px narrow widths.
- Color contrast hardened: health/banner/error text uses darker ratios (e.g. `#166534`, `#78350f`, `#7f1d1d`) on light backgrounds.
- Displays: damaging score + per-feature explanation (RSA, conservation, Δphyschem) with `aria-label` on score circle and `aria-valuenow` on progress bar.
- 503 abstention handling: frontend detects 503 / "not yet released" message and shows dedicated abstention error rather than generic failure.
- Honest banner (`role="alert"`) when `model_loaded=false`.
- Clean scientific-dashboard design (no generic boilerplate).

### 4. `tests/` — pytest

- Synthetic tiny PDB fixture (hand-written, hand-verifiable geometry) covering multiple header variants, missing optional columns, comment lines, insertion codes, empty file, hetero-only edge cases.
- Synthetic FASTA MSA fixture (conserved vs variable columns) covering comment lines (`;`, `#`), leading-whitespace headers, empty headers, spaced/ lower-case sequences, differing-length error, empty-file error, ambiguous AA (`B/Z/X/*`) handling, gap chars (`-` and `.`), multi-line sequences.
- Checks: RSA correctness, SS correctness (including single-residue and multi-chain), physicochemical table correctness (including Grantham ordering), conservation scoring (conserved > variable, gap/ambiguous handling), fusion model recovers injected signal, backend release gate (health/readiness/503, CORS, lifespan, 422/400 validation, corrupted artifact handling, `/model/info`), auth stub, CLI variant parsing edge cases.

## Data flow

```
PDB/AlphaFold .pdb/.cif  ──┐
                            ├─► structure.py ─┐
MSA FASTA  ─────────────────┘                  ├─► model.py ─► score [0,1]
Physicochemical tables ────────────────────────┘
```

Training: synthetic labeled set with INJECTED real signal
(buried+conserved+large Δ → damaging) or real ClinVar if available offline.

## Evaluation

- Metrics: AUROC, AUPRC (sklearn).
- Baselines: conservation-only, hydrophobicity-change-only.
- Fusion model must beat baselines on synthetic injected-signal set (verified in tests).

## Real-run procedure (Kaggle/Modal)

```bash
# 1. Download structure (example: AlphaFold DB per UniProt ID)
curl -o AF-P53.pdb https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v4.pdb

# 2. Download MSA (example: UniProt orthologs via Clustal/MAFFT, or ConSurf)
# 3. (Optional) Install DSSP: conda install -c salilab dssp  or  apt-get install dssp

# 4. Train / evaluate
python -m data_pipeline.cli --structure-path AF-P53.pdb --msa-path orthologs.fasta --train --labels clinvar_labels.csv

# 5. Predict single variant
python -m data_pipeline.cli --structure-path AF-P53.pdb --msa-path orthologs.fasta --predict --variant P53_R175H
```

## Security & correctness

- No fabricated datasets/APIs/metrics. All sources cited in `docs/data_sources.md`.
- No heavy downloads in sandbox; tests use synthetic fixtures.
- Backend never fabricates loaded state.
