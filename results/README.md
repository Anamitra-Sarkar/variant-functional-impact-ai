# Real data run — status (2026-09-04)

## Real sources confirmed reachable and correct

- Structure: `AlphaFold DB` `https://alphafold.ebi.ac.uk/files/AF-P04637-F1-model_v4.pdb`
  (human p53, matching the repo's own worked example) — real, verified 200.
- Labels: NCBI ClinVar bulk `variant_summary.txt.gz` (real, 442MB, filterable
  to `GeneSymbol == TP53` for real pathogenic/benign missense labels via the
  `Name` column's `p.ArgXXXHis`-style HGVS protein notation).
- A real orthologous-sequence alignment source exists too: Pfam's curated
  seed alignment for the P53 DNA-binding domain (`PF00870`), fetched via
  InterPro's REST API (`.../entry/pfam/PF00870/?annotation=alignment:seed`,
  real, verified 200, Stockholm format).

## Real blocker found — not attempted, to avoid a silently-wrong result

`data_pipeline/conservation.py`'s `conservation_scores()` treats **MSA
column index == protein residue position directly** (1-indexed, no
handling for internal gaps in the reference sequence, and `cli.py` calls
it with no offset). That is a correct assumption only if the MSA's human
sequence has **zero internal gap characters relative to its own full-length
numbering** end to end.

The only real alignment source found (Pfam `PF00870`) covers only the p53
DNA-binding domain (roughly residues 94-312 of the 393-residue protein), not
the full protein, and — like any real multiple sequence alignment — it may
insert gap columns into the human row itself wherever other species have
insertions relative to human. Using it directly, or naively splicing it into
a "full-length" file, would silently shift ClinVar variant positions
(e.g. R175H, R248Q, R273H) onto the wrong alignment column with no error
raised — the code has no way to detect a wrong-but-plausible offset, and
neither did checking here.

Building a genuinely correct full-length ortholog MSA for TP53 (real
sequences, real alignment, and a runtime check that the human row is
gap-free at every ClinVar-labelled position before trusting it) is real,
necessary work that was not attempted rather than risk exactly the kind of
plausible-but-wrong result this task explicitly warned against.

## What would close this out

1. Fetch real full-length TP53 ortholog sequences (UniProt REST, e.g. mouse
   `P02340`, rat, a handful of mammals — reachable, confirmed above for the
   human sequence itself).
2. Run a real alignment tool (e.g. MAFFT) with the human sequence pinned as
   the reference row (`--keeplength` or equivalent), so alignment columns
   stay 1:1 with full-length human numbering by construction.
3. Verify programmatically that the human row has no `-`/`.` characters
   before treating any ClinVar position as reachable, rather than trusting
   `chain_offset` blindly.
4. Then run `data_pipeline.cli --structure-path AF-P04637.pdb --msa-path
   <verified MSA> --labels <ClinVar-derived CSV> --train`.
