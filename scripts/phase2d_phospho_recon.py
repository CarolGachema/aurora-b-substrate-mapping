"""
Aurora B Substrate Mapping — Phase 2d: Phosphoproteomics Recon
------------------------------------------------------------------
Two things are still unknown:

  1. ID format: Phase 1's candidates are UniProt ID + position.
     CPTAC phosphoproteomics tables are usually indexed by GENE SYMBOL
     (sometimes gene + site as a MultiIndex).

  2. Sample composition: CPTAC tables often include BOTH tumor and
     adjacent-normal columns per patient. We only want tumor samples,
     and only the ones flagged TNBC-positive in tnbc_status_by_patient.csv.

This script only inspects.
"""

import re
from pathlib import Path

import cptac
import pandas as pd
from Bio import SeqIO

DATA_DIR = Path("data")
PSSM_CANDIDATES_CSV = DATA_DIR / "aurkb_candidate_sites_pssm.csv"
TNBC_STATUS_CSV = DATA_DIR / "tnbc_status_by_patient.csv"
PROTEOME_FASTA = DATA_DIR / "human_reviewed_proteome.fasta"



#  Load phosphoproteomics 

def load_phosphoproteomics():
    print("Loading CPTAC BRCA phosphoproteomics...")
    brca = cptac.Brca()
    for src in ("umich", "bcm", "broad", "washu"):
        try:
            phospho = brca.get_phosphoproteomics(source=src)
            print(f"Loaded phosphoproteomics from source: {src}")
            return phospho, src
        except Exception as e:
            print(f"  Source '{src}' failed ({type(e).__name__}), trying next...")
    raise RuntimeError("No phosphoproteomics source worked — paste this output and we'll adjust.")



#  Inspect column structure (gene symbol? multi-index? peptide?)

def inspect_columns(phospho: pd.DataFrame):
    print("\n" + "=" * 60)
    print("COLUMN STRUCTURE")
    print("=" * 60)
    print("Column index type:", type(phospho.columns))
    if isinstance(phospho.columns, pd.MultiIndex):
        print("Column level names:", phospho.columns.names)
        for level in range(phospho.columns.nlevels):
            print(f"\nLevel {level} sample values:")
            print(phospho.columns.get_level_values(level)[:10].tolist())
    else:
        print("First 20 columns:")
        print(phospho.columns[:20].tolist())



# Inspect sample/patient ID overlap with the TNBC list
def inspect_samples(phospho: pd.DataFrame, tnbc_status: pd.DataFrame):
    print("\n" + "=" * 60)
    print("SAMPLE / PATIENT ID OVERLAP")
    print("=" * 60)
    print("First 10 row index values (samples):", phospho.index[:10].tolist())
    print(f"Total samples in phospho table: {len(phospho.index)}")

    
    suffix_like = [i for i in phospho.index if not str(i)[-1].isdigit()]
    print(f"\nSample IDs with a non-numeric ending (possible normal-tissue "
          f"or replicate flag): {len(suffix_like)}")
    if suffix_like:
        print("Examples:", suffix_like[:10])

    tnbc_ids = set(tnbc_status["patient_id"])
    matched = [i for i in phospho.index
               if i in tnbc_ids or str(i).split(".")[0] in tnbc_ids]
    unmatched = [i for i in phospho.index if i not in matched]
    print(f"\nPhospho sample IDs matching a known TNBC-status patient_id: "
          f"{len(matched)} / {len(phospho.index)}")
    if unmatched:
        print("Examples of non-matching phospho sample IDs:", unmatched[:10])



# Build uniprot_id -> gene_symbol map from the FASTA headers

def load_gene_symbol_map(fasta_path: Path) -> dict:
    print("\nBuilding uniprot_id -> gene_symbol map from proteome FASTA (GN= field)...")
    mapping = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        uniprot_id = record.id.split("|")[1] if "|" in record.id else record.id
        m = re.search(r"GN=(\S+)", record.description)
        if m:
            mapping[uniprot_id] = m.group(1)
    print(f"Mapped {len(mapping)} UniProt IDs to gene symbols "
          f"(some entries lack a GN= field and will be unmapped).")
    return mapping


def inspect_candidate_overlap(phospho: pd.DataFrame, candidates: pd.DataFrame, gene_map: dict):
    print("\n" + "=" * 60)
    print("CANDIDATE SITE <-> PHOSPHO COLUMN OVERLAP (rough check)")
    print("=" * 60)
    candidates = candidates.copy()
    candidates["gene_symbol"] = candidates["uniprot_id"].map(gene_map)
    n_unmapped = candidates["gene_symbol"].isna().sum()
    print(f"Candidates with no gene symbol resolved: {n_unmapped} / {len(candidates)}")

    if isinstance(phospho.columns, pd.MultiIndex):
        col_genes = set(phospho.columns.get_level_values(0))
    else:
        col_genes = set(str(c).split("-")[0].split("_")[0] for c in phospho.columns)

    candidate_genes = set(candidates["gene_symbol"].dropna())
    gene_overlap = candidate_genes & col_genes
    print(f"\nCandidate gene symbols found among phospho columns: "
          f"{len(gene_overlap)} / {len(candidate_genes)}")
    print("Sample overlap:", list(gene_overlap)[:10])
    print("Sample phospho column 'genes' for comparison:", list(col_genes)[:10])


def main():
    phospho, source = load_phosphoproteomics()
    print("\nPhosphoproteomics table shape:", phospho.shape)

    inspect_columns(phospho)

    tnbc_status = pd.read_csv(TNBC_STATUS_CSV)
    inspect_samples(phospho, tnbc_status)

    candidates = pd.read_csv(PSSM_CANDIDATES_CSV)
    gene_map = load_gene_symbol_map(PROTEOME_FASTA)
    inspect_candidate_overlap(phospho, candidates, gene_map)

    print("\n" + "=" * 60)
    print("Paste this whole output back and I'll write the real Phase 2d filter.")
    print("=" * 60)


if __name__ == "__main__":
    main()
