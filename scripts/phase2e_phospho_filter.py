"""
Aurora B Substrate Mapping — Phase 2e: TNBC Phosphoproteomic Filter
------------------------------------------------------------------
Phase 2d's reports.
  - Phospho columns are a MultiIndex: (Name=gene symbol, Site, Peptide,
    Database_ID=Ensembl protein ID). 
  - 'Site' can pack multiple residues into one string for multiply-
    phosphorylated peptides (e.g. "S19S22"), needs splitting into
    individual sites before it's comparable to Phase 1's per-residue
    candidates.
  - Sample IDs already match tnbc_status_by_patient.csv's patient_id
    format directly (no suffix stripping needed)  152/158 matched,
    the rest lack a clinical TNBC call and are correctly excluded.
  - No tumor/normal suffix was found, so no extra tissue filtering step
    is needed beyond the TNBC-positive split itself.

This script:
  1. Restricts the phosphoproteomics table to TNBC-positive tumor samples.
  2. For every site column, computes what fraction of those TNBC samples
     actually have a detected (non-null) value there.
  3. Explodes multi-residue Site strings into individual (gene, residue,
     position) sites.
  4. Matches Phase 1's PSSM candidates (UniProt ID + position) against
     the observed sites via a UniProt->gene symbol map built from the
     FASTA headers.


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
OUTPUT_CSV = DATA_DIR / "aurkb_candidates_tnbc_filtered.csv"

# Minimum fraction of TNBC tumor samples a site must be detected
# (non-null) in to count as "actually observed phosphorylated in TNBC
# tissue." The sweep printed below shows the tradeoff before this gets
# used — adjust and rerun if the default doesn't fit.
MIN_TNBC_DETECTION_FRAC = 0.2

SITE_TOKEN_RE = re.compile(r"[ST](\d+)")  # Aurora B only phosphorylates S/T


def load_phosphoproteomics() -> pd.DataFrame:
    print("Loading CPTAC BRCA phosphoproteomics (source=umich)...")
    brca = cptac.Brca()
    phospho = brca.get_phosphoproteomics(source="umich")
    print(f"Loaded phosphoproteomics: {phospho.shape}")
    return phospho


def get_tnbc_tumor_samples(phospho: pd.DataFrame, tnbc_status: pd.DataFrame) -> list:
    tnbc_ids = set(tnbc_status.loc[tnbc_status["is_tnbc"], "patient_id"])
    samples = [i for i in phospho.index if i in tnbc_ids]
    print(f"TNBC-positive samples present in phospho table: {len(samples)} / {len(phospho.index)}")
    return samples


def compute_detection_fractions(phospho: pd.DataFrame, tnbc_samples: list) -> pd.Series:
    phospho_tnbc = phospho.loc[tnbc_samples]
    n = len(tnbc_samples)
    detection_counts = phospho_tnbc.notna().sum(axis=0)
    return detection_counts / n


def explode_sites(detection_frac: pd.Series) -> pd.DataFrame:
    """Turn (Name, Site, Peptide, Database_ID)-indexed fractions into one
    row per individual (gene, residue, position), keeping the max
    detection fraction across any column reporting that exact site."""
    rows = []
    for (name, site, peptide, db_id), frac in detection_frac.items():
        for m in SITE_TOKEN_RE.finditer(site):
            residue = site[m.start()]
            position = int(m.group(1))
            rows.append({"gene": name, "residue": residue, "position": position, "detection_frac": frac})
    exploded = pd.DataFrame(rows)
    exploded = exploded.groupby(["gene", "residue", "position"], as_index=False)["detection_frac"].max()
    print(f"Exploded to {len(exploded)} unique (gene, residue, position) S/T sites.")
    return exploded


def load_gene_symbol_map(fasta_path: Path) -> dict:
    print("Building uniprot_id -> gene_symbol map from proteome FASTA...")
    mapping = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        uniprot_id = record.id.split("|")[1] if "|" in record.id else record.id
        m = re.search(r"GN=(\S+)", record.description)
        if m:
            mapping[uniprot_id] = m.group(1)
    return mapping


def main():
    DATA_DIR.mkdir(exist_ok=True)

    phospho = load_phosphoproteomics()
    tnbc_status = pd.read_csv(TNBC_STATUS_CSV)
    tnbc_samples = get_tnbc_tumor_samples(phospho, tnbc_status)
    if not tnbc_samples:
        raise RuntimeError("No TNBC-positive samples matched — check tnbc_status_by_patient.csv.")

    detection_frac = compute_detection_fractions(phospho, tnbc_samples)
    observed = explode_sites(detection_frac)

    print("\nDetection-fraction sweep (min frac | sites passing | as % of exploded sites):")
    for frac in (0.0, 0.1, 0.2, 0.3, 0.5):
        n_pass = (observed["detection_frac"] >= frac).sum()
        marker = "  <- selected" if frac == MIN_TNBC_DETECTION_FRAC else ""
        print(f"  >= {frac:.1f} | {n_pass:>7,} | {n_pass/len(observed):.1%}{marker}")

    observed_pass = observed[observed["detection_frac"] >= MIN_TNBC_DETECTION_FRAC]
    residues_by_gene_pos = {}
    for _, row in observed_pass.iterrows():
        residues_by_gene_pos.setdefault((row["gene"], row["position"]), set()).add(row["residue"])

    gene_map = load_gene_symbol_map(PROTEOME_FASTA)
    candidates = pd.read_csv(PSSM_CANDIDATES_CSV)
    candidates["gene_symbol"] = candidates["uniprot_id"].map(gene_map)

    def check_site(row):
        residues_here = residues_by_gene_pos.get((row["gene_symbol"], row["position"]))
        if not residues_here:
            return pd.Series({"observed_in_tnbc": False, "position_residue_mismatch": False})
        if row["residue"] in residues_here:
            return pd.Series({"observed_in_tnbc": True, "position_residue_mismatch": False})
        return pd.Series({"observed_in_tnbc": False, "position_residue_mismatch": True})

    checked = candidates.apply(check_site, axis=1)
    candidates = pd.concat([candidates, checked], axis=1)
    candidates = candidates.sort_values("pssm_score", ascending=False)
    candidates.to_csv(OUTPUT_CSV, index=False)

    n_observed = int(candidates["observed_in_tnbc"].sum())
    n_mismatch = int(candidates["position_residue_mismatch"].sum())
    n_unmapped_gene = int(candidates["gene_symbol"].isna().sum())

    print(f"\nSaved {len(candidates)} candidates (all, with flags) to {OUTPUT_CSV}")
    print(f"  -> {n_observed} ARE observed phosphorylated in TNBC tissue "
          f"(>= {MIN_TNBC_DETECTION_FRAC:.0%} of {len(tnbc_samples)} TNBC samples)")
    print(f"  -> {n_mismatch} had a gene+position hit but a residue-letter mismatch — "
          f"likely UniProt/Ensembl numbering conflicts, worth spot-checking a few manually")
    print(f"  -> {n_unmapped_gene} had no gene symbol at all and couldn't be checked")
    print(f"\nYour real Phase 2 shortlist: the {n_observed} candidates with "
          f"observed_in_tnbc == True, already sorted by pssm_score.")


if __name__ == "__main__":
    main()
