"""
Aurora B Substrate Mapping — Phase 3a: ColabFold Input Prep
------------------------------------------------------------------

Steps:
  1. Fetch AURKB's canonical sequence + domain annotation from UniProt
     live (no hardcoded residue numbers — printed so you can sanity-check).
  2. Load the top N orphan, TNBC-observed candidates from Phase 2e.
  3. Pull a wider peptide window (21-mer, site-centered) directly from
     the proteome FASTA — Phase 1's 10-mer window was sized for motif
     scoring, not structural docking.
  4. Write one multimer FASTA (kinase:peptide) per candidate, plus a
     manifest CSV, ready to zip and upload to the ColabFold batch notebook.
"""

from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO

DATA_DIR = Path("data")
CANDIDATES_CSV = DATA_DIR / "aurkb_candidates_tnbc_filtered.csv"
PROTEOME_FASTA = DATA_DIR / "human_reviewed_proteome.fasta"
OUT_DIR = Path("phase3_colabfold_inputs")
MANIFEST_CSV = OUT_DIR / "manifest.csv"

TOP_N = 25
PEPTIDE_FLANK = 10  # residues each side of the site -> 21-mer total


def fetch_aurkb_kinase_domain() -> str:
    print("Fetching AURKB (human) entry from UniProt...")
    resp = requests.get(
        "https://rest.uniprot.org/uniprotkb/search",
        params={"query": "gene:AURKB AND organism_id:9606 AND reviewed:true", "format": "json"},
        timeout=60,
    )
    resp.raise_for_status()
    results = resp.json()["results"]
    if not results:
        raise RuntimeError("No reviewed human AURKB entry found on UniProt.")
    entry = results[0]
    accession = entry["primaryAccession"]
    full_seq = entry["sequence"]["value"]
    print(f"Using UniProt accession {accession}, {len(full_seq)} aa.")

    kinase_domain = None
    for feature in entry.get("features", []):
        if feature.get("type") == "Domain" and "protein kinase" in feature.get("description", "").lower():
            start = feature["location"]["start"]["value"]
            end = feature["location"]["end"]["value"]
            kinase_domain = full_seq[start - 1:end]
            print(f"Found kinase domain annotation: residues {start}-{end} ({len(kinase_domain)} aa).")
            break

    if kinase_domain is None:
        print("WARNING: no 'Protein kinase' domain feature found — using the full-length "
              "sequence instead. Check the UniProt entry manually (accession above) before trusting this.")
        kinase_domain = full_seq

    return kinase_domain


def load_proteome(fasta_path: Path) -> dict:
    print("Loading proteome for peptide extraction...")
    seqs = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        uniprot_id = record.id.split("|")[1] if "|" in record.id else record.id
        seqs[uniprot_id] = str(record.seq)
    return seqs


def extract_peptide(seq: str, position_1based: int, flank: int = PEPTIDE_FLANK) -> str:
    idx = position_1based - 1
    start = max(0, idx - flank)
    end = min(len(seq), idx + flank + 1)
    return seq[start:end]


def main():
    OUT_DIR.mkdir(exist_ok=True)

    kinase_domain = fetch_aurkb_kinase_domain()
    proteome = load_proteome(PROTEOME_FASTA)

    candidates = pd.read_csv(CANDIDATES_CSV)
    shortlist = (
        candidates[(~candidates["known_in_omnipath"]) & (candidates["observed_in_tnbc"])]
        .sort_values("pssm_score", ascending=False)
        .head(TOP_N)
        .copy()
    )
    print(f"\nBuilding ColabFold inputs for top {len(shortlist)} orphan, TNBC-observed candidates...")

    manifest_rows = []
    n_skipped = 0
    for _, row in shortlist.iterrows():
        seq = proteome.get(row["uniprot_id"])
        if seq is None or len(seq) < 5:
            n_skipped += 1
            continue
        peptide = extract_peptide(seq, int(row["position"]))
        if len(peptide) < 5:
            n_skipped += 1
            continue

        gene = row["gene_symbol"] if pd.notna(row["gene_symbol"]) else row["uniprot_id"]
        job_id = f"{gene}_{row['residue']}{int(row['position'])}"

        fasta_path = OUT_DIR / f"{job_id}.fasta"
        fasta_path.write_text(f">{job_id}\n{kinase_domain}:{peptide}\n")

        manifest_rows.append({
            "job_id": job_id,
            "gene_symbol": gene,
            "uniprot_id": row["uniprot_id"],
            "position": int(row["position"]),
            "residue": row["residue"],
            "pssm_score": row["pssm_score"],
            "peptide": peptide,
            "peptide_site_index_0based": min(PEPTIDE_FLANK, int(row["position"]) - 1),
            "fasta_file": fasta_path.name,
        })

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(MANIFEST_CSV, index=False)

    print(f"\nWrote {len(manifest_rows)} per-candidate FASTA files to {OUT_DIR.resolve()}")
    if n_skipped:
        print(f"Skipped {n_skipped} candidates (protein not found or peptide too short at a sequence boundary).")
    print(f"Manifest saved to {MANIFEST_CSV}")
    print("\nNext: zip the phase3_colabfold_inputs/ folder, upload it to the ColabFold "
          "batch notebook, and run. Bring back the output PDBs + confidence JSONs and "
          "I'll help you evaluate physical plausibility.")


if __name__ == "__main__":
    main()
