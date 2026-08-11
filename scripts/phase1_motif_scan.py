"""
Aurora B Substrate Mapping — Phase 1: Motif Scanning
------------------------------------------------------
Scans the human reviewed proteome for the Aurora B kinase consensus
phosphorylation motif, and cross-references hits against known AURKB
substrates pulled live from OmniPath to flag which candidates are
already known vs. "orphan" leads worth investigating further.


import re
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO
from omnipath.requests import Enzsub

# ---------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------
DATA_DIR = Path("data")
DATA_DIR.mkdir(exist_ok=True)

PROTEOME_FASTA = DATA_DIR / "human_reviewed_proteome.fasta"
OUTPUT_CSV = DATA_DIR / "aurkb_candidate_sites.csv"


# ---------------------------------------------------------------
# Download the human reviewed (Swiss-Prot) proteome
# ---------------------------------------------------------------
def download_human_proteome(output_path: Path) -> None:
    if output_path.exists():
        print(f"Proteome already downloaded at {output_path}, skipping.")
        return

    print("Downloading reviewed human proteome from UniProt...")
    url = (
        "https://rest.uniprot.org/uniprotkb/stream"
        "?query=organism_id:9606+AND+reviewed:true"
        "&format=fasta"
    )
    response = requests.get(url, timeout=180)
    response.raise_for_status()
    output_path.write_text(response.text)
    print(f"Saved proteome to {output_path}")


# ---------------------------------------------------------------
# Pull known AURKB substrates live from OmniPath
# ---------------------------------------------------------------
def load_known_aurkb_sites() -> pd.DataFrame:
    """
    Queries OmniPath's enzyme-substrate (Enzsub) endpoint for every
    known AURKB (human) substrate site. This replaces the old manual
    PhosphoSitePlus download — no account or license needed.
    """
    print("Querying OmniPath for known AURKB substrate sites...")
    df = Enzsub.get(enzymes="AURKB", genesymbols=True, organisms="human")
    print("Columns returned by OmniPath:", list(df.columns))

    # Column names occasionally shift slightly between omnipath versions
    # if this KeyErrors, print(df.head()) and adjust the names below to match.
    
    known = df[["substrate", "substrate_genesymbol", "residue_offset", "residue_type"]].copy()
    known.columns = ["uniprot_id", "substrate_name", "position", "residue_type"]
    known = known.dropna(subset=["position"])
    known["position"] = known["position"].astype(float)

    print(f"Found {len(known)} known AURKB substrate sites via OmniPath.")
    return known


# ---------------------------------------------------------------
# Define the Aurora B consensus motif
# ---------------------------------------------------------------
# Literature consensus: [R/K]-x-[S/T]-[hydrophobic]
#   position -2 : R or K
#   position -1 : any residue
#   position  0 : S or T   <- the phosphorylated residue
#   position +1 : hydrophobic (I, L, V, F, M, A)
#
# This is a simplification. Step 5 below checks
# how many of the KNOWN sites this regex actually catches, if that
# recovery rate is low, tighten/loosen the pattern before trusting it
# on the full proteome scan.
AURKB_MOTIF = re.compile(r"[RK].[ST][ILVFMA]")


# ---------------------------------------------------------------
# Scan the proteome for motif matches
# ---------------------------------------------------------------
def scan_proteome_for_motif(fasta_path: Path) -> pd.DataFrame:
    hits = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        seq = str(record.seq)
        uniprot_id = record.id.split("|")[1] if "|" in record.id else record.id
        protein_name = record.description

        for match in AURKB_MOTIF.finditer(seq):
            # match.start() is the 0-based index of the R/K residue.
            # The phospho-residue (S/T) sits 2 positions later, and
            # biological residue numbering is 1-based -> +3 total.
            phospho_position = match.start() + 3
            hits.append({
                "uniprot_id": uniprot_id,
                "protein_name": protein_name,
                "position": phospho_position,
                "matched_motif": match.group(),
                "phospho_residue": match.group()[2],
            })

    hits_df = pd.DataFrame(hits)
    print(f"Found {len(hits_df)} raw motif matches across the proteome.")
    return hits_df


# ---------------------------------------------------------------
# Cross-reference against known AURKB sites
# ---------------------------------------------------------------
def flag_known_sites(hits_df: pd.DataFrame, known_df: pd.DataFrame) -> pd.DataFrame:
    known_pairs = set(zip(known_df["uniprot_id"], known_df["position"]))

    def is_known(row):
        return (row["uniprot_id"], float(row["position"])) in known_pairs

    hits_df["known_in_omnipath"] = hits_df.apply(is_known, axis=1)
    return hits_df


# ---------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------
def main():
    download_human_proteome(PROTEOME_FASTA)
    known_sites = load_known_aurkb_sites()

    hits_df = scan_proteome_for_motif(PROTEOME_FASTA)
    hits_df = flag_known_sites(hits_df, known_sites)

    hits_df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSaved {len(hits_df)} candidate sites to {OUTPUT_CSV}")

    n_known = int(hits_df["known_in_omnipath"].sum())
    n_orphan = len(hits_df) - n_known
    recovery_rate = n_known / len(known_sites) if len(known_sites) else 0

    print(f"  -> {n_known} matches overlap already-known AURKB sites")
    print(f"  -> Recovery rate: {recovery_rate:.1%} of known sites caught by this motif")
    print(f"  -> {n_orphan} are 'orphan' candidates worth a closer look")

    if recovery_rate < 0.5:
        print(
            "\n  Recovery rate is under 50% - the regex motif may be too strict. "
            "Consider loosening it (e.g. widening the hydrophobic residue set) "
            "before treating the orphan list as meaningful."
        )


if __name__ == "__main__":
    main()
