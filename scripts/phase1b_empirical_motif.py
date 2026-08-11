"""
Aurora B Substrate Mapping 
Phase 1b: Empirical Motif Building
------------------------------------------------------------------

This phase builds the motif FROM DATA instead of
guessing it: it pulls the real sequence window around each of the
281 known AURKB sites (from OmniPath), and uses those to build a
Position-Specific Scoring Matrix (PSSM), a per-position log-odds
score for every amino acid, relative to how common that amino acid
is in the proteome generally. Every Ser/Thr in the human proteome is
then scored against this matrix, giving a ranked, continuous score
instead of a crude yes/no regex match.

Reuses the proteome FASTA already downloaded in Phase 1 

"""

import math
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from omnipath.requests import Enzsub

# ---------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------
DATA_DIR = Path("data")
PROTEOME_FASTA = DATA_DIR / "human_reviewed_proteome.fasta"
OUTPUT_CSV = DATA_DIR / "aurkb_candidate_sites_pssm.csv"

FLANK = 5              # 5 residues before the site, 4 after (10 total, site at index 5)
WINDOW_LEN = FLANK * 2  # = 10
PSEUDOCOUNT = 0.5
AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

# Which percentile of the KNOWN site score distribution to use as the final
# cutoff. Lower = keep more known sites (higher recall) but far more noise.
# Higher = stricter, fewer false leads.

CUTOFF_PERCENTILE = 50


# ---------------------------------------------------------------
#  Load the proteome into a lookup dict {uniprot_id: sequence}
# ---------------------------------------------------------------
def load_proteome(fasta_path: Path) -> dict:
    print("Loading proteome into memory...")
    seqs = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        uniprot_id = record.id.split("|")[1] if "|" in record.id else record.id
        seqs[uniprot_id] = str(record.seq)
    print(f"Loaded {len(seqs)} protein sequences.")
    return seqs


# ---------------------------------------------------------------
#  Pull known AURKB substrate sites from OmniPath
# ---------------------------------------------------------------
def load_known_aurkb_sites() -> pd.DataFrame:
    print("Querying OmniPath for known AURKB substrate sites...")
    df = Enzsub.get(enzymes="AURKB", genesymbols=True, organisms="human")
    known = df[["substrate", "substrate_genesymbol", "residue_offset", "residue_type"]].copy()
    known.columns = ["uniprot_id", "substrate_name", "position", "residue_type"]
    known = known.dropna(subset=["position"])
    known["position"] = known["position"].astype(int)

    # OmniPath aggregates multiple source databases, so the same site can
    # appear more than once (once per source that reported it). Dedupe so
    # every downstream count reflects unique biological sites.
    n_before = len(known)
    known = known.drop_duplicates(subset=["uniprot_id", "position"])
    n_duplicates = n_before - len(known)
    if n_duplicates:
        print(f"Removed {n_duplicates} duplicate rows (same site reported by multiple sources).")

    print(f"Found {len(known)} unique known AURKB substrate sites.")
    return known


# -----------------------------------------------------------------
# Extract the real sequence window around each known site
# -----------------------------------------------------------------

def extract_window(seq: str, position_1based: int, flank: int = FLANK) -> str | None:
    site_idx = position_1based - 1  # convert to 0-based
    start = site_idx - flank
    end = site_idx + flank  # exclusive - gives `flank*2` residues total
    if start < 0 or end > len(seq):
        return None
    return seq[start:end]


def build_known_windows(known_df: pd.DataFrame, proteome: dict) -> list[str]:
    windows = []
    n_missing_protein = 0
    n_boundary_skip = 0
    n_residue_mismatch = 0

    for _, row in known_df.iterrows():
        seq = proteome.get(row["uniprot_id"])
        if seq is None:
            n_missing_protein += 1
            continue

        window = extract_window(seq, row["position"])
        if window is None:
            n_boundary_skip += 1
            continue

        #  Mismatches usually mean
        # a UniProt version/isoform difference between OmniPath and the
        # downloaded proteome.
        site_residue = window[FLANK]
        expected = row["residue_type"]
        if expected and site_residue.upper() != str(expected).upper():
            n_residue_mismatch += 1
            continue

        windows.append(window)

    print(f"Built {len(windows)} usable known-site windows out of {len(known_df)} total known sites.")
    print(f"  Skipped: {n_missing_protein} (protein not found), "
          f"{n_boundary_skip} (too close to sequence end), "
          f"{n_residue_mismatch} (residue mismatch vs OmniPath)")
    return windows


# -----------------------------------------------------------------
#  Background amino acid frequency across the whole proteome
# --------------------------------------------------------------------

def compute_background_freq(proteome: dict) -> dict:
    counts = Counter()
    total = 0
    for seq in proteome.values():
        counts.update(seq)
        total += len(seq)
    return {aa: counts.get(aa, 0) / total for aa in AMINO_ACIDS}


# ---------------------------------------------------------------
# Build the PSSM from known-site windows
# ---------------------------------------------------------------
def build_pssm(known_windows: list[str], background: dict) -> list[dict]:
    position_counts = [Counter() for _ in range(WINDOW_LEN)]
    for window in known_windows:
        for i, aa in enumerate(window):
            position_counts[i][aa] += 1

    n = len(known_windows)
    pssm = []
    for i in range(WINDOW_LEN):
        row = {}
        for aa in AMINO_ACIDS:
            freq_known = (position_counts[i].get(aa, 0) + PSEUDOCOUNT) / (n + PSEUDOCOUNT * 20)
            bg_freq = background.get(aa, 1e-6)
            row[aa] = math.log2(freq_known / bg_freq)
        pssm.append(row)
    return pssm


def score_window(window: str, pssm: list[dict]) -> float:
    return sum(pssm[i].get(aa, 0.0) for i, aa in enumerate(window))


# --------------------------------------------------------------
#  Score every S/T in the proteome against the PSSM
# ---------------------------------------------------------------
def scan_proteome_with_pssm(proteome: dict, pssm: list[dict]) -> pd.DataFrame:
    rows = []
    for uniprot_id, seq in proteome.items():
        for i, residue in enumerate(seq):
            if residue not in ("S", "T"):
                continue
            window = extract_window(seq, i + 1)  # i is 0-based -> convert to 1-based
            if window is None:
                continue
            score = score_window(window, pssm)
            rows.append({
                "uniprot_id": uniprot_id,
                "position": i + 1,
                "residue": residue,
                "window": window,
                "pssm_score": score,
            })
    print(f"Scored {len(rows)} S/T sites across the proteome.")
    return pd.DataFrame(rows)


# ---------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------
def main():
    proteome = load_proteome(PROTEOME_FASTA)
    known_df = load_known_aurkb_sites()
    known_windows = build_known_windows(known_df, proteome)

    background = compute_background_freq(proteome)
    pssm = build_pssm(known_windows, background)

    known_scores = sorted(score_window(w, pssm) for w in known_windows)
    all_scores_df = scan_proteome_with_pssm(proteome, pssm)

    #  Cutoff sweep: show the actual tradeoff
    print("\nCutoff sweep (percentile of known-site scores | recovery | total candidates):")
    print(f"{'Percentile':>10} | {'Cutoff':>8} | {'Known recovered':>16} | {'Total candidates':>17}")
    for pct in (10, 25, 50, 75, 90):
        pct_cutoff = known_scores[int(len(known_scores) * pct / 100)]
        pct_candidates = all_scores_df[all_scores_df["pssm_score"] >= pct_cutoff]
        pct_known_recovered = sum(1 for s in known_scores if s >= pct_cutoff)
        marker = "  <- selected" if pct == CUTOFF_PERCENTILE else ""
        print(f"{pct:>9}% | {pct_cutoff:>8.2f} | "
              f"{pct_known_recovered:>6}/{len(known_scores)} ({pct_known_recovered/len(known_scores):.0%})"
              f"{'':>2} | {len(pct_candidates):>17,}{marker}")

    cutoff = known_scores[int(len(known_scores) * CUTOFF_PERCENTILE / 100)]
    print(f"\nUsing {CUTOFF_PERCENTILE}th percentile as the final cutoff: {cutoff:.2f}")
    print("(Change CUTOFF_PERCENTILE near the top of the script and rerun if you want a different tradeoff.)")

    candidates_df = all_scores_df[all_scores_df["pssm_score"] >= cutoff].copy()
    candidates_df = candidates_df.sort_values("pssm_score", ascending=False)

    known_pairs = set(zip(known_df["uniprot_id"], known_df["position"]))
    candidates_df["known_in_omnipath"] = candidates_df.apply(
        lambda r: (r["uniprot_id"], r["position"]) in known_pairs, axis=1
    )

    candidates_df.to_csv(OUTPUT_CSV, index=False)

    n_known_recovered = candidates_df["known_in_omnipath"].sum()
    n_orphan = len(candidates_df) - n_known_recovered
    recovery_rate = n_known_recovered / len(known_windows) if known_windows else 0

    print(f"\nSaved {len(candidates_df)} candidates above cutoff to {OUTPUT_CSV}")
    print(f"  -> {n_known_recovered} known sites recovered ({recovery_rate:.1%} of {len(known_windows)} usable known sites)")
    print(f"  -> {n_orphan} orphan candidates — your actual shortlist to explore further")
    print("\nReminder: recovery here is measured against the same sites used to build")
    print("the PSSM, so treat it as an internal consistency check, not proof the")
    print("model generalizes. Worth mentioning explicitly in your writeup.")


if __name__ == "__main__":
    main()
