"""
Aurora B Substrate Mapping 
Phase 6: Composite Ranking
------------------------------------------------------------------
Combines the independent evidence layers into ONE actual composite
score per candidate site, instead of the AND-filter approach used in
every figure up to this point (motif AND observed-in-TNBC AND network
-proximal, treated as separate yes/no checks.

"""

from pathlib import Path

import pandas as pd

NETWORK_CSV = Path("data") / "aurkb_candidates_network_context.csv"
CIN_CSV = Path("data") / "phase5_cin_correlation.csv"
OUTPUT_CSV = Path("data") / "phase6_composite_ranking.csv"

WEIGHT_MOTIF = 0.50
WEIGHT_NETWORK = 0.25
WEIGHT_CIN = 0.25

# Manual snapshot from Phase 3's Mol* checks -- only 3 genes have any
# structural data at all (see caveat above for why this isn't scored).
STRUCTURAL_NOTES = {
    "RGL2": "Confirmed surface-exposed (RGL2 S736, PDB 8B69) -- structurally validated",
    "PBRM1": "Structure exists (7Y8R) but this residue unresolved -- inconclusive",
    "SRRM2": "Structure exists (8C6J) but this residue unresolved -- inconclusive",
}


def main():
    candidates = pd.read_csv(NETWORK_CSV)
    cin = pd.read_csv(CIN_CSV)
    print(f"Loaded {len(candidates)} candidate sites, {len(cin)} gene-level CIN correlations.")

    merged = candidates.merge(cin, on="gene_symbol", how="left")
    n_with_cin = merged["ANEUPLOIDY_SCORE_rho"].notna().sum()
    print(f"{n_with_cin} / {len(merged)} candidate sites have CIN correlation data "
          f"(rest get a neutral score for this component, not a penalty).")

    # --- Motif component: percentile rank of pssm_score ---
    merged["motif_percentile"] = merged["pssm_score"].rank(pct=True)

    # --- Network component: binary, hub-excluded direct interactor ---
    merged["network_component"] = merged["direct_mitotic_interactor_specific"].fillna(False).astype(float)

    # --- CIN component: significance-weighted rho, soft not hard-thresholded ---
    rho = merged["ANEUPLOIDY_SCORE_rho"]
    pval = merged["ANEUPLOIDY_SCORE_pval"]
    cin_weighted = rho * (1 - pval)
    cin_filled = cin_weighted.fillna(cin_weighted.mean())

    # Combining raw percentile (0-1 continuous), raw binary (0/1), and raw
    # min-max scaled values directly is a bug, not a design choice: those
    # three scales don't have comparable spread, so whichever one happens to
    # jump the most (the binary network flag, here) silently dominates the
    # "weighted" sum regardless of the documented weights. Standardizing
    # (z-scoring) each component first makes the stated 50/25/25 weights
    # actually true in practice, not just in the docstring.
    def zscore(s: pd.Series) -> pd.Series:
        std = s.std()
        return (s - s.mean()) / std if std > 0 else s * 0

    merged["motif_z"] = zscore(merged["motif_percentile"])
    merged["network_z"] = zscore(merged["network_component"])
    merged["cin_z"] = zscore(cin_filled)

    merged["composite_index"] = (
        WEIGHT_MOTIF * merged["motif_z"]
        + WEIGHT_NETWORK * merged["network_z"]
        + WEIGHT_CIN * merged["cin_z"]
    )
    # Report as a percentile too -- easier to read than a raw z-score sum.
    merged["composite_score"] = merged["composite_index"].rank(pct=True)

    merged["structural_note"] = merged["gene_symbol"].map(STRUCTURAL_NOTES).fillna("No structural data yet")

    merged = merged.sort_values("composite_score", ascending=False)
    OUTPUT_CSV.parent.mkdir(exist_ok=True)
    merged.to_csv(OUTPUT_CSV, index=False)

    print(f"\nSaved composite ranking for {len(merged)} candidate sites to {OUTPUT_CSV}")

    top20 = merged.head(20)
    n_network_hits_in_top20 = int(top20["direct_mitotic_interactor_specific"].sum())
    print(f"\nSanity check: {n_network_hits_in_top20} / 20 top candidates are network hits "
          f"(should be a mix, not all 20 -- if it's still all/none, the balance needs another look).")

    print("\nTop 15 by composite score:")
    display_cols = ["gene_symbol", "residue", "position", "pssm_score", "motif_z",
                     "direct_mitotic_interactor_specific", "network_z",
                     "ANEUPLOIDY_SCORE_rho", "ANEUPLOIDY_SCORE_pval", "cin_z",
                     "composite_score", "structural_note"]
    print(merged[display_cols].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
