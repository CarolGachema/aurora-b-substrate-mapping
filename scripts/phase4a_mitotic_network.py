"""
Aurora B Substrate Mapping — Phase 4a: Mitotic Network Context
------------------------------------------------------------------
Checks whether each of the TNBC-validated candidates (Phase 2e's
1,482-candidate shortlist) sits in the immediate interaction
neighborhood of known mitotic fidelity machinery -- the kinetochore,
the spindle assembly checkpoint (SAC), and Aurora B's own activating
complex (the Chromosomal Passenger Complex, CPC).

Two live-pulled sources, no hand-typed gene lists:
  1. EBI's QuickGO API for the reference "known mitotic machinery" gene
     set, built from three GO terms:
       - GO:0000776  kinetochore
       - GO:0007094  mitotic spindle assembly checkpoint signaling
       - GO:0032133  chromosome passenger complex
     The GO term name QuickGO actually returns for each ID is printed
     -- check it matches what's listed above before trusting the run;
     if a term ID has been retired/merged since this was written, the
     printed name will make that obvious immediately.
  2. OmniPath (same package Phase 1 used for known AURKB sites) for
     real protein-protein interaction data, to check whether each
     candidate directly interacts with anything in that reference set.

This is a DIRECT (1-hop) interactor check -- deliberately conservative.
No hit doesn't disqualify a candidate, it just means this particular
piece of evidence doesn't support it (yet); extending to 2-hop
neighbors is a natural next iteration if the direct-hit list turns out
too small to be useful on its own.
"""

from pathlib import Path

import pandas as pd
import requests
from omnipath.interactions import AllInteractions

CANDIDATES_CSV = Path("data") / "aurkb_candidates_tnbc_filtered.csv"
OUTPUT_CSV = Path("data") / "aurkb_candidates_network_context.csv"

GO_TERMS = {
    "GO:0000776": "kinetochore",
    "GO:0007094": "mitotic spindle assembly checkpoint signaling",
    "GO:0032133": "chromosome passenger complex",
}

QUICKGO_URL = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"


def fetch_go_gene_set(go_id: str, expected_name: str) -> set:
    genes = set()
    page = 1
    while True:
        resp = requests.get(
            QUICKGO_URL,
            params={"goId": go_id, "taxonId": 9606, "geneProductType": "protein",
                    "limit": 100, "page": page},
            headers={"Accept": "application/json"},
            timeout=60,
        )
        resp.raise_for_status()
        data = resp.json()
        results = data.get("results", [])
        if page == 1:
            actual_name = (results[0].get("goName") if results else None) or "(no goName returned)"
            match_note = "" if actual_name.lower() == expected_name.lower() else "  <-- MISMATCH, check this term ID"
            print(f"{go_id}: expected '{expected_name}', QuickGO returned '{actual_name}'{match_note}")
        for r in results:
            qualifier = (r.get("qualifier") or "")
            if qualifier.upper().startswith("NOT"):
                continue  # skip negated annotations
            symbol = r.get("symbol")
            if symbol:
                genes.add(symbol.upper())
        total_pages = data.get("pageInfo", {}).get("total", 1)
        if page >= total_pages:
            break
        page += 1
    return genes


def build_reference_gene_set() -> dict:
    ref = {}
    for go_id, name in GO_TERMS.items():
        try:
            genes = fetch_go_gene_set(go_id, name)
            print(f"  -> {len(genes)} genes annotated to {go_id} ({name})")
        except Exception as e:
            print(f"  -> FAILED to fetch {go_id} ({name}): {type(e).__name__}: {e}")
            print(f"     Continuing with the other terms; this category will be under-represented.")
            genes = set()
        for g in genes:
            ref.setdefault(g, set()).add(name)
    return ref


def main():
    print("=" * 60)
    print("BUILDING REFERENCE MITOTIC MACHINERY GENE SET")
    print("=" * 60)
    reference_genes = build_reference_gene_set()
    print(f"\nTotal unique reference genes across all 3 categories: {len(reference_genes)}")

    print("\n" + "=" * 60)
    print("LOADING TNBC-VALIDATED CANDIDATES")
    print("=" * 60)
    candidates = pd.read_csv(CANDIDATES_CSV)
    candidates = candidates[candidates["observed_in_tnbc"] == True].copy()
    print(f"Working from the {len(candidates)} TNBC-validated candidates (Phase 2e's shortlist).")
    candidates["gene_symbol_upper"] = candidates["gene_symbol"].str.upper()
    candidate_genes = set(candidates["gene_symbol_upper"].dropna())
    print(f"{len(candidate_genes)} unique candidate genes to check.")

    print("\n" + "=" * 60)
    print("FETCHING HUMAN PROTEIN-PROTEIN INTERACTIONS FROM OMNIPATH")
    print("=" * 60)
    interactions = AllInteractions.get(genesymbols=True, organisms="human")
    print(f"Loaded {len(interactions)} interactions. Columns: {list(interactions.columns)}")

    partners = {}
    for _, row in interactions.iterrows():
        a, b = row.get("source_genesymbol"), row.get("target_genesymbol")
        if pd.isna(a) or pd.isna(b):
            continue
        a, b = str(a).upper(), str(b).upper()
        partners.setdefault(a, set()).add(b)
        partners.setdefault(b, set()).add(a)

    # Reference-set genes with an extremely high number of TOTAL interaction
    # partners (well-studied hub proteins like ATM) will connect to almost
    # anything by chance, not because of real mitotic-specific proximity.
    # Report the raw connectivity so this is visible, then compute BOTH an
    # "any" count (includes hubs) and a "specific" count (hubs excluded) --
    # don't silently pick one, since the difference between them IS the story.
    ref_degrees = {g: len(partners.get(g, set())) for g in reference_genes}
    print("\nReference gene connectivity (top 15 by degree -- watch for hub outliers):")
    for g, d in sorted(ref_degrees.items(), key=lambda x: -x[1])[:15]:
        print(f"  {g}: {d} total interaction partners")

    HUB_DEGREE_THRESHOLD = 300
    hub_genes = {g for g, d in ref_degrees.items() if d > HUB_DEGREE_THRESHOLD}
    if hub_genes:
        print(f"\nExcluding {len(hub_genes)} hub gene(s) from the 'specific' interactor check "
              f"(>{HUB_DEGREE_THRESHOLD} total partners -- too promiscuous to count as meaningful "
              f"evidence of mitotic-specific proximity): {sorted(hub_genes)}")
    reference_genes_specific = {g: cats for g, cats in reference_genes.items() if g not in hub_genes}

    results = []
    for gene in sorted(candidate_genes):
        direct_partners = partners.get(gene, set())
        hits_any = direct_partners & reference_genes.keys()
        hits_specific = direct_partners & reference_genes_specific.keys()
        results.append({
            "gene_symbol_upper": gene,
            "direct_mitotic_interactor_any": len(hits_any) > 0,
            "direct_mitotic_interactor_specific": len(hits_specific) > 0,
            "interacting_partners_any": ";".join(sorted(hits_any)) if hits_any else "",
            "interacting_partners_specific": ";".join(sorted(hits_specific)) if hits_specific else "",
            "partner_categories_specific": ";".join(sorted({c for h in hits_specific for c in reference_genes[h]})) if hits_specific else "",
        })
    network_df = pd.DataFrame(results)

    merged = candidates.merge(network_df, on="gene_symbol_upper", how="left")
    merged["direct_mitotic_interactor_any"] = merged["direct_mitotic_interactor_any"].fillna(False)
    merged["direct_mitotic_interactor_specific"] = merged["direct_mitotic_interactor_specific"].fillna(False)
    merged = merged.sort_values(["direct_mitotic_interactor_specific", "pssm_score"], ascending=[False, False])
    OUTPUT_CSV.parent.mkdir(exist_ok=True)
    merged.to_csv(OUTPUT_CSV, index=False)

    n_hits_any = int(network_df["direct_mitotic_interactor_any"].sum())
    n_hits_specific = int(network_df["direct_mitotic_interactor_specific"].sum())
    print(f"\nSaved {len(merged)} candidates with network context to {OUTPUT_CSV}")
    print(f"  -> {n_hits_any} / {len(network_df)} candidate genes interact with SOME reference gene (hubs included)")
    print(f"  -> {n_hits_specific} / {len(network_df)} candidate genes interact with a NON-hub reference gene "
          f"-- this is the number worth trusting")
    if n_hits_specific:
        print("\nSpecific direct hits (hub-excluded):")
        print(merged[merged["direct_mitotic_interactor_specific"]][
            ["gene_symbol", "residue", "position", "pssm_score", "interacting_partners_specific", "partner_categories_specific"]
        ].to_string(index=False))


if __name__ == "__main__":
    main()
