"""
Aurora B Substrate Mapping — Phase 5b: TCGA CIN Correlation
------------------------------------------------------------------
Correlates candidate gene expression against real chromosomal
instability metrics (Aneuploidy Score, Fraction Genome Altered) across
TCGA-BRCA's Basal-like ("BRCA_Basal") subtype samples.

CAVEAT worth stating explicitly in any writeup: this Pan-Cancer Atlas
clinical release doesn't include IHC-based ER/PR/HER2 status directly,
so Basal-like PAM50 subtype is used as the TNBC proxy here instead of
the same-definition TNBC status used in Phase 2. Basal-like and TNBC
overlap heavily in the literature (most TNBC tumors are Basal-like)
but aren't a perfect 1:1 match -- some Basal-like tumors aren't
clinically TNBC and vice versa.

For each of the TNBC-validated candidates (Phase 2e), tests whether
higher mRNA expression associates with higher genomic instability
across real tumors -- a cohort-scale signal independent of everything
computed in Phases 1-4.

Steps:
  1. Pull SAMPLE-level clinical data (the CIN metrics).
  2. Pull PATIENT-level clinical data (SUBTYPE), filter to Basal-like.
  3. Map candidate gene symbols -> Entrez IDs in chunks (cBioPortal's
     expression endpoint keys on Entrez, not symbols) -- dropped genes
     are printed, not silently lost.
  4. Fetch mRNA expression (RNA Seq V2 RSEM) for those genes across the
     Basal-like samples, also chunked (large single requests have been
     known to time out server-side).
  5. Spearman-correlate each gene's expression against both CIN metrics.
  6. Save one clean results CSV, ranked by correlation strength.

Requires: pip install scipy (if not already installed)
"""

from pathlib import Path
import time

import pandas as pd
import requests
from scipy.stats import spearmanr

BASE_URL = "https://www.cbioportal.org/api"
STUDY_ID = "brca_tcga_pan_can_atlas_2018"
EXPRESSION_PROFILE = f"{STUDY_ID}_rna_seq_v2_mrna"
CIN_ATTRIBUTES = ["ANEUPLOIDY_SCORE", "FRACTION_GENOME_ALTERED"]

CANDIDATES_CSV = Path("data") / "aurkb_candidates_tnbc_filtered.csv"
OUTPUT_CSV = Path("data") / "phase5_cin_correlation.csv"

CHUNK_SIZE = 200
EXPRESSION_MAX_RETRIES = 3


def chunked(items: list, size: int):
    for i in range(0, len(items), size):
        yield items[i:i + size]


def fetch_sample_clinical() -> pd.DataFrame:
    print("Fetching SAMPLE-level clinical data (CIN metrics)...")
    resp = requests.get(
        f"{BASE_URL}/studies/{STUDY_ID}/clinical-data",
        params={"clinicalDataType": "SAMPLE", "projection": "SUMMARY"},
        timeout=120,
    )
    resp.raise_for_status()
    data = pd.DataFrame(resp.json())
    print(f"  {len(data)} sample-level clinical records")
    wide = data.pivot_table(index=["sampleId", "patientId"], columns="clinicalAttributeId",
                             values="value", aggfunc="first").reset_index()
    print("  Columns available:", [c for c in wide.columns if c not in ("sampleId", "patientId")])
    return wide


def fetch_patient_clinical() -> pd.DataFrame:
    print("Fetching PATIENT-level clinical data (subtype)...")
    resp = requests.get(
        f"{BASE_URL}/studies/{STUDY_ID}/clinical-data",
        params={"clinicalDataType": "PATIENT", "projection": "SUMMARY"},
        timeout=120,
    )
    resp.raise_for_status()
    data = pd.DataFrame(resp.json())
    wide = data.pivot_table(index="patientId", columns="clinicalAttributeId",
                             values="value", aggfunc="first").reset_index()
    if "SUBTYPE" in wide.columns:
        print("  SUBTYPE value counts:")
        print(wide["SUBTYPE"].value_counts(dropna=False))
    else:
        print("  WARNING: no SUBTYPE column found. Available columns:", list(wide.columns))
    return wide


def map_symbols_to_entrez(symbols: list) -> dict:
    print(f"\nMapping {len(symbols)} gene symbols to Entrez IDs (chunked)...")
    mapping = {}
    for chunk in chunked(symbols, CHUNK_SIZE):
        try:
            resp = requests.post(
                f"{BASE_URL}/genes/fetch",
                params={"geneIdType": "HUGO_GENE_SYMBOL"},
                json=chunk,
                headers={"Content-Type": "application/json"},
                timeout=120,
            )
            resp.raise_for_status()
            for g in resp.json():
                if "hugoGeneSymbol" in g and "entrezGeneId" in g:
                    mapping[g["hugoGeneSymbol"]] = g["entrezGeneId"]
        except Exception as e:
            print(f"  Chunk failed ({type(e).__name__}: {e}), continuing with remaining chunks...")
    print(f"  Mapped {len(mapping)} / {len(symbols)} symbols "
          f"({len(symbols) - len(mapping)} not found -- often outdated/alias symbols)")
    return mapping


def fetch_expression(entrez_ids: list, sample_ids: list) -> pd.DataFrame:
    print(f"\nFetching expression for {len(entrez_ids)} genes across {len(sample_ids)} samples (chunked)...")
    all_data = []
    for chunk in chunked(entrez_ids, CHUNK_SIZE):
        for attempt in range(1, EXPRESSION_MAX_RETRIES + 1):
            try:
                resp = requests.post(
                    f"{BASE_URL}/molecular-profiles/{EXPRESSION_PROFILE}/molecular-data/fetch",
                    params={"projection": "SUMMARY"},
                    json={"entrezGeneIds": chunk, "sampleIds": sample_ids},
                    headers={"Content-Type": "application/json"},
                    timeout=180,
                )
                resp.raise_for_status()
                all_data.extend(resp.json())
                break
            except Exception as e:
                if attempt < EXPRESSION_MAX_RETRIES:
                    print(f"  Chunk attempt {attempt} failed ({type(e).__name__}: {e}), retrying...")
                    time.sleep(3)
                else:
                    print(f"  Chunk failed after {EXPRESSION_MAX_RETRIES} attempts "
                          f"({type(e).__name__}: {e}), giving up on this chunk "
                          f"({len(chunk)} genes will be missing from results).")
    data = pd.DataFrame(all_data)
    print(f"  Retrieved {len(data)} expression data points")
    return data


def main():
    sample_clinical = fetch_sample_clinical()
    patient_clinical = fetch_patient_clinical()

    if "SUBTYPE" not in patient_clinical.columns:
        raise RuntimeError("Can't find SUBTYPE column -- paste this output back to debug.")

    basal_patients = patient_clinical[patient_clinical["SUBTYPE"].astype(str).str.contains("Basal", na=False)]
    print(f"\n{len(basal_patients)} patients with a Basal-like SUBTYPE call "
          f"(TNBC proxy -- see caveat in the docstring).")

    merged_clinical = sample_clinical.merge(basal_patients[["patientId"]], on="patientId", how="inner")
    print(f"{len(merged_clinical)} Basal-like samples with clinical data.")

    cin_cols = [c for c in CIN_ATTRIBUTES if c in merged_clinical.columns]
    if not cin_cols:
        raise RuntimeError(f"None of {CIN_ATTRIBUTES} found in sample clinical data -- "
                            f"paste this output back to debug.")
    print(f"CIN metrics available: {cin_cols}")
    for c in cin_cols:
        merged_clinical[c] = pd.to_numeric(merged_clinical[c], errors="coerce")

    candidates = pd.read_csv(CANDIDATES_CSV)
    candidates = candidates[candidates["observed_in_tnbc"] == True].copy()
    candidate_genes = sorted(candidates["gene_symbol"].dropna().unique().tolist())
    print(f"\n{len(candidate_genes)} unique TNBC-validated candidate genes to test.")

    symbol_to_entrez = map_symbols_to_entrez(candidate_genes)
    if not symbol_to_entrez:
        raise RuntimeError("No genes mapped to Entrez IDs -- paste this output back to debug.")
    entrez_ids = list(symbol_to_entrez.values())
    entrez_to_symbol = {v: k for k, v in symbol_to_entrez.items()}

    sample_ids = merged_clinical["sampleId"].tolist()
    expression = fetch_expression(entrez_ids, sample_ids)
    if expression.empty:
        raise RuntimeError("No expression data returned -- paste this output back to debug.")

    expr_wide = expression.pivot_table(index="sampleId", columns="entrezGeneId", values="value", aggfunc="first")
    merged = merged_clinical.set_index("sampleId").join(expr_wide, how="inner")
    print(f"\n{len(merged)} samples with both CIN metrics and expression data -- running correlations...")

    results = []
    for entrez_id in entrez_ids:
        if entrez_id not in merged.columns:
            continue
        gene_symbol = entrez_to_symbol.get(entrez_id, str(entrez_id))
        try:
            expr_values = pd.to_numeric(merged[entrez_id], errors="coerce")
            row = {"gene_symbol": gene_symbol}
            n_used = 0
            for c in cin_cols:
                valid = merged[[c]].join(expr_values.rename("expr")).dropna()
                n_used = max(n_used, len(valid))
                if len(valid) < 10:
                    row[f"{c}_rho"] = None
                    row[f"{c}_pval"] = None
                    continue
                rho, pval = spearmanr(valid["expr"], valid[c])
                row[f"{c}_rho"] = rho
                row[f"{c}_pval"] = pval
            row["n_samples"] = n_used
            results.append(row)
        except Exception as e:
            print(f"  Skipping {gene_symbol}: {type(e).__name__}: {e}")

    results_df = pd.DataFrame(results)
    primary_col = f"{cin_cols[0]}_rho"
    results_df = results_df.sort_values(primary_col, ascending=False)
    OUTPUT_CSV.parent.mkdir(exist_ok=True)
    results_df.to_csv(OUTPUT_CSV, index=False)

    n_sig = int((results_df[f"{cin_cols[0]}_pval"] < 0.05).sum())
    print(f"\nSaved {len(results_df)} gene correlations to {OUTPUT_CSV}")
    print(f"  -> {n_sig} genes show a significant correlation (p<0.05) with {cin_cols[0]}")
    print("\nTop 10 positively correlated:")
    print(results_df.head(10).to_string(index=False))
    print("\nBottom 10 (most negatively correlated):")
    print(results_df.tail(10).to_string(index=False))


if __name__ == "__main__":
    main()
