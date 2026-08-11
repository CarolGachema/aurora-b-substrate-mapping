"""
Aurora B Substrate Mapping 
Phase 2b: cBioPortal Clinical Recon
--------------------------------------------------------------------

This script:
  1. Lists all available clinical attributes for the study, so we can
     see the REAL attribute names rather than guessing.
  2. Pulls the full patient-level clinical data table.
  3. Checks whether patient IDs match the "01BR001"-style IDs used by
     the phosphoproteomics table from Phase 2a (needed to join later).

"""

import requests
import pandas as pd

BASE_URL = "https://www.cbioportal.org/api"
STUDY_ID = "brca_cptac_2020"


def main():
    # list clinical attributes available for this study 
    print("=" * 60)
    print("AVAILABLE CLINICAL ATTRIBUTES")
    print("=" * 60)
    resp = requests.get(f"{BASE_URL}/studies/{STUDY_ID}/clinical-attributes", timeout=60)
    resp.raise_for_status()
    attrs = pd.DataFrame(resp.json())
    print(f"Total clinical attributes available: {len(attrs)}")
    print("\nAll attribute IDs:")
    print(attrs[["clinicalAttributeId", "displayName"]].to_string())

    keywords = ["er_updated", "pr_clinical", "erbb2", "tnbc", "pam50"]
    relevant = attrs[attrs["clinicalAttributeId"].str.lower().str.contains("|".join(keywords))]
    print("\nAttributes that look receptor/subtype-related:")
    print(relevant[["clinicalAttributeId", "displayName", "description"]].to_string())

    # pull full patient-level clinical data 
    print("\n" + "=" * 60)
    print("PATIENT CLINICAL DATA")
    print("=" * 60)
    resp = requests.get(
        f"{BASE_URL}/studies/{STUDY_ID}/clinical-data",
        params={"clinicalDataType": "PATIENT", "projection": "SUMMARY"},
        timeout=120,
    )
    resp.raise_for_status()
    data = pd.DataFrame(resp.json())
    print("Raw shape (long format):", data.shape)
    print("\nColumns:", list(data.columns))

    # Pivot from long format (one row per patient+attribute) to wide
    # (one row per patient, one column per attribute) for easier use.
    wide = data.pivot(index="patientId", columns="clinicalAttributeId", values="value")
    print("\nPivoted shape (wide format):", wide.shape)

    print("\nFirst 10 patient IDs (checking ID format(s) present):")
    print(list(wide.index[:10]))

    # Two ID formats seem to be mixed in this study ("X01BR001" vs
    # "CPT000814"). This splits them apart so we can see how many of each
    # exist, and whether the "X01BR..." ones (after stripping the "X")
    # actually match our phosphoproteomics table's "01BR..." IDs.
    x_prefixed = [i for i in wide.index if i.upper().startswith("X0") or i.upper().startswith("X1")]
    other_format = [i for i in wide.index if i not in x_prefixed]
    print(f"\n'X0...'-style IDs: {len(x_prefixed)}")
    print(f"Other-format IDs (e.g. 'CPT...'): {len(other_format)}")
    if other_format:
        print("Sample of the other format:", other_format[:5])

    relevant_ids = relevant["clinicalAttributeId"].tolist()
    existing_relevant = [c for c in relevant_ids if c in wide.columns]
    print("\nReceptor/subtype columns actually present in the pivoted data:", existing_relevant)
    if existing_relevant:
        print("\nValue counts for each:")
        for col in existing_relevant:
            print(f"\n--- {col} ---")
            print(wide[col].value_counts(dropna=False))


if __name__ == "__main__":
    main()
