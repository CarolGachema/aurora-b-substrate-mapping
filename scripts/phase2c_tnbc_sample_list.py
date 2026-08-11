"""
Aurora B Substrate Mapping — Phase 2c: TNBC Sample List
------------------------------------------------------------
Builds a clean, patient-level TNBC status table from cBioPortal's
brca_cptac_2020 study, ready to join against the phosphoproteomics
data pulled via the cptac package.

Two cleanup steps are required.
  1. Status values are inconsistently capitalized ("negative" vs
     "Negative"), normalized here before use.
  2. cBioPortal prefixes numeric-style patient IDs with "X" (e.g.
     "X01BR001"), while the cptac package's phosphoproteomics table
     uses "01BR001" directly — stripped here so the two can be joined
     later. A handful of patients ("CPT000814"-style IDs) use a
     completely different naming scheme and won't match anything in
     the phosphoproteomics table, that's expected, and
     they're just left as-is (they'll simply just fail to join later).
"""

from pathlib import Path

import pandas as pd
import requests

BASE_URL = "https://www.cbioportal.org/api"
STUDY_ID = "brca_cptac_2020"
DATA_DIR = Path("data")
OUTPUT_CSV = DATA_DIR / "tnbc_status_by_patient.csv"


def fetch_clinical_data() -> pd.DataFrame:
    print("Fetching clinical data from cBioPortal...")
    resp = requests.get(
        f"{BASE_URL}/studies/{STUDY_ID}/clinical-data",
        params={"clinicalDataType": "PATIENT", "projection": "SUMMARY"},
        timeout=120,
    )
    resp.raise_for_status()
    data = pd.DataFrame(resp.json())
    wide = data.pivot(index="patientId", columns="clinicalAttributeId", values="value")
    print(f"Pulled clinical data for {len(wide)} patients.")
    return wide


def normalize_patient_id(raw_id: str) -> str:
    """'X01BR001' -> '01BR001'; 'CPT000814' is left untouched (won't join later)."""
    if raw_id.upper().startswith("X") and raw_id[1:3].isdigit():
        return raw_id[1:]
    return raw_id


def main():
    DATA_DIR.mkdir(exist_ok=True)
    wide = fetch_clinical_data()

    if "TNBC_UPDATED_CLINICAL_STATUS" not in wide.columns:
        raise RuntimeError(
            "TNBC_UPDATED_CLINICAL_STATUS column is missing — the cBioPortal "
            "API response may have changed. Rerun phase2b_cbioportal_clinical.py "
            "and check the attribute list again."
        )

    status = wide["TNBC_UPDATED_CLINICAL_STATUS"].dropna()
    status_clean = status.str.strip().str.capitalize()  # "negative"/"Negative" -> "Negative"

    print("\nNormalized TNBC status counts:")
    print(status_clean.value_counts())

    result = pd.DataFrame({
        "patient_id_raw": status_clean.index,
        "tnbc_status": status_clean.values,
    })
    result["patient_id"] = result["patient_id_raw"].apply(normalize_patient_id)
    result["is_tnbc"] = result["tnbc_status"] == "Positive"

    # Flag which IDs are in the "won't join later" bucket, so it's visible
    # now instead of silently losing patients further down the pipeline.
    unmatched_format = result[result["patient_id"] == result["patient_id_raw"]]
    unmatched_format = unmatched_format[~unmatched_format["patient_id_raw"].str.match(r"^\d")]
    if len(unmatched_format):
        print(f"\n{len(unmatched_format)} patients have an ID format that won't match "
              f"the phosphoproteomics table (e.g. {unmatched_format['patient_id_raw'].iloc[0]}) "
              f"— they'll be excluded once we join, not an error.")

    result.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSaved {len(result)} patients' TNBC status to {OUTPUT_CSV}")
    print(f"  -> {result['is_tnbc'].sum()} TNBC-positive")
    print(f"  -> {(~result['is_tnbc']).sum()} TNBC-negative")


if __name__ == "__main__":
    main()
