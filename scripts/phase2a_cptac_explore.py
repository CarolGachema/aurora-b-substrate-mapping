"""
Aurora B Substrate Mapping — Phase 2a: CPTAC Recon
------------------------------------------------------
Before writing the real filtering logic, we need to see the ACTUAL
shape of CPTAC's breast cancer data — column naming, clinical fields,
which data source has phosphoproteomics — rather than guessing and
risking a wasted download.

This will trigger a real download the first time you run it (CPTAC
data isn't tiny — could take a while depending on your connection).
That's expected, it's a one-time cost; it caches locally after that.

Just run this, then paste me everything it prints.
"""

import cptac

def main():
    print("Initializing CPTAC breast cancer (BRCA) dataset...")
    print("(First run downloads data — this may take a few minutes.)\n")
    brca = cptac.Brca()

    print("=" * 60)
    print("AVAILABLE DATA SOURCES / DATA TYPES")
    print("=" * 60)
    sources_df = brca.list_data_sources()
    print(sources_df)

    print("\n" + "=" * 60)
    print("CLINICAL DATA — looking for ER/PR/HER2 or subtype columns")
    print("=" * 60)
    clinical = brca.get_clinical(source="mssm")
    print("Clinical table shape:", clinical.shape)
    print("\nAll clinical columns:")
    print(list(clinical.columns))

    # No clean ER/PR/HER2 column exists in this table — receptor status is
    # likely buried in the free-text pathology fields instead. Show us what's
    # actually in there so we know how to parse it.
    ihc_cols = [c for c in clinical.columns if "immunohistochemistry" in c.lower()]
    print("\nImmunohistochemistry-related columns:", ihc_cols)
    if ihc_cols:
        for col in ihc_cols:
            non_null = clinical[col].dropna()
            print(f"\n--- Sample non-null values from '{col}' ({len(non_null)} total) ---")
            print(non_null.head(10).to_list())

    # Also check medical_history (another mssm-only data type) in case
    # subtype/receptor info lives there instead.
    print("\n" + "=" * 60)
    print("MEDICAL HISTORY DATA — checking for subtype/receptor info")
    print("=" * 60)
    try:
        med_history = brca.get_medical_history(source="mssm")
        print("Medical history columns:", list(med_history.columns))
        print(med_history.head(5))
    except Exception as e:
        print(f"Couldn't load medical_history ({type(e).__name__}): {e}")

    print("\n" + "=" * 60)
    print("PHOSPHOPROTEOMICS DATA — checking structure")
    print("=" * 60)
    # Try each source until one works, since not every source provides
    # phosphoproteomics for BRCA.
    phospho = None
    used_source = None
    for src in ("umich", "bcm", "broad", "washu"):
        try:
            phospho = brca.get_phosphoproteomics(source=src)
            used_source = src
            break
        except Exception as e:
            print(f"  Source '{src}' didn't work for phosphoproteomics ({type(e).__name__}), trying next...")

    if phospho is None:
        print("\nCouldn't load phosphoproteomics from any known source automatically.")
        print("Paste this whole output to me and I'll adjust the source list.")
        return

    print(f"\nUsing source: {used_source}")
    print("Phosphoproteomics table shape:", phospho.shape)
    print("\nColumn type:", type(phospho.columns))
    print("First 10 columns:")
    print(phospho.columns[:10])
    print("\nFirst 5 row index values (samples):")
    print(phospho.index[:5])
    print("\nA small data preview:")
    print(phospho.iloc[:5, :5])


if __name__ == "__main__":
    main()
