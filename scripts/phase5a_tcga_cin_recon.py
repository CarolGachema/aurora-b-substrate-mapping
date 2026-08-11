"""
Aurora B Substrate Mapping — Phase 5a: TCGA CIN Signature Recon
------------------------------------------------------------------
Phases 1-4 used CPTAC (protein-level, ~150 patients) and OmniPath
(interaction data). Neither touches chromosomal instability directly.
This phase needs TCGA's much larger BRCA cohort (~1,000 patients),
which has real genomic instability metrics (aneuploidy score, fraction
genome altered) and mRNA expression data -- letting us test whether
candidate genes' expression actually correlates with genomic
instability across real tumors, not just whether they're phosphorylated.

Same "list what's actually there before building on it" approach as
Phase 2a/2b -- don't guess attribute names, print them and check.

Target study: brca_tcga_pan_can_atlas_2018 (TCGA PanCancer Atlas) --
the standard, most complete TCGA breast cancer resource on cBioPortal.
If that exact ID has changed/moved, this script searches by keyword
instead of failing silently.
"""

import requests
import pandas as pd

BASE_URL = "https://www.cbioportal.org/api"
STUDY_ID = "brca_tcga_pan_can_atlas_2018"


def main():
    print("=" * 60)
    print(f"CHECKING STUDY EXISTS: {STUDY_ID}")
    print("=" * 60)
    resp = requests.get(f"{BASE_URL}/studies/{STUDY_ID}", timeout=30)
    if resp.status_code != 200:
        print(f"Study ID not found (status {resp.status_code}). "
              f"Searching cBioPortal for the right TCGA breast cancer study instead...")
        resp2 = requests.get(f"{BASE_URL}/studies", params={"keyword": "breast", "projection": "SUMMARY"}, timeout=30)
        resp2.raise_for_status()
        studies = pd.DataFrame(resp2.json())
        print(studies[["studyId", "name"]].to_string())
        print("\nPick the right studyId from above, update STUDY_ID at the top of this "
              "script, and rerun.")
        return
    print("Study found:", resp.json().get("name"))

    print("\n" + "=" * 60)
    print("AVAILABLE CLINICAL ATTRIBUTES")
    print("=" * 60)
    resp = requests.get(f"{BASE_URL}/studies/{STUDY_ID}/clinical-attributes", timeout=60)
    resp.raise_for_status()
    attrs = pd.DataFrame(resp.json())
    print(f"Total clinical attributes available: {len(attrs)}")
    print(attrs[["clinicalAttributeId", "displayName"]].to_string())

    keywords = ["aneuploidy", "cin", "instability", "fraction_genome_altered",
                "msi", "subtype", "receptor", "er_status", "pr_status", "her2"]
    relevant = attrs[attrs["clinicalAttributeId"].str.lower().str.contains("|".join(keywords))]
    print("\nCIN / subtype-relevant attributes found:")
    print(relevant[["clinicalAttributeId", "displayName", "description"]].to_string())

    print("\n" + "=" * 60)
    print("MOLECULAR PROFILES (checking for mRNA expression data)")
    print("=" * 60)
    resp = requests.get(f"{BASE_URL}/studies/{STUDY_ID}/molecular-profiles", timeout=60)
    resp.raise_for_status()
    profiles = pd.DataFrame(resp.json())
    print(profiles[["molecularProfileId", "name", "molecularAlterationType"]].to_string())

    print("\nPaste this whole output back and I'll write the real Phase 5 pull "
          "using whatever attribute names and profile IDs actually exist here.")


if __name__ == "__main__":
    main()
