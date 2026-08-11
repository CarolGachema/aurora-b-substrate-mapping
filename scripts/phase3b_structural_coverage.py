"""
Aurora B Substrate Mapping — Phase 3b: Structural Coverage & Accessibility Check
------------------------------------------------------------------------------
Alternative/complement to full AlphaFold-Multimer docking: rather than
predicting a kinase-substrate complex from scratch, this checks whether
each candidate phosphosite already has EXPERIMENTAL structural evidence
(X-ray/cryo-EM/NMR), and if so, whether that residue is solvent-
accessible -- a real, established, purely structural plausibility filter
(buried residues are very unlikely to be phosphorylatable without local
unfolding).

IMPORTANT -- read before trusting the output:
  Many of your candidates sit in large, intrinsically disordered
  nuclear/splicing-associated proteins (SRRM2, ZFC3H1, TRA2A, etc.)
  that are UNLIKELY to have full-length experimental structures at
  all. "No structure found" here does NOT mean "implausible" --
  disordered regions are actually enriched for real kinase substrates
  precisely because they're accessible. Treat "no structure" as
  inconclusive-leaning-plausible, not a strike against a candidate.

  For candidates WITH structural coverage, this also makes a
  simplifying assumption when mapping UniProt position -> PDB residue
  numbering (linear interpolation from the SIFTS start/end anchors,
  not a full per-residue alignment) -- fine for first-pass triage, but
  worth spot-checking any borderline "buried" call manually in Mol*
  before ruling a candidate out.

Steps:
  1. Load the same 25-candidate shortlist used for Phase 3a.
  2. Query PDBe's SIFTS "best_structures" API per UniProt ID to check
     for experimental structural coverage of the phosphosite position.
  3. Where coverage exists, download the structure and compute relative
     solvent accessibility for that residue with freesasa.
  4. Save one summary CSV -- use it to pick which 3-5 candidates
     actually deserve a manual look (and screenshot) in Mol*.

Requires: pip install freesasa requests pandas
"""

import time
from pathlib import Path

import freesasa
import pandas as pd
import requests

MANIFEST_CSV = Path("phase3_colabfold_inputs") / "manifest.csv"
OUTPUT_CSV = Path("data") / "phase3_structural_coverage.csv"
PDB_CACHE_DIR = Path("data") / "pdb_cache"

BURIED_THRESHOLD = 0.20   # relative SASA below this = buried
EXPOSED_THRESHOLD = 0.50  # relative SASA above this = clearly exposed


def get_best_structures(uniprot_id: str) -> list:
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_id}"
    resp = requests.get(url, timeout=30)
    if resp.status_code == 404:
        return []
    resp.raise_for_status()
    data = resp.json()
    return data.get(uniprot_id, [])


def find_covering_structure(structures: list, position: int):
    for s in structures:
        unp_start = s.get("unp_start")
        unp_end = s.get("unp_end")
        if unp_start is not None and unp_end is not None and unp_start <= position <= unp_end:
            return s
    return None


def download_pdb(pdb_id: str):
    PDB_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    out_path = PDB_CACHE_DIR / f"{pdb_id.lower()}.pdb"
    if out_path.exists():
        return out_path
    resp = requests.get(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb", timeout=60)
    if resp.status_code != 200:
        return None
    out_path.write_text(resp.text)
    return out_path


def _residue_number(value):
    """PDBe's API has been inconsistent about whether start/end come back as a
    plain int or a {'author_residue_number': ...} dict -- handle both."""
    if isinstance(value, dict):
        return value.get("author_residue_number", value.get("residue_number"))
    return value


def compute_relative_sasa(pdb_path: Path, chain_id: str, pdb_residue_number: int):
    structure = freesasa.Structure(str(pdb_path))
    result = freesasa.calc(structure)
    residue_areas = result.residueAreas()
    chain_areas = residue_areas.get(chain_id)
    if chain_areas is None:
        return None
    area = chain_areas.get(str(pdb_residue_number))
    if area is None:
        return None
    return area.relativeTotal


def main():
    manifest = pd.read_csv(MANIFEST_CSV)
    rows = []

    for _, row in manifest.iterrows():
        uniprot_id = row["uniprot_id"]
        position = int(row["position"])
        job_id = row["job_id"]
        print(f"\n{job_id} ({uniprot_id}, position {position})")

        try:
            structures = get_best_structures(uniprot_id)
        except Exception as e:
            print(f"  API error: {type(e).__name__}: {e}")
            rows.append({"job_id": job_id, "uniprot_id": uniprot_id, "position": position,
                         "has_structure": False, "note": f"API error: {e}"})
            continue

        match = find_covering_structure(structures, position)
        if match is None:
            print(f"  No experimental structure covers position {position} "
                  f"({len(structures)} structures found for this protein, none reach that residue).")
            rows.append({"job_id": job_id, "uniprot_id": uniprot_id, "position": position,
                         "has_structure": False,
                         "note": "No structural coverage - presumed disordered/inconclusive"})
            time.sleep(0.3)
            continue

        try:
            pdb_id = match["pdb_id"]
            chain_id = match["chain_id"]
            unp_start = match["unp_start"]
            pdb_start = _residue_number(match["start"])
            pdb_position = pdb_start + (position - unp_start)

            print(f"  Covered by {pdb_id.upper()} chain {chain_id} "
                  f"({match.get('experimental_method', '?')}, "
                  f"resolution {match.get('resolution', '?')}) -> PDB residue ~{pdb_position}")

            pdb_path = download_pdb(pdb_id)
            rel_sasa = None
            note = ""
            if pdb_path is None:
                note = "Structure found but download failed"
            else:
                try:
                    rel_sasa = compute_relative_sasa(pdb_path, chain_id, pdb_position)
                    if rel_sasa is None:
                        note = "Residue not resolved at this exact number in the PDB file (gap/insertion) - check manually"
                except Exception as e:
                    note = f"SASA calc failed: {type(e).__name__}: {e}"

            accessibility = None
            if rel_sasa is not None:
                if rel_sasa < BURIED_THRESHOLD:
                    accessibility = "buried"
                elif rel_sasa > EXPOSED_THRESHOLD:
                    accessibility = "exposed"
                else:
                    accessibility = "partially exposed"
                print(f"  Relative SASA: {rel_sasa:.2f} -> {accessibility}")

            rows.append({
                "job_id": job_id, "uniprot_id": uniprot_id, "position": position,
                "has_structure": True, "pdb_id": pdb_id, "chain_id": chain_id,
                "pdb_residue_number": pdb_position,
                "resolution": match.get("resolution"),
                "experimental_method": match.get("experimental_method"),
                "relative_sasa": rel_sasa, "accessibility": accessibility, "note": note,
            })
        except Exception as e:
            # Never let one candidate's weird API response kill the whole run --
            # log it and keep going, same spirit as the source-fallback loops
            # elsewhere in this pipeline.
            print(f"  Unexpected error processing this candidate: {type(e).__name__}: {e}")
            rows.append({"job_id": job_id, "uniprot_id": uniprot_id, "position": position,
                         "has_structure": True,
                         "note": f"Structure match found but processing failed: {type(e).__name__}: {e}"})

        time.sleep(0.3)  # be polite to the API

    out_df = pd.DataFrame(rows)
    OUTPUT_CSV.parent.mkdir(exist_ok=True)
    out_df.to_csv(OUTPUT_CSV, index=False)

    n_covered = int(out_df["has_structure"].sum())
    n_exposed = int((out_df.get("accessibility") == "exposed").sum()) if "accessibility" in out_df else 0
    print(f"\nSaved {len(out_df)} candidates to {OUTPUT_CSV}")
    print(f"  -> {n_covered} have experimental structural coverage of the phosphosite")
    print(f"  -> {len(out_df) - n_covered} have no structure (disordered/inconclusive, not necessarily bad)")
    print(f"  -> {n_exposed} confirmed surface-exposed among the covered ones")
    print("\nUse this to decide which 3-5 candidates are worth opening manually in Mol* "
          "for a publication screenshot -- no need to browse all 25 blind.")


if __name__ == "__main__":
    main()
