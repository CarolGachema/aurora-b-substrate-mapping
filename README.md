# Aurora B Substrate Mapping

**In silico identification of Aurora B kinase (AURKB) substrates driving chromosomal instability in triple-negative breast cancer (TNBC).**

> Status: active / work-in-progress research pipeline. This README reflects what has actually been run and validated so far — see [Planned / Not Yet Implemented](#planned--not-yet-implemented) for what's still open.

---

## Overview

Aurora B kinase is a core regulator of mitotic fidelity, and its overexpression in TNBC has been linked to chromosomal instability (CIN) — a hallmark that drives tumor heterogeneity and therapy resistance. Only a handful of AURKB substrates are well characterized. This project builds a computational pipeline to identify and prioritize candidate AURKB substrates by integrating:

- Consensus motif prediction
- TNBC-specific phosphoproteomic evidence (CPTAC)
- Structural context from experimental structures (PDB)
- Functional network proximity to known mitotic fidelity machinery (OmniPath)

The goal is a prioritized, testable substrate list as a starting point for future wet-lab validation. CIN/TCGA correlation and true kinase–substrate structural modeling are planned extensions, not yet implemented (see below).

---

## Pipeline

### Phase 1 — Motif Scanning
Scans the human proteome for AURKB consensus phosphorylation motifs and scores candidate S/T sites with a PSSM (position-specific scoring matrix). Produces the full candidate pool (21,293 sites) that every downstream phase filters from.

### Phase 2 — TNBC Phosphoproteomic Evidence (CPTAC)
- **2a**: Explored the CPTAC BRCA dataset (`umich` phosphoproteomics source, 158 patients) to understand available data types and structure.
- **2b/2c**: Since CPTAC's own clinical table lacks a clean receptor-subtype field, pulled TNBC status separately from cBioPortal (`brca_cptac_2020` study). Normalized inconsistent status capitalization and patient ID formats (cBioPortal's `X`-prefixed IDs vs. CPTAC's raw IDs) so the two sources can be joined. Result: 119 patients with usable TNBC status (28 positive / 91 negative).
- **2d**: Reconciliation pass — built a UniProt→gene-symbol map from the proteome FASTA and checked how well Phase 1 candidate genes overlap with phosphoproteomics columns before committing to a filter.
- **2e**: Joined the phosphosite table to TNBC-positive samples (29/158 matched), exploded to unique (gene, residue, position) sites, and applied a ≥20% detection-fraction filter within TNBC samples. **Output: 1,482 TNBC-observed candidate substrates**, ranked by PSSM score — the core shortlist used by every later phase.

### Phase 3 — Structural Context
- **3a**: Selected the top 25 orphan (not previously known), TNBC-observed candidates and prepared per-candidate FASTA inputs for ColabFold/AlphaFold batch folding.
- **3b**: ColabFold batch runs hit repeated tooling friction (Drive folder handling, Colab session limits) and were not completed for the full candidate set. Pivoted to a **structure-database lookup** instead: for each candidate, checked whether an existing experimental PDB structure covers that residue, and if so, whether the site is solvent-accessible.
  - **3/25 candidates have experimental structural coverage of the actual phosphosite**: PBRM1_S39, SRRM2_S1923, RGL2_S736.
  - **0/3 covered sites are confirmed surface-exposed** — worth treating as a real (if inconvenient) result, but also worth a manual SASA sanity-check given some `UNK`-residue warnings during the run.
  - ⚠️ RGL2_S736's structural "coverage" needs a manual check — the reported PDB residue (~96) doesn't line up with the full-length position (736), which usually signals a fragment/construct numbering offset rather than genuine coverage of that exact residue.

  *Note: this phase currently answers "is this site structurally observable/accessible?" 

### Phase 4 — Network Proximity to Mitotic Machinery
- **4a**: Built a 220-gene reference set from GO annotations (kinetochore, spindle assembly checkpoint, chromosomal passenger complex), pulled the human interactome from OmniPath (290,512 interactions), and excluded 5 highly promiscuous hub genes (>300 partners: ATM, SMARCA4, GTF2B, SIN3A, PLK1) from the "specific interactor" check. Of 1,070 unique candidate genes, **190 have a specific (non-hub) interaction with a known mitotic-machinery gene** — the network-validated subset.
- **4b**: Figures summarizing Phases 1–4a (see `figures/`). ⚠️ One known rendering bug: the "Mitotic partner(s)" column in the top-15 combined-evidence table overflows for rows with many partners (ATR S435, RCC1 S11) — needs a text-wrap or truncation fix before this figure is final.

---

## Key Results So Far

| Metric | Value |
|---|---|
| Total candidate sites (motif scan) | 21,293 |
| TNBC-observed candidates | 1,482 |
| — previously known AURKB sites, TNBC-observed | 31 |
| — orphan (novel) sites, TNBC-observed | 1,451 |
| Candidates with experimental structural coverage | 3 / 25 (top orphan subset checked) |
| Candidate genes with specific mitotic-network interaction | 190 / 1,070 |

Top prioritized candidates (motif + TNBC phosphoproteomic evidence + network proximity) include LBR S97, PRKD1 S205, ATR S435, PRKCD S302, and TP53BP1 S1686 — see `figures/fig4_top_combined_evidence.png` for the full top-15.

---

## Repository Structure

```
aurora-b-substrate-mapping/
├── README.md
├── scripts/
│   ├── phase1_motif_scan/
│   ├── phase2_cptac_tnbc/
│   │   ├── phase2c_tnbc_sample_list.py
│   │   └── ...
│   ├── phase3_structural_context/
│   └── phase4_network_proximity/
├── data/                      
├── outputs/                  
├── figures/                   

```


---

## Data Sources

- **CPTAC** (Clinical Proteomic Tumor Analysis Consortium) — BRCA phosphoproteomics, `umich` source, via the [`cptac`](https://pypi.org/project/cptac/) Python package
- **cBioPortal** — clinical/TNBC status, `brca_cptac_2020` study, via the cBioPortal REST API
- **UniProt** — AURKB reference sequence/domain annotation, proteome FASTA for gene-symbol mapping
- **OmniPath** — human protein–protein interaction network
- **PDB** — experimental structure lookup for Phase 3b
- *(Planned) TCGA* — CIN/aneuploidy signature correlation, not yet integrated

---

## Environment

Developed in a `conda` environment (`aurora-b`) on Python 3.x. Key dependencies observed across phase scripts:

```
pandas
requests
cptac
omnipath (or the omnipath REST endpoints directly)
freesasa
biopython          # if used for PDB parsing
```


---

## Known Issues / Open Questions

- RGL2_S736 structural coverage claim needs manual verification (residue-numbering mismatch, see Phase 3).
- 0/3 structurally-covered sites are surface-exposed — confirm this isn't a SASA artifact from incomplete/`UNK`-residue structures before treating it as a finding.
- Minor TNBC sample-count discrepancy: Phase 2c reports 28 TNBC-positive patients (clinical), Phase 2e reports 29 TNBC-positive samples matched in the phospho table — worth a quick trace.
- `fig4_top_combined_evidence` table overflow bug (see Phase 4b).
- No composite confidence score yet — current ranking is sequential filtering (motif AND TNBC-observed AND network-proximal), not a combined weighted score.

## Planned / Not Yet Implemented

- **True kinase–substrate structural modeling** (AURKB kinase domain + substrate peptide co-folding, e.g. via AlphaFold3/AlphaFold-Multimer) — current Phase 3 only checks experimental structure coverage and accessibility, not engagement.
- **TCGA CIN/aneuploidy signature correlation** — no TCGA data has been pulled yet.
- **Composite confidence score** combining motif, phosphoproteomic, structural, and network evidence into a single ranked metric.

---
