"""
Aurora B Substrate Mapping 
Phase 4b: Network Context Figures
------------------------------------------------------------------
Builds figures from Phase 4a's output (aurkb_candidates_network_context.csv):
which TNBC-validated candidates sit directly next to known mitotic
fidelity machinery, after excluding non-specific "hub" proteins (ATM
and friends) that would otherwise inflate the count.

Figures saved to figures/:
  1. fig1_hub_filtering_impact.png       -- why hub exclusion mattered
  2. fig2_reference_gene_connectivity.png -- which reference genes got excluded, and why
  3. fig3_partner_category_breakdown.png -- kinetochore vs SAC vs CPC, among specific hits
  4. fig4_top_combined_evidence.png      -- headline table: orphan + TNBC-observed +
                                             network-validated, your strongest candidates overall


"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

DATA_DIR = Path("data")
CANDIDATES_CSV = DATA_DIR / "aurkb_candidates_network_context.csv"
FIG_DIR = Path("figures")

# Snapshot from Phase 4a's last printed run 
REFERENCE_GENE_DEGREES = {
    "ATM": 1178, "SMARCA4": 428, "GTF2B": 381, "SIN3A": 368, "PLK1": 334,
    "HSF1": 273, "SMARCB1": 259, "AURKB": 258, "AURKA": 204, "CCNB1": 180,
    "ACTB": 163, "BIRC5": 158, "APC": 142, "KAT5": 126, "KAT2B": 125,
}
HUB_DEGREE_THRESHOLD = 300

COLOR_SPECIFIC = "#3C8067"       # green -trustworthy specific hit
COLOR_HUB_ONLY = "#F2A541"       # amber -hub-driven only, questionable
COLOR_NONE = "#C4C4C4"           # grey -no interactor evidence
COLOR_HUB_EXCLUDED = "#C0453D"   # red -excluded hub gene
COLOR_KEPT = "#2E5266"           # blue-grey - kept reference gene
COLOR_KNOWN = "#2E5266"

plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "font.size": 11,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def load_candidates() -> pd.DataFrame:
    df = pd.read_csv(CANDIDATES_CSV)
    print(f"Loaded {len(df)} candidates from {CANDIDATES_CSV}")
    return df


def fig1_hub_filtering_impact(df: pd.DataFrame):
    per_gene = df.drop_duplicates("gene_symbol_upper")
    n_total = len(per_gene)
    n_specific = int(per_gene["direct_mitotic_interactor_specific"].sum())
    n_any = int(per_gene["direct_mitotic_interactor_any"].sum())
    n_hub_only = n_any - n_specific
    n_none = n_total - n_any

    categories = ["No interactor\nevidence", "Hub-driven only\n(e.g. via ATM)", "Specific mitotic\ninteractor"]
    counts = [n_none, n_hub_only, n_specific]
    colors = [COLOR_NONE, COLOR_HUB_ONLY, COLOR_SPECIFIC]

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    bars = ax.bar(categories, counts, color=colors)
    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + n_total * 0.01,
                f"{count:,}", ha="center", fontsize=10)
    ax.set_ylabel("Number of candidate genes")
    ax.set_title(f"Why hub filtering matters (n={n_total:,} candidate genes)")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig1_hub_filtering_impact.png")
    plt.close(fig)
    print("Saved fig1_hub_filtering_impact.png")


def fig2_reference_gene_connectivity():
    genes = list(REFERENCE_GENE_DEGREES.keys())
    degrees = list(REFERENCE_GENE_DEGREES.values())
    colors = [COLOR_HUB_EXCLUDED if d > HUB_DEGREE_THRESHOLD else COLOR_KEPT for d in degrees]

    fig, ax = plt.subplots(figsize=(7, 5.5))
    y_pos = range(len(genes))
    ax.barh(y_pos, degrees, color=colors)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(genes)
    ax.invert_yaxis()
    ax.axvline(HUB_DEGREE_THRESHOLD, color="black", linestyle="--", linewidth=1)
    ax.text(HUB_DEGREE_THRESHOLD + 15, len(genes) - 0.5, f"exclusion cutoff\n({HUB_DEGREE_THRESHOLD} partners)",
            fontsize=9, va="top")
    ax.set_xlabel("Total interaction partners in OmniPath")
    ax.set_title("Reference gene connectivity: hub proteins skew the count")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig2_reference_gene_connectivity.png")
    plt.close(fig)
    print("Saved fig2_reference_gene_connectivity.png")


def fig3_partner_category_breakdown(df: pd.DataFrame):
    hits = df[df["direct_mitotic_interactor_specific"] == True].drop_duplicates("gene_symbol_upper")
    label_map = {
        "kinetochore": "Kinetochore",
        "mitotic spindle assembly checkpoint signaling": "Spindle checkpoint\n(SAC)",
        "chromosome passenger complex": "Chromosomal passenger\ncomplex (CPC)",
    }
    counts = {v: 0 for v in label_map.values()}
    for cats in hits["partner_categories_specific"].dropna():
        for c in str(cats).split(";"):
            if c in label_map:
                counts[label_map[c]] += 1

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    labels = list(counts.keys())
    values = list(counts.values())
    bars = ax.bar(labels, values, color=COLOR_SPECIFIC)
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(values) * 0.01,
                f"{v}", ha="center", fontsize=10)
    ax.set_ylabel(f"Candidate genes (of {len(hits)} specific hits)")
    ax.set_title("Which mitotic machinery category candidates connect to")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig3_partner_category_breakdown.png")
    plt.close(fig)
    print("Saved fig3_partner_category_breakdown.png")


def fig4_top_combined_evidence(df: pd.DataFrame, n: int = 15):
    top = (
        df[(~df["known_in_omnipath"]) & (df["observed_in_tnbc"]) & (df["direct_mitotic_interactor_specific"])]
        .sort_values("pssm_score", ascending=False)
        .head(n)
        .copy()
    )
    if top.empty:
        print("No candidates satisfy orphan + TNBC-observed + specific-network-hit simultaneously "
              "-- skipping fig4. Worth loosening one criterion if this list matters to you.")
        return

    top["site"] = top["gene_symbol"] + " " + top["residue"] + top["position"].astype(int).astype(str)

    def truncate_partners(s: str, max_items: int = 3) -> str:
        items = s.split(";")
        if len(items) <= max_items:
            return s
        return ";".join(items[:max_items]) + f" +{len(items) - max_items} more"

    top["interacting_partners_specific"] = top["interacting_partners_specific"].apply(truncate_partners)
    top_display = top[["site", "pssm_score", "interacting_partners_specific"]].reset_index(drop=True)
    top_display["pssm_score"] = top_display["pssm_score"].round(2)
    top_display.columns = ["Site", "PSSM score", "Mitotic partner(s)"]
    top_display.index = top_display.index + 1

    fig, ax = plt.subplots(figsize=(8.5, 0.4 * len(top_display) + 1))
    ax.axis("off")
    table = ax.table(
        cellText=top_display.astype(str).values,
        colLabels=list(top_display.columns),
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.4)
    table.auto_set_column_width(col=list(range(len(top_display.columns))))
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(weight="bold", color="white")
            cell.set_facecolor(COLOR_KNOWN)
        else:
            cell.set_facecolor("#F7F7F7" if row % 2 == 0 else "white")

    ax.set_title(f"Top {len(top_display)} candidates: novel + TNBC-observed + "
                 f"network-validated\n(motif, real tumor phosphorylation, AND mitotic-machinery proximity)",
                 fontsize=11, pad=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig4_top_combined_evidence.png", bbox_inches="tight")
    plt.close(fig)
    print("Saved fig4_top_combined_evidence.png")


def main():
    FIG_DIR.mkdir(exist_ok=True)
    df = load_candidates()

    fig1_hub_filtering_impact(df)
    fig2_reference_gene_connectivity()
    fig3_partner_category_breakdown(df)
    fig4_top_combined_evidence(df)

    print(f"\nAll figures saved to {FIG_DIR.resolve()}")


if __name__ == "__main__":
    main()
