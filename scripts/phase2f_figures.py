"""
Aurora B Substrate Mapping 
Phase 2f: Summary Figures
------------------------------------------------------------------

Figures saved to figures/:
  1. fig1_pssm_score_distribution.png   — score distribution, known vs orphan
  2. fig2_known_vs_orphan_tnbc.png      — 2x2: known/orphan x observed/not
  3. fig3_pssm_score_by_tnbc_evidence.png — does TNBC evidence track with score?
  4. fig4_top20_orphan_candidates.png   — headline table: novel, TNBC-observed

Also writes top20_orphan_tnbc_candidates.csv alongside the PNGs.

Requires matplotlib 
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

DATA_DIR = Path("data")
CANDIDATES_CSV = DATA_DIR / "aurkb_candidates_tnbc_filtered.csv"
FIG_DIR = Path("figures")

# Pallete
COLOR_KNOWN = "#2E5266"          # deep blue-grey — known AURKB sites
COLOR_ORPHAN = "#F2A541"         # amber — orphan/candidate sites
COLOR_OBSERVED = "#3C8067"       # green — observed in TNBC
COLOR_NOT_OBSERVED = "#C4C4C4"   # grey — not observed

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
    print(f"  Known AURKB sites: {df['known_in_omnipath'].sum():,}")
    print(f"  Orphan candidates: {(~df['known_in_omnipath']).sum():,}")
    print(f"  Observed in TNBC:  {df['observed_in_tnbc'].sum():,}")
    return df


def fig1_pssm_distribution(df: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    known = df[df["known_in_omnipath"]]
    orphan = df[~df["known_in_omnipath"]]

    bins = 40
    ax.hist(orphan["pssm_score"], bins=bins, alpha=0.7,
            label=f"Orphan (n={len(orphan):,})", color=COLOR_ORPHAN,
            edgecolor="white", linewidth=0.3)
    ax.hist(known["pssm_score"], bins=bins, alpha=0.85,
            label=f"Known AURKB site (n={len(known):,})", color=COLOR_KNOWN,
            edgecolor="white", linewidth=0.3)

    ax.set_xlabel("PSSM score")
    ax.set_ylabel("Number of candidate sites")
    ax.set_title("Motif score distribution: known vs. orphan candidates")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig1_pssm_score_distribution.png")
    plt.close(fig)
    print("Saved fig1_pssm_score_distribution.png")


def fig2_known_vs_orphan_tnbc(df: pd.DataFrame):
    df = df.copy()
    df["status"] = df["known_in_omnipath"].map({True: "Known", False: "Orphan"})
    df["evidence"] = df["observed_in_tnbc"].map({True: "Observed in TNBC", False: "Not observed"})

    counts = df.groupby(["status", "evidence"]).size().unstack(fill_value=0)
    counts = counts.reindex(["Known", "Orphan"])
    counts = counts[["Observed in TNBC", "Not observed"]]

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    x = range(len(counts))
    width = 0.35
    ax.bar([i - width / 2 for i in x], counts["Observed in TNBC"], width,
           label="Observed in TNBC", color=COLOR_OBSERVED)
    ax.bar([i + width / 2 for i in x], counts["Not observed"], width,
           label="Not observed", color=COLOR_NOT_OBSERVED)

    for i, status in enumerate(counts.index):
        ax.text(i - width / 2, counts.loc[status, "Observed in TNBC"] + max(counts.values.max() * 0.01, 1),
                f"{counts.loc[status, 'Observed in TNBC']:,}", ha="center", fontsize=9)
        ax.text(i + width / 2, counts.loc[status, "Not observed"] + max(counts.values.max() * 0.01, 1),
                f"{counts.loc[status, 'Not observed']:,}", ha="center", fontsize=9)

    ax.set_xticks(list(x))
    ax.set_xticklabels(counts.index)
    ax.set_ylabel("Number of candidate sites")
    ax.set_title("TNBC phosphoproteomic evidence, by known/orphan status")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig2_known_vs_orphan_tnbc.png")
    plt.close(fig)
    print("Saved fig2_known_vs_orphan_tnbc.png")


def fig3_score_by_evidence(df: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    observed = df[df["observed_in_tnbc"]]["pssm_score"]
    not_observed = df[~df["observed_in_tnbc"]]["pssm_score"]

    ax.hist(not_observed, bins=40, density=True, alpha=0.6,
            label=f"Not observed (n={len(not_observed):,})", color=COLOR_NOT_OBSERVED,
            edgecolor="white", linewidth=0.3)
    ax.hist(observed, bins=40, density=True, alpha=0.75,
            label=f"Observed in TNBC (n={len(observed):,})", color=COLOR_OBSERVED,
            edgecolor="white", linewidth=0.3)

    ax.set_xlabel("PSSM score")
    ax.set_ylabel("Density")
    ax.set_title("Does TNBC phosphoproteomic evidence track with motif score?")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig3_pssm_score_by_tnbc_evidence.png")
    plt.close(fig)
    print("Saved fig3_pssm_score_by_tnbc_evidence.png")


def fig4_top_orphan_table(df: pd.DataFrame, n: int = 20):
    top = (
        df[(~df["known_in_omnipath"]) & (df["observed_in_tnbc"])]
        .sort_values("pssm_score", ascending=False)
        .head(n)
        .copy()
    )
    top["site"] = top["gene_symbol"].fillna(top["uniprot_id"]) + " " + top["residue"] + top["position"].astype(str)
    top_display = top[["site", "uniprot_id", "position", "pssm_score"]].reset_index(drop=True)
    top_display["pssm_score"] = top_display["pssm_score"].round(2)
    top_display.index = top_display.index + 1

    top_display.to_csv("top20_orphan_tnbc_candidates.csv")
    print(f"Saved top20_orphan_tnbc_candidates.csv ({len(top_display)} rows)")

    fig, ax = plt.subplots(figsize=(7, 0.4 * len(top_display) + 1))
    ax.axis("off")
    table = ax.table(
        cellText=top_display.astype(str).values,
        colLabels=["Site", "UniProt ID", "Position", "PSSM score"],
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.4)
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(weight="bold", color="white")
            cell.set_facecolor(COLOR_KNOWN)
        else:
            cell.set_facecolor("#F7F7F7" if row % 2 == 0 else "white")

    ax.set_title(f"Top {n} novel candidates: not previously known, "
                 f"observed phosphorylated in TNBC tissue", fontsize=11, pad=10)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "fig4_top20_orphan_candidates.png", bbox_inches="tight")
    plt.close(fig)
    print("Saved fig4_top20_orphan_candidates.png")


def main():
    FIG_DIR.mkdir(exist_ok=True)
    df = load_candidates()

    fig1_pssm_distribution(df)
    fig2_known_vs_orphan_tnbc(df)
    fig3_score_by_evidence(df)
    fig4_top_orphan_table(df)

    print(f"\nAll figures saved to {FIG_DIR.resolve()}")


if __name__ == "__main__":
    main()
