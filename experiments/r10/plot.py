"""Figures from results/bins.tsv (pooled over strains). No mappings are run."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import pandas as pd

import r10lib as r
from report import CONFIGS, TECH, pooled_bins

STYLE = {("pb", "default"): ("tab:blue", "--", "o"), ("pb", "best_ont"): ("tab:blue", "-", "s"),
         ("ont", "default"): ("tab:red", "--", "o"), ("ont", "best_ont"): ("tab:red", "-", "s")}
MIN_N = 30  # do not draw bins with fewer truth reads


def bin_figure(bins, kind, labels, xlabel, fname):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)
    for ax, (col, title) in zip(axes, [("sens", "Sensitivity (any MAPQ)"), ("sens60", "Sensitivity (mapquik MAPQ 60 only)")]):
        for tech, setting in CONFIGS:
            p = pooled_bins(bins, kind, r.STRAINS, tech, setting)
            ok = p.n_truth >= MIN_N
            color, ls, mk = STYLE[(tech, setting)]
            ax.plot([i for i, o in enumerate(ok) if o], 100 * p[col][ok], color=color, ls=ls, marker=mk, label=f"{TECH[tech]} {setting}")
        ax.set_xticks(range(len(labels)))
        n = pooled_bins(bins, kind, r.STRAINS, "ont", "default").n_truth
        ax.set_xticklabels([f"{l}\nn={n[l]:,.0f}" for l in labels], fontsize=8)
        ax.set_xlabel(xlabel + "  (n = ONT truth reads)")
        ax.set_title(title)
        ax.set_ylim(40 if col == "sens" else 0, 101)
        ax.grid(alpha=.3)
    axes[0].set_ylabel("% of minimap2 MAPQ>=30 reads (9 strains pooled)")
    axes[0].legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(r.RESULTS / fname, dpi=130)
    plt.close(fig)


def main():
    bins = pd.read_csv(r.RESULTS / "bins.tsv", sep="\t")
    bin_figure(bins, "id_bin", r.ID_LABELS, "minimap2 read identity bin (%)", "sens_by_identity.png")
    bin_figure(bins, "len_bin", r.LEN_LABELS, "read length bin", "sens_by_length.png")
    print("wrote", r.RESULTS / "sens_by_identity.png", r.RESULTS / "sens_by_length.png")


if __name__ == "__main__":
    main()
