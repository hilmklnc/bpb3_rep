#!/usr/bin/env python3
"""Plot the reconstructed four-panel structure of NAR 2015 Figure 6."""

import csv
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "benchmark/data"
OUTPUT = ROOT / "benchmark/reports/figure6_reconstruction.png"
COLORS = {"BayesPI-BAR": "#4c72b0", "sTRAP": "#55a868", "is-rSNP": "#c44e52"}


def read(name):
    with (DATA / name).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    curves = read("figure6_cumulative_accuracy.tsv")
    pairs = read("figure6_pairwise_ranks.tsv")
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))

    for ax, panel, title in (
        (axes[0, 0], "A_direction_not_considered", "Direction not considered"),
        (axes[0, 1], "B_considering_direction", "Considering direction"),
    ):
        for method in COLORS:
            rows = [x for x in curves if x["panel"] == panel and x["method"] == method]
            ax.step([int(x["rank_cutoff"]) for x in rows], [float(x["cumulative_accuracy"]) for x in rows],
                    where="post", label=method, color=COLORS[method])
        ax.set(xlim=(1, 50), ylim=(0, 1), xlabel="Rank", ylabel="Cumulative accuracy", title=title)
        ax.legend(frameon=False)

    bayes = [int(x["bayespibar_rank_capped_51"]) for x in pairs]
    comparisons = (
        (axes[1, 0], "sTRAP", [int(x["strap_rank_capped_51"]) for x in pairs], "0.015"),
        (axes[1, 1], "is-rSNP", [int(x["is_rsnp_rank_capped_51"]) for x in pairs], "0.00078"),
    )
    for ax, method, other, published_p in comparisons:
        bayes_better = sum(x < y for x, y in zip(bayes, other)) / len(bayes)
        other_better = sum(y < x for x, y in zip(bayes, other)) / len(bayes)
        ax.scatter(bayes, other, facecolors="none", edgecolors="#4c72b0")
        ax.plot([0, 51], [0, 51], color="#c44e52", linewidth=1)
        ax.text(4, 39, f"BayesPI-BAR better: {bayes_better:.0%}\npublished p={published_p}")
        ax.text(28, 8, f"{method} better: {other_better:.0%}")
        ax.set(xlim=(0, 52), ylim=(0, 52), xlabel="Rank by BayesPI-BAR", ylabel=f"Rank by {method}",
               title=f"BayesPI-BAR vs. {method}")

    for label, ax in zip("ABCD", axes.flat):
        ax.text(-0.15, 1.06, label, transform=ax.transAxes, fontweight="bold", fontsize=14)
    fig.tight_layout()
    fig.savefig(OUTPUT, dpi=200)
    print(OUTPUT)


if __name__ == "__main__":
    main()
