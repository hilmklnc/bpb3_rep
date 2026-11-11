#!/usr/bin/env python3
"""Reconstruct NAR 2015 Figure 6 ranks from the archived prediction tables.

The matching aliases and reduction logic are transcribed from the original
BayesPI-BAR demo's code/r/accuracy.R.
"""

import csv
import math
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "benchmark/data"
REPORTS = ROOT / "benchmark/reports"
METHODS = ("BayesPI-BAR", "sTRAP", "is-rSNP")

# Original 2015 accuracy.R aliases, restricted to labels present in the 67 set.
ALIASES = {
    "NFKB": ["NFKB", "NF-kappaB"],
    "NRF2": ["NRF2", "NFE2L2", "NFE2", "NRF-2"],
    "NFY": ["NFY", "NF-Y"], "YY1": ["YY1"],
    "USF1": ["USF"], "USF-1": ["USF"], "USF2": ["USF"], "USF": ["USF"],
    "C-MYB": ["MYB"], "IRF1": ["IRF1", "IRF-1"], "HLF": ["HLF"],
    "P53": ["P53"], "HNF-1": ["HNF1"], "SPI1": ["SPI1"],
    "AP-2ALPHA": ["AP-2", "AP2", "TFAP2", "AP2TF"],
    "AP-2-ALPHA": ["AP-2", "AP2", "TFAP2", "AP2TF"],
    "CCAAT": ["CRP1", "CCAAT", "CEBPZ", "CEBPA", "CEBPG", "CEBPB", "CEBP", "C/EBP"],
    "CACCC": ["CACCC", "EKLF", "KLF1", "KLF-1", "KLF"],
    "SP1": ["SP1", "SP-1", "TFSP1"], "SP-1": ["SP1", "SP-1", "TFSP1"],
    "MSX1": ["MSX1", "MSX-1", "HOX7", "ECTD3", "HYD1", "STHAG1"],
    "GATA-1": ["GATA"], "GATA-A": ["GATA"],
    "HNF-4": ["HNF4a", "HNF4g", "HNF4", "TCF14", "TCF-14", "MODY", "MODY1"],
    "CEBP": ["CRP1", "CEBP", "CCAAT", "C/EBP"],
    "CEBP-DELTA": ["CRP1", "CEBP", "CCAAT", "C/EBP"],
    "AR": ["_AR_", "^ar$"], "C-FOS": ["AP1_", "AP-1"],
    "C-REL": ["REL_", "c-Rel"], "EGR1": ["EGR1", "Egr-1"],
    "ETS": ["ETS"], "FOXA1": ["FOXA"], "MYF-5": ["MYF"],
    "OCT-1": ["POU2F2", "OCT"], "SOX10": ["SOX10"], "SRY": ["SRY"],
    "STAT1": ["STAT"], "TBP": ["TATA", "TBP"],
    "TCF4": ["TCF4", "TCF-4"], "TEL2": ["ETV7"],
}


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_targets(repository_label):
    # The published archive sanitized filenames to uppercase and removed '/'.
    label = repository_label.upper().replace("/", "")
    fields = label.split("_")[1:]
    targets, directions = [], []
    for field in fields:
        parts = field.split("|")
        targets.append(parts[0])
        directions.append(parts[1] if len(parts) > 1 else "")
    for index, direction in enumerate(directions):
        if not direction:
            directions[index] = directions[-1]
    keep = [
        index for index, target in enumerate(targets)
        if target not in {str(x) for x in range(1, 11)} and target in ALIASES
    ]
    return [(targets[index], directions[index]) for index in keep]


def matches_model(target, model):
    return any(re.search(pattern, model, flags=re.IGNORECASE) for pattern in ALIASES[target])


def best_match(rows, targets, directed):
    matches = []
    for row in rows:
        for target, expected_direction in targets:
            if matches_model(target, row["model_or_tf"]):
                direction_ok = (
                    not directed or
                    expected_direction != "UNKNOWN" and
                    ((expected_direction == "GAIN") == (row["predicted_change"] == "Positive"))
                )
                if direction_ok:
                    matches.append(row)
                # Original R code stops at the first implicated TF matching a model.
                break
    if not matches:
        return None
    best_rank = min(int(row["rank_within_change"]) for row in matches)
    return next(row for row in matches if int(row["rank_within_change"]) == best_rank)


def write_tsv(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def main():
    variants = read_tsv(DATA / "variants_67.tsv")
    predictions = read_tsv(DATA / "published_predictions.tsv")
    by_variant_method = defaultdict(list)
    for row in predictions:
        by_variant_method[(row["benchmark_id"], row["method"])].append(row)

    rank_rows = []
    curve_rows = []
    ranks_by_panel_method = defaultdict(list)
    for panel, directed in (("A_direction_not_considered", False), ("B_considering_direction", True)):
        for variant in variants:
            targets = parse_targets(variant["repository_label"])
            if directed and not any(direction in {"GAIN", "LOSS"} for _, direction in targets):
                continue
            for method in METHODS:
                best = best_match(by_variant_method[(variant["benchmark_id"], method)], targets, directed)
                rank = int(best["rank_within_change"]) if best else None
                capped = rank if rank is not None and rank <= 50 else 51
                ranks_by_panel_method[(panel, method)].append(capped)
                rank_rows.append({
                    "panel": panel, "benchmark_id": variant["benchmark_id"],
                    "repository_label": variant["repository_label"], "method": method,
                    "experimental_targets": ";".join(target for target, _ in targets),
                    "experimental_directions": ";".join(direction for _, direction in targets),
                    "best_rank": "" if rank is None else rank, "rank_capped_51": capped,
                    "matched_model_or_tf": "NA" if best is None else best["model_or_tf"],
                    "matched_predicted_change": "NA" if best is None else best["predicted_change"],
                })

        for method in METHODS:
            ranks = ranks_by_panel_method[(panel, method)]
            for cutoff in range(1, 51):
                correct = sum(rank <= cutoff for rank in ranks)
                curve_rows.append({
                    "panel": panel, "method": method, "rank_cutoff": cutoff,
                    "correct_count": correct, "denominator": len(ranks),
                    "cumulative_accuracy": correct / len(ranks),
                })

    write_tsv(DATA / "figure6_best_ranks.tsv", list(rank_rows[0]), rank_rows)
    write_tsv(DATA / "figure6_cumulative_accuracy.tsv", list(curve_rows[0]), curve_rows)

    panel_a = {method: ranks_by_panel_method[("A_direction_not_considered", method)] for method in METHODS}
    pair_rows = []
    for index, variant in enumerate(variants):
        pair_rows.append({
            "benchmark_id": variant["benchmark_id"], "repository_label": variant["repository_label"],
            "bayespibar_rank_capped_51": panel_a["BayesPI-BAR"][index],
            "strap_rank_capped_51": panel_a["sTRAP"][index],
            "is_rsnp_rank_capped_51": panel_a["is-rSNP"][index],
        })
    write_tsv(DATA / "figure6_pairwise_ranks.tsv", list(pair_rows[0]), pair_rows)

    for panel in ("A_direction_not_considered", "B_considering_direction"):
        print(panel)
        for method in METHODS:
            ranks = ranks_by_panel_method[(panel, method)]
            print(method, "denominator", len(ranks), "top10", sum(x <= 10 for x in ranks), "top50", sum(x <= 50 for x in ranks))
    for other in ("sTRAP", "is-rSNP"):
        bayes = panel_a["BayesPI-BAR"]; comparison = panel_a[other]
        print(other, "BayesPI-BAR_better", sum(x < y for x, y in zip(bayes, comparison)), "/", len(bayes),
              "other_better", sum(y < x for x, y in zip(bayes, comparison)), "/", len(bayes))


if __name__ == "__main__":
    main()
