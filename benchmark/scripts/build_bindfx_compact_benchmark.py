#!/usr/bin/env python3
"""Build one compact row per strict bindFx-overlap SNP."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "benchmark/data"
OUTPUT = DATA / "bindfx_benchmark_strict_38.tsv"


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def unique_join(values):
    return ";".join(dict.fromkeys(value for value in values if value))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("model_summary", type=Path)
    args = parser.parse_args()

    variants = {row["benchmark_id"]: row for row in read_tsv(DATA / "variants_67.tsv")}
    overlap = [
        row for row in read_tsv(DATA / "bindfx_experimental_tf_overlap.tsv")
        if row["primary_benchmark_eligible"] == "yes"
    ]
    with args.model_summary.open(encoding="utf-8-sig", newline="") as handle:
        model_rows = list(csv.DictReader(handle))
    model_file = {row["TF_Name"]: row["TF Files (ChIP+PBM)"] for row in model_rows}

    grouped = defaultdict(list)
    for row in overlap:
        grouped[row["benchmark_id"]].append(row)

    output = []
    for benchmark_id in sorted(grouped):
        variant = variants[benchmark_id]
        matches = grouped[benchmark_id]
        models = unique_join(
            model for row in matches for model in row["bindfx_model_names"].split(";") if model
        )
        model_files = unique_join(model_file[model] for model in models.split(";") if model)
        effects = unique_join(row["experimental_effect"] for row in matches)
        output.append({
            "benchmark_id": benchmark_id,
            "repository_label": variant["repository_label"],
            "source_group": variant["source_group"],
            "gene": variant["gene"],
            "rsid": variant["rsid"],
            "hgmd_id": variant["hgmd_id"],
            "ref_allele": variant["ref_allele"],
            "alt_allele": variant["alt_allele"],
            "variant_offset_0based": variant["variant_offset_0based"],
            "ref_sequence": variant["ref_sequence"],
            "alt_sequence": variant["alt_sequence"],
            "experimental_tf": unique_join(row["reported_tf_or_motif"] for row in matches),
            "bindfx_model": models,
            "bindfx_model_file": model_files,
            "experimental_effect": effects,
            "direction_evaluable": "yes" if effects and "Unknown" not in effects.split(";") else "no",
            "allele_order_match": variant["allele_order_match"],
            "reference_pubmed_id": variant["reference_pubmed_id"],
        })

    if len(output) != 38:
        raise ValueError(f"Expected 38 strict-overlap SNPs, found {len(output)}")
    fields = list(output[0])
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    print(f"wrote {len(output)} SNP rows to {OUTPUT}")
    print("direction_evaluable", sum(row["direction_evaluable"] == "yes" for row in output))
    print("allele_orientation_flags", sum(row["allele_order_match"] != "exact" for row in output))


if __name__ == "__main__":
    main()
