#!/usr/bin/env python3
"""Validate benchmark keys, coverage, evidence labels, ranks, and sequence invariants."""

import csv
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "benchmark/data"


def read(name):
    with (DATA / name).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    variants = read("variants_67.tsv")
    effects = read("experimental_tf_effects.tsv")
    predictions = read("published_predictions.tsv")
    assert len(variants) == 67
    ids = {x["benchmark_id"] for x in variants}
    assert len(ids) == 67
    assert all(len(x["ref_sequence"]) == len(x["alt_sequence"]) == 61 for x in variants)
    assert all(sum(a != b for a, b in zip(x["ref_sequence"], x["alt_sequence"])) == 1 for x in variants)
    assert all(x["benchmark_id"] in ids for x in effects + predictions)
    assert {x["evidence_class"] for x in effects} == {"Experimental observation"}
    assert {x["evidence_class"] for x in predictions} == {"Original published prediction"}
    coverage = Counter((x["method"], x["benchmark_id"]) for x in predictions)
    assert all(coverage[(method, variant)] > 0 for method in ("BayesPI-BAR", "sTRAP", "is-rSNP") for variant in ids)
    ranks = defaultdict(list)
    for row in predictions:
        ranks[(row["method"], row["benchmark_id"], row["predicted_change"])].append(int(row["rank_within_change"]))
    assert all(sorted(values) == list(range(1, max(values) + 1)) for values in ranks.values())
    print("validation=passed")
    print("variants", len(variants), Counter(x["source_group"] for x in variants))
    print("experimental_effect_rows", len(effects), Counter(x["effect_reported"] for x in effects))
    print("allele_order", Counter(x["allele_order_match"] for x in variants))
    print("prediction_rows", len(predictions), Counter(x["method"] for x in predictions))


if __name__ == "__main__":
    main()
