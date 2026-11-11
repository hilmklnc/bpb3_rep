#!/usr/bin/env python3
"""Verify that committed lightweight assets reproduce established results."""

from __future__ import annotations

import csv
import hashlib
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "benchmark/data"
REQUIRED_PATHS = (
    "docs/BPB3_HANDOFF.md",
    "benchmark/README.md",
    "benchmark/requirements.txt",
    "benchmark/schemas/bindfx_predictions.tsv",
    "benchmark/data/variants_67.tsv",
    "benchmark/data/experimental_tf_effects.tsv",
    "benchmark/data/published_predictions.tsv",
    "benchmark/data/repository_67_manifest.tsv",
    "benchmark/data/figure6_best_ranks.tsv",
    "benchmark/data/figure6_cumulative_accuracy.tsv",
    "benchmark/data/figure6_pairwise_ranks.tsv",
    "benchmark/data/bindfx_experimental_tf_overlap.tsv",
    "benchmark/data/bindfx_benchmark_strict_38.tsv",
    "benchmark/reports/source_inventory.md",
    "benchmark/reports/coverage_audit.md",
    "benchmark/reports/bindfx_overlap_audit.md",
    "benchmark/reports/figure6_reconstruction.md",
    "benchmark/reports/figure6_reconstruction.png",
    "benchmark/reports/original_workflow_verification.md",
    "benchmark/reports/handoff_file_inventory.tsv",
    "benchmark/scripts/build_benchmark_tables.py",
    "benchmark/scripts/evaluate_bindfx_predictions.py",
    "benchmark/scripts/reproduce_figure6.py",
    "benchmark/scripts/plot_figure6_reconstruction.py",
    "final_demo_data/demo1_67snps/in/snps_ref.fasta",
    "final_demo_data/demo1_67snps/in/snps_alt.fasta",
    "bpb3/bin/Linux/bayesPI_affinity",
    "bpb3/bin/Mac/bayesPI_affinity",
    "bpb3/script/script_high/other/stored_tfs.py",
)
EXPECTED_HASHES = {
    "final_demo_data/demo1_67snps/in/snps_ref.fasta": "4efd7bf4b7bb0b12d077e397d08354fbfeaff6b4b1bd39f64aaff79f86b090db",
    "final_demo_data/demo1_67snps/in/snps_alt.fasta": "65de7f50b4d932282d9c13f9facf1ea53d47668bba72974818b88dedcfc17874",
    "bpb3/script/script_high/other/stored_tfs.py": "ee031d9d60ef3c878467a82fd37a0edc336b66b52d52fee33c9f90782f19e43a",
}


def rows(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    missing = [path for path in REQUIRED_PATHS if not (ROOT / path).is_file()]
    assert not missing, f"Missing handoff assets: {missing}"
    for path, expected in EXPECTED_HASHES.items():
        assert sha256(ROOT / path) == expected, f"Checksum changed: {path}"

    variants = rows("variants_67.tsv")
    effects = rows("experimental_tf_effects.tsv")
    predictions = rows("published_predictions.tsv")
    ranks = rows("figure6_best_ranks.tsv")
    curves = rows("figure6_cumulative_accuracy.tsv")
    pairs = rows("figure6_pairwise_ranks.tsv")
    overlap = rows("bindfx_experimental_tf_overlap.tsv")
    compact = rows("bindfx_benchmark_strict_38.tsv")
    with (ROOT / "benchmark/reports/handoff_file_inventory.tsv").open(newline="") as handle:
        inventory = list(csv.DictReader(handle, delimiter="\t"))
    missing_inventory_assets = [
        row["path"] for row in inventory
        if row["included_in_git"] == "yes" and not (ROOT / row["path"]).exists()
    ]
    assert not missing_inventory_assets, f"Missing inventoried Git assets: {missing_inventory_assets}"
    assert len(list((ROOT / "final_demo_data/pwm").glob("*.mlp"))) == 1772

    assert len(variants) == 67 and len({row["benchmark_id"] for row in variants}) == 67
    assert len(effects) == 70
    assert Counter(row["effect_reported"] for row in effects) == {"Loss": 31, "Gain": 23, "Unknown": 16}
    assert len(predictions) == 16000
    assert Counter(row["method"] for row in predictions) == {"sTRAP": 6700, "BayesPI-BAR": 4743, "is-rSNP": 4557}
    assert len(ranks) == 354 and len(curves) == 300 and len(pairs) == 67
    assert len(overlap) == 70 and len(compact) == 38
    assert len({row["canonical_tf_concept"] for row in overlap if row["primary_benchmark_eligible"] == "yes"}) == 18
    assert sum(row["primary_benchmark_eligible"] == "yes" for row in overlap) == 39

    expected = {
        ("A_direction_not_considered", "BayesPI-BAR", "10"): ("47", "67"),
        ("A_direction_not_considered", "BayesPI-BAR", "50"): ("54", "67"),
        ("B_considering_direction", "BayesPI-BAR", "10"): ("29", "51"),
        ("B_considering_direction", "BayesPI-BAR", "50"): ("35", "51"),
    }
    by_key = {(row["panel"], row["method"], row["rank_cutoff"]): row for row in curves}
    for key, values in expected.items():
        row = by_key[key]
        assert (row["correct_count"], row["denominator"]) == values, f"Unexpected metric for {key}"

    png = (ROOT / "benchmark/reports/figure6_reconstruction.png").read_bytes()
    assert png.startswith(b"\x89PNG\r\n\x1a\n")
    print("handoff_validation=passed")
    print("assets", len(REQUIRED_PATHS), "inventory_rows", len(inventory),
          "variants", len(variants), "observations", len(effects),
          "published_predictions", len(predictions))
    print("figure6 BayesPI-BAR top10", "47/67 direction-agnostic", "29/51 direction-aware")


if __name__ == "__main__":
    main()
