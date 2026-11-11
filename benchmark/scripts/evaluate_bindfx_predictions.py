#!/usr/bin/env python3
"""Validate and evaluate externally generated bindFx predictions.

This script never runs bindFx and never assumes model coverage. It joins a
user-supplied prediction export to all 70 experimental TF-effect observations,
records unavailable/missing matches explicitly, and emits both observation-
and variant-level cumulative accuracy curves.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "benchmark/data"
REQUIRED_COLUMNS = {
    "benchmark_id", "chrom", "position_1based", "ref", "alt", "tf",
    "delta_b", "rank_score",
}
SUMMARY_CUTOFFS = (1, 5, 10, 50)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header in {path}")
        return list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def norm_tf(value: str) -> str:
    """Normalize only case and surrounding whitespace; do not infer families."""
    return value.strip().upper()


def finite_float(value: str, field: str, row_number: int) -> float:
    try:
        result = float(value)
    except ValueError as exc:
        raise ValueError(f"Row {row_number}: {field} is not numeric: {value!r}") from exc
    if not math.isfinite(result):
        raise ValueError(f"Row {row_number}: {field} must be finite")
    return result


def competition_ranks(rows: list[dict[str, object]]) -> None:
    rows.sort(key=lambda row: (-float(row["rank_score"]), str(row["tf_normalized"])))
    previous = None
    previous_rank = 0
    for index, row in enumerate(rows, 1):
        score = float(row["rank_score"])
        if previous is None or score != previous:
            previous_rank = index
            previous = score
        row["rank_within_variant"] = previous_rank


def prediction_direction(delta_b: float, positive_means: str) -> str:
    if delta_b == 0:
        return "Neutral"
    positive = delta_b > 0
    if positive_means == "gain":
        return "Gain" if positive else "Loss"
    return "Loss" if positive else "Gain"


def aggregate_duplicates(
    rows: list[dict[str, object]], policy: str
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["benchmark_id"]), str(row["tf_normalized"]))].append(row)

    output = []
    duplicates = {key: values for key, values in grouped.items() if len(values) > 1}
    if duplicates and policy == "error":
        examples = ", ".join(f"{key[0]}:{key[1]}({len(value)})" for key, value in list(duplicates.items())[:5])
        raise ValueError(
            "Duplicate variant-TF predictions found. Select --duplicate-policy explicitly "
            f"or export one row per TF. Examples: {examples}"
        )

    for key, values in grouped.items():
        if len(values) == 1:
            selected = dict(values[0])
        elif policy == "best-rank-score":
            selected = dict(max(values, key=lambda row: (float(row["rank_score"]), abs(float(row["delta_b"])))))
        elif policy == "mean":
            selected = dict(values[0])
            selected["delta_b"] = sum(float(row["delta_b"]) for row in values) / len(values)
            selected["rank_score"] = sum(float(row["rank_score"]) for row in values) / len(values)
            selected["model_id"] = ";".join(dict.fromkeys(str(row.get("model_id", "")) for row in values if row.get("model_id")))
        else:
            selected = dict(values[0])
        selected["source_row_count"] = len(values)
        output.append(selected)
    return output


def make_curves(
    evaluation_rows: list[dict[str, object]], variants: list[dict[str, str]], max_rank: int
) -> list[dict[str, object]]:
    curves: list[dict[str, object]] = []

    def add_curve(panel: str, unit: str, scope: str, ranks: list[int | None]) -> None:
        for cutoff in range(1, max_rank + 1):
            correct = sum(rank is not None and rank <= cutoff for rank in ranks)
            curves.append({
                "panel": panel, "evaluation_unit": unit, "denominator_scope": scope,
                "rank_cutoff": cutoff, "correct_count": correct,
                "denominator": len(ranks),
                "cumulative_accuracy": "" if not ranks else correct / len(ranks),
            })

    all_a = [int(row["rank"]) if row["rank"] != "" else None for row in evaluation_rows]
    covered_a = [int(row["rank"]) for row in evaluation_rows if row["rank"] != ""]
    known = [row for row in evaluation_rows if row["experimental_effect"] in {"Gain", "Loss"}]
    all_b = [int(row["rank"]) if row["direction_correct"] == "yes" else None for row in known]
    covered_b = [int(row["rank"]) if row["direction_correct"] == "yes" else None for row in known if row["rank"] != ""]
    add_curve("direction_not_considered", "experimental_tf_effect", "all", all_a)
    add_curve("direction_not_considered", "experimental_tf_effect", "covered", covered_a)
    add_curve("considering_direction", "experimental_tf_effect", "all_known_direction", all_b)
    add_curve("considering_direction", "experimental_tf_effect", "covered_known_direction", covered_b)

    by_variant: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in evaluation_rows:
        by_variant[str(row["benchmark_id"])].append(row)
    variant_ids = [row["benchmark_id"] for row in variants]
    variant_a_all: list[int | None] = []
    variant_a_covered: list[int] = []
    variant_b_all: list[int | None] = []
    variant_b_covered: list[int | None] = []
    for benchmark_id in variant_ids:
        observations = by_variant[benchmark_id]
        ranks = [int(row["rank"]) for row in observations if row["rank"] != ""]
        best = min(ranks) if ranks else None
        variant_a_all.append(best)
        if best is not None:
            variant_a_covered.append(best)

        known_rows = [row for row in observations if row["experimental_effect"] in {"Gain", "Loss"}]
        if known_rows:
            correct_ranks = [int(row["rank"]) for row in known_rows if row["direction_correct"] == "yes"]
            correct_best = min(correct_ranks) if correct_ranks else None
            variant_b_all.append(correct_best)
            if any(row["rank"] != "" for row in known_rows):
                variant_b_covered.append(correct_best)
    add_curve("direction_not_considered", "variant", "all", variant_a_all)
    add_curve("direction_not_considered", "variant", "covered", variant_a_covered)
    add_curve("considering_direction", "variant", "all_known_direction", variant_b_all)
    add_curve("considering_direction", "variant", "covered_known_direction", variant_b_covered)
    return curves


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "benchmark/results/bindfx")
    parser.add_argument("--positive-means", required=True, choices=("gain", "loss"),
                        help="Explicit biological interpretation of a positive delta_b")
    parser.add_argument("--mapping-mode", choices=("primary", "primary-and-family"), default="primary")
    parser.add_argument("--duplicate-policy", choices=("error", "best-rank-score", "mean"), default="error")
    parser.add_argument("--max-rank", type=int, default=50)
    args = parser.parse_args()
    if args.max_rank < 1:
        parser.error("--max-rank must be positive")

    with args.predictions.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required prediction columns: {sorted(missing)}")
        raw_predictions = list(reader)

    variants = read_tsv(DATA / "variants_67.tsv")
    effects = read_tsv(DATA / "experimental_tf_effects.tsv")
    mappings = read_tsv(DATA / "bindfx_experimental_tf_overlap.tsv")
    variant_by_id = {row["benchmark_id"]: row for row in variants}
    mapping_by_key = {(row["benchmark_id"], row["reported_tf_or_motif"]): row for row in mappings}

    input_audit: list[dict[str, object]] = []
    valid: list[dict[str, object]] = []
    for row_number, row in enumerate(raw_predictions, 2):
        benchmark_id = row["benchmark_id"].strip()
        tf = norm_tf(row["tf"])
        delta_b = finite_float(row["delta_b"], "delta_b", row_number)
        rank_score = finite_float(row["rank_score"], "rank_score", row_number)
        status = "valid"
        allele_status = "not_checked"
        variant = variant_by_id.get(benchmark_id)
        if variant is None:
            status = "unknown_variant"
        elif row["ref"].upper() == variant["ref_allele"] and row["alt"].upper() == variant["alt_allele"]:
            allele_status = "exact_fasta_orientation"
        elif row["ref"].upper() == variant["alt_allele"] and row["alt"].upper() == variant["ref_allele"]:
            status, allele_status = "allele_orientation_reversed", "reversed"
        else:
            status, allele_status = "allele_mismatch", "mismatch"
        if not tf:
            status = "missing_tf"
        audited = dict(row)
        audited.update({
            "input_row": row_number, "tf_normalized": tf, "input_status": status,
            "allele_status": allele_status,
        })
        input_audit.append(audited)
        if status == "valid":
            audited["delta_b"] = delta_b
            audited["rank_score"] = rank_score
            valid.append(audited)

    normalized = aggregate_duplicates(valid, args.duplicate_policy)
    by_variant_predictions: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in normalized:
        by_variant_predictions[str(row["benchmark_id"])].append(row)
    for rows in by_variant_predictions.values():
        competition_ranks(rows)
    prediction_lookup = {
        (str(row["benchmark_id"]), str(row["tf_normalized"])): row
        for rows in by_variant_predictions.values() for row in rows
    }

    evaluation_rows: list[dict[str, object]] = []
    for observation in effects:
        key = (observation["benchmark_id"], observation["tf_or_motif_reported"])
        mapping = mapping_by_key.get(key)
        if mapping is None:
            raise ValueError(f"Missing reviewed mapping for {key}")
        primary = mapping["primary_benchmark_eligible"] == "yes"
        family = mapping["secondary_family_analysis_eligible"] == "yes"
        allowed = primary or (args.mapping_mode == "primary-and-family" and family)
        candidates = [norm_tf(value) for value in mapping["bindfx_model_names"].split(";") if value] if allowed else []
        matches = [prediction_lookup[(observation["benchmark_id"], candidate)] for candidate in candidates
                   if (observation["benchmark_id"], candidate) in prediction_lookup]
        selected = min(matches, key=lambda row: (int(row["rank_within_variant"]), -float(row["rank_score"]))) if matches else None
        if not candidates:
            status = "no_eligible_model_under_mapping_mode"
        elif selected is None:
            status = "eligible_model_not_predicted"
        else:
            status = "matched"
        predicted = "" if selected is None else prediction_direction(float(selected["delta_b"]), args.positive_means)
        experimental = observation["effect_reported"]
        direction_correct = ""
        if selected is not None and experimental in {"Gain", "Loss"}:
            direction_correct = "yes" if predicted == experimental else "no"
        evaluation_rows.append({
            "benchmark_id": observation["benchmark_id"],
            "repository_label": observation["repository_label"],
            "observation_index": observation["observation_index"],
            "experimental_tf_or_motif": observation["tf_or_motif_reported"],
            "canonical_tf_concept": mapping["canonical_tf_concept"],
            "match_tier": mapping["match_tier"],
            "candidate_bindfx_tfs": ";".join(candidates),
            "experimental_effect": experimental,
            "match_status": status,
            "matched_bindfx_tf": "" if selected is None else selected["tf_normalized"],
            "delta_b": "" if selected is None else selected["delta_b"],
            "rank_score": "" if selected is None else selected["rank_score"],
            "rank": "" if selected is None else selected["rank_within_variant"],
            "predicted_direction": predicted,
            "direction_correct": direction_correct,
        })

    curves = make_curves(evaluation_rows, variants, args.max_rank)
    summary = [row for row in curves if int(row["rank_cutoff"]) in SUMMARY_CUTOFFS]
    output = args.output_dir
    audit_fields = list(input_audit[0]) if input_audit else list(raw_predictions[0]) + ["input_row", "tf_normalized", "input_status", "allele_status"] if raw_predictions else sorted(REQUIRED_COLUMNS) + ["input_row", "tf_normalized", "input_status", "allele_status"]
    normalized_fields = [
        "benchmark_id", "chrom", "position_1based", "ref", "alt", "tf", "tf_normalized",
        "delta_b", "rank_score", "rank_within_variant", "model_id", "p_value", "q_value",
        "source_row_count",
    ]
    write_tsv(output / "bindfx_input_audit.tsv", audit_fields, input_audit)
    write_tsv(output / "bindfx_predictions_normalized.tsv", normalized_fields,
              [row for rows in by_variant_predictions.values() for row in rows])
    write_tsv(output / "bindfx_observation_evaluation.tsv", list(evaluation_rows[0]), evaluation_rows)
    write_tsv(output / "bindfx_cumulative_accuracy.tsv", list(curves[0]), curves)
    write_tsv(output / "bindfx_evaluation_summary.tsv", list(summary[0]), summary)
    print(f"input_rows={len(raw_predictions)} valid_rows={len(valid)} normalized_predictions={len(normalized)}")
    print(f"experimental_observations={len(evaluation_rows)} matched={sum(row['match_status'] == 'matched' for row in evaluation_rows)}")
    print(f"outputs={output}")


if __name__ == "__main__":
    main()
