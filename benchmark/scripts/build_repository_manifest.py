#!/usr/bin/env python3
"""Build a provenance-preserving inventory from the committed bpb3 67-SNP demo."""

from __future__ import annotations

import ast
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REF_FASTA = ROOT / "final_demo_data/demo1_67snps/in/snps_ref.fasta"
ALT_FASTA = ROOT / "final_demo_data/demo1_67snps/in/snps_alt.fasta"
ANNOTATIONS = ROOT / "bpb3/script/script_high/other/stored_tfs.py"
OUTPUT = ROOT / "benchmark/data/repository_67_manifest.tsv"

GROUP_VARIABLES = {
    "fredrikssonMuts": "fredriksson",
    "hgmdMuts": "hgmd",
    "andersonEpsteinMuts": "andersen_epstein",
}


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if line.startswith(">"):
            name = line[1:]
            if not name or name in records:
                raise ValueError(f"Invalid or duplicate FASTA ID in {path}: {name!r}")
            records[name] = ""
        elif line:
            if name is None:
                raise ValueError(f"Sequence before first FASTA ID in {path}")
            records[name] += line.upper()
    return records


def read_groups(path: Path) -> dict[str, str]:
    tree = ast.parse(path.read_text(), filename=str(path))
    groups: dict[str, str] = {}
    for node in tree.body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Name) or target.id not in GROUP_VARIABLES:
            continue
        for variant_id in ast.literal_eval(node.value):
            if variant_id in groups:
                raise ValueError(f"Variant assigned to multiple groups: {variant_id}")
            groups[variant_id] = GROUP_VARIABLES[target.id]
    return groups


def label_direction(label: str) -> str:
    effects = []
    for token in label.split("|")[1:]:
        head = token.split("_", 1)[0].upper()
        if head in {"GAIN", "LOSS", "UNKNOWN"}:
            effects.append(head.lower())
    unique = sorted(set(effects))
    return ";".join(unique) if unique else "unparsed"


def main() -> None:
    ref = read_fasta(REF_FASTA)
    alt = read_fasta(ALT_FASTA)
    groups = read_groups(ANNOTATIONS)
    if set(ref) != set(alt) or set(ref) != set(groups):
        raise ValueError("Reference FASTA, alternate FASTA, and annotation IDs differ")

    fields = [
        "benchmark_id", "source_group", "repository_label", "label_direction",
        "ref_sequence", "alt_sequence", "sequence_length", "variant_offset_0based",
        "ref_allele", "alt_allele", "difference_count", "rsid", "genome_build",
        "chrom", "position_1based", "experimental_tf_curated", "experimental_effect_curated",
        "assay", "primary_reference", "evidence_class", "curation_status",
        "ref_fasta_sha256", "alt_fasta_sha256",
    ]
    import hashlib
    ref_hash = hashlib.sha256(REF_FASTA.read_bytes()).hexdigest()
    alt_hash = hashlib.sha256(ALT_FASTA.read_bytes()).hexdigest()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for index, variant_id in enumerate(ref, start=1):
            differences = [i for i, (r, a) in enumerate(zip(ref[variant_id], alt[variant_id])) if r != a]
            if len(ref[variant_id]) != len(alt[variant_id]):
                raise ValueError(f"Sequence lengths differ for {variant_id}")
            writer.writerow({
                "benchmark_id": f"BPB3_67_{index:03d}",
                "source_group": groups[variant_id],
                "repository_label": variant_id,
                "label_direction": label_direction(variant_id),
                "ref_sequence": ref[variant_id],
                "alt_sequence": alt[variant_id],
                "sequence_length": len(ref[variant_id]),
                "variant_offset_0based": differences[0] if len(differences) == 1 else "",
                "ref_allele": ref[variant_id][differences[0]] if len(differences) == 1 else "",
                "alt_allele": alt[variant_id][differences[0]] if len(differences) == 1 else "",
                "difference_count": len(differences),
                "evidence_class": "repository_annotation",
                "curation_status": "primary_source_not_yet_verified",
                "ref_fasta_sha256": ref_hash,
                "alt_fasta_sha256": alt_hash,
            })
    print(f"wrote {len(ref)} records to {OUTPUT}")


if __name__ == "__main__":
    main()
