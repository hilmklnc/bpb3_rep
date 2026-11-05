#!/usr/bin/env python3
"""Parse the 2015 BayesPI-BAR author archives into benchmark TSV tables."""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import io
import re
import tarfile
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
REF_FASTA = ROOT / "final_demo_data/demo1_67snps/in/snps_ref.fasta"
ALT_FASTA = ROOT / "final_demo_data/demo1_67snps/in/snps_alt.fasta"
ANNOTATIONS = ROOT / "bpb3/script/script_high/other/stored_tfs.py"
DATA = ROOT / "benchmark/data"
GROUPS = {"fredrikssonMuts": "fredriksson", "hgmdMuts": "hgmd", "andersonEpsteinMuts": "andersen_epstein"}
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main", "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
EXPECTED = {
    "Mutations.tar.gz": "9108a14efab606898e0c9c8d94fac1f842b717d79ac9294c919c83ac1805158d",
    "rankings.tar.gz": "cf8d614dcd6e98be10af737a09820c6e18def0dd0be253d362d58755a226ae0a",
}
ARCHIVE_COMMIT = "316bf2a6bcd072336f4ed698776cad6b87f64db6"
ARCHIVE_URL = "https://github.com/junbaiw/BayesPI-BAR"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:]
            records[name] = ""
        elif line.strip():
            if name is None:
                raise ValueError(f"Sequence before identifier in {path}")
            records[name] += line.strip().upper()
    return records


def read_groups() -> dict[str, str]:
    groups = {}
    for node in ast.parse(ANNOTATIONS.read_text()).body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            variable = node.targets[0].id
            if variable in GROUPS:
                for variant in ast.literal_eval(node.value):
                    groups[variant] = GROUPS[variable]
    return groups


def norm(value: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", value.upper())


def column_index(cell: str) -> int:
    letters = re.match(r"[A-Z]+", cell).group(0)
    result = 0
    for char in letters:
        result = result * 26 + ord(char) - 64
    return result - 1


def xlsx_rows(data: bytes) -> list[tuple[str, list[dict[str, str]]]]:
    output = []
    with zipfile.ZipFile(io.BytesIO(data)) as archive:
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            shared = ["".join(t.text or "" for t in item.iter(f"{{{NS['m']}}}t")) for item in root]
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rel_targets = {x.attrib["Id"]: x.attrib["Target"] for x in relationships}
        for sheet in workbook.find("m:sheets", NS):
            target = rel_targets[sheet.attrib[f"{{{NS['r']}}}id"]].lstrip("/")
            if not target.startswith("xl/"):
                target = "xl/" + target
            root = ET.fromstring(archive.read(target))
            rows = []
            for row in root.findall(".//m:row", NS):
                values = {}
                for cell in row.findall("m:c", NS):
                    value_node = cell.find("m:v", NS)
                    value = "" if value_node is None else value_node.text or ""
                    if cell.attrib.get("t") == "s" and value:
                        value = shared[int(value)]
                    values[column_index(cell.attrib["r"])] = value.strip()
                if values:
                    width = max(values) + 1
                    rows.append([values.get(i, "") for i in range(width)])
            header_index = next(i for i, row in enumerate(rows) if "Name" in row)
            header = rows[header_index]
            parsed = []
            for row in rows[header_index + 1:]:
                padded = row + [""] * (len(header) - len(row))
                record = dict(zip(header, padded))
                if record.get("Name"):
                    parsed.append(record)
            output.append((sheet.attrib["name"], parsed))
    return output


def member_bytes(archive_path: Path, suffix: str) -> bytes:
    with tarfile.open(archive_path, "r:gz") as archive:
        matches = [m for m in archive.getmembers() if m.name.endswith(suffix) and not Path(m.name).name.startswith("._")]
        if len(matches) != 1:
            raise ValueError(f"Expected one {suffix} in {archive_path}, found {len(matches)}")
        return archive.extractfile(matches[0]).read()


def verify_archive(path: Path) -> None:
    observed = sha256_bytes(path.read_bytes())
    expected = EXPECTED[path.name]
    if observed != expected:
        raise ValueError(f"SHA-256 mismatch for {path}: {observed} != {expected}")


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mutations-archive", required=True, type=Path)
    parser.add_argument("--rankings-archive", required=True, type=Path)
    args = parser.parse_args()
    verify_archive(args.mutations_archive)
    verify_archive(args.rankings_archive)

    ref, alt, groups = read_fasta(REF_FASTA), read_fasta(ALT_FASTA), read_groups()
    if set(ref) != set(alt) or set(ref) != set(groups):
        raise ValueError("FASTA and group identifier sets differ")
    canonical = {norm(x): x for x in ref}

    workbook = member_bytes(args.mutations_archive, "All Mutations.xlsx")
    annotations = {}
    sheet_by_variant = {}
    for sheet, rows in xlsx_rows(workbook):
        for row in rows:
            key = norm(row["Name"])
            if key not in canonical:
                continue
            variant = canonical[key]
            if variant in annotations:
                raise ValueError(f"Duplicate workbook row: {variant}")
            annotations[variant] = row
            sheet_by_variant[variant] = sheet
    if set(annotations) != set(ref):
        raise ValueError(f"Workbook coverage differs: missing={set(ref)-set(annotations)}")

    variants = []
    observations = []
    for index, variant in enumerate(ref, 1):
        row = annotations[variant]
        differences = [i for i, pair in enumerate(zip(ref[variant], alt[variant])) if pair[0] != pair[1]]
        if len(differences) != 1:
            raise ValueError(f"Expected one difference for {variant}")
        offset = differences[0]
        seq_change = f"{ref[variant][offset]} -> {alt[variant][offset]}"
        published_change = row.get("Nucleotide change", "").replace("→", "->")
        cleaned_pub = re.sub(r"\s+", " ", published_change).strip()
        reverse_change = f"{alt[variant][offset]} -> {ref[variant][offset]}"
        match = "exact" if cleaned_pub == seq_change else "reversed" if cleaned_pub == reverse_change else "mismatch"
        benchmark_id = f"BPB3_67_{index:03d}"
        variants.append({
            "benchmark_id": benchmark_id, "repository_label": variant, "source_group": groups[variant],
            "workbook_sheet": sheet_by_variant[variant], "gene": row.get("Gene", ""), "rsid": row.get("dbSNP ID", "").strip(),
            "hgmd_id": row.get("HGMD ID", ""), "reference_pubmed_id": row.get("Reference PubMed ID", ""),
            "published_nucleotide_change": cleaned_pub, "sequence_nucleotide_change": seq_change,
            "allele_order_match": match, "location": row.get("Location", ""), "disease_or_function": row.get("Disease / function", ""),
            "ref_sequence": ref[variant], "alt_sequence": alt[variant], "sequence_length": len(ref[variant]),
            "variant_offset_0based": offset, "ref_allele": ref[variant][offset], "alt_allele": alt[variant][offset],
            "genome_build": "", "chrom": "", "position_1based": "",
            "annotation_evidence_class": "Experimental observation",
            "annotation_verification": "author_supporting_workbook;primary_article_not_independently_checked",
        })
        tf_raw = row.get("TF or motif", "")
        for observation_index, tf in enumerate((x.strip() for x in tf_raw.split(",") if x.strip()), 1):
            observations.append({
                "benchmark_id": benchmark_id, "repository_label": variant, "observation_index": observation_index,
                "tf_or_motif_reported": tf, "tf_or_motif_normalized": tf.rstrip("*").strip(),
                "effect_reported": row.get("TF-binding site affected", ""), "gene": row.get("Gene", ""),
                "reference_pubmed_id": row.get("Reference PubMed ID", ""), "assay": "",
                "evidence_class": "Experimental observation",
                "verification_status": "author_supporting_workbook;primary_article_and_assay_not_independently_checked",
                "source_archive_sha256": EXPECTED["Mutations.tar.gz"],
            })

    prediction_rows = []
    with tarfile.open(args.rankings_archive, "r:gz") as archive:
        members = [m for m in archive.getmembers() if m.isfile() and m.name.endswith(".tsv") and "/._" not in m.name]
        method_map = {"BayesPI-BAR": "BayesPI-BAR", "sTRAP": "sTRAP", "is-rSNP": "is-rSNP"}
        by_norm = {norm(v["repository_label"]): v for v in variants}
        for member in sorted(members, key=lambda x: x.name):
            parts = Path(member.name).parts
            method_folder = next((x for x in parts if x in method_map), None)
            if method_folder is None:
                continue
            filename_variant = Path(member.name).stem
            variant_row = by_norm.get(norm(filename_variant))
            if variant_row is None:
                raise ValueError(f"Cannot map prediction file: {member.name}")
            text = io.TextIOWrapper(archive.extractfile(member), encoding="utf-8-sig")
            for prediction in csv.DictReader(text, delimiter="\t"):
                prediction_rows.append({
                    "benchmark_id": variant_row["benchmark_id"], "repository_label": variant_row["repository_label"],
                    "method": method_map[method_folder], "predicted_change": prediction.get("change", ""),
                    "rank_within_change": prediction.get("rank", ""), "model_or_tf": prediction.get("fileName", ""),
                    "score": prediction.get("score", ""), "mean_ddba_or_method_value": prediction.get("mean.ddba", ""),
                    "median_p_reference": prediction.get("med.p.ref", ""), "median_p_alternate": prediction.get("med.p.alt", ""),
                    "evidence_class": "Original published prediction",
                    "source_description": "BayesPI-BAR study author-archived comparative ranking",
                    "source_repository": ARCHIVE_URL, "source_commit": ARCHIVE_COMMIT,
                    "source_archive_sha256": EXPECTED["rankings.tar.gz"], "availability_status": "recovered",
                })

    write_tsv(DATA / "variants_67.tsv", list(variants[0]), variants)
    write_tsv(DATA / "experimental_tf_effects.tsv", list(observations[0]), observations)
    write_tsv(DATA / "published_predictions.tsv", list(prediction_rows[0]), prediction_rows)
    print(f"variants={len(variants)} observations={len(observations)} predictions={len(prediction_rows)}")


if __name__ == "__main__":
    main()

