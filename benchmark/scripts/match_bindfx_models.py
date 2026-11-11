#!/usr/bin/env python3
"""Match reported experimental TF/motif labels to the bindFx model inventory."""

import argparse
import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
EFFECTS = ROOT / "benchmark/data/experimental_tf_effects.tsv"
OUTPUT = ROOT / "benchmark/data/bindfx_experimental_tf_overlap.tsv"

# Primary matches require the same gene or a well-established gene alias.
# Secondary matches are deliberately excluded from the primary denominator.
MAPPING = {
    "AP-2-ALPHA": ("TFAP2A", "established_alias", ["TFAP2A"], "AP-2alpha is an established name for TFAP2A"),
    "AP-2ALPHA": ("TFAP2A", "established_alias", ["TFAP2A"], "AP-2alpha is an established name for TFAP2A"),
    "AR": ("AR", "exact_gene", ["AR"], "exact gene symbol"),
    "C/EBP": ("C/EBP family", "family_or_motif", ["CEBPA", "CEBPB", "CEBPD", "CEBPG"], "reported family does not identify one CEBP gene"),
    "C/EBP-DELTA": ("CEBPD", "established_alias", ["CEBPD"], "C/EBP-delta is CEBPD"),
    "CACCC": ("CACCC-binding KLF family", "family_or_motif", ["KLF1", "KLF9", "KLF11", "KLF12"], "motif label does not identify one TF; available KLF models are candidates only"),
    "CCAAT": ("CCAAT/NF-Y motif", "family_or_motif", ["NFYA", "NFYB", "NFYC"], "motif is commonly bound by NF-Y but is not a unique TF label"),
    "ETS": ("ETS family", "family_or_motif", ["ETS1"], "reported ETS family label is not gene-specific"),
    "EGR1": ("EGR1", "established_alias", ["EGR1"], "case-normalized gene symbol"),
    "FOXA1": ("FOXA1", "exact_gene", ["FOXA1"], "exact gene symbol"),
    "GATA-1": ("GATA1", "established_alias", ["GATA1"], "GATA-1 is GATA1"),
    "HLF": ("HLF", "exact_gene", ["HLF"], "exact gene symbol"),
    "HNF-1": ("HNF1 family", "unavailable", [], "neither HNF1A nor HNF1B is available"),
    "HNF-4": ("HNF4 family", "family_or_motif", ["HNF4A", "HNF4G"], "reported HNF-4 label does not distinguish HNF4A and HNF4G"),
    "IRF1": ("IRF1", "exact_gene", ["IRF1"], "asterisk stripped where present; exact gene model"),
    "MSX1": ("MSX1", "unavailable", [], "MSX1 is absent; MSX2/MSX3 are not treated as equivalent"),
    "MYF-5": ("MYF5", "unavailable", [], "MYF5 model is absent"),
    "NF-Y": ("NF-Y complex", "family_or_motif", ["NFYA", "NFYB", "NFYC"], "heterotrimeric NF-Y complex; component models are available"),
    "NFKB": ("NF-kappaB family", "unavailable", [], "no NFKB1/NFKB2/RELA/RELB/REL model is available"),
    "NRF-2": ("NFE2L2", "established_alias", ["NFE2L2"], "NRF2 is the established alias of NFE2L2"),
    "OCT-1": ("POU2F1", "established_alias", ["POU2F1"], "OCT-1 is the established alias of POU2F1"),
    "SOX10": ("SOX10", "unavailable", [], "SOX10 model is absent"),
    "SP1": ("SP1", "exact_gene", ["SP1"], "case-normalized gene symbol; asterisk stripped where present"),
    "SPI1": ("SPI1", "exact_gene", ["SPI1"], "exact gene symbol"),
    "SRY": ("SRY", "unavailable", [], "SRY model is absent"),
    "STAT1": ("STAT1", "unavailable", [], "STAT1 model is absent"),
    "TBP": ("TBP", "exact_gene", ["TBP"], "exact gene symbol"),
    "TCF4": ("TCF7L2", "established_alias", ["TCF7L2"], "TCF4 in beta-catenin signaling is the historical alias of TCF7L2"),
    "TEL2": ("ETV7", "unavailable", [], "TEL2 is an alias of ETV7; ETV7 model is absent"),
    "USF": ("USF family", "family_or_motif", ["USF1", "USF2"], "generic USF label does not distinguish USF1 and USF2"),
    "USF-1": ("USF1", "established_alias", ["USF1"], "USF-1 is USF1; asterisk stripped"),
    "USF1": ("USF1", "exact_gene", ["USF1"], "exact gene symbol"),
    "USF2": ("USF2", "exact_gene", ["USF2"], "exact gene symbol"),
    "YY1": ("YY1", "exact_gene", ["YY1"], "exact gene symbol"),
    "BETA-CATENIN": ("CTNNB1", "unavailable", [], "beta-catenin is CTNNB1; model is absent and it is not a sequence-specific TF"),
    "C-FOS": ("FOS", "established_alias", ["FOS"], "c-FOS is FOS"),
    "C-MYB": ("MYB", "unavailable", [], "MYB is absent; MYBL2 is not treated as equivalent"),
    "C-REL": ("REL", "unavailable", [], "REL model is absent"),
    "P53": ("TP53", "unavailable", [], "TP53 model is absent"),
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("model_summary", type=Path)
    args = parser.parse_args()
    with args.model_summary.open(encoding="utf-8-sig", newline="") as handle:
        models = list(csv.DictReader(handle))
    if len(models) != 403:
        raise ValueError(f"Expected 403 bindFx models, found {len(models)}")
    by_name = {row["TF_Name"].upper(): row for row in models}
    if len(by_name) != 403:
        raise ValueError("bindFx TF_Name values are not unique")

    with EFFECTS.open(newline="") as handle:
        effects = list(csv.DictReader(handle, delimiter="\t"))
    output = []
    for row in effects:
        key = row["tf_or_motif_normalized"].upper()
        if key not in MAPPING:
            raise ValueError(f"Unreviewed experimental label: {row['tf_or_motif_normalized']}")
        concept, tier, candidates, rationale = MAPPING[key]
        missing = [x for x in candidates if x.upper() not in by_name]
        if missing:
            raise ValueError(f"Mapped bindFx model(s) absent: {missing}")
        selected = [by_name[x.upper()] for x in candidates]
        output.append({
            "benchmark_id": row["benchmark_id"], "repository_label": row["repository_label"],
            "reported_tf_or_motif": row["tf_or_motif_reported"], "experimental_effect": row["effect_reported"],
            "canonical_tf_concept": concept, "match_tier": tier,
            "bindfx_model_names": ";".join(x["TF_Name"] for x in selected),
            "bindfx_ensembl_ids": ";".join(x["Ensembl_Gene_ID"] for x in selected),
            "bindfx_model_count": len(selected),
            "primary_benchmark_eligible": "yes" if tier in {"exact_gene", "established_alias"} else "no",
            "secondary_family_analysis_eligible": "yes" if tier == "family_or_motif" else "no",
            "mapping_rationale": rationale,
        })

    fields = list(output[0])
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(output)

    concepts = {}
    for row in output:
        concepts[row["canonical_tf_concept"]] = row["match_tier"]
    primary_concepts = {x for x, tier in concepts.items() if tier in {"exact_gene", "established_alias"}}
    broad_concepts = {x for x, tier in concepts.items() if tier != "unavailable"}
    primary_variants = {x["benchmark_id"] for x in output if x["primary_benchmark_eligible"] == "yes"}
    broad_variants = {x["benchmark_id"] for x in output if x["match_tier"] != "unavailable"}
    print("rows", len(output), Counter(x["match_tier"] for x in output))
    print("unique_experimental_concepts", len(concepts))
    print("primary_comparable_concepts", len(primary_concepts), "variants", len(primary_variants))
    print("broad_comparable_concepts", len(broad_concepts), "variants", len(broad_variants))


if __name__ == "__main__":
    main()
