# bindFx / BayesPI-BAR benchmark reconstruction

This directory is additive: it does not modify the original prediction code or demo inputs.

Evidence classes must remain separate:

1. `experimental`: observations and metadata verified against the primary paper or supplement.
2. `published_prediction`: method outputs copied from a paper, supplement, or archived author result.
3. `reproduced_prediction`: outputs generated locally from documented code, inputs, and parameters.

## Main tables

- `data/variants_67.tsv`: 67 sequence pairs joined to the BayesPI-BAR author workbook metadata.
- `data/experimental_tf_effects.tsv`: one row per reported variant–TF/motif effect (70 rows; multi-TF variants are not collapsed).
- `data/published_predictions.tsv`: recovered author-archived BayesPI-BAR, sTRAP, and is-rSNP rankings (16,000 rows).
- `data/repository_67_manifest.tsv`: preliminary repository-only inventory retained for provenance.

Every table carries an evidence or verification field. An `Experimental observation` label means the observation is reported in the study's supporting workbook with a primary PubMed citation; it does not mean this reconstruction independently reviewed every primary assay.

Planned inputs not currently available in this workspace:

- bindFx 403-model inventory and aliases
- bindFx variant-level delta-B predictions
- manuscript/Mproject files
- primary-paper supplementary tables and archived method predictions

## Reproduce

```bash
bash benchmark/scripts/fetch_published_sources.sh
python benchmark/scripts/build_benchmark_tables.py \
  --mutations-archive benchmark/cache/BayesPI-BAR/Mutations.tar.gz \
  --rankings-archive benchmark/cache/BayesPI-BAR/rankings.tar.gz
python benchmark/scripts/validate_benchmark.py
```

The parser checks pinned SHA-256 hashes before reading either archive. It uses only the Python standard library.

The published ranking files rank positive and negative effects separately. `rank_within_change` must not be interpreted as a single global TF rank without an explicitly documented transformation.

## Evidence categories

- **Original published prediction**: recovered from the BayesPI-BAR study author archive.
- **Reproduced computational prediction**: generated locally; kept outside `published_predictions.tsv` unless a later table explicitly identifies it as reproduced.
- **Experimental observation**: reported TF/motif and effect from the author workbook, linked to a PubMed ID.
- **Unavailable or unverified**: searched but not recovered, or not checked against the underlying primary experiment.

See `reports/source_inventory.md` and `reports/coverage_audit.md` before analysis.


