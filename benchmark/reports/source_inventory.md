# Source inventory

Inventory date: 2026-11-05 (Europe/Istanbul).

## Recovered sources

| Resource | Classification | Location / identifier | Contents and use |
|---|---|---|---|
| BayesPI-BAR supporting mutation archive | Experimental observation | `junbaiw/BayesPI-BAR`, commit `316bf2a6bcd072336f4ed698776cad6b87f64db6`, `download/results/Mutations.tar.gz`; SHA-256 `9108a14efab606898e0c9c8d94fac1f842b717d79ac9294c919c83ac1805158d` | `All Mutations.xlsx`: 67 labels, genes, TFs/motifs, PubMed IDs, rsIDs, HGMD IDs, nucleotide changes, reported effects, locations, and disease/function annotations. The workbook is an author-curated secondary compilation; primary assays remain to be checked individually. |
| BayesPI-BAR comparative ranking archive | Original published prediction | `junbaiw/BayesPI-BAR`, same commit, `download/results/rankings.tar.gz`; SHA-256 `cf8d614dcd6e98be10af737a09820c6e18def0dd0be253d362d58755a226ae0a` | Variant-level rankings for BayesPI-BAR, sTRAP, and is-rSNP for all 67 variants. The repository page explicitly describes these as supporting ranking results from the study. Parsed into `data/published_predictions.tsv`. |
| Original 67 paired FASTAs | Experimental dataset representation / repository annotation | `final_demo_data/demo1_67snps/in/snps_ref.fasta`, `snps_alt.fasta` in `bpb3/bpb3` | 67 paired 61-bp sequences, one central substitution each. Joined to the mutation workbook in `data/variants_67.tsv`. |
| bpb3 TF/effect group lists | Repository annotation | `bpb3/script/script_high/other/stored_tfs.py` | Defines 4 Fredriksson, 29 HGMD, and 34 Andersen/Epstein identifiers plus broad TF aliases used by the demo accuracy checker. Not treated as primary experimental evidence. |
| Historical bpb3 full-run logs | Unavailable or unverified (outputs absent) | `final_demo/demo1_67snps/slurm-5779204.out`, `slurm-6020446.out` | Show a successful 67-sequence × 1,772-PWM, 10,000-shuffle bpb3 execution, but variant-level result tables are not committed. These logs do not supply predictions for this benchmark. |
| Local minimal bpb3 run | Reproduced computational prediction | `benchmark/reports/original_workflow_verification.md` | One variant × three unmodified PWMs at 100 shuffles. Confirms the pipeline path only; excluded from the published prediction table and from benchmarking. |

## Publications and method context

- BayesPI-BAR study: Nucleic Acids Research 43(21):e147 (2015), PubMed `26202972`. The recovered author archive is treated as the source of the original comparative predictions.
- BayesPI-BAR2: Scientific Reports (2017), DOI `10.1038/s41598-017-07226-4`. The public `junbaiw/BayesPI-BAR2` repository was searched.
- bpb3/BayesPI-BAR3: Frontiers in Genetics (2019), DOI `10.3389/fgene.2019.00282`, PubMed `31001324`.
- The 67-variant demo comments cite PubMed `19641089`, `25383969`, and `24077912` as collection sources; the recovered workbook provides a per-variant PubMed reference instead of assigning one citation to every dataset group.

## Search outcomes

### BayesPI-BAR

**Original published prediction — recovered.** All 67 variant files are present in the author ranking archive. The archive includes signed-side (`Positive`/`Negative`) rankings, model names, scores, method values, and reference/alternate median p-value fields.

### sTRAP

**Original published comparative prediction — recovered.** All 67 variant files are present in the same BayesPI-BAR study archive. This label means the predictions were published as the study's sTRAP comparison; it does not assert they came from an official sTRAP software repository.

`Wittelab/sTRAP_SNP_TF_Data` was also inspected. It contains a separate prostate-cancer GWAS sTRAP dataset and is not the source of these 67 predictions.

### is-rSNP

**Original published comparative prediction — recovered.** All 67 variant files are present in the BayesPI-BAR study archive. As for sTRAP, provenance is the comparative archive, not an independently recovered official is-rSNP run directory.

### BayesPI-BAR2 / bpb3

**Unavailable or unverified for this 67-variant dataset.** The public `junbaiw/BayesPI-BAR2` repository contains documentation and a skin-cancer cohort archive, but no 67-SNP ranking files were found. The `bpb3/bpb3` repository contains the paired inputs, scripts, group annotations, and successful-run logs, but not the resulting rankings. No rows from these methods were inserted into `published_predictions.tsv`.

## Access limitations

The managed environment permits GitHub access but not general publisher or NCBI browsing, and no connected browser-control surface was available. Therefore the author archives and repository histories were inspected directly, while journal-hosted supplementary pages and each cited primary experimental article were not independently reviewed. Those checks remain **Unavailable or unverified** in this reconstruction.

