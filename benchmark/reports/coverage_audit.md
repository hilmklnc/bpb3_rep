# Coverage and validation audit

## Dataset coverage

| Item | Validated | Missing / unresolved |
|---|---:|---:|
| Variant sequence pairs | 67 | 0 |
| Single central substitutions | 67 | 0 |
| Per-variant PubMed IDs | 67 | 0 |
| dbSNP rsIDs | 47 | 20 |
| HGMD IDs | 57 | 10 |
| Genome build | 0 | 67 |
| Chromosome and genomic coordinate | 0 | 67 |
| Explicit assay type | 0 | 67 |
| Author-workbook TF/motif-effect rows | 70 | primary assay verification pending for all 70 |

The 70 observations represent 67 variants because three labels report two implicated factors: `CHGB_SRY|LOSS_YY1|LOSS`, `FOXE1_USF1|GAIN_USF2|GAIN`, and `MYC_TCF4|GAIN_beta-catenin|GAIN`.

Reported effect rows: 23 gain, 31 loss, and 16 unknown. Direction-agreement analyses therefore have a maximum author-workbook denominator of 54 observations before TF-model coverage or other exclusions; the 16 unknown observations must not be counted as direction successes or failures.

## Allele orientation audit

The workbook nucleotide change matches the reference-to-alternate FASTA orientation for 65/67 variants. Two are reversed:

- `ITGA2_SP1|LOSS`: workbook `C -> T`; FASTA reference-to-alternate `T -> C`.
- `F9_C/EBP|LOSS`: workbook `G -> C`; FASTA reference-to-alternate `C -> G`.

`variants_67.tsv` preserves both representations and records `allele_order_match=reversed`. Downstream direction comparisons must establish which allele is experimentally defined as reference/baseline before scoring these two records.

## Published prediction coverage

| Method | Evidence class | Variants | Prediction rows | Distinct model labels | Rows per variant |
|---|---|---:|---:|---:|---:|
| BayesPI-BAR | Original published prediction | 67/67 | 4,743 | 1,374 | 25–100 |
| sTRAP | Original published comparative prediction | 67/67 | 6,700 | 1,673 | 100 |
| is-rSNP | Original published comparative prediction | 67/67 | 4,557 | 644 | 41–122 |
| BayesPI-BAR2 | Unavailable or unverified | 0/67 | 0 | — | — |
| bpb3/BayesPI-BAR3 | Unavailable or unverified | 0/67 | 0 | — | — |

The archived tables rank `Positive` and `Negative` predictions separately. `rank_within_change=1` is the best model on that sign-specific side, not necessarily rank 1 in a combined TF universe. Method candidate universes differ substantially, so raw Top-K rates are not directly comparable without reporting the method-specific denominator and a harmonized-universe sensitivity analysis.

## TF annotation ambiguities

- The 70 rows contain 41 literal normalized TF/motif labels before biological alias harmonization. Case and spelling variants remain (`SP1`/`Sp1`, `OCT-1`/`Oct-1`, `AP-2-alpha`/`AP-2alpha`).
- Some labels are motif families or elements rather than unique human TF genes: `CACCC`, `CCAAT`, `GATA-A`, `NF-Y`, and `C/EBP`.
- `beta-catenin` is a co-regulator (`CTNNB1`) rather than a sequence-specific TF model in many libraries.
- Three labels carry an asterisk in the workbook (`IRF1*`, `Sp1*`, `USF-1*`). The workbook image/footnote defining the asterisk has not yet been independently interpreted, so the raw and stripped forms are both retained.
- The repository accuracy checker uses deliberately broad aliases (for example, one implicated TF may match several family members). Those aliases are appropriate for reproducing the repository's reported accuracy logic but must not be silently reused for bindFx gene-level evaluation.
- A curated TF mapping table should record exact gene, family-level acceptable matches, bindFx-model availability, and the rationale for every non-exact match.

## Validation performed

`scripts/validate_benchmark.py` verifies:

- 67 unique benchmark identifiers and paired 61-bp sequences.
- Exactly one substitution per pair.
- Referential integrity of observation and prediction rows.
- Required evidence labels.
- Prediction coverage for all 67 variants in each recovered method.
- Consecutive ranks within each method, variant, and predicted-change side.

The validation does not establish primary experimental assay validity, genomic coordinates/build, TF alias equivalence, or cross-method score calibration.
