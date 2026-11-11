# bindFx overlap with experimental TF/motif observations

Input: user-supplied `403TFmodel_summary.csv` containing 403 unique `TF_Name` models.

## Recommended primary benchmark

Use only exact gene-symbol matches and established aliases:

- **18 unique experimental TF concepts / 18 bindFx models**
- **39 experimental TF-effect rows**
- **38 of 67 variants** have at least one primary-eligible TF match
- **34 TF-effect rows across 33 variants** have both a primary match and a known gain/loss direction

Primary bindFx models:

`AR, CEBPD, EGR1, FOS, FOXA1, GATA1, HLF, IRF1, NFE2L2, POU2F1, SP1, SPI1, TBP, TCF7L2, TFAP2A, USF1, USF2, YY1`

Established alias examples include NRF2 → NFE2L2, OCT-1 → POU2F1, AP-2alpha → TFAP2A, c-FOS → FOS, and the beta-catenin-pathway name TCF4 → TCF7L2. Each mapping and rationale is recorded in `data/bindfx_experimental_tf_overlap.tsv`.

## Secondary family/motif-assisted analysis

If explicitly reported as a separate sensitivity analysis, family or motif labels can be linked to candidate bindFx models:

- **25 comparable experimental concepts**
- **31 distinct bindFx models**
- **56 TF-effect rows across 55 variants**
- **42 known-direction rows across 41 variants**

The seven additional concepts are C/EBP family, CACCC-binding KLF family, CCAAT/NF-Y motif, ETS family, HNF4 family, NF-Y complex, and generic USF family. These are not one-to-one TF matches and must not be included in the primary denominator.

## Unavailable concepts

Twelve experimental concepts lack an acceptable bindFx model under this mapping:

`CTNNB1 (beta-catenin), ETV7 (TEL2), HNF1 family, MSX1, MYB, MYF5, NF-kappaB family, REL (c-Rel), SOX10, SRY, STAT1, TP53 (p53)`

Available paralogs were not substituted silently: for example, MSX2/MSX3 do not replace MSX1, MYBL2 does not replace MYB, and other SOX proteins do not replace SOX10 or SRY.

## Benchmark denominators

- TF-ranking benchmark, strict: 39 experimental TF-effect rows in 38 variants, subject to how multiple TF observations are scored.
- Direction benchmark, strict: 34 known-direction TF-effect rows in 33 variants.
- Do not use all 67 variants as the direction denominator because 16 experimental observations have unknown direction and several implicated TFs lack bindFx models.
- Report both observation-level and unique-variant denominators; multi-TF variants otherwise receive implicit extra weight.
