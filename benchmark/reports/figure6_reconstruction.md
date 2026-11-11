# Reconstruction of BayesPI-BAR 2015 Figure 6

## Can `published_predictions.tsv` reproduce Figure 6?

Yes. It contains the complete variant-level ranking tables used for BayesPI-BAR, sTRAP, and is-rSNP. The exact aggregation logic was recovered from `code/r/accuracy.R` inside the original `bayesPI-BAR_GLibc2.5.tar.gz` demo archive on the `junbaiw/BayesPI-BAR` repository's historical branch.

The reconstructed inputs are:

- `data/figure6_best_ranks.tsv`: one selected mutation-level rank per method and panel.
- `data/figure6_cumulative_accuracy.tsv`: cumulative counts and fractions at cutoffs 1–50.
- `data/figure6_pairwise_ranks.tsv`: direction-agnostic ranks capped at 51 for panels C–D.
- `scripts/reproduce_figure6.py`: standard-library implementation of the original matching logic.

## Why are there multiple rows for the same TF?

Each archived row is a TF-binding model/PWM, not necessarily a unique TF gene. One TF can have several models from TRANSFAC, JASPAR, protein-binding microarrays, or different family members accepted by the paper's alias rules. Positive and negative affinity changes are ranked independently, so rank 1 can occur once on each side.

The original analysis reduces these rows as follows:

1. Parse the experimentally implicated TF or motif from the mutation label.
2. Expand it with the paper's curated alias rules. These rules can be broad; for example, `GATA-1` matches any model name containing `GATA`, and `STAT1` matches `STAT`.
3. Select every PWM/model row whose name matches an accepted alias.
4. Take the minimum rank among matching rows. This is the single rank assigned to that mutation and method.
5. If no accepted model appears in the retained rankings, treat it as missing; panels C–D cap it at 51.

Thus, multiple PWMs are treated as alternative representations of the experimentally implicated TF. They are not averaged, and the number of PWM rows is not the benchmark denominator.

## Panel definitions

### A: direction not considered

- Denominator: all 67 mutations.
- A matching TF model can be on either the positive or negative side.
- At cutoff K, accuracy is the fraction whose minimum matching rank is ≤K.
- Reconstructed Top-10: BayesPI-BAR 47/67, sTRAP 41/67, is-rSNP 34/67.

### B: considering direction

- The 16 mutations annotated `UNKNOWN` are excluded; denominator: 51.
- `GAIN` requires a `Positive` prediction and `LOSS` requires a `Negative` prediction.
- Reconstructed Top-10: BayesPI-BAR 29/51, sTRAP 27/51, is-rSNP 22/51.

### C and D: pairwise ranks

- Use the direction-agnostic ranks from panel A.
- Missing or rank >50 is plotted as 51.
- BayesPI-BAR has a lower rank than sTRAP for 28/67 = 41.8% (reported as 42%); sTRAP is lower for 14/67 = 20.9% (21%).
- BayesPI-BAR has a lower rank than is-rSNP for 30/67 = 44.8% (45%); is-rSNP is lower for 11/67 = 16.4% (16%).
- The paper reports one-sided paired Wilcoxon signed-rank p-values of 0.015 and 0.00078, respectively.

## Important implications for bindFx

Figure 6 can be reconstructed as a historical-method replication, but its broad alias rules should not define the primary bindFx comparison. They allow family-level matches and select the best of multiple PWMs, which can improve apparent rank performance for methods with larger or more redundant model libraries.

For bindFx, report at least two analyses:

1. **Strict common-TF analysis:** the 18 exact/established-alias bindFx models already curated, with one bindFx model per TF concept wherever possible.
2. **Historical-rule sensitivity analysis:** reproduce the paper's broad alias/minimum-rank convention, clearly noting that model-universe sizes and redundancy differ.

Direction concordance should remain separate from TF retrieval. The original panel B combines both requirements: a matching TF must also occur on the correct sign side.
