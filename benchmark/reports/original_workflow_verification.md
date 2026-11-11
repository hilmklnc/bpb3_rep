# Original BayesPI-BAR workflow verification

Verification date: 2026-11-05. This report records a completed smoke test; it
is not a full 67-variant benchmark.

## Successfully executed

The original paired-FASTA command path completed end to end for
`AGTRL1_SP1|GAIN` and three unmodified SP1 PWM files. It used all original
chemical potentials (`none -10 -13 -15 -18 -20`), seed 1, one core, and 100
shuffles rather than the production 10,000. Exit status was 0; elapsed time was
0.508 seconds (0.178 seconds user, 0.147 seconds system).

```bash
python -m bpb3.bpb3 bayespi_bar \
  --reference_sequences inputs/ref.fasta \
  --alternate_sequences inputs/alt.fasta \
  --pwm_folder pwm \
  --chemical_potentials none -10 -13 -15 -18 -20 \
  --iterations 100 \
  --result_folder results/three-pwm \
  --max_rank 50 --seed 1 --use_cores 1 --max_nodes 1
```

Actual integrated mean delta-dbA values were M00008 = 0.0, M00196 =
0.01704463738, and M00665 = 0.0. PCA scores were -0.028752158373098705,
0.05750431674619741, and -0.028752158373098705, respectively. The positive
rank-1 result was `SP1_known2_from_SP1_2_from_Sp1_transfac_M00196.mlp`, with
PCA score 0.05750431674619741 and mean delta-dbA 0.01704463738.

Output checksums:

```text
2e792a45a704a345cb1cd6b7ec2b29d658d5c7169d1419f4ce3da057b0c3ee09  ddba_integrated_using_mean_ddba.tsv
2c79e567fee0d2f358a283b33f6c9b94aff4e2a77ecaf7e9e9bc3b5b471f604d  ddba_integrated_using_pca.tsv
25423afb2cde90d88c55eded69bc1155f1b06e0f5581e36708ec1f6c30b6fd77  AGTRL1_SP1_GAIN_ranking.tsv
```

Runtime: Python 3.9.25, NumPy 1.24.3, pandas 1.4.4, SciPy 1.10.1, numba
0.60.0, matplotlib 3.5.3, seaborn 0.11.2, multiprocess 0.70.14, and setuptools
65.6.3. The binary executed here was the bundled x86-64 Linux ELF. The
repository also contains an x86-64 Mach-O `bpb3/bin/Mac/bayesPI_affinity`,
which bpb3 selects on Darwin; that binary was not tested here and may require
Rosetta on Apple Silicon.

## Observed warnings and errors

- A one-PWM run exited successfully but produced NumPy PCA
  degrees-of-freedom/division warnings; one-feature PCA is degenerate.
- `/usr/bin/time` was absent, so shell built-in timing was used.
- Matplotlib warned about unwritable cloud config/font caches until
  `MPLCONFIGDIR=/tmp/bpb3-mpl` was set.
- Bedtools was not on PATH; the verified cloud executable was bedtools 2.31.1
  in a local prefix.
- A bedtools 2.30.0 package download failed when the environment proxy was
  unavailable, so binary-level equality was not tested.

Bedtools is not called by paired-FASTA mutation scoring. Version differences
can matter only in upstream sequence/background construction. Bedtools 2.31.1
passed representative `getfasta`, `intersect`, `sort`, `merge`, and `slop`
coordinate microtests.

## Explicitly untested

- Full 67 variants x 1,772 PWMs x 10,000 shuffles
- Production runtime and all-PWM numerical outputs
- Bedtools 2.30.0 vs 2.31.1 byte-for-byte comparison
- MuSSD, background modeling, significance testing, gene-expression filtering,
  plotting, and cleanup
- Genome FASTA/annotation download and validation

The original demo command and inputs remain in
`final_demo/demo1_67snps/bpb3_on_67snps_may22.py` and
`final_demo_data/demo1_67snps/in/`.
