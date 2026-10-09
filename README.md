# liver-transcriptome-study

Code and results for the paper *Expected Expression Trends in Healthy vs Diseased Liver: Expression Distributions in 300 Transcriptomes and Transcriptome Prediction* (RECOMB 2027 submission).

The study simulates 300 liver RNA-seq profiles (150 healthy, 150 MASLD/MASH with fibrosis stages F0--F4) from an explicit negative-binomial model, so that every expected trend is known exactly. It then asks what expression distribution and trends to expect, and how well whole transcriptomes can be predicted.

<!-- TODO before publishing: add the author list, a citation line or DOI, and a LICENSE file (then name the license here). -->

## What `liver_transcriptome_study.R` does

One self-contained script, seeded with `set.seed(20261006)`:

1. **Simulates the cohort.** 6,000 genes: 670 in six disease modules (hepatocyte metabolism/CYP, lipid handling, inflammation, fibrosis/ECM, ductular reaction, plasma-protein synthesis), 6 sex-linked genes and 5,324 background genes. Each diseased liver has a latent severity; module effects follow known shapes of that severity. Counts are drawn from a negative binomial with shared latent factors, module activity, batch and sex effects, and 20% extra dispersion in diseased livers.
2. **Describes the distributions.** Library-wide profile, biological coefficient of variation (BCV) by module and expression level, per-gene variance tests, and the variance explained by fibrosis stage.
3. **Differential expression and trends.** DESeq2 (`~ batch + sex + group`), scored against the ground truth (empirical FDR, power, fold-change bias). limma-trend with orthogonal linear and quadratic stage terms for trend and curvature, and module scores by stage.
4. **Prediction** (stratified 5-fold cross-validation, every selection step inside the training fold). Disease status (10 principal components with LDA; a 20-gene signature with logistic regression; COL1A1 alone, 10 repeats). Fibrosis stage (ridge regression). Whole transcriptomes from K landmark genes chosen by k-means (K = 10 to 400), compared with label-only baselines and with an oracle ceiling.

## Requirements

R 4.5.2 was used for the reported results (Windows). Packages: `MASS`, `matrixStats`, and the Bioconductor packages `limma` and `DESeq2`.

```r
install.packages(c("MASS", "matrixStats", "BiocManager"))
BiocManager::install(c("limma", "DESeq2"))
```

The figures need only Python 3.9 or later and matplotlib (3.11 was used): `pip install matplotlib`.

## Run

```bash
Rscript liver_transcriptome_study.R
```

It takes about 8 minutes (7.8 minutes for the reference run). Outputs go to `./results/` and **overwrite the files shipped here**; the console output is also written to `results/summary.txt`. DESeq2 prints "converting counts to integer mode" and, in one step, a note that a local dispersion fit was substituted; both are expected.

## Figures

```bash
python3 make_figures.py
```

reads `results/*.csv` and writes the paper's three figures to `figures/`: `fig_variance.pdf`, `fig_trends.pdf` and `fig_prediction.pdf`. They are vector PDFs drawn at print size (6.4 in wide), with every label at 10 pt in an embedded STIX (Times-like) font. The script takes a few seconds, uses no random numbers and writes no creation date, so rerunning it with the same matplotlib version gives identical files. It also checks the numbers that the paper's text quotes against `results/` and stops if they differ. The figures shipped here are the ones in the paper.

## Results files and where the paper uses them

| Paper item | File(s) in `results/` |
|---|---|
| Section 3.1 numbers (library profile, dispersion, variance tests) | `summary.txt` (section "Distribution") |
| Figure 1 (a) BCV by gene class | `module_variance.csv` |
| Figure 1 (b) BCV by expression level | `dispersion_by_expression.csv` |
| Section 3.2 DESeq2 numbers | `summary.txt`, `deseq2_diseased_vs_healthy.csv` |
| Table 1 gene modules (the simulation's inputs) | defined in `liver_transcriptome_study.R`, section 2 (`mods` list) |
| Table 2 canonical marker genes (14 of the 16 rows) | `marker_genes.csv` |
| Figure 2 module trends by stage | `module_stage_trend.csv` |
| Table 3 trend and curvature by module | `module_curvature.csv` |
| Table 4 disease classification (mean of 10 repeats, also in `summary.txt`) | `classification_cv.csv` |
| Section 3.3 fibrosis-stage regression | `summary.txt` |
| Figure 3 (a) per-liver r by number of landmarks | `transcriptome_prediction.csv` |
| Figure 3 (b) per-gene R² by module | `transcriptome_prediction_by_module.csv` |
| Simulated data and ground truth | `liver_counts.csv.gz`, `sample_metadata.csv`, `gene_truth.csv` |

The script writes the numbers behind every table and figure; the figures in the paper are drawn from these CSV files by `make_figures.py` (see Figures), not by the R script.

Module keys used in the files: `metabolic` (hepatocyte metabolism, CYP), `lipid`, `inflam`, `ecm` (fibrosis / extracellular matrix), `ductular`, `synthetic` (plasma-protein synthesis), `sex` (six sex-linked genes, `gene_truth.csv` only), `background`.

### Columns

- `sample_metadata.csv`: `sample`, `group` (Healthy/Diseased), `stage` (Healthy, F0-F1, F2, F3, F4), `severity` (latent severity s), `sex`, `batch`, `library_size`.
- `liver_counts.csv.gz`: integer count matrix, genes in rows (first column is the gene name) and the 300 samples in columns.
- `gene_truth.csv`: `gene`, `module`, `base_log2` (baseline log2 abundance), `true_mean_log2FC` (mean true effect over the 150 diseased livers), `nb_dispersion`.
- `deseq2_diseased_vs_healthy.csv`: `gene`, `module`, `baseMean`, `log2FC`, `lfcSE`, `pvalue`, `padj`, `true_log2FC`; sorted by `padj`.
- `marker_genes.csv`: `gene`, `module`, `healthy_log2CPM`, `true_log2FC`, `est_log2FC`, `padj`.
- `module_variance.csv`: per module, expressed genes, median BCV in healthy and diseased livers, median variance ratio (diseased/healthy) of log2 CPM, share of genes significantly more variable in disease, and median share of diseased-group variance explained by stage.
- `dispersion_by_expression.csv`: `bin` (mean log2 CPM), `genes`, median BCV in healthy and diseased livers.
- `module_stage_trend.csv`: `module`, `label`, `stage`, then median, quartiles (`q1`, `q3`) of the module score across livers, and the `expected` score from the simulation truth.
- `module_curvature.csv`: per module, percentage of expressed genes with a significant linear trend (up or down) and curvature (convex or concave), FDR 5%.
- `classification_cv.csv`: `model` (`pca` = 10 PCs + LDA, `sig` = 20-gene signature, `one` = COL1A1), `rep`, `auc`, `acc`, `sens`, `spec`.
- `transcriptome_prediction.csv`: per model (`Global mean`, `Disease label`, `Fibrosis stage`, `LM<K>`, `LM<K>+stage`), median per-gene cross-validated R², percentage of genes with R² > 0.5, and median per-liver centered correlation (overall, healthy, diseased).
- `transcriptome_prediction_by_module.csv`: per module, median per-gene R² for the disease label, the stage, 100 landmarks, 100 landmarks plus stage, and the oracle ceiling.

## Reproducibility

- The simulation and the analysis are fully seeded. For the same R and package versions the results are reproducible; other versions of DESeq2/limma or a different BLAS can change the last digits and, rarely, a borderline significance call.
- The files in `results/` are the output of one run (R 4.5.2, Windows). The script has not yet been re-run on another platform for this repository.
- The simulated count matrix is included, so the analyses can be inspected without re-simulating.
