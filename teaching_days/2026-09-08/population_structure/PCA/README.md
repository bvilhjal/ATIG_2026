# ATiG 2026 PCA exercise

A runnable adaptation of the 2025 PCA exercise, retaining questions Q1–Q14 and
using the same 1000 Genomes dataset and relatedness/outlier thresholds.
The verified full run took 20 minutes on two cores; see [the run report](RUN_REPORT.md).

- `PCA_tutorial.Rmd`: student exercise.
- `PCA_tutorial.html`: rendered student worksheet, without worked answers.
- `PCA_tutorial_answers.Rmd`: worked answers, with calculated results and plots.
- `PCA_tutorial_answers.html`: updated worked answers using the recorded PCA fit.
- `results/`: PC scores, variance explained, related pairs, outlier scores,
  fitted PCA objects, numerical checks, and software/input provenance.

Open either R Markdown file in RStudio and choose **Knit to HTML**. To run
chunks interactively, use **Session > Set Working Directory > To Source File
Location**. Install the packages listed at the start of the exercise if needed.
The first run downloads the input files into `data/`; later runs reuse them.
PLINK 2 is downloaded for the detected architecture if it is absent. The data,
PLINK executable and `results/pca_results.rds` are local caches excluded from
Git. When transferring an existing cache to another platform, remove the cached
PLINK executable so the appropriate build is downloaded.

From this folder, run the complete answer sheet and validation with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 Rscript run_exercise.R
```

The run uses two cores, seed 2026 for the initial PCA and seed 2027 for the
refitted PCA. Set `ATIG_PCA_NCORES` to change the core count. RStudio supplies
Pandoc for HTML rendering; command-line installations need Pandoc available.
PDF output additionally requires LaTeX and was not part of this walkthrough.

After a full run, inspect the fitted results without recomputing PCA using
`readRDS("results/pca_results.rds")` in R.

## Teaching update

The handouts retain all 14 questions in four stages: relatedness, fitting and
interpretation, outlier assessment, and projection. Students consider sample
selection before fitting PCA and assess outliers before refitting. Prediction
prompts, dimension checks and four coding starters support the practical;
the starters use `eval=FALSE` so students can adapt and run them interactively.

The questions address population labels, LD weighting, variance calculations,
cutoff sensitivity, projection and changes in component order. The final plot
marks reference individuals and both kinds of projected individuals separately.
A small sign-flip experiment and a three-sentence analysis note finish the
exercise. The worked answers contain eight numbered figures and two tables.

The original full computation passed 13 checks. The teaching refresh reused
those fits, recalculated the variance and summaries, executed all four student
starters, and passed 28 checks.
The following optional instructor command replays the recorded fit. It requires
the exact saved inputs and PCA cache from the original course folder; those
large files are not included in Git. Use the full-run command above to carry
out a new analysis from a fresh clone.

```sh
Rscript results/teaching_update/refresh.R
```

This took about 30 seconds. It checks the fitting code and input hashes against
the recorded run; use the full-run command above when changing the analysis.
The original executed sources and HTML are preserved in
`results/full_run_sources/`, the earlier teaching version in
`results/teaching_update_v1/`, and current refresh outputs in
`results/teaching_update/`.

## Corrections carried into 2026

Removed machine-specific working directories, a stray markup fragment, and
unused package dependencies. PLINK architecture selection is automatic; ancestry
labels are matched by sample ID; empty related/projected sets are handled.

The worked answer now uses squared singular values divided by the total sum of
squares of the actual standardized PCA matrix. The first 20 PCs therefore need
not sum to 100%. Q8 consistently asks about 90% of total variance and recognizes
that a truncated PCA may not determine how many PCs are needed. PCA comparisons
and removed-sample counts are calculated rather than copied from an old run.

The original threshold of 0.5 for the sample outlier statistic is retained.
This is a teaching threshold, and an unusual ancestry pattern is not itself
proof of a genotyping problem. Related individuals are excluded from fitting
and projected back into the final output.

The method follows the [bigsnpr PCA vignette](https://privefl.github.io/bigsnpr/articles/bedpca.html).
Input hashes and the local 2025 source state are recorded in
`results/source_manifest.json`; see `results/run_info.txt`,
`results/sessionInfo.txt`, `results/plink_version.txt`, and
`results/validation.tsv` for the executed run. The recorded walkthrough reused
the existing data and PLINK cache; a fresh download was not required.

[Population-structure materials](../README.md) · [Course index](../../../../README.md)
