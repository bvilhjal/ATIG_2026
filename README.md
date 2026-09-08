# ATIG 2026 — reproducible teaching code

Simulation and figure-generation code, organized by teaching date.
Version **0.2.0** contains the population-structure lecture and produces
**23 numbered figure families** in PNG, PDF and SVG, with calculated data and
numerical checks.

## Teaching days

| Teaching date | Topic | Material |
|---|---|---|
| [8 September 2026](teaching_days/2026-09-08/README.md) | [Population structure](teaching_days/2026-09-08/population_structure/README.md) | PCA, ancestry, stratification, genomic control, LD scores, simulations and figure sources. |

**Table 1.** Teaching days currently represented in this repository.

## Reproduce the figures

Install R and the packages `bigsnpr`, `bigstatsr` and `jsonlite`. On macOS, install
Poppler (`pdftocairo`) for SVG conversion; for example, `brew install poppler`.
The other platforms use R's Cairo graphics devices. A working Cairo-capable R
installation is required there.

From this repository's root:

```sh
Rscript scripts/install_dependencies.R   # once, if the packages are missing
Rscript scripts/reproduce.R
```

The second command runs all calculations, writes figures to
`outputs/2026-09-08/figures/`, checks numerical invariants, and records package
versions and elapsed times in `outputs/2026-09-08/run_manifest.json`. It needs
no network access. An optional output root can be supplied:

```sh
Rscript scripts/reproduce.R outputs_classroom
```

The teaching date is retained inside that output root. The shared entrypoint
and dependency installer are in `scripts/`; source code, data, figures, methods
and provenance belong to their teaching-day folder. For the spacious LD heatmap
with actual LaTeX equations, use the
[optional LaTeX workflow](teaching_days/2026-09-08/population_structure/optional/ld_score_latex/README.md).
The normal R workflow already reproduces its data and scientific content.

![Figure 18. LD-score heatmap.](teaching_days/2026-09-08/population_structure/figures/F18_ld_score_heatmap.png)

**Figure 18.** Exact illustrative population LD. The highlighted row has score
3.35, including self-correlation 1. Four other SNPs exceed an illustrative
`r² ≥ 0.20` cutoff. That cutoff is not used to calculate the LD score. Repository
figure IDs are independent of slide and equation numbers in an editable lecture.

## Reuse and limitations

The genotype and phenotype examples are synthetic. Two small published numerical
summaries are explicitly labelled and cited. Lecture files, participant data,
downloaded paper panels and private course documents are not stored here.

The code is teaching material, not a production GWAS pipeline. A single selected
simulation does not establish general calibration, power or robustness. The
LDSC examples simulate a working model for test statistics rather than real
genotype GWAS data. See the [methods](teaching_days/2026-09-08/population_structure/docs/METHODS.md)
and [validation](teaching_days/2026-09-08/population_structure/docs/VALIDATION.md).

This repository is publicly available. No open-source license is currently
specified. Third-party software and published material retain their own licenses.
