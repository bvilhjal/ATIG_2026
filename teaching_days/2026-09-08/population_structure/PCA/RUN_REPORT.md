# ATiG 2026 PCA walkthrough — 9 September 2026

The complete worked answer sheet ran on the cached 1000 Genomes dataset in
20 minutes and produced finite scores for every individual. All 13 numerical
and output checks passed. A subsequent teaching update reused these fits,
recalculated the variance and new summaries, rendered both handouts, and passed
28 checks in about 31 seconds. Both HTML documents retain all 14 questions;
the updated worked answers contain eight numbered figures and two tables.

| Quantity | Result |
|---|---:|
| Input individuals | 2,490 |
| Input variants | 1,664,852 |
| Pairs above kinship 2^-4.5 | 31 |
| Distinct individuals in these pairs | 58 |
| Pairs above kinship 2^-3.5 | 0 |
| Individuals in the initial PCA | 2,432 |
| Variants in the initial PCA | 322,376 |
| Samples excluded at S >= 0.5 | 2 |
| Individuals in the final reference PCA | 2,430 |
| Variants in the final PCA | 322,263 |
| Individuals projected back | 60 |
| PCs for every individual | 20 |

**Table 1.** Counts from the executed walkthrough. Sample exclusions apply to
fitting the reference PCs; all individuals receive final scores.

In the initial PCA, PC1 explains **7.67%** and the first
20 PCs together explain **13.86%** of the total
variance of the retained, standardized genotype matrix. They do not reach 90%;
the exact number of PCs required is not determined by this truncated fit.

Both LD-filtering loops converged without reaching their iteration limits.
The original executed HTML in `results/full_run_sources/` preserves the
per-iteration filtering logs. The current HTML identifies its use of saved
fits. A successful teaching walkthrough is not general validation of a
population-QC pipeline.

## Teaching changes and validation

The handouts now follow four stages with a table of contents. Questions about
sample selection and outlier assessment precede the corresponding fits.
Prediction prompts and dimension checks help students inspect their results;
four coding starters cover kinship plotting, sample-label alignment and PC
plotting, comparison of fitted PCs, and invariance to a sign flip. The final
task is a three-sentence analysis note supported by the computed results.

The update clarifies related pairs versus distinct people, singular vectors
versus scores, total variance versus variance within a truncated PCA, and
reference scaling during projection. It adds population mean plots and makes
the outlier threshold a sensitivity exercise: cutoffs 0.4, 0.5 and 0.6 flag
4, 2 and 1 individuals, respectively, on the same initial fit.

The first eight PCs have absolute correlations above 0.997 across fits. The
later-PC comparison shows why matching only the diagonal can mislead: initial
PC19 matches refitted PC18 almost perfectly. Plot signs are aligned for visual
comparison; the exported PC scores retain their original signs. Shapes in
Figure 7 of the worked answers distinguish the 2,430 reference individuals,
58 related individuals projected back, and two projected outliers.

The refresh checked that the eight fitting-related blocks of R statements and
the input hashes match the original run. Sample selection was split from the
initial fitting chunk to place Q4–Q5 between those steps; concatenating the two
chunks gives the original code. It reused the saved PCA models and outlier
scores, recalculated the total-variance denominator, and executed the new
summaries and plotting code. The original 13 checks and eight additional
consistency checks passed, including unchanged raw scores and variance.
Seven further checks cover the new teaching code: score definitions,
projection roles, sign invariance, aligned student labels, construction of
the student plot, and agreement of both comparison starters with the answers.
All four student starters executed successfully in the verified fit context.
See `results/teaching_update/validation.tsv` and
`results/teaching_update/refresh_info.json`.
The student HTML is a code-only worksheet; its shared fitting code matches
the worked answers. All new or changed plots were inspected directly.

The run used R 4.6.0 on macOS arm64, bigsnpr 1.12.21, and two cores. Exact
versions, seeds, elapsed time, and the execution log are recorded under
`results/`. The existing data and PLINK cache were reused; fresh downloads
were not exercised. The browser preview of local HTML was unavailable; the
HTML structure and generated plot files were inspected directly.

See [the worked answers](PCA_tutorial_answers.html),
[student worksheet](PCA_tutorial.html), and [run instructions](README.md).
