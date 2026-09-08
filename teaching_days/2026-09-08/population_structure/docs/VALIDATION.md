# Validation of version 0.2.0

Version 0.2.0 groups sources, documentation, provenance and previews by the
teaching date, 8 September 2026. The root reproduction command is
unchanged; its outputs now include the teaching date. The numerical models,
seeds and stable figure IDs are preserved.

The complete R workflow ran on 8 September 2026 with R 4.6.0, bigsnpr 1.12.21,
bigstatsr 1.6.2 and jsonlite 2.0.0 on an arm64 macOS host. It generated all 23
figure families in three formats and passed 30 numerical/output checks.
The committed machine-readable records are [the run manifest](../provenance/validated_run.json) and
[the check results](../provenance/validation.json). Installed direct and transitive package versions
are listed in [the dependency record](../provenance/dependencies.tsv); hashes of the validated R sources
are in [the source hashes](../provenance/validated_code_hashes.json).

| Calculation | Observed result |
|---|---:|
| Unadjusted null-SNP estimate | 0.374658 |
| PC1-adjusted null-SNP estimate | −0.083134 |
| Squared correlation of PC1 with group | 0.991238 |
| SNPs retained after simulated LD filtering | 2,468 of 3,000 |
| GC inflation factor before correction | 1.390651 |
| GC inflation factor after correction | 1.000000, by construction |
| Nominal 0.05 null rejection rate before/after GC | 0.09739 / 0.05107 |
| Four-SNP likelihood optimum | 0.217655 |
| Exact population-LD score for SNP 4 | 3.35 |
| Other friends under the illustrative r² cutoff | 4 |

**Table 5.** Results of the checked teaching calculations. The null-SNP
coefficient differs from the saved lecture snapshot by less than $4\times10^{-12}$;
the LD-filtered SNP count agrees exactly. This establishes reproduction of the
displayed realization, not general method performance.

The checks cover centered/scaled matrix identities, symmetry and positive
semidefiniteness, eigendecomposition and score scaling, regression against an
independent `lm` fit, Frisch–Waugh–Lovell residualization, correct confidence
interval degrees of freedom, disjoint projection samples, ancestry components
summing to one, GC rank preservation, identical GC/LDSC normal draws and exact
LD-score arithmetic. The corrected small PCA kernel was also run independently
in Python, with rational arithmetic verifying the GRM entries.

The version 0.2.0 simulation-and-figure run took about 34 seconds on this host:
14.6 seconds for the main population examples, 4.9 seconds for GC, 2.0 seconds
for LD, 1.2 seconds for the small calculations and 11.3 seconds for rendering.
These measured phase times exclude validation and are not a portable
performance guarantee.

The reorganized workflow reproduced all 23 PNG previews pixel for pixel.
All 23 rendered figure families were reviewed on a contact sheet. The LD and GC
plots and the LD-filtering example were also inspected at large size. The
optional LaTeX LD heatmap rebuilt to a 3200 × 2000 PNG and vector SVG using
Tectonic, Poppler, lxml and sharp 0.35.4.

The `bigsnpr` random stream and eigenvector signs can differ with versions and
platforms. Exact PNG bytes are not a numerical reproducibility criterion. Linux
and Windows were not available for a live end-to-end test; their standard R
Cairo path is implemented, while the validated macOS path uses Quartz and
Poppler. The code makes no claim of empirical GWAS reproduction, repeated-run
calibration or biobank-scale performance.
