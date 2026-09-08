# Figure and provenance inventory

Repository figure IDs remain stable when lecture slides are reordered. The
current deck's figure numbers can differ. The R plotting script produces the
PNG, PDF and SVG files named in Table 2; committed PNGs are convenient previews.
Their captions distinguish original computations from published summaries.

| ID | Filename stem | Data or calculation |
|---|---|---|
| F01 | `standardized_pca` | Corrected seven-SNP GRM, eigenvectors, PC geometry and projections. |
| F02 | `ancestry_likelihood` | Four-SNP binomial likelihood and mixture-proportion optimum. |
| F03 | `reference_panels` | Seed 905; 30 synthetic genomes and two reference panels. |
| F04 | `local_global_ancestry` | Deterministic chromosome-copy schematic and length-weighted ancestry. |
| F05 | `stratification_diamond` | Original diamond schematic for the running null-SNP example. |
| F06 | `pooled_association` | Running simulation, before PC adjustment. |
| F07 | `background_pca` | The same 600 individuals and 3,000 background SNPs. |
| F08 | `residual_association` | Both phenotype and dosage adjusted for intercept and PC1. |
| F09 | `known_effects` | Null and causal SNP estimates with 95% t intervals. |
| F10 | `pc_sensitivity` | Null-SNP estimate for 0, 1, 2, 5, 10 and 20 PCs. |
| F11 | `ld_pca` | Seed 901; PCA before/after `snp_autoSVD` filtering. |
| F12 | `ld_loadings` | All 3,000 loadings before filtering in F11. |
| F13 | `batch_pca` | Seed 902; the same PCs colored by group and by batch. |
| F14 | `oadp_projection` | Seed 903; 480 reference and 120 held-out samples. |
| F15 | `genomic_control` | Seed 904; uniform null inflation before/after genomic control. |
| F16 | `matched_qq` | Same normal draws as F15; two inflation scenarios. |
| F17 | `ldsc_scenarios` | Same statistics as F16; binned LD-score conditional means. |
| F18 | `ld_score_heatmap` | Exact 12-SNP population LD; weighted score and example friend count. |
| F19 | `pedigree_relationships` | Original pedigree and expected additive-relationship matrix. |
| F20 | `parental_transmission` | Original transmitted-genotype and family-environment schematic. |
| F21 | `local_ancestry_model` | Original ancestry-specific association-model schematic. |
| F22 | `rare_variant_summary` | Five published values from Hanson et al. (2026). |
| F23 | `moba_summary` | Four published estimates, sample sizes and intervals from Corfield et al. (2026). |

**Table 2.** Portable figure families. The deterministic four-SNP raw-correlation
and two-individual PCA appendix calculations are saved as JSON/RDS by
`small_calculations.R`; they are worked tables rather than extra plots. The five
progressive PCA teaching states are represented by F01 and preserved in the
historical layout source. Numerical outputs match the final standardized example.

## Saved generated-image prompts

| Source-version label | File | Status |
|---|---|---|
| Figure 13, admixture v1 | `Figure_13_admixture_20260908_v1.prompt.txt` | Original generated teaching schematic, not empirical data. |
| Figure 21, stratification v1 | `Figure_21_population_stratification_20260908_v1.prompt.txt` | First redesigned raster illustration. |
| Figure 21, stratification v2 | `Figure_21_population_stratification_20260908_v2.prompt.txt` | Revised diamond layout requested by the lecturer. |

**Table 3.** The saved prompts are in `population_structure/prompts/`. They record
the requested geometry and semantics, not a random seed or a reproducible image
model version. Generating them again may change the pixels; those outputs require
visual checking. F04/F05/F21 supply deterministic vector alternatives for the
main scientific relationships. The original generated raster images are not
needed by the numerical pipeline.

## Published panels used in the lecture

These are source figures, not figures simulated by this repository. No pixels or
participant-level coordinates have been copied here. The historical crop
inventory preserves panel framing and source attribution, and the corresponding
layout definitions are in `provenance/layout_sources/`.

| Lecture subject | Source recorded in the lecture | Treatment here |
|---|---|---|
| Height-PGS accuracy across genetic distance | Ding et al. (2023), [doi:10.1038/s41586-023-06079-4](https://doi.org/10.1038/s41586-023-06079-4), Figure 3b. | Published panel excluded; crop provenance retained. |
| Fine-scale ancestry components | Hu et al. (2025), [doi:10.1038/s41588-024-02035-8](https://doi.org/10.1038/s41588-024-02035-8), Figure 3c and its legend. | Published panel excluded; crop provenance retained. |
| All of Us PCA-density groups | Sharma et al. (2025), [doi:10.1038/s41467-025-59351-8](https://doi.org/10.1038/s41467-025-59351-8), Figure 1c. | Published panel excluded; the lecture distinguishes seven PCA-density groups from thirteen separately defined UMAP groups. |
| UK Biobank PCA | Privé et al. (2022), [doi:10.1016/j.ajhg.2021.11.008](https://doi.org/10.1016/j.ajhg.2021.11.008). | Supplied lecture image excluded; no point coordinates reconstructed. |
| 1000 Genomes PCA | Original supplied 2025 lecture panels. | Source coordinates were not available; image excluded. |
| Population trees | Jakobsson et al. (2008), [doi:10.1038/nature06742](https://doi.org/10.1038/nature06742); Li et al. (2008), [doi:10.1126/science.1153717](https://doi.org/10.1126/science.1153717). | Published/source-lecture trees excluded; topology and branch lengths not reconstructed. |
| HGDP FRAPPE ancestry proportions | Tang et al. (2005), [doi:10.1002/gepi.20064](https://doi.org/10.1002/gepi.20064); Li et al. (2008), above. | Empirical ancestry panel excluded. F03 is a separate synthetic likelihood example. |
| Microsatellite STRUCTURE proportions | Rosenberg et al. (2002), as cited by the supplied lecture. | Source plot excluded; not regenerated from synthetic data. |
| Historical height and type 1 diabetes maps | Source references retained in the supplied lecture. | Geographical images excluded. |
| Historical diabetes/rheumatoid-arthritis Q–Q panels | Source references retained in the supplied lecture. | Empirical panels excluded; F15/F16 are separate original simulations. |
| Rare-variant inflation and MoBa attenuation | Hanson et al. (2026) and Corfield et al. (2026), linked in [Methods](METHODS.md). | Numerical summaries checked against the primary papers and redrawn in F22/F23. |

**Table 4.** Inventory of published/source-derived panels. References in this
table document the lecture provenance; they are not claims that the code
reproduces the underlying empirical analyses.
