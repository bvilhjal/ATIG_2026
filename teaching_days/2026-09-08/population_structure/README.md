# Population structure — 8 September 2026

Start with the [methods](docs/METHODS.md) for the notation, models and seeds.
The [figure inventory](docs/FIGURES.md) identifies all 23 figure families and
separates original simulations from published summaries and source panels.
The [validation record](docs/VALIDATION.md) reports the reproduction checks.

- `R/`: simulations, deterministic worked examples, plotting and numerical checks.
- `PCA/`: the 1000 Genomes [PCA practical](PCA/README.md), with student and worked-answer handouts and a separate full-analysis runner.
- `data/`: two source-cited published summary tables; no individual data.
- `figures/`: PNG previews from the validated workflow and the original LaTeX LD heatmap.
- `prompts/`: saved prompts for the generated admixture and stratification illustrations.
- `optional/`: the standalone workflow for the LD heatmap with actual LaTeX equations.
- `provenance/`: lecture baselines, source hashes, validation records and historical figure/math definitions.

From the repository root, run `Rscript scripts/reproduce.R`. It writes the
calculated data and PNG/PDF/SVG figures under `outputs/2026-09-08/`. Simulated
examples are synthetic; the original model seeds and repository figure IDs are
preserved when the lecture is rearranged.

[Teaching day](../README.md) · [Course index](../../../README.md)
