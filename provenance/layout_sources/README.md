# Historical figure and mathematical sources

These files preserve source definitions from the lecture revision process. Use
`Rscript scripts/reproduce.R` for the supported scientific reproduction workflow.

1. `redesign_figures.mjs` retains the original native figure and schematic
   definitions, exposed as a function accepting a presentation, calculated data,
   equation assets and an asset directory. Creating, importing, merging and
   exporting PPTX files has been removed. Some text refers to historical
   centered-only examples and old figure numbers.
2. `genomic_control_figure.mjs` preserves the v5 native Q–Q chart layout. Its
   caller supplies calculated simulation data and rendered equation assets.
3. `standardized_pca.py` is the corrected numerical PCA kernel, including an
   independent rational-arithmetic check. It requires Python and NumPy and
   writes `small.json`; it contains no presentation-building machinery.
4. The JSON files preserve the LaTeX expressions by source version. Old equation
   numbers are historical metadata. The standardized and v6 LMM registries
   supersede their earlier counterparts where keys overlap. They are not a
   single current-deck equation list.

The JavaScript definitions use the `@oai/artifact-tool` presentation object API
available in the authoring environment. That runtime and published image assets
are not distributed here. The functions are archival, not a promise of a
self-contained PPTX rebuild. The R workflow reproduces the scientific plots
without this runtime. Source identities are in `../source_hashes.json`.

The current September 8 v6 lecture has its own sequential figure and equation
numbering. Repository figures use F01–F23; the optional LaTeX heatmap preserves
the earlier source-version Figure 28a label. None of these labels should be used
as an implicit cross-reference to a newly edited deck.
