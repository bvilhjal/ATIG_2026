# Optional LaTeX LD-score figure

The main R workflow already creates this scientific example in PNG, PDF and SVG.
This optional workflow preserves the spacious original Figure 28a layout and
uses actual LaTeX for the two equations.

Requirements: R with `jsonlite`, Python 3 with `lxml`, Node.js with `sharp`,
[Tectonic](https://tectonic-typesetting.github.io/) and Poppler's `pdftocairo`.
All executables are discovered on `PATH`. For a custom installation, the
`RSCRIPT`, `TECTONIC`, `PDFTOCAIRO`, `NODE` and `SHARP_MODULE` environment variables
accept explicit paths. The latter points to the installed `sharp` module.

From the repository root:

```sh
python3 -m pip install -r population_structure/optional/ld_score_latex/requirements.txt
npm --prefix population_structure/optional/ld_score_latex install
python3 population_structure/optional/ld_score_latex/build_figure.py
```

Outputs go to `outputs/ld_score_latex/`. Use `--output DIRECTORY` to select
another folder. `--offline` asks Tectonic to use only cached LaTeX packages.
`--equation-number 22` adds optional lecture-specific equation tags. Numbering
is omitted by default because the current deck can be reordered independently.

The SVG contains accessible title/description text and vector equations. Its
PNG is 3200 × 2000 pixels. Minor font metrics can vary across operating systems.
The original scientific model, focal score 3.35 and thresholded count 4 are
checked before drawing. Neither the plot nor its color patterns are empirical
LD measurements.
