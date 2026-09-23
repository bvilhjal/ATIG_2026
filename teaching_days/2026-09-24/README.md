# 24 September 2026

[LT-FH exercise](https://bvilhjal.github.io/ATIG_2026/teaching_days/2026-09-24/LTFH/LTFH_exercise.html): a practical on family-history genetic liability scores
with the Python package [ltpred](https://github.com/bvilhjal/ltpred), following
the liability-model lecture of 22 September.

**[Open the exercise in Google Colab](https://colab.research.google.com/github/bvilhjal/ATIG_2026/blob/main/teaching_days/2026-09-24/LTFH/LTFH_exercise.ipynb)**: nothing to install.
Choose *File → Save a copy in Drive*, then run the cells in order.

- [Student exercise](LTFH/LTFH_exercise.ipynb) and [view exercise in browser](https://bvilhjal.github.io/ATIG_2026/teaching_days/2026-09-24/LTFH/LTFH_exercise.html).
- [Worked answers](LTFH/LTFH_exercise_answers.ipynb) and [view answers in browser](https://bvilhjal.github.io/ATIG_2026/teaching_days/2026-09-24/LTFH/LTFH_exercise_answers.html).
- [`run_exercise.py`](LTFH/run_exercise.py) prints every answer without Jupyter.

Written for master's students with mixed biology and bioinformatics backgrounds,
working in pairs. A five-sentence primer and a glossary cover the liability model, and a
short Python cheat sheet and a hint under every coding task cover the code. The 11
questions use one simulated register of 5,075 people whose genetic liability is known;
students make six plots, each from a ready template with one `...` to fill: the score
against the truth, calibration under three assumed heritabilities, score histograms for
later cases and non-cases, predicted against observed risk by fifth of the score, R² of
own status, a yes/no family history and the score, and observed-scale heritability
against the case fraction. In class (Part 0 and Parts A and B, about an hour): what a
wrong heritability does, and how well the score predicts later diagnosis. Homework
(Parts C to E): the family-history baseline, heritability from parent–offspring pairs by
tetrachoric correlation, and observed- versus liability-scale heritability (Lee et al.
2011). In Colab the first cell installs ltpred 0.7.1; to run it on
your own computer instead, follow the
[setup page](https://bvilhjal.github.io/ATIG_2026/setup.html). No data files and no R;
the notebook runs in under a minute. The handouts were executed on 23 September with ltpred 0.7.1.

[Course index](../../README.md)
