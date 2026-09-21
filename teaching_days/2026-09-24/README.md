# 24 September 2026

[LT-FH exercise](https://bvilhjal.github.io/ATIG_2026/teaching_days/2026-09-24/LTFH/LTFH_exercise.html): a practical on family-history genetic liability scores
with the Python package [ltpred](https://github.com/bvilhjal/ltpred), following
the liability-model lecture of 22 September.

- [Student exercise](LTFH/LTFH_exercise.ipynb) and [view exercise in browser](https://bvilhjal.github.io/ATIG_2026/teaching_days/2026-09-24/LTFH/LTFH_exercise.html).
- [Worked answers](LTFH/LTFH_exercise_answers.ipynb) and [view answers in browser](https://bvilhjal.github.io/ATIG_2026/teaching_days/2026-09-24/LTFH/LTFH_exercise_answers.html).
- [`run_exercise.py`](LTFH/run_exercise.py) prints every answer without Jupyter.

Its 13 questions in four parts cover what a wrong heritability does to the
score, how far a posterior mean liability is from a disease risk, estimating
heritability from the register itself with a tetrachoric correlation, and
moving an observed-scale GWAS h² (Lee et al. 2011; 2012), with effective
sample sizes, onto the liability scale. Needs Python 3.9+ with NumPy, SciPy
and ltpred installed from source (not on PyPI; see the notebook's Setup
section) — no data files, no R, and Windows, macOS and Linux behave the same.
Part 0 rebuilds the
[ltpred tutorial](https://bvilhjal.github.io/ltpred/tutorial/) register from
fixed seeds, so no prior session is needed. Parts A and B take one hour;
Parts C and D are homework. The handouts were executed on 21 September with
ltpred 0.7.0 on a simulated register whose genetic liability is known.

[Course index](../../README.md)
