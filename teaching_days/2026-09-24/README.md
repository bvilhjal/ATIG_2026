# 24 September 2026

[LT-FH exercise](https://bvilhjal.github.io/ATIG_2026/teaching_days/2026-09-24/LTFH/LTFH_exercise.html): a practical on family-history genetic liability scores
with the Python package [ltpred](https://github.com/bvilhjal/ltpred), following
the liability-model lecture of 22 September.

- [Student exercise](LTFH/LTFH_exercise.ipynb) and [view exercise in browser](https://bvilhjal.github.io/ATIG_2026/teaching_days/2026-09-24/LTFH/LTFH_exercise.html).
- [Worked answers](LTFH/LTFH_exercise_answers.ipynb) and [view answers in browser](https://bvilhjal.github.io/ATIG_2026/teaching_days/2026-09-24/LTFH/LTFH_exercise_answers.html).
- [`run_exercise.py`](LTFH/run_exercise.py) prints every answer without Jupyter.

Its 15 questions in five parts use one simulated register of 5,075 people whose
genetic liability is known (the model of the
[ltpred tutorial](https://bvilhjal.github.io/ltpred/tutorial/), five times larger).
Part 0 simulates it and shows the scoring call once, argument by argument; two
helpers, `score(h2)` and `score_at_40(h2)`, then rescore it for every question, and
each student cell asks for one or two expressions marked `...`. In class (Parts A
and B, one hour): what a wrong heritability does to the score, and how far a
posterior mean liability is from a disease risk, measured as AUC and as calibrated
absolute risk by rank fifth. Homework (Parts C to E): the score against a 0/1
family-history indicator; heritability from parent–offspring pairs by tetrachoric
correlation, and what an iPSYCH-style case-cohort does to it; and the
observed-to-liability-scale transformation of Lee et al. (2011; 2012) with
effective sample sizes. Install Python, NumPy, SciPy and ltpred from the
[setup page](https://bvilhjal.github.io/ATIG_2026/setup.html) (macOS, Windows, or
Anaconda). No data files and no R; the notebook runs in under a minute with less
than 1 GB of memory. The handouts were executed on 23 September with ltpred 0.7.1.

[Course index](../../README.md)
