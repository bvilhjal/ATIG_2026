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
working in pairs. A five-sentence primer and a glossary cover the liability model;
a four-line Python cheat sheet and a collapsible hint under every coding task cover
the code, where each student cell asks for one or two expressions marked `...`.
Questions marked **Discuss** are for the pair, and Q11, Q13 and Q14 are optional
extensions. The 15 questions in five parts use one simulated register of 5,075
people whose genetic liability is known (the model of the
[ltpred tutorial](https://bvilhjal.github.io/ltpred/tutorial/), five times larger).
In class (Part 0 and Parts A and B, about an hour): what a wrong heritability does
to the score, and how far the score is from a disease risk, measured as AUC and as
calibrated absolute risk by rank fifth. Homework (Parts C to E): the score against
a yes/no family-history indicator; heritability from parent–offspring pairs by
tetrachoric correlation, and what an iPSYCH-style case-cohort does to it; and the
observed-to-liability-scale transformation of Lee et al. (2011) with
effective sample sizes. In Colab the first cell installs ltpred 0.7.1; to run it on
your own computer instead, follow the
[setup page](https://bvilhjal.github.io/ATIG_2026/setup.html). No data files and no R;
the notebook runs in under a minute. The handouts were executed on 23 September with ltpred 0.7.1.

[Course index](../../README.md)
