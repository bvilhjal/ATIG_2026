# 29 October 2026

LD score regression practical. Students read published GWAS results for BMI, triglycerides
and diastolic blood pressure, match them to a HapMap3 LD reference, and estimate SNP
heritability and genetic correlations with [PLDSC](https://github.com/bvilhjal/PLDSC).
Seven questions, about 75 minutes. All code is given; students read the output, interpret
it and change the settings marked `# TRY`. An optional last section adds more traits and
asks students to have an AI model implement LDSC and compare it with PLDSC.

[Open the exercise in Google Colab](https://colab.research.google.com/github/bvilhjal/ATIG_2026/blob/main/teaching_days/2026-10-29/LDSC/LDSC_practical.ipynb)
(the first cell installs PLDSC; choose *File → Save a copy in Drive*).

- [Exercise](LDSC/LDSC_practical.ipynb) ([HTML](https://bvilhjal.github.io/ATIG_2026/teaching_days/2026-10-29/LDSC/LDSC_practical.html))
- [Answers](LDSC/LDSC_practical_answers.ipynb) ([HTML](https://bvilhjal.github.io/ATIG_2026/teaching_days/2026-10-29/LDSC/LDSC_practical_answers.html))
- [`run_exercise.py`](LDSC/run_exercise.py) runs the calculations without Jupyter.

The GWAS files (about 240 MB) are downloaded when the notebook runs and are not stored
here. Sources and checksums are in [`LDSC/data/sources.json`](LDSC/data/sources.json).

[Course index](../../README.md)
