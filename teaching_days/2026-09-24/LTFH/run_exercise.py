"""ATiG 2026, 24 September: complete LT-FH exercise calculations."""
import matplotlib
matplotlib.use("Agg")

import numpy as np
from scipy.stats import norm, rankdata
import matplotlib.pyplot as plt

import ltpred
from ltpred import (simulate_pedigree, simulate_register_liabilities, estimate_liabilities,
                    kaplan_meier_cip)
from dataclasses import replace
from collections import Counter
print("ltpred", ltpred.__version__)


def corr(x, y):
    """Pearson correlation of two arrays."""
    return np.corrcoef(x, y)[0, 1]


def auc(score, case):
    """AUC: the chance that a random case scores above a random non-case.

    Computed from ranks (the Mann-Whitney statistic); `case` is a boolean array."""
    r = rankdata(score)
    n1, n0 = case.sum(), (~case).sum()
    return (r[case].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

H2 = 0.5                                        # true liability-scale heritability
AGES = np.arange(0, 121.0)                      # ages 0, 1, ..., 120
CIP_M = 0.12 / (1 + np.exp((58 - AGES) / 8))    # men: lifetime 12%, half of it by age 58
CIP_F = 0.08 / (1 + np.exp((62 - AGES) / 8))    # women: lifetime 8%, half of it by age 62
K = 0.10                                        # simplified lifetime prevalence for Part E only

ids, father, mother = simulate_pedigree(np.random.default_rng(1), n_founder_pairs=500, gens=2)
fathers, mothers = set(father), set(mother)
coin = np.random.default_rng(2).random(len(ids)) < 0.5  # a random sex for people who never became parents
male = np.array([p in fathers or (p not in mothers and c) for p, c in zip(ids, coin)])
sex = np.where(male, "M", "F")


def simulate(cip):
    """The register with one incidence curve for everyone (same seed = same liabilities)."""
    return simulate_register_liabilities(np.random.default_rng(1), ids, father, mother,
                                         h2=H2, cip_ages=AGES, cip_values=cip, eval_age=70)


as_men, as_women = simulate(CIP_M), simulate(CIP_F)
# each person's records follow their own sex's curve: sex-specific thresholds, as in LT-FH++
reg = replace(as_men, status=np.where(male, as_men.status, as_women.status),
              age=np.where(male, as_men.age, as_women.age),
              onset=np.where(male, as_men.onset, as_women.onset))
g = reg.genetic                                 # the truth, known only because we simulated it

print(f"{len(reg.ids)} people ({male.sum()} men), {reg.status.sum()} diagnosed by age 70")
print(f"diagnosed: men {reg.status[male].mean():.1%}, women {reg.status[~male].mean():.1%}")

row = {p: i for i, p in enumerate(reg.ids)}                      # id -> row number
years, n_born = np.unique(reg.birth_time, return_counts=True)    # people born in each year
couples = Counter((f, m) for f, m in zip(reg.father, reg.mother) if f in row and m in row)
sizes, n_couples = np.unique(list(couples.values()), return_counts=True)   # children per couple

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.5))
ax1.bar(years, n_born, width=6)       # how many people were born in each year
ax1.set_xlabel("birth year")
ax1.set_ylabel("people")
ax2.bar(sizes, n_couples)             # how many couples have each number of children
ax2.set_xlabel("children per couple")
ax2.set_ylabel("couples")
plt.show()
print("born per year:", dict(zip(years.astype(int).tolist(), n_born.tolist())))
print("couples by number of children:", dict(zip(sizes.tolist(), n_couples.tolist())))

bins = np.arange(0, 75, 5)                                      # 5-year age bins
fig, ax = plt.subplots(figsize=(6, 3.5))
ax.hist(reg.onset[male & reg.status], bins=bins, alpha=0.6, label="men")
ax.hist(reg.onset[~male & reg.status], bins=bins, alpha=0.6, label="women")   # the women's ages at diagnosis
ax.set_xlabel("age at diagnosis")
ax.set_ylabel("cases")
ax.legend()
plt.show()
print(f"cases: {(male & reg.status).sum()} men, {(~male & reg.status).sum()} women; median age at "
      f"diagnosis {np.median(reg.onset[male & reg.status]):.0f} and {np.median(reg.onset[~male & reg.status]):.0f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
for label, keep, cip, colour in (("men", male, CIP_M, "C0"), ("women", ~male, CIP_F, "C1")):
    # Kaplan-Meier: everyone enters at birth and leaves at diagnosis or at 70
    km = kaplan_meier_cip(np.zeros(keep.sum()), reg.age[keep], reg.status[keep])
    print(f"{label}: {km.values[km.ages <= 50][-1]:.1%} diagnosed by 50, {km.values[-1]:.1%} by 70")
    ages, inc = np.r_[0, km.ages], np.r_[0, km.values]      # start both curves at age 0
    ax1.step(ages, 100 * (1 - inc), where="post", color=colour, label=label)   # survival S(t)
    ax2.step(ages, 100 * inc, where="post", color=colour, label=f"{label}: observed")
    ax2.plot(AGES[:71], 100 * cip[:71], ls="--", color=colour, label=f"{label}: true curve")
ax1.set_title("Kaplan-Meier survival S(t)")
ax1.set_ylabel("% still undiagnosed")
ax2.set_title("cumulative incidence = 1 - S(t)")
ax2.set_ylabel("% diagnosed by this age")
for ax in (ax1, ax2):
    ax.set_xlabel("age")
    ax.legend()
plt.show()

row = {p: i for i, p in enumerate(reg.ids)}                  # id -> row number
has_parents = np.array([f in row for f in reg.father])       # parents recorded in the register
parent_dx = np.array([any(reg.status[row[p]] for p in (f, m) if p in row)
                      for f, m in zip(reg.father, reg.mother)])   # a parent diagnosed by 70

fig, ax = plt.subplots(figsize=(6.5, 4))
for label, keep in (("a diagnosed parent", parent_dx),
                    ("no diagnosed parent", has_parents & ~parent_dx)):   # parents recorded, neither diagnosed
    km = kaplan_meier_cip(np.zeros(keep.sum()), reg.age[keep], reg.status[keep])
    ax.step(km.ages, 100 * km.values, where="post", label=f"{label} (n = {keep.sum()})")
    print(f"{label}: n = {keep.sum()}, {km.values[-1]:.1%} diagnosed by 70")
ax.set_xlabel("age")
ax.set_ylabel("% diagnosed by this age")
ax.legend()
plt.show()

fig, ax = plt.subplots(figsize=(6, 3.5))
ax.hist(g[~reg.status], bins=50, density=True, alpha=0.5, label="not diagnosed by 70")
ax.hist(g[reg.status], bins=50, density=True, alpha=0.5, label="diagnosed by 70")    # g of the cases
ax.set_xlabel("true genetic liability g")
ax.legend()
plt.show()
print(f"mean g: cases {g[reg.status].mean():.2f}, non-cases {g[~reg.status].mean():.2f}; "
      f"of the people with g > 1, {reg.status[g > 1].mean():.0%} were diagnosed")

curves = {"M": (AGES, CIP_M, 0.12), "F": (AGES, CIP_F, 0.08)}   # sex -> (ages, curve, lifetime K)

scores = estimate_liabilities(
    reg.ids, reg.father, reg.mother,          # the pedigree: who is whose parent
    probands=reg.ids,                         # whom to score: everyone
    status=reg.status, age=reg.age,           # diagnosed by 70? age at diagnosis or at 70
    strata=sex, cip_by_stratum=curves,        # each sex has its own incidence curve, so its own thresholds
    h2=H2,                                    # heritability: you supply it, the scorer never estimates it
    use="gwas")                               # also use each person's own diagnosis

est, var = scores.est, scores.var
print(f"corr(score, g) = {corr(est, g):.3f}")
print(f"var(score) = {est.var():.3f}  +  mean posterior variance = {var.mean():.3f}"
      f"  =  {est.var() + var.mean():.3f}")

free40 = reg.onset > 40                     # still undiagnosed on their 40th birthday


def score(h2):
    """Everyone's score from all records up to 70 (use="gwas"): (mean, variance)."""
    s = estimate_liabilities(reg.ids, reg.father, reg.mother, probands=reg.ids,
                             status=reg.status, age=reg.age, strata=sex, cip_by_stratum=curves,
                             h2=h2, use="gwas")
    return s.est, s.var


def score_at_40(h2):
    """Scores as known on each person's 40th birthday (use="prediction"), for the
    people still undiagnosed then: their own status and every later record are hidden."""
    s = estimate_liabilities(reg.ids, reg.father, reg.mother,
                             probands=[p for p, keep in zip(reg.ids, free40) if keep],
                             status=reg.status, age=reg.age, strata=sex, cip_by_stratum=curves,
                             h2=h2, use="prediction",
                             birth_time=reg.birth_time,                  # everyone's birth date
                             index_time=(reg.birth_time + 40)[free40])   # each proband's 40th birthday
    return s.est, s.var


est40, var40 = score_at_40(H2)
print(f"{free40.sum()} people undiagnosed at 40; corr(score at 40, g) = {corr(est40, g[free40]):.3f}")

fig, ax = plt.subplots(figsize=(5, 4))
ax.scatter(est, g, s=3, alpha=0.4)             # x: the score est, y: the truth g
ax.set_xlabel("score (posterior mean)")
ax.set_ylabel("true genetic liability g")
ax.set_title(f"corr = {corr(est, g):.2f}")
plt.show()

# people still undiagnosed at 70 with both parents recorded; their status and age match,
# but their own sex-specific thresholds and other family records can still differ
undiagnosed = ~reg.status & has_parents
n_dx = np.array([sum(reg.status[row[p]] for p in (f, m) if p in row)
                 for f, m in zip(reg.father, reg.mother)])            # diagnosed parents: 0, 1 or 2
parent_onset = np.array([min([reg.onset[row[p]] for p in (f, m) if p in row and reg.status[row[p]]], default=np.nan)
                         for f, m in zip(reg.father, reg.mother)])    # age of the (first) diagnosed parent
one = undiagnosed & (n_dx == 1)                                              # exactly one diagnosed parent

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
ax1.boxplot([est[undiagnosed & (n_dx == k)] for k in (0, 1, 2)])     # the scores of the people in undiagnosed with k diagnosed parents
ax1.set_xticks([1, 2, 3], ["0", "1", "2"])
ax1.set_xlabel("diagnosed parents")
ax1.set_ylabel("score")
ax2.boxplot([est[one & (parent_onset < 50)], est[one & (parent_onset >= 50) & (parent_onset < 60)],
             est[one & (parent_onset >= 60)]])
ax2.set_xticks([1, 2, 3], ["before 50", "50 to 59", "60 to 70"])
ax2.set_xlabel("age at which the parent was diagnosed")
plt.show()
for k in (0, 1, 2):
    print(f"{k} diagnosed parents: n = {(undiagnosed & (n_dx == k)).sum()}, mean score {est[undiagnosed & (n_dx == k)].mean():.2f}")
for label, lo, hi in (("before 50", 0, 50), ("50 to 59", 50, 60), ("60 to 70", 60, 71)):
    sel = one & (parent_onset >= lo) & (parent_onset < hi)
    print(f"parent diagnosed {label}: n = {sel.sum()}, mean score {est[sel].mean():.2f}")

y = reg.status[free40]      # diagnosed between 40 and 70, one entry per person in est40
print(f"{y.sum()} of {len(y)} people undiagnosed at 40 were diagnosed by 70")

fig, ax = plt.subplots(figsize=(6, 3.5))
ax.hist(est40[~y], bins=40, density=True, alpha=0.5, label="not diagnosed")
ax.hist(est40[y], bins=40, density=True, alpha=0.5, label="diagnosed 40 to 70")   # the cases' scores
ax.set_xlabel("score at 40")
ax.legend()
plt.show()

print(f"AUC score at 40        {auc(est40, y):.3f}")
print(f"AUC true g             {auc(g[free40], y):.3f}")
print(f"AUC use='gwas' score   {auc(est[free40], y):.3f}")

cip40 = np.where(male, np.interp(40, AGES, CIP_M), np.interp(40, AGES, CIP_F))[free40]
cip70 = np.where(male, np.interp(70, AGES, CIP_M), np.interp(70, AGES, CIP_F))[free40]
T40, T70 = norm.isf(cip40), norm.isf(cip70)            # each person's LT-FH++ thresholds at 40 and 70
def risk_40_to_70(mean_g, var_g, h2):
    """Normal approximation to P(diagnosed 40–70 | relatives, undiagnosed at 40)."""
    sd = np.sqrt(var_g + 1 - h2)          # genetic uncertainty + environmental variance
    below40 = norm.cdf((T40 - mean_g) / sd)
    below70 = norm.cdf((T70 - mean_g) / sd)
    return (below40 - below70) / below40  # condition on being undiagnosed at 40


risk = risk_40_to_70(est40, var40, H2)

fifths = np.array_split(np.argsort(est40, kind="stable"), 5)   # lowest to highest score
observed = [y[idx].mean() for idx in fifths]                   # share diagnosed in each fifth
predicted = [risk[idx].mean() for idx in fifths]               # mean risk in each fifth

fig, ax = plt.subplots(figsize=(6, 3.5))
ax.bar(np.arange(5) - 0.2, observed, width=0.4, label="observed")
ax.bar(np.arange(5) + 0.2, predicted, width=0.4, label="predicted risk")
ax.set_xticks(range(5), ["lowest", "2nd", "middle", "4th", "highest"])
ax.set_xlabel("fifth of the score at 40")
ax.set_ylabel("share diagnosed 40 to 70")
ax.legend()
plt.show()
print(f"mean predicted risk {risk.mean():.3f}   observed {y.mean():.3f}")
for name, o, p in zip(("lowest", "2nd", "middle", "4th", "highest"), observed, predicted):
    print(f"{name:>8} fifth: observed {o:.3f}, predicted {p:.3f}")

row = {p: i for i, p in enumerate(reg.ids)}          # id -> row number
children = {}                                        # (father, mother) -> rows of their children
for i, (f, m) in enumerate(zip(reg.father, reg.mother)):
    if f in row and m in row:
        children.setdefault((f, m), []).append(i)


def first_degree(i):
    """Rows of person i's parents and full siblings."""
    f, m = reg.father[i], reg.mother[i]
    parents = [row[p] for p in (f, m) if p in row]
    siblings = [j for j in children.get((f, m), []) if j != i]
    return np.array(parents + siblings, dtype=int)


fdr = [first_degree(i) for i in range(len(reg.ids))]
diag_time = reg.birth_time + reg.onset               # calendar time of each diagnosis
birth40 = reg.birth_time + 40                        # calendar time of each 40th birthday

# yes/no family history: any affected parent or sibling, by 70 and by one's own 40th birthday
fh = np.array([reg.status[r].any() for r in fdr])
fh40 = np.array([(reg.status[r] & (diag_time[r] <= birth40[i])).any() for i, r in enumerate(fdr)])
print(f"family-history positive: {fh.sum()} by age 70, {fh40.sum()} at their 40th birthday")

names = ["own status", "FH indicator", "score (use='gwas')"]
r2 = [corr(x, g) ** 2 for x in (reg.status, fh, est)]      # R²: the squared correlation of each with g

fig, ax = plt.subplots(figsize=(6, 3.2))
ax.bar(names, r2)
ax.set_ylabel("R² with the true g")
plt.show()
print(f"prediction at 40: AUC FH indicator {auc(fh40[free40], y):.3f}, AUC score {auc(est40, y):.3f}")
print("R2:", ", ".join(f"{n} {v:.3f}" for n, v in zip(names, r2)))

from ltpred import tetrachoric

# every parent-child pair in the register, as row numbers
pairs = [(row[p], i) for i, (f, m) in enumerate(zip(reg.father, reg.mother))
         for p in (f, m) if p in row]
par, kid = np.array(pairs).T


def table(a, b):
    """2x2 table of two True/False arrays: both, first only, second only, neither."""
    return np.array([(a & b).sum(), (a & ~b).sum(), (~a & b).sum(), (~a & ~b).sum()])

p_status, k_status = reg.status[par], reg.status[kid]
t = tetrachoric(p_status, k_status)
print(f"{len(par)} pairs, table {table(p_status, k_status)}")
print(f"h2 = 2 rho = {2 * t.rho:.2f}   (generating value {H2})")
print(f"naive SE = {2 * t.se:.2f}; treats pairs as independent, not family-adjusted")

from ltpred import liability_to_observed_h2

Ps = np.linspace(0.02, 0.6, 50)                # case fraction in the sample
h2_obs = [float(liability_to_observed_h2(H2, K, P)) for P in Ps]   # observed-scale h² for each P

fig, ax = plt.subplots(figsize=(6, 3.5))
ax.plot(Ps, h2_obs)
ax.axhline(H2, ls="--", color="grey", label="liability scale: 0.5")
ax.axvline(K, ls=":", color="C1", label="population sample: P = K")
ax.set_xlabel("case fraction in the sample, P")
ax.set_ylabel("observed-scale h²")
ax.legend()
plt.show()
print(f"population sample {float(liability_to_observed_h2(H2, K, None)):.3f}, "
      f"1:1 study {float(liability_to_observed_h2(H2, K, 0.5)):.3f}")

h2_wrong = float(liability_to_observed_h2(H2, K, None))  # illustrative wrong-scale input
est_wrong, var_wrong = score(h2_wrong)
print(f"corr(score, g): correct h² {corr(est, g):.3f}; wrong h² {corr(est_wrong, g):.3f}")
print(f"slope of g on score: correct {np.polyfit(est, g, 1)[0]:.2f}; "
      f"wrong {np.polyfit(est_wrong, g, 1)[0]:.2f} (1 means correct scale)")

est40_wrong, var40_wrong = score_at_40(h2_wrong)
risk_wrong = risk_40_to_70(est40_wrong, var40_wrong, h2_wrong)
print(f"AUC at 40: correct h² {auc(est40, y):.4f}; wrong h² {auc(est40_wrong, y):.4f}")
print(f"mean risk: correct {risk.mean():.3f}; wrong {risk_wrong.mean():.3f}; observed {y.mean():.3f}")
# Keep the original Q9 groups so we compare risks for the same people.
for label, idx in (("lowest", fifths[0]), ("highest", fifths[-1])):
    print(f"{label} fifth: correct risk {risk[idx].mean():.3f}; "
          f"wrong risk {risk_wrong[idx].mean():.3f}; observed {y[idx].mean():.3f}")


import platform, scipy
print(f"Python {platform.python_version()}, NumPy {np.__version__}, "
      f"SciPy {scipy.__version__}, ltpred {ltpred.__version__}")
