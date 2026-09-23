"""ATiG 2026, 24 September: LT-FH exercise with ltpred -- every answer as one script.

Prints the answers to LTFH_exercise.ipynb in order (the figures are only in the
notebook). Needs Python >= 3.9 with NumPy, SciPy and ltpred; runs in under a minute.
"""
import numpy as np
from scipy.stats import norm, rankdata

import ltpred
from ltpred import simulate_pedigree, simulate_register_liabilities, estimate_liabilities
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

H2 = 0.5                                    # true liability-scale heritability
K = 0.10                                    # lifetime prevalence
AGES = np.arange(0, 121.0)                  # ages 0, 1, ..., 120
CIP = K / (1 + np.exp((60 - AGES) / 8))     # cumulative incidence by age: half of K by age 60

ids, father, mother = simulate_pedigree(np.random.default_rng(1), n_founder_pairs=500, gens=2)
reg = simulate_register_liabilities(np.random.default_rng(1), ids, father, mother,
                                    h2=H2, cip_ages=AGES, cip_values=CIP, eval_age=70)
g = reg.genetic                             # the truth, known only because we simulated it

print(f"{len(reg.ids)} people, {reg.status.sum()} diagnosed by age 70 ({reg.status.mean():.1%})")
print(f"var(g) = {g.var():.3f}   (target {H2})")

scores = estimate_liabilities(
    reg.ids, reg.father, reg.mother,          # the pedigree: who is whose parent
    probands=reg.ids,                         # whom to score: everyone
    status=reg.status, age=reg.age,           # diagnosed by 70? age at diagnosis or at 70
    cip_ages=AGES, cip_values=CIP, k_pop=K,   # incidence curve: turns an age into a liability threshold
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
                             status=reg.status, age=reg.age,
                             cip_ages=AGES, cip_values=CIP, k_pop=K, h2=h2, use="gwas")
    return s.est, s.var


def score_at_40(h2):
    """Scores as known on each person's 40th birthday (use="prediction"), for the
    people still undiagnosed then: their own status and every later record are hidden."""
    s = estimate_liabilities(reg.ids, reg.father, reg.mother,
                             probands=[p for p, keep in zip(reg.ids, free40) if keep],
                             status=reg.status, age=reg.age,
                             cip_ages=AGES, cip_values=CIP, k_pop=K, h2=h2, use="prediction",
                             birth_time=reg.birth_time,                  # everyone's birth date
                             index_time=(reg.birth_time + 40)[free40])   # each proband's 40th birthday
    return s.est, s.var


est40, var40 = score_at_40(H2)
print(f"{free40.sum()} people undiagnosed at 40; corr(score at 40, g) = {corr(est40, g[free40]):.3f}")

print("\n== Q2")
by_h2 = {}                                   # kept for the figure below
print("assumed h2   corr  var(score)  mean var    sum  slope")
for h2 in (0.2, 0.5, 0.8):
    e, v = score(h2)
    c = corr(e, g)
    b = np.polyfit(e, g, 1)[0]
    by_h2[h2] = e
    print(f"{h2:10.1f} {c:6.3f} {e.var():10.3f} {v.mean():9.3f} {e.var() + v.mean():6.3f} {b:6.2f}")

print("\n== Q4")
y = reg.status[free40]       # everyone in est40 was undiagnosed at 40, so "by 70" means "40 to 70"
print(f"{y.sum()} incident cases among {len(y)} people ({y.mean():.1%})")

print("\n== Q5")
print(f"AUC score at 40             {auc(est40, y):.3f}")
print(f"AUC true g                  {auc(g[free40], y):.3f}")
print(f"AUC use='gwas' score        {auc(est[free40], y):.3f}")

print("\n== Q6")
print(f"var(g)                        {g.var():.3f}")
print(f"var(score), use='gwas'        {est.var():.3f}")
print(f"var(score at 40)              {est40.var():.3f}")

print("\n== Q7")
T40, T70 = norm.isf(np.interp([40, 70], AGES, CIP))   # liability thresholds at ages 40 and 70
sd = np.sqrt(var40 + 1 - H2)                          # spread of full liability given the relatives
below_T40 = norm.cdf((T40 - est40) / sd)              # P(liability < T40): undiagnosed at 40
below_T70 = norm.cdf((T70 - est40) / sd)              # P(liability < T70): undiagnosed at 70
risk = (below_T40 - below_T70) / below_T40

print(f"mean risk {risk.mean():.3f}   observed incidence {y.mean():.3f}")
fifths = np.array_split(np.argsort(est40, kind="stable"), 5)    # lowest to highest score
for name, idx in zip(("lowest", "second", "middle", "fourth", "highest"), fifths):
    print(f"{name:>8} fifth: {y[idx].sum():3d} cases, observed {y[idx].mean():.3f}, "
          f"risk {risk[idx].mean():.3f}")

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

print("\n== Q8")
for name, x in (("own status", reg.status), ("FH indicator 0/1", fh), ("use='gwas' score", est)):
    print(f"R2  {name:18s} {corr(x, g) ** 2:.3f}")
for name, x in (("FH indicator at 40", fh40[free40]), ("score at 40", est40)):
    print(f"AUC {name:18s} {auc(x, y):.3f}")

from ltpred import tetrachoric

# every parent-child pair in the register, as row numbers
pairs = [(row[p], i) for i, (f, m) in enumerate(zip(reg.father, reg.mother))
         for p in (f, m) if p in row]
par, kid = np.array(pairs).T


def table(a, b):
    """2x2 table of two True/False arrays: both, first only, second only, neither."""
    return np.array([(a & b).sum(), (a & ~b).sum(), (~a & b).sum(), (~a & ~b).sum()])

print("\n== Q9")
p_status, k_status = reg.status[par], reg.status[kid]
t = tetrachoric(p_status, k_status)
print(f"{len(par)} pairs, table {table(p_status, k_status)}")
print(f"h2 = 2 rho = {2 * t.rho:.2f} +/- {2 * t.se:.2f}   (truth {H2})")

print("\n== Q11")
from ltpred import liability_to_observed_h2, observed_to_liability_h2

for name, P in (("population sample", None), ("1:1 case/control", 0.5), ("1:4 case/control", 0.2)):
    obs = liability_to_observed_h2(H2, K, P)
    back = observed_to_liability_h2(obs, K, P)
    print(f"{name:18s} h2_obs = {float(obs):.3f}   converted back = {float(back):.3f}")

import platform, scipy
print(f"Python {platform.python_version()}, NumPy {np.__version__}, "
      f"SciPy {scipy.__version__}, ltpred {ltpred.__version__}")
