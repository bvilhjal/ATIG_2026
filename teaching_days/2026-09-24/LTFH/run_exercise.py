"""ATiG 2026, 24 September: LT-FH exercise with ltpred -- every answer, as a script.

Rebuilds tutorial steps 1, 3 and 4 with the tutorial's seeds, adds the larger
register, then prints the answers to Parts A-E of ``LTFH_exercise.ipynb`` in
order. Needs Python >= 3.9 with NumPy, SciPy and ltpred >= 0.7.0 installed from
source (https://github.com/bvilhjal/ltpred); no other files, same on Windows,
macOS and Linux. Runs from any directory in about two minutes and needs under
1 GB of memory. The figures are only in the notebook.
"""
import numpy as np
from scipy.stats import norm, rankdata

import ltpred
print("ltpred", ltpred.__version__)
from ltpred import simulate_pedigree, simulate_register_liabilities

H2 = 0.5                                      # liability-scale heritability
CIP_K, CIP_MID, CIP_SLOPE = 0.10, 60.0, 1.0 / 8.0   # logistic onset curve: 10% lifetime, half of it by age 60
EVAL_AGE = 70.0                               # everyone is followed to here
INDEX_AGE = 40.0                              # the prospective cut in step 4

# the generating cumulative-incidence curve, on a 1-year age grid
AGE_GRID = np.arange(0, 121, 1.0)
TRUE_CIP = CIP_K / (1.0 + np.exp((CIP_MID - AGE_GRID) * CIP_SLOPE))

ids, father, mother = simulate_pedigree(
    np.random.default_rng(20260921), n_founder_pairs=100, gens=2)

cohort = simulate_register_liabilities(
    np.random.default_rng(20260921), ids, father, mother,
    h2=H2, cip_ages=AGE_GRID, cip_values=TRUE_CIP, eval_age=EVAL_AGE)

print(f"{len(ids)} people, {int(cohort.status.sum())} diagnosed by age "
      f"{EVAL_AGE:.0f} ({cohort.status.mean():.1%})")
print(f"var(true genetic liability) = {cohort.genetic.var():.4f}  (target {H2})")
from ltpred import estimate_liabilities

scores = estimate_liabilities(
    cohort.ids, cohort.father, cohort.mother,
    probands=cohort.ids,
    status=cohort.status.astype(int), age=cohort.age,
    use="gwas",                          # the proband's own diagnosis is used
    cip_ages=AGE_GRID, cip_values=TRUE_CIP, k_pop=CIP_K,
    h2=H2)

est = np.asarray(scores.est)
print(f"corr(score, true genetic liability) = {np.corrcoef(est, cohort.genetic)[0, 1]:.4f}")
print(f"score sd = {est.std():.4f}   mean posterior variance = {np.asarray(scores.var).mean():.4f}")
at_risk = cohort.onset > INDEX_AGE          # still undiagnosed at the cut
index_time = cohort.birth_time + INDEX_AGE
probands_at_risk = [p for p, keep in zip(cohort.ids, at_risk) if keep]

predicted = estimate_liabilities(
    cohort.ids, cohort.father, cohort.mother,
    probands=probands_at_risk,
    status=cohort.status.astype(int), age=cohort.age,
    use="prediction",                    # the proband's own diagnosis is hidden
    cip_ages=AGE_GRID, cip_values=TRUE_CIP, k_pop=CIP_K, h2=H2,
    birth_time=cohort.birth_time, index_time=index_time[at_risk])

pred_est = np.asarray(predicted.est)
print(f"{int(at_risk.sum())} of {len(ids)} probands are disease-free at age {INDEX_AGE:.0f}")
print(f"corr(score, truth) = {np.corrcoef(pred_est, cohort.genetic[at_risk])[0, 1]:.4f}"
      f"   score sd = {pred_est.std():.4f}")

print("\n== Part 0 -- the larger register")
big_ids, big_father, big_mother = simulate_pedigree(
    np.random.default_rng(1), n_founder_pairs=500, gens=2)
big = simulate_register_liabilities(
    np.random.default_rng(1), big_ids, big_father, big_mother,
    h2=H2, cip_ages=AGE_GRID, cip_values=TRUE_CIP, eval_age=EVAL_AGE)


def score_big(use, h2=H2):
    """Score `big` as tutorial step 3 (use="gwas") or step 4 (use="prediction")."""
    prospective = use == "prediction"
    keep = big.onset > INDEX_AGE if prospective else np.ones(len(big.ids), bool)
    return estimate_liabilities(
        big.ids, big.father, big.mother,
        probands=[p for p, k in zip(big.ids, keep) if k],
        status=big.status.astype(int), age=big.age, use=use,
        cip_ages=AGE_GRID, cip_values=TRUE_CIP, k_pop=CIP_K, h2=h2,
        birth_time=big.birth_time if prospective else None,
        index_time=(big.birth_time + INDEX_AGE)[keep] if prospective else None)


big_gwas = score_big("gwas")                 # every person, own diagnosis used
big_pred = score_big("prediction")           # disease-free at 40, own diagnosis hidden
big_est = np.asarray(big_gwas.est)
big_at_risk = big.onset > INDEX_AGE
big_m, big_v = np.asarray(big_pred.est), np.asarray(big_pred.var)

print(f"{len(big.ids)} people, {int(big.status.sum())} diagnosed by age 70 "
      f"({big.status.mean():.1%}); var(genetic) = {big.genetic.var():.4f}")
print(f"use='gwas':       corr(score, truth) = {np.corrcoef(big_est, big.genetic)[0, 1]:.4f}")
print(f"use='prediction': corr(score, truth) = "
      f"{np.corrcoef(big_m, big.genetic[big_at_risk])[0, 1]:.4f} on "
      f"{int(big_at_risk.sum())} probands disease-free at 40")

print("\n== Part A -- the same register scored under three assumed heritabilities (Q2)")
by_h2 = {}
print(f"{'assumed h2':>10} {'corr':>7} {'var(score)':>11} {'mean var':>9} {'sum':>7} {'slope':>7}")
for h2_assumed in (0.2, 0.5, 0.8):
    s = estimate_liabilities(
        cohort.ids, cohort.father, cohort.mother, probands=cohort.ids,
        status=cohort.status.astype(int), age=cohort.age, use="gwas",
        cip_ages=AGE_GRID, cip_values=TRUE_CIP, k_pop=CIP_K,
        h2=h2_assumed)
    e, v = np.asarray(s.est), np.asarray(s.var)
    corr = np.corrcoef(e, cohort.genetic)[0, 1]
    slope = np.polyfit(e, cohort.genetic, 1)[0]
    by_h2[h2_assumed] = (e, slope)
    print(f"{h2_assumed:>10.1f} {corr:>7.4f} {e.var():>11.4f} {v.mean():>9.4f} "
          f"{e.var() + v.mean():>7.4f} {slope:>7.3f}")

print("\n== Part B -- who is diagnosed between 40 and 70? (Q4)")
y = (big_at_risk & big.status)[big_at_risk]
print(f"{int(y.sum())} incident cases among {len(y)} probands ({y.mean():.1%})")

print("\n== Q5")
def auc(x, y):
    r = rankdata(x)
    n1, n0 = y.sum(), (~y).sum()
    return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)

def auc_se(x, y):
    # Hanley & McNeil (1982) standard error
    a = auc(x, y)
    n1, n0 = y.sum(), (~y).sum()
    q1, q2 = a / (2 - a), 2 * a**2 / (1 + a)
    return np.sqrt((a * (1 - a) + (n1 - 1) * (q1 - a**2) + (n0 - 1) * (q2 - a**2)) / (n1 * n0))

for label, x in (("prospective score", big_m),
                 ("true genetic liability", big.genetic[big_at_risk]),
                 ("use='gwas' score, same people", big_est[big_at_risk])):
    print(f"AUC {label:30s} = {auc(x, y):.3f} +/- {auc_se(x, y):.3f}")

print("\n== Q7")
m, v = big_m, big_v
T40 = norm.ppf(1.0 - np.interp(INDEX_AGE, AGE_GRID, TRUE_CIP))
T70 = norm.ppf(1.0 - np.interp(EVAL_AGE, AGE_GRID, TRUE_CIP))
s = np.sqrt(v + 1.0 - H2)
p_free_at_40 = norm.cdf((T40 - m) / s)
risk = (p_free_at_40 - norm.cdf((T70 - m) / s)) / p_free_at_40

print(f"mean model risk {risk.mean():.4f}   observed incidence {y.mean():.4f}")

fifths = np.array_split(np.argsort(m, kind="stable"), 5)
for label, g in zip(("lowest", "second", "middle", "fourth", "highest"), fifths):
    print(f"{label:>8} fifth  n={len(g)}  events={int(y[g].sum()):3d}  "
          f"observed {y[g].mean():.3f}  model risk {risk[g].mean():.3f}")

print("\n== Part C -- the score against a 0/1 indicator (Q8)")
big_idx = {p: i for i, p in enumerate(big.ids)}
sibs = {}                                     # children grouped by (father, mother)
for i, (fa, mo) in enumerate(zip(big.father, big.mother)):
    if fa in big_idx and mo in big_idx:
        sibs.setdefault((fa, mo), []).append(i)


def first_degree(i):
    """Row indices of person i's parents and full siblings present in the register."""
    fa, mo = big.father[i], big.mother[i]
    rels = [big_idx[p] for p in (fa, mo) if p in big_idx]
    if fa in big_idx and mo in big_idx:
        rels += [j for j in sibs[(fa, mo)] if j != i]
    return np.asarray(rels, int)


relatives = [first_degree(i) for i in range(len(big.ids))]
big_index_time = big.birth_time + INDEX_AGE

n_fdr = np.array([big.status[r].sum() for r in relatives], float)
n_fdr_40 = np.array([(big.birth_time[r] + big.onset[r] <= big_index_time[i]).sum()
                     for i, r in enumerate(relatives)], float)
fh01, fh01_40 = (n_fdr > 0).astype(float), (n_fdr_40 > 0).astype(float)
print(f"FH-positive: {int(fh01.sum())} of {len(fh01)} with follow-up to 70, "
      f"{int(fh01_40.sum())} at the proband's 40th birthday")

for name, x in (("own status", big.status.astype(float)), ("FH indicator 0/1", fh01),
                ("# affected FDRs", n_fdr), ("use='gwas' score", big_est)):
    print(f"R2  {name:18s} {np.corrcoef(x, big.genetic)[0, 1] ** 2:.3f}")
for name, x in (("FH indicator 0/1 at 40", fh01_40[big_at_risk]),
                ("# affected FDRs at 40", n_fdr_40[big_at_risk]),
                ("prospective score", big_m)):
    print(f"AUC {name:22s} {auc(x, y):.3f}")

print("\n== Part D -- h2 from parent-offspring pairs (Q9)")
from ltpred.tetrachoric import tetrachoric, tetrachoric_table


def parent_offspring_status(reg):
    idx = {p: i for i, p in enumerate(reg.ids)}
    parent, child = [], []
    for kid, fa, mo in zip(reg.ids, reg.father, reg.mother):
        for par in (fa, mo):
            if par in idx:
                parent.append(reg.status[idx[par]])
                child.append(reg.status[idx[kid]])

    return np.asarray(parent), np.asarray(child)


def table(parent, child):
    """Counts: both diagnosed, parent only, child only, neither."""
    return np.array([(parent & child).sum(), (parent & ~child).sum(),
                     (~parent & child).sum(), (~parent & ~child).sum()])


for label, reg in (("tutorial register", cohort), ("larger register", big)):
    par, kid = parent_offspring_status(reg)
    t = tetrachoric(par, kid)
    print(f"{label:17s}: {len(par):5d} pairs, table {table(par, kid)}  ->  "
          f"h2 = {2 * t.rho:+.2f} +/- {2 * t.se:.2f}   (truth {H2})")

print("\n== Q11 -- iPSYCH-style case-cohort")
Q_KEEP = 0.2
keep = big.status | (np.random.default_rng(20260924).random(len(big.ids)) < Q_KEEP)
pi = np.where(big.status, 1.0, Q_KEEP)                  # inclusion probability

par, kid, w = [], [], []
for i in np.flatnonzero(keep):
    for p in (big.father[i], big.mother[i]):
        if p in big_idx and keep[big_idx[p]]:
            par.append(big.status[big_idx[p]])
            kid.append(big.status[i])
            w.append(1.0 / (pi[big_idx[p]] * pi[i]))
par, kid, w = np.asarray(par), np.asarray(kid), np.asarray(w)

t = tetrachoric(par, kid)
print(f"kept {int(keep.sum())} of {len(big.ids)}: case rate "
      f"{big.status[keep].mean():.1%} vs {big.status.mean():.1%} in the population")
print(f"{len(par)} pairs, table {table(par, kid)}  ->  naive h2 = {2 * t.rho:+.2f} +/- {2 * t.se:.2f}")

weighted = np.array([w[par & kid].sum(), w[par & ~kid].sum(),
                     w[~par & kid].sum(), w[~par & ~kid].sum()])
tw = tetrachoric_table(*np.round(weighted).astype(int))
t_pop = tetrachoric(*parent_offspring_status(big))
print(f"weighted table {np.round(weighted).astype(int)}  ->  h2 = {2 * tw.rho:+.2f}")
print(f"big's own table {table(*parent_offspring_status(big))}  ->  h2 = {2 * t_pop.rho:+.2f}")
# the weighted table has the population's size, so its own SE is far too small;
# rescale the cells to the number of pairs actually observed for an honest one
scaled = np.round(weighted / weighted.sum() * len(par)).astype(int)
print(f"SE at the {len(par)} pairs actually kept: +/- {2 * tetrachoric_table(*scaled).se:.2f}")

print("\n== Part E -- observed and liability scale (Q12)")
from ltpred import observed_to_liability_h2, liability_to_observed_h2

for name, P in (("population sample", None), ("1:1 case/control", 0.5), ("1:4 case/control", 0.2)):
    obs = liability_to_observed_h2(H2, CIP_K, P)
    back = observed_to_liability_h2(obs, CIP_K, P)
    print(f"{name:18s} P={CIP_K if P is None else P:.2f}  h2_obs={float(obs):.4f}  round trip={float(back):.4f}")

print("\n== Q13")
def score_cohort(h2):
    s = estimate_liabilities(
        cohort.ids, cohort.father, cohort.mother, probands=cohort.ids,
        status=cohort.status.astype(int), age=cohort.age, use="gwas",
        cip_ages=AGE_GRID, cip_values=TRUE_CIP, k_pop=CIP_K, h2=h2)
    e = np.asarray(s.est)
    return np.corrcoef(e, cohort.genetic)[0, 1], np.polyfit(e, cohort.genetic, 1)[0]


obs_population = float(liability_to_observed_h2(H2, CIP_K, None))
corr, slope = score_cohort(obs_population)
print(f"h2 = {obs_population:.4f} unconverted: corr = {corr:.4f}  slope = {slope:.3f}")

for P in (0.5, 0.2):
    forgotten = float(observed_to_liability_h2(liability_to_observed_h2(H2, CIP_K, P), CIP_K))
    try:
        corr, slope = score_cohort(forgotten)
        print(f"P = {P}: forgetting the P factor gives h2 = {forgotten:.2f}, accepted: "
              f"corr = {corr:.4f}  slope = {slope:.3f}")
    except ValueError as err:
        print(f"P = {P}: forgetting the P factor gives h2 = {forgotten:.2f}, rejected: "
              f"{str(err).split(' -- ')[0]}")

print("\n== Q14 -- effective sample sizes")
for name, n_case, n_ctrl in (("big, population", int(big.status.sum()), len(big.ids) - int(big.status.sum())),
                             ("1:1 study", 2_500, 2_500),
                             ("biobank 1:9", 5_000, 45_000),
                             ("consortium 1:3", 50_000, 150_000)):
    n = n_case + n_ctrl
    neff = 4 * n_case * n_ctrl / n
    print(f"{name:16s} N={n:7d}  P={n_case / n:.3f}  Neff={neff:8.0f}  Neff/N={neff / n:.3f}")
import platform, scipy
print(f"Python {platform.python_version()}  NumPy {np.__version__}  "
      f"SciPy {scipy.__version__}  ltpred {ltpred.__version__}")
