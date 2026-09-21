"""ATiG 2026, 24 September: LT-FH exercise with ltpred -- every answer, as a script.

Rebuilds tutorial steps 1, 3 and 4 with the tutorial's seeds, then prints the
answers to Parts A-E of ``LTFH_exercise.ipynb``. Needs Python >= 3.9
with NumPy, SciPy and ltpred >= 0.7.0 installed from source
(https://github.com/bvilhjal/ltpred); no other files, same on Windows, macOS
and Linux. Runs from any directory in about a minute; Part C factorises a
10,000-person kinship matrix and needs a few GB of RAM.
"""
import platform

import numpy as np
import scipy
from scipy.stats import norm, rankdata

import ltpred
from ltpred import (estimate_liabilities, liability_to_observed_h2,
                    observed_to_liability_h2, simulate_pedigree,
                    simulate_register_liabilities)
from ltpred.tetrachoric import tetrachoric

# ------------------------------------------------------------ Part 0: the tutorial
H2 = 0.5
CIP_K, CIP_MID, CIP_SLOPE = 0.10, 60.0, 1.0 / 8.0
EVAL_AGE, INDEX_AGE = 70.0, 40.0
AGE_GRID = np.arange(0, 121, 1.0)
TRUE_CIP = CIP_K / (1.0 + np.exp((CIP_MID - AGE_GRID) * CIP_SLOPE))

ids, father, mother = simulate_pedigree(
    np.random.default_rng(20260921), n_founder_pairs=100, gens=2)
cohort = simulate_register_liabilities(
    np.random.default_rng(20260921), ids, father, mother,
    h2=H2, cip_ages=AGE_GRID, cip_values=TRUE_CIP, eval_age=EVAL_AGE)
truth = cohort.genetic
print(f"Part 0 -- {len(ids)} people, {int(cohort.status.sum())} diagnosed by "
      f"{EVAL_AGE:.0f}; var(genetic) = {truth.var():.4f}")


def score(h2, use="gwas", probands=None, index_time=None):
    return estimate_liabilities(
        cohort.ids, cohort.father, cohort.mother,
        probands=cohort.ids if probands is None else probands,
        status=cohort.status.astype(int), age=cohort.age, use=use,
        cip_ages=AGE_GRID, cip_values=TRUE_CIP, k_pop=CIP_K, h2=h2,
        birth_time=None if use == "gwas" else cohort.birth_time,
        index_time=index_time)


scores = score(H2)
est = np.asarray(scores.est)
print(f"step 3: corr(score, truth) = {np.corrcoef(est, truth)[0, 1]:.4f}")

at_risk = cohort.onset > INDEX_AGE
probands_at_risk = [p for p, keep in zip(cohort.ids, at_risk) if keep]
predicted = score(H2, use="prediction", probands=probands_at_risk,
                  index_time=(cohort.birth_time + INDEX_AGE)[at_risk])
pred_est = np.asarray(predicted.est)
print(f"step 4: corr(score, truth) = "
      f"{np.corrcoef(pred_est, truth[at_risk])[0, 1]:.4f}\n")

# ------------------------------------------------------------ Part A: wrong h2
print("Part A -- the same register scored under three assumed heritabilities")
print(f"{'assumed h2':>10} {'corr':>7} {'var(score)':>11} {'mean var':>9} "
      f"{'sum':>7} {'slope':>7}")
for h2_assumed in (0.2, 0.5, 0.8):
    s = score(h2_assumed)
    e, v = np.asarray(s.est), np.asarray(s.var)
    print(f"{h2_assumed:>10.1f} {np.corrcoef(e, truth)[0, 1]:>7.4f} "
          f"{e.var():>11.4f} {v.mean():>9.4f} {e.var() + v.mean():>7.4f} "
          f"{np.polyfit(e, truth, 1)[0]:>7.3f}")
print()

# ------------------------------------------------ Part B: from liability to risk
y = (at_risk & cohort.status)[at_risk]          # diagnosed in (40, 70]


def auc(x, y):
    r = rankdata(x)
    n1, n0 = y.sum(), (~y).sum()
    return (r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


print("Part B -- who is diagnosed between 40 and 70?")
print(f"Q4: {int(y.sum())} incident cases among {len(y)} probands ({y.mean():.1%})")
print(f"Q5: AUC prospective score         = {auc(pred_est, y):.3f}")
print(f"    AUC true genetic liability    = {auc(truth[at_risk], y):.3f}")
print(f"    AUC step-3 score, same people = {auc(est[at_risk], y):.3f}")

m, v = pred_est, np.asarray(predicted.var)
T40 = norm.ppf(1.0 - np.interp(INDEX_AGE, AGE_GRID, TRUE_CIP))
T70 = norm.ppf(1.0 - np.interp(EVAL_AGE, AGE_GRID, TRUE_CIP))
s = np.sqrt(v + 1.0 - H2)
p_free_at_40 = norm.cdf((T40 - m) / s)
risk = (p_free_at_40 - norm.cdf((T70 - m) / s)) / p_free_at_40
print(f"Q7: mean model risk {risk.mean():.4f}   observed incidence {y.mean():.4f}")
for label, g in zip(("low", "middle", "high"),
                    np.array_split(np.argsort(m, kind="stable"), 3)):
    print(f"    {label:>7}  n={len(g)}  observed {y[g].mean():.3f}  "
          f"model risk {risk[g].mean():.3f}")
print()

# ------------------------------------------------ Part C: h2 from the register
print("Part C -- h2 from parent-offspring pairs, twice the tetrachoric")


def parent_offspring_status(reg):
    idx = {p: i for i, p in enumerate(reg.ids)}
    parent, child = [], []
    for kid, fa, mo in zip(reg.ids, reg.father, reg.mother):
        for par in (fa, mo):
            if par in idx:
                parent.append(reg.status[idx[par]])
                child.append(reg.status[idx[kid]])
    return np.asarray(parent), np.asarray(child)


def report(label, reg):
    par, kid = parent_offspring_status(reg)
    t = tetrachoric(par, kid)
    print(f"{label}: {len(reg.ids)} people, {int(reg.status.sum())} diagnosed; "
          f"{len(par)} pairs, {int((par & kid).sum())} both diagnosed")
    print(f"    rho = {t.rho:+.3f} +/- {t.se:.3f}  ->  h2 = {2 * t.rho:+.2f} "
          f"+/- {2 * t.se:.2f}   (truth {H2})")


report("tutorial register", cohort)
big_ids, big_fa, big_mo = simulate_pedigree(
    np.random.default_rng(1), n_founder_pairs=1000, gens=2)
big = simulate_register_liabilities(
    np.random.default_rng(1), big_ids, big_fa, big_mo,
    h2=H2, cip_ages=AGE_GRID, cip_values=TRUE_CIP, eval_age=EVAL_AGE)
report("larger register  ", big)

# --------------------------- Part D: scales and effective sample sizes
print("\nPart D -- observed scale, liability scale, effective sample size")
n_cases = int(cohort.status.sum())
for name, P in (("register as sampled", n_cases / len(cohort.ids)),
                ("balanced case/control", 0.5), ("population sample", None)):
    obs = liability_to_observed_h2(H2, CIP_K, P)
    back = observed_to_liability_h2(obs, CIP_K, P)
    print(f"Q10: {name:20s} P={CIP_K if P is None else P:.3f}  "
          f"h2_obs={float(obs):.4f}  round trip={float(back):.4f}")
obs_register = float(liability_to_observed_h2(
    H2, CIP_K, n_cases / len(cohort.ids)))
e_raw = np.asarray(score(obs_register).est)
print(f"Q11: unconverted h2={obs_register:.4f}: "
      f"corr={np.corrcoef(e_raw, truth)[0, 1]:.4f} "
      f"slope={np.polyfit(e_raw, truth, 1)[0]:.3f}")
print(f"    balanced h2 without the P factor -> {float(observed_to_liability_h2(liability_to_observed_h2(H2, CIP_K, 0.5), CIP_K)):.2f}"
      "  (>1, rejected by the scorer)")
for name, n_case, n_ctrl in (("register as sampled", n_cases, len(cohort.ids) - n_cases),
                             ("balanced 492/492", 492, 492),
                             ("biobank 5k/45k", 5_000, 45_000),
                             ("consortium 50k/150k", 50_000, 150_000)):
    n = n_case + n_ctrl
    print(f"Q12: {name:20s} N={n:6d}  P={n_case / n:.3f}  "
          f"Neff={4 * n_case * n_ctrl / n:8.0f}  "
          f"Neff/N={4 * (n_case / n) * (n_ctrl / n):.3f}")

# ------------------------- Part E: 0/1 baseline and iPSYCH-style ascertainment
print("\nPart E -- the score against a 0/1 indicator, and a case-cohort sample")
idx = {p: i for i, p in enumerate(cohort.ids)}
sibs = {}
for i, (f_, m_) in enumerate(zip(cohort.father, cohort.mother)):
    if f_ in idx and m_ in idx:
        sibs.setdefault((f_, m_), []).append(i)


def fdr_count(i):
    rels = [idx[p_] for p_ in (cohort.father[i], cohort.mother[i]) if p_ in idx]
    f_, m_ = cohort.father[i], cohort.mother[i]
    if f_ in idx and m_ in idx:
        rels += [j for j in sibs[(f_, m_)] if j != i]
    return int(cohort.status[rels].sum())


n_fdr = np.array([fdr_count(i) for i in range(len(cohort.ids))], float)
fh01 = (n_fdr > 0).astype(float)
print(f"Q13: {int(fh01.sum())} of {len(fh01)} FH-positive; R2 with truth:")
for name, x in (("own status", cohort.status.astype(float)),
                ("FH indicator 0/1", fh01), ("# affected FDRs", n_fdr),
                ("ltpred score", est)):
    print(f"     {name:18s} {np.corrcoef(x, truth)[0, 1] ** 2:.3f}")
for name, x in (("FH indicator 0/1", fh01[at_risk]),
                ("# affected FDRs", n_fdr[at_risk])):
    print(f"     AUC {name:18s} {auc(x, y):.3f}   (score 0.612)")

Q_KEEP = 0.2
keep = big.status | (np.random.default_rng(20260924).random(len(big.ids)) < Q_KEEP)
kept = set(np.flatnonzero(keep).tolist())
bidx = {p: i for i, p in enumerate(big.ids)}
pi = np.where(big.status, 1.0, Q_KEEP)
par, kid, w = [], [], []
for i in kept:
    for p_ in (big.father[i], big.mother[i]):
        if p_ in bidx and bidx[p_] in kept:
            par.append(big.status[bidx[p_]]); kid.append(big.status[i])
            w.append(1.0 / (pi[bidx[p_]] * pi[i]))
par, kid, w = np.asarray(par), np.asarray(kid), np.asarray(w)
t_keep = tetrachoric(par, kid)
print(f"Q14: {int(keep.sum())} of {len(big.ids)} kept, case rate "
      f"{big.status[np.flatnonzero(keep)].mean():.1%} vs {big.status.mean():.1%}; "
      f"naive h2 = {2 * t_keep.rho:+.2f} +/- {2 * t_keep.se:.2f}")
p_w = w / w.sum()
rhos = []
for rep in range(20):
    d = np.random.default_rng(rep).choice(len(par), size=len(par), p=p_w)
    rhos.append(2 * tetrachoric(par[d], kid[d]).rho)
print(f"     20 inverse-probability remixes: mean h2 = {np.mean(rhos):+.2f} "
      f"(spread {np.std(rhos):.2f})")

print(f"\nPython {platform.python_version()}  NumPy {np.__version__}  "
      f"SciPy {scipy.__version__}  ltpred {ltpred.__version__}")
