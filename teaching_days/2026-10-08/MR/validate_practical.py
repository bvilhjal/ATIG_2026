"""Check the prepared data and the notebook's direct calculations against references."""
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
from scipy.linalg import lstsq

ROOT = Path(__file__).resolve().parent
student = json.loads((ROOT / 'MR_practical.ipynb').read_text())
answers = json.loads((ROOT / 'MR_practical_answers.ipynb').read_text())
code_cells = [c for c in student['cells'] if c['cell_type'] == 'code']
answer_code = [c['source'] for c in answers['cells'] if c['cell_type'] == 'code']
checks = []


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    checks.append(name)


check('Identical student and answer code', [c['source'] for c in code_cells] == answer_code)
check('No exercise placeholders', not any('= ...' in ''.join(c['source']) for c in code_cells))
scope = {'display': lambda *_: None}
with redirect_stdout(StringIO()):
    for cell in code_cells:
        exec(''.join(cell['source']), scope)
scope['plt'].close('all')
check('Embedded prepared data match CSV', scope['CSV_DATA'] == (ROOT / 'data/bmi_chd_prepared.csv').read_text())
provenance = json.loads((ROOT / 'data/sources.json').read_text())
for filename, record in provenance['files'].items():
    check('Checksum ' + filename,
          hashlib.sha256((ROOT / 'data' / filename).read_bytes()).hexdigest() == record['sha256'])

bmi = pd.read_csv(ROOT / 'data/bmi_gwas.csv')
chd = pd.read_csv(ROOT / 'data/chd_gwas.csv')
joined = bmi.merge(chd, on='snp', suffixes=('_bmi', '_chd'), validate='one_to_one')
check('All source rows match', len(joined) == len(bmi) == len(chd) == 79)
check('Coordinates agree', (joined.chr_bmi == joined.chr_chd).all()
      and (joined.pos_bmi == joined.pos_chd).all())
check('Original alleles agree', (joined.effect_allele_bmi == joined.effect_allele_chd).all()
      and (joined.other_allele_bmi == joined.other_allele_chd).all())
check('Exposure variants are genome-wide significant', (joined.p_bmi < 5e-8).all())
palindrome = (joined.effect_allele_bmi + joined.other_allele_bmi).isin(['AT', 'TA', 'CG', 'GC'])
raw = joined.loc[~palindrome].reset_index(drop=True)
data = scope['data']
check('Four palindromes removed, 75 retained', palindrome.sum() == 4 and len(data) == 75)
check('Variant order matches', data.snp.equals(raw.snp))
check('Original allele labels retained', data.effect_allele.equals(raw.effect_allele_bmi)
      and data.other_allele.equals(raw.other_allele_bmi))
check('Standard errors unchanged', np.allclose(data.bmi_se, raw.se_bmi)
      and np.allclose(data.chd_se, raw.se_chd))
check('Original signed effects retained', np.allclose(data.bmi_effect, raw.beta_bmi)
      and np.allclose(data.chd_effect, raw.beta_chd))
numeric_columns = ['bmi_effect', 'bmi_se', 'chd_effect', 'chd_se']
check('Finite effects and positive standard errors', np.isfinite(data[numeric_columns]).all().all()
      and (data[['bmi_se', 'chd_se']] > 0).all().all())
check('Ratio standard errors stay positive with signed BMI effects',
      (scope['snp_standard_errors'] > 0).all()
      and (data.bmi_effect < 0).any() and (data.bmi_effect > 0).any())

# Independent weighted-regression reference, using original allele orientation.
X = (raw.beta_bmi / raw.se_chd).to_numpy()[:, None]
y = (raw.beta_chd / raw.se_chd).to_numpy()
reference_beta = lstsq(X, y)[0][0]
reference_q = np.sum((y - X[:, 0] * reference_beta)**2)
reference_se = np.sqrt(max(1, reference_q / (len(data)-1)) / np.sum(X**2))
check('IVW equals independent least squares', np.isclose(scope['mr_estimate'], reference_beta, rtol=1e-12))
check('Ratio and regression Q agree', np.isclose(scope['Q'], reference_q, rtol=1e-12))
check('Adjusted SE agrees', np.isclose(scope['mr_standard_error'], reference_se, rtol=1e-12))
check('Interval agrees', np.isclose(scope['lower_limit'], reference_beta - 1.96*reference_se)
      and np.isclose(scope['upper_limit'], reference_beta + 1.96*reference_se))

# Changing the effect allele for both traits must leave MR unchanged.
flipped = dict(scope)
flipped['bmi_effect'] = -scope['bmi_effect']
flipped['chd_effect'] = -scope['chd_effect']
ivw_code = next(''.join(c['source']) for c in code_cells if 'ivw' in c['metadata']['tags'])
with redirect_stdout(StringIO()):
    exec(ivw_code, flipped)
check('MR is invariant to consistent allele reversal',
      np.allclose(flipped['snp_estimates'], scope['snp_estimates'])
      and np.allclose(flipped['snp_standard_errors'], scope['snp_standard_errors'])
      and np.isclose(flipped['mr_estimate'], scope['mr_estimate']))

# Exercise every requested simulation setting using the notebook's actual code.
simulation_code = next(''.join(c['source']) for c in code_cells if 'simulation' in c['metadata']['tags'])
result_code = next(''.join(c['source']) for c in code_cells if 'simulation-result' in c['metadata']['tags'])
results = {}
for label, causal, pleiotropic in [('A', .3, 0.), ('B', .3, .2), ('C', 0., .3)]:
    changed = simulation_code.replace('causal_effect = 0.3', f'causal_effect = {causal}')
    changed = changed.replace('pleiotropy = 0.0', f'pleiotropy = {pleiotropic}')
    with redirect_stdout(StringIO()):
        exec(changed, scope)
        exec(result_code, scope)
    results[label] = (scope['outcome_effect'].copy(), scope['simulated_estimate'])
scope['plt'].close('all')
check('A and C give identical data', np.array_equal(results['A'][0], results['C'][0]))
check('A and C give identical estimates', results['A'][1] == results['C'][1])
check('Added proportional pleiotropy shifts the estimate', np.isclose(results['B'][1]-results['A'][1], .2))
print(f'Passed {len(checks)} checks.')
print(f'Real-data OR {np.exp(reference_beta):.6f}; adjusted 95% interval '
      f'{np.exp(reference_beta-1.96*reference_se):.6f}–{np.exp(reference_beta+1.96*reference_se):.6f}.')
print('Simulation estimates:', {name: round(value[1], 6) for name, value in results.items()})
