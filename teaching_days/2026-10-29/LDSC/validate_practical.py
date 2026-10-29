"""Check the LDSC practical against the raw files and against invariances of the method.

Usage: python validate_practical.py [WORK_DIR]
WORK_DIR holds (or receives) gwas_downloads/; default is the current directory.
Writes validation.json next to this script.
"""
from contextlib import redirect_stdout
from dataclasses import replace
from io import StringIO
from pathlib import Path
import hashlib
import json
import os
import platform
import sys
import time

import allelign
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import pldsc

ROOT = Path(__file__).resolve().parent
student = json.loads((ROOT / 'LDSC_practical.ipynb').read_text())
answers = json.loads((ROOT / 'LDSC_practical_answers.ipynb').read_text())
code_cells = [c for c in student['cells'] if c['cell_type'] == 'code']
answer_cells = [c for c in answers['cells'] if c['cell_type'] == 'code']
sources = json.loads((ROOT / 'data/sources.json').read_text())
checks = []


def check(name, condition):
    if not condition:
        raise AssertionError(name)
    checks.append(name)


def tagged(tag):
    return next(''.join(c['source']) for c in code_cells if tag in c['metadata']['tags'])


check('Identical student and answer code', [c['source'] for c in code_cells]
      == [c['source'] for c in answer_cells])
student_text = [c['source'] for c in student['cells'][1:] if c['cell_type'] == 'markdown']
answer_text = [c['source'] for c in answers['cells'][1:] if c['cell_type'] == 'markdown'
               and not ''.join(c['source']).startswith('**Answer.**')]
check('Identical question text', student_text == answer_text)
check('No exercise placeholders', not any('= ...' in ''.join(c['source']) for c in code_cells))
check('Student notebook has no outputs', not any(c['outputs'] for c in code_cells))
run = [c for c in answer_cells if 'skip-execution' not in c['metadata']['tags']]
check('Answers notebook executed without errors',
      all(c['execution_count'] for c in run)
      and not any(o['output_type'] == 'error' for c in run for o in c['outputs']))

os.chdir(Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd())
# Jupyter's table layout, so the re-run can be compared with the stored outputs.
pd.set_option('display.max_columns', 20)
pd.set_option('display.width', 80)
shown = []
scope = {'display': lambda value: shown.append(repr(value))}
start = time.perf_counter()
reruns = []
for cell in run:
    shown = []
    printed = StringIO()
    with redirect_stdout(printed):
        exec(''.join(cell['source']), scope)
    reruns.append((printed.getvalue(), shown))
elapsed = time.perf_counter() - start
scope['plt'].close('all')
for cell, (printed, displayed) in zip(run, reruns):
    stored_print = ''.join(''.join(o['text']) for o in cell['outputs']
                           if o['output_type'] == 'stream' and o['name'] == 'stdout')
    stored_display = [''.join(o['data']['text/plain']) for o in cell['outputs']
                      if o['output_type'] == 'display_data' and 'image/png' not in o['data']]
    check(f"Stored outputs reproduce: {cell['metadata']['tags'][0]}",
          printed == stored_print and displayed == stored_display)
for trait, study in scope['STUDIES'].items():
    record = sources['gwas'][trait]
    check(f'{trait} source matches sources.json',
          (record['accession'], record['file']) == (study, scope['paths'][trait].name))
    check(f'{trait} download matches its checksum',
          hashlib.sha256(scope['paths'][trait].read_bytes()).hexdigest() == record['sha256'])

ref, gwas, aligned = scope['ref'], scope['gwas'], scope['aligned']
ldsc_h2_sumstats, ldsc_rg_sumstats = scope['ldsc_h2_sumstats'], scope['ldsc_rg_sumstats']
check('Reference matches its provenance',
      ref.provenance['bundle_sha256'] == sources['ld_reference']['bundle_sha256']
      and len(ref.variants) + ref.long_range_ld['n_dropped'] == sources['ld_reference']['n_variants']
      and ref.long_range_ld['n_dropped'] == sources['ld_reference']['long_range_ld_dropped']
      and ref.m_snps_5_50 == sources['ld_reference']['m_snps_5_50'])

# Raw files read independently with pandas.
columns = {'BMI': ('SNP', 'Tested_Allele', 'Other_Allele', 'BETA', 'SE'),
           'TG': ('hm_rsid', 'hm_effect_allele', 'hm_other_allele', 'hm_beta', 'standard_error'),
           'DBP': ('MarkerName', 'Allele1', 'Allele2', 'Effect', 'StdErr')}
separators = {'BMI': '\t', 'TG': '\t', 'DBP': r'\s+'}
reference = pd.DataFrame({'row': np.arange(len(ref.variants)), 'id': ref.variants.id,
                          'chrom': ref.variants.chrom, 'pos': ref.variants.pos,
                          'a1': ref.variants.a1, 'a2': ref.variants.a2})
for trait, (snp, ea, oa, beta, se) in columns.items():
    raw = pd.read_csv(scope['paths'][trait], sep=separators[trait], usecols=[snp, ea, oa, beta, se])
    check(f'{trait}: every row read', len(raw) == len(gwas[trait]) == sources['gwas'][trait]['rows'])
    raw['ea'], raw['oa'] = raw[ea].str.upper(), raw[oa].str.upper()
    raw['z'] = raw[beta] / raw[se]
    if trait == 'DBP':
        marker = raw[snp].str.split(':', expand=True)
        check('DBP positions come from MarkerName',
              np.array_equal(gwas['DBP'].pos, marker[1].astype(int).to_numpy())
              and np.array_equal(gwas['DBP'].chrom, marker[0].to_numpy()))
        raw['chrom'], raw['pos'] = marker[0], marker[1].astype(int)
        matched = raw.merge(reference, on=['chrom', 'pos'])
    else:
        matched = raw.merge(reference, left_on=snp, right_on='id')
    same = (matched.ea == matched.a1) & (matched.oa == matched.a2)
    swapped = (matched.ea == matched.a2) & (matched.oa == matched.a1)
    matched = matched[same | swapped].drop_duplicates('row', keep=False)
    expected = pd.Series(np.where(same[matched.index], 1, -1) * matched.z.to_numpy(),
                         index=matched.row.to_numpy())
    h = aligned[trait]
    observed = pd.Series(h.beta / h.se, index=h.var_index)
    common = observed.index.intersection(expected.index)
    check(f'{trait}: allelign signs agree with an independent merge',
          len(common) >= 0.999 * len(observed)
          and np.allclose(observed[common], expected[common], rtol=1e-6))

# Invariances: M rescales h2 only; allele coding and trait order do not matter.
bmi = gwas['BMI']
M, M_all = ref.m_snps_5_50, ref.m_snps
fit = ldsc_h2_sumstats(bmi, ref.variants, ref.scores, m_snps=M).estimate
fit_all = ldsc_h2_sumstats(bmi, ref.variants, ref.scores, m_snps=M_all).estimate
check('h2 scales with M', np.isclose(fit_all.h2 / fit.h2, M_all / M, rtol=1e-9))
check('Intercept does not depend on M', np.isclose(fit_all.intercept, fit.intercept, rtol=1e-9))
reversed_bmi = replace(bmi, ea=bmi.oa, oa=bmi.ea, beta=-bmi.beta)
fit_reversed = ldsc_h2_sumstats(reversed_bmi, ref.variants, ref.scores, m_snps=M).estimate
check('h2 does not depend on allele coding', np.isclose(fit_reversed.h2, fit.h2, rtol=1e-9))

rg = ldsc_rg_sumstats(bmi, gwas['TG'], ref.variants, ref.scores, m_snps=M, ncores=2).estimate
swapped_order = ldsc_rg_sumstats(gwas['TG'], bmi, ref.variants, ref.scores, m_snps=M, ncores=2).estimate
rg_all = ldsc_rg_sumstats(bmi, gwas['TG'], ref.variants, ref.scores, m_snps=M_all, ncores=2).estimate
rg_reversed = ldsc_rg_sumstats(reversed_bmi, gwas['TG'], ref.variants, ref.scores, m_snps=M,
                               ncores=2).estimate
check('rg does not depend on trait order', np.isclose(swapped_order.rg, rg.rg, rtol=1e-8)
      and np.isclose(swapped_order.rg_se, rg.rg_se, rtol=1e-8))
check('rg does not depend on M', np.isclose(rg_all.rg, rg.rg, rtol=1e-8))
check('rg does not depend on allele coding', np.isclose(rg_reversed.rg, rg.rg, rtol=1e-8))

h2_table = {t: {'h2': f.h2, 'h2_se': f.h2_se, 'intercept': f.intercept,
                'intercept_se': f.intercept_se, 'mean_chi2': f.mean_chisq}
            for t, f in scope['fits'].items()}
rg_table = {k: {'rg': v['rg'], 'rg_se': v['SE'], 'cross_intercept': v['intercept']}
            for k, v in scope['rg_table'].items()}
check('BMI-TG rg in the BMI-TG row of Table 5', np.isclose(rg_table['BMI–TG']['rg'], rg.rg, rtol=1e-12))
check('TG intercept below one', h2_table['TG']['intercept'] < 1)
answer_markdown = '\n'.join(''.join(c['source']) for c in answers['cells']
                             if ''.join(c['source']).startswith('**Answer.**'))
quoted = [f"{h2_table[t]['h2']:.3f}" for t in ('BMI', 'TG', 'DBP')]
quoted += [f"{h2_table[t]['intercept']:.3f}" for t in ('BMI', 'TG', 'DBP')]
quoted += [f"{rg_table[p]['rg']:.3f}".replace('-', '−') for p in rg_table]
quoted += [f"{rg_table[p]['rg_se']:.3f}" for p in rg_table]
quoted += [f"{rg_table['BMI–TG']['cross_intercept']:.3f}"]
missing = [value for value in quoted if value not in answer_markdown]
check('Answer text quotes the computed estimates', not missing)
check('BMI-DBP cross-trait intercept close to zero', abs(rg_table['BMI–DBP']['cross_intercept']) < 0.02)

record = {
    'teaching_date': student['metadata']['atig']['teaching_date'],
    'course_version': student['metadata']['atig']['version'],
    'pldsc_ref': student['metadata']['atig']['pldsc_ref'],
    'versions': {'python': platform.python_version(), 'pldsc': pldsc.__version__,
                 'allelign': allelign.__version__, 'numpy': np.__version__,
                 'pandas': pd.__version__, 'matplotlib': matplotlib.__version__},
    'executed_code_cells': len(run),
    'elapsed_seconds': round(elapsed, 1),
    'scientific_checks': len(checks),
    'h2': {t: {k: round(v, 4) for k, v in r.items()} for t, r in h2_table.items()},
    'rg': {p: {k: round(v, 4) for k, v in r.items()} for p, r in rg_table.items()},
    'colab': 'Live Colab execution not verified.',
    'notebooks': {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                  for name in ('LDSC_practical.ipynb', 'LDSC_practical_answers.ipynb')},
}
(ROOT / 'validation.json').write_text(json.dumps(record, indent=2, ensure_ascii=False) + '\n')
print(f'Passed {len(checks)} checks in {elapsed:.0f} s of notebook code.')
print(json.dumps({'h2': record['h2'], 'rg': record['rg']}, indent=1, ensure_ascii=False))
