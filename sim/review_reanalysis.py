"""Reanalyse retained observations for the revised manuscript; no simulation is run.

All resampling uses independent substrate seeds as units. Coverage tests permute
whole groups within missingness arm, testing exchangeability across coverage under
the null. Detector contrasts use paired group differences and a two-sided normal
approximation. The post-hoc detector family includes every coverage test, pooled
contrast and distinct per-coverage contrast generated here; overlapping coarse-grid
cells are counted once. Holm controls familywise error without a dependence
assumption. Benjamini-Hochberg is a descriptive sensitivity under dependent tests.

Run from any directory: python sim/review_reanalysis.py
"""
from pathlib import Path
import collections
import hashlib
import json
import math
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
AGG = HERE / 'aggregates'
SEED = 20260911
DRAWS = 19999
DETECTORS = ('latency', 'obs_giant_frac', 'obs_mean_degree', 'obs_susceptibility')
SOURCES = {}


def load(name):
    data = (AGG / name).read_bytes()
    SOURCES[name] = hashlib.sha256(data).hexdigest()
    return json.loads(data)


def summary(values, seed=SEED, draws=4000):
    x = np.asarray(values, float)
    assert len(x) and np.isfinite(x).all()
    rng = np.random.default_rng(seed)
    means = np.concatenate([x[rng.integers(len(x), size=(min(250, draws-i), len(x)))].mean(1)
                            for i in range(0, draws, 250)])
    return {'n': len(x), 'mean': float(x.mean()), 'se': float(means.std(ddof=1)),
            'ci95': np.quantile(means, [.025, .975]).tolist()}


def paired_test(x):
    x = np.asarray(x, float)
    se = float(x.std(ddof=1) / np.sqrt(len(x)))
    p = math.erfc(abs(float(x.mean())) / se / math.sqrt(2)) if se else (1.0 if x.mean() == 0 else 0.0)
    return {'n': len(x), 'difference': float(x.mean()), 'se': se, 'p': p,
            'ci95': [float(x.mean()-1.95996398454*se), float(x.mean()+1.95996398454*se)]}


def adjust(tests):
    n = len(tests)
    order = sorted(range(n), key=lambda i: tests[i]['p'])
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (n-rank)*tests[i]['p'])
        tests[i]['holm_p'] = min(1.0, running)
    running = 1.0
    for rank in range(n-1, -1, -1):
        i = order[rank]
        running = min(running, n*tests[i]['p']/(rank+1))
        tests[i]['bh_p'] = running
    return {'n': n, 'holm_survivors': sum(t['holm_p'] <= .05 for t in tests),
            'bh_survivors': sum(t['bh_p'] <= .05 for t in tests), 'tests': tests}


def grouped(rows, nested=True):
    g = collections.defaultdict(list)
    for r in rows:
        if r['scen'] == 'control':
            g[(r['c'], r['arm'], r['rep'])].append(r)
    keys = sorted(g)
    if nested:
        vals = np.array([[np.mean([r[d]['lead_fired'] for r in g[k]]) for d in DETECTORS] for k in keys])
    else:
        vals = np.array([[np.mean([r[d] for r in g[k]]) for d in ('lead_fired', 'lag_fired')] for k in keys])
    return keys, vals, g


def heterogeneity(keys, vals):
    """Permutation of independent groups, stratified by missingness arm.

    The statistic is the coverage-cell sum of squared mean deviations, weighted
    by cell group counts. Its scale is immaterial to permutation inference.
    """
    cs = np.array([k[0] for k in keys])
    cells = [np.flatnonzero(cs == c) for c in sorted(set(cs))]
    arms = [np.array([i for i,k in enumerate(keys) if k[1] == a]) for a in sorted({k[1] for k in keys})]
    def stat(v):
        return sum(len(ix)*(v[ix].mean(0)-v.mean(0))**2 for ix in cells)
    obs = stat(vals)
    rng = np.random.default_rng(SEED)
    exceed = np.zeros(vals.shape[1], int)
    for _ in range(DRAWS):
        shuffled = vals.copy()
        for ix in arms:
            shuffled[ix] = vals[rng.permutation(ix)]
        exceed += stat(shuffled) >= obs-1e-14
    return [{'statistic': float(s), 'p': float((k+1)/(DRAWS+1)),
             'exceedances': int(k), 'draws': DRAWS} for s,k in zip(obs, exceed)]


def graph_analysis():
    tests, grids = [], {}
    previous = {}
    for grid, name in [('coarse', 'three_observable_block_b_results_20260903.json'),
                       ('fine', 'coverage_finegrid_block_b_20260904.json')]:
        rows = load(name)['per_run_rows']
        keys, vals, groups = grouped(rows)
        degradation = collections.defaultdict(list)
        for row in rows:
            if row['scen'] == 'degrade':
                degradation[(row['c'], row['arm'], row['rep'])].append(row)
        hs = heterogeneity(keys, vals)
        cells = sorted({k[0] for k in keys})
        out = {'groups': len(keys), 'coverage_cells': cells, 'detectors': {}}
        for j,d in enumerate(DETECTORS):
            tests.append({'comparison': f'{grid}: coverage heterogeneity, {d}', **hs[j]})
            item = {'rate': summary(vals[:,j]), 'heterogeneity': hs[j],
                    'invariant_groups': sum(len({r[d]['lead_fired'] for r in groups[k]}) == 1 for k in keys)}
            if j:
                def alarm(row, detector):
                    return min(320, row[detector]['t_lead']) if row[detector]['lead_fired'] else 320
                gains = [np.mean([alarm(row, 'latency')-alarm(row, d) for row in degradation[k]])
                         for k in sorted(degradation)]
                item['restricted_alarm_gain_vs_latency'] = summary(gains)
                item['paired_comparator_difference'] = paired_test(vals[:,j]-vals[:,0])
                tests.append({'comparison': f'{grid}: pooled {d} minus comparator', **item['paired_comparator_difference']})
                item['per_coverage'] = {}
                for c in cells:
                    ix = [i for i,k in enumerate(keys) if k[0] == c]
                    diff = vals[ix,j]-vals[ix,0]
                    pair = paired_test(diff)
                    item['per_coverage'][str(c)] = pair
                    label = (d,c)
                    if label in previous:
                        assert np.array_equal(diff, previous[label]), 'Overlapping coverage cell changed'
                    else:
                        tests.append({'comparison': f'coverage {c}: {d} minus comparator', **pair})
                        previous[label] = diff
            out['detectors'][d] = item
        grids[grid] = out
    keys, vals, groups = grouped(load('ea_headline_recal.json')['per_run_rows'], False)
    hs = heterogeneity(keys, vals)
    for j,d in enumerate(('formation latency', 'lagging output')):
        tests.append({'comparison': f'primary campaign: coverage heterogeneity, {d}', **hs[j]})
    primary = {'groups': len(keys), 'weighting': 'equal weight per independent group; equal weight per noise copy within group',
               'group_sizes': dict(collections.Counter(len(v) for v in groups.values())),
               'leading': summary(vals[:,0]), 'lagging': summary(vals[:,1]), 'heterogeneity': hs}
    return grids, primary, adjust(tests)


def sensitivity():
    d = load('remediation_runs.json')
    rows = [r for r in d['per_run_rows'] if r['scen'] == 'degrade' and r['c'] == 1]
    base_name = d['design']['base_variant']
    base = {r['seed']: r for r in rows if r['variant'] == base_name}
    out = {}
    for variant in sorted({r['variant'] for r in rows if r['variant'].startswith('knee')}):
        select = {r['seed']:r for r in rows if r['variant'] == variant}
        assert set(select) == set(base)
        diff = [select[s]['dt_censored']-base[s]['dt_censored'] for s in sorted(base)]
        out[variant] = {'alarm_gain': summary([select[s]['dt_censored'] for s in sorted(base)]),
                        'paired_difference_from_base': summary(diff)}
    return {'base': base_name, 'variants': out, 'inference': 'Unadjusted paired percentile intervals; exploratory sensitivity, not equivalence tests.'}


def recovery():
    rows = load('tranche2_runs.json')['a23_small_shock']['per_run_rows']
    cells = [[r for r in rows if r['frac'] == .35 and r['a'] == a] for a in (.48,.39)]
    rng = np.random.default_rng(20260904)
    rates = [np.array([r['recovery_rate'] for r in rs]) for rs in cells]
    assert all(len(x) == 12 for x in rates)
    samples = [x[rng.integers(0,12,size=(4000,12))].mean(1) for x in rates]
    ratio = samples[0]/samples[1]
    assert np.isfinite(ratio).all() and (ratio > 0).all()
    return {'times': [float(1/x.mean()) for x in rates], 'ratio': float(rates[0].mean()/rates[1].mean()),
            'ratio_ci95': np.quantile(ratio,[.025,.975]).tolist(), 'draws':4000, 'seed':20260904,
            'method':'Independent within-cell resampling of 12 rates, ratio of means; no fitted rates or bootstrap draws discarded.',
            'negative_fits': [sum(r['recovery_rate'] < 0 for r in rs) for rs in cells],
            'poor_fits': [sum(r['fit_r2'] < .12 for r in rs) for rs in cells]}


def auxiliary():
    d = load('aux_arms_recal.json')
    out = {}
    for condition in ('static_full','revocation'):
        rows = [r for r in d['ensign_rows'] if r['condition'] == condition]
        raw = np.array([r['abs_bias_raw'] for r in rows])
        corrected = np.array([r['abs_bias_ipw'] for r in rows])
        rng = np.random.default_rng(SEED)
        ix = rng.integers(len(rows),size=(4000,len(rows)))
        reduced = 1-corrected[ix].sum(1)/raw[ix].sum(1)
        out[condition] = {'n':len(rows), 'pooled_bias_reduction':float(1-corrected.sum()/raw.sum()),
                          'ci95':np.quantile(reduced,[.025,.975]).tolist(),
                          'mean_estimated_consent':float(np.mean([r['p_hat'] for r in rows]))}
    return out


def main():
    grids, primary, family = graph_analysis()
    out = {'method':__doc__, 'seed':SEED, 'draws':DRAWS, 'graph':grids,
           'primary_false_alarm':primary, 'detector_multiplicity':family,
           'comparator_sensitivity':sensitivity(), 'recovery':recovery(), 'bias_correction':auxiliary(),
           'coverage_bootstrap':load('remediation_analysis.json')['floor_bootstrap'],
           'exclusions':[
               'Withdrawn run-level tests are replaced by group-level tests, not retained alongside them.',
               'The coarse graph grid is a subset of the fine grid; identical per-cell contrasts appear once.',
               'The bootstrap fraction of positive coverage-floor gaps is not a null p-value and is excluded.',
               'Coupling, recovery, propensity, comparator sensitivity and threshold intervals are exploratory, unadjusted summaries outside this detector family.',
               'The historical statement of 60 gaming tests has no enumerated test registry; no claim of comprehensive correction across the study is made.'
           ], 'input_sha256':SOURCES}
    target = AGG/'review_reanalysis.json'
    encoded = json.dumps(out,indent=2,allow_nan=False)+'\n'
    if '--check' in sys.argv:
        assert target.read_text(encoding='utf-8') == encoded, 'Retained reanalysis differs from recomputation'
    else:
        target.write_bytes(encoded.encode('utf-8'))
    print(json.dumps({'multiplicity':{k:v for k,v in family.items() if k!='tests'},
                      'survivors':[t for t in family['tests'] if t['holm_p'] <= .05],
                      'recovery':out['recovery'], 'bias':out['bias_correction']},indent=2))


if __name__ == '__main__':
    main()
