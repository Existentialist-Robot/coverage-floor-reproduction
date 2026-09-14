#!/usr/bin/env python3
"""emit_missingness_arm_interval.py -- persist the MNAR-minus-MCAR achieved-edge interval that
main.tex prints, as an aggregate cell, so it can be bound like every other reported number.

WHY THIS EXISTS. An audit on 2026-09-12 enumerated every numeral the manuscript prints and tried to
bind each to a cell in sim/aggregates/. Three of 302 matched nothing at all, and one of the three was
the lower endpoint of a published 95% interval:

    "The mean MNAR-minus-MCAR achieved-edge difference at c=0.5 is 0.0012 [-0.0060, +0.0084]
     (Wald 95% interval, 100 independent replicates per arm; noise colours are copies)."

A direct search of all 106 committed aggregates found ZERO cells within half a printed unit of
-0.0060, exactly ONE within half a unit of +0.0084 (a standard deviation in an unrelated bifurcation
per-seed row, a coincidence rather than a source), and 25 unrelated correlation values near 0.0012.

⚠ THE NUMBERS ARE CORRECT AND THEY ALWAYS WERE. That is the finding, and the first framing of it
overstated the problem: the interval REPRODUCES EXACTLY from committed data. What did not exist was a
persisted DERIVED CELL, so an auditor binding printed numbers to aggregates could not reach it. This
script writes that cell. It changes no reported value.

METHOD, stated so the reproduction is checkable rather than asserted:
  * Source is ea_headline_recal.json's per_run_rows. That is the RECALIBRATED campaign, which
    RESULTS.md's own correction names as superseding campaign 1 for every number the manuscript
    reports. Running this against ea_headline.json instead gives +0.0031 [-0.0110, +0.0172] over 30
    replicates -- a different campaign, not a discrepancy.
  * Rows are filtered to c == 0.5 exactly.
  * "Noise colours are copies" is VERIFIED rather than assumed, and it HOLDS: all nine colour-by-
    severity rows sharing a (seed, rep, scen) carry an identical achieved_edge_frac. This script
    raises if any group disagrees, because collapsing a group that genuinely varies would silently
    discard replicates and shrink the interval.
  * ⚠ THE VERIFICATION CAUGHT SOMETHING THE FIRST DRAFT OF THIS SCRIPT HAD WRONG. Keying on
    (seed, rep) alone raises: 100 of 100 such groups carry TWO distinct values, because `scen` takes
    control and degrade and the achieved edge fraction differs between them. The independent unit is
    (seed, rep, scen), which gives 200 per arm rather than 100.
  * ⚠ SO THE MANUSCRIPT'S "100 independent replicates per arm" IS THE DEGRADE SCENARIO ONLY, and the
    sentence does not say so. All three readings are emitted below and all three straddle zero, so
    the paper's inferential claim does not depend on the choice:
        degrade only   n=100/arm   +0.0012 [-0.0060, +0.0084]   <- what the manuscript prints
        control only   n=100/arm   +0.0025 [-0.0047, +0.0096]
        pooled         n=200/arm   +0.0018 [-0.0033, +0.0069]
    The printed interval is CORRECT and reproducible at its stated denominator. What is missing from
    the sentence is one clause naming the scenario.
  * The interval is an UNPAIRED two-sample Wald 95% interval, z = 1.96, with the sample variance
    using the n-1 denominator.
  * ⚠ NO PAIRED FORM EXISTS, and this was also found by the code rather than assumed: the two arms'
    (seed, rep, scen) sets intersect in ZERO identifiers, so there is nothing to pair. "We report the
    unpaired interval" and "no paired interval is computable from these rows" are different
    statements, and only the second one is true here.

Usage:   python sim/emit_missingness_arm_interval.py            # writes the aggregate
         python sim/emit_missingness_arm_interval.py --check    # recomputes, byte-diffs, writes
                                                                # nothing, exits 1 on any difference
"""
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE / 'aggregates' / 'ea_headline_recal.json'
OUT = HERE / 'aggregates' / 'missingness_arm_interval.json'

COVERAGE = 0.5
Z_95 = 1.96
ARMS = ('mnar', 'mcar')
FIELD = 'achieved_edge_frac'


def _mean(xs):
    return sum(xs) / len(xs)


def _var(xs):
    m = _mean(xs)
    return sum((x - m) ** 2 for x in xs) / (len(xs) - 1)


UNIT = ('seed', 'rep', 'scen')


def collapse(rows):
    """One value per independent unit, raising if a unit disagrees -- see METHOD above."""
    groups = {}
    for r in rows:
        key = tuple(r[k] for k in UNIT)
        groups.setdefault(key, []).append(r[FIELD])
    for key, vals in groups.items():
        if max(vals) != min(vals):
            raise ValueError(
                f'{UNIT} {key} carries {len(set(vals))} distinct {FIELD} values '
                f'({min(vals)} .. {max(vals)}); "noise colours are copies" does not hold here and '
                f'collapsing would discard real replicates')
    return {k: v[0] for k, v in groups.items()}


def interval(a_vals, b_vals):
    diff = _mean(a_vals) - _mean(b_vals)
    se = math.sqrt(_var(a_vals) / len(a_vals) + _var(b_vals) / len(b_vals))
    return {'n_per_arm': len(a_vals), 'difference_mnar_minus_mcar': diff,
            'standard_error': se, 'ci95_low': diff - Z_95 * se, 'ci95_high': diff + Z_95 * se}


def render():
    src = json.loads(SOURCE.read_text(encoding='utf-8'))
    rows = [r for r in src['per_run_rows'] if r['c'] == COVERAGE]
    per_arm = {a: collapse([r for r in rows if r['arm'] == a]) for a in ARMS}
    scenarios = sorted({k[2] for k in per_arm['mnar']})

    def vals(arm, scen=None):
        return [v for k, v in sorted(per_arm[arm].items()) if scen is None or k[2] == scen]

    readings = {s: interval(vals('mnar', s), vals('mcar', s)) for s in scenarios}
    readings['pooled'] = interval(vals('mnar'), vals('mcar'))

    # ⚠ THE ARMS SHARE NO IDENTIFIERS, SO NO PAIRED FORM EXISTS. The first draft of this script
    # emitted one and divided by zero, which is how this was found: the (seed, rep, scen) sets of the
    # two arms intersect in 0 ids. That makes the UNPAIRED interval the only available reading rather
    # than a choice between two, and it is recorded here because "we report the unpaired interval"
    # and "no paired interval is computable from these rows" are different statements and only the
    # second one is true.
    shared_ids = len(set(per_arm['mnar']) & set(per_arm['mcar']))

    out = {
        '_source': {
            'result_file': 'sim/aggregates/ea_headline_recal.json',
            'key': 'per_run_rows',
            'filter': f'c == {COVERAGE}',
            'independent_unit': list(UNIT),
            '_note': 'RECOMPUTED from committed per-run rows, not copied from any stored cell. This '
                     'aggregate exists because the manuscript printed this interval while no '
                     'aggregate cell held it, so an audit binding printed numbers to aggregates '
                     'could not reach it. No reported value changes.',
        },
        'coverage': COVERAGE,
        'field': FIELD,
        'arm_mean_by_scenario': {s: {a: _mean(vals(a, s)) for a in ARMS} for s in scenarios},
        'manuscript_prints': 'degrade',
        '_manuscript_note': 'The manuscript states "100 independent replicates per arm" without '
                            'naming the scenario. That denominator is the degrade scenario alone; '
                            'the independent unit across the file is (seed, rep, scen) at 200 per '
                            'arm. All readings below straddle zero, so the inferential claim does '
                            'not depend on the choice.',
        'unpaired_by_scenario': readings,
        'paired_form_unavailable': {
            '_what': 'No paired interval is computable from these rows. The two arms share '
                     f'{shared_ids} identifiers, so every unit in one arm has no counterpart in the '
                     'other and the unpaired interval is the only available reading.',
            'shared_identifiers': shared_ids,
        },
    }
    return json.dumps(out, indent=2, ensure_ascii=True) + '\n'


def main():
    text = render()
    if '--check' in sys.argv[1:]:
        if not OUT.exists():
            print(f'emit_missingness_arm_interval --check: {OUT.name} does not exist')
            sys.exit(1)
        if OUT.read_bytes() != text.encode('utf-8'):
            print(f'emit_missingness_arm_interval --check: FAILED -- {OUT.name} does not match '
                  f'what {SOURCE.name} renders')
            sys.exit(1)
        print(f'emit_missingness_arm_interval --check: {OUT.name} matches {SOURCE.name}')
        sys.exit(0)
    OUT.write_bytes(text.encode('utf-8'))
    doc = json.loads(text)
    print(f'emit_missingness_arm_interval: wrote {OUT.name}')
    for name, r in doc['unpaired_by_scenario'].items():
        flag = '  <- the manuscript prints this one' if name == doc['manuscript_prints'] else ''
        print(f"  {name:8} n={r['n_per_arm']:3}/arm  {r['difference_mnar_minus_mcar']:+.4f} "
              f"[{r['ci95_low']:+.4f}, {r['ci95_high']:+.4f}]{flag}")


if __name__ == '__main__':
    main()
