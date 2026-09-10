"""Re-score the three observed-subgraph detectors on a 16-cell coverage grid.

WHY THIS EXISTS.  The 2026-09-03 campaign established that graph-valued detectors have
false-alarm rates depending on observation coverage, on a four-cell grid c in
{0.25, 0.50, 0.75, 1.00}.  That is a homogeneity test on 3 degrees of freedom, and 3 df
is thin for a claim that a rate varies with a continuous quantity: it can say the four
cells differ, and it cannot say whether the rate falls, rises, or moves without pattern.
This run puts the same detectors on a 0.0625-step grid, 16 cells, 15 degrees of freedom.

TWO DESIGN COMMITMENTS, BOTH DELIBERATE.

1. THE DETECTION RULES ARE NOT RECALIBRATED.  The thresholds frozen on block A controls
   on 2026-09-03 -- 13.00 for observed giant-component fraction, 7.50 for observed mean
   repeat-tie degree, 14.75 for observed susceptibility -- are loaded unchanged from
   configs/.  Recalibrating on the new grid would produce a "pre-declared" rule written
   after the first campaign's outcome was known, which is the exact error this campaign
   already made once with its pooled criterion and recorded rather than smoothed over.
   The new cells are therefore an out-of-design extension evaluated under the original
   frozen rule, and that is the only reading of them that is honest.

2. THE GRID IS A STRICT SUPERSET OF THE OLD ONE.  0.0625-steps contain 0.25, 0.50, 0.75
   and 1.00 exactly.  Seeds are drawn per cell from seed_for("ea|{c}|{arm}", rep), so
   those four cells receive the same seeds and the same rules as on 2026-09-03 and MUST
   reproduce their committed false-alarm counts to the integer.  The script asserts this
   against aggregates/three_observable_block_b_results_20260903.json and fails if any of
   the four disagrees.  A finer grid whose shared cells do not reproduce is measuring
   something other than what it claims to measure.

Block B (odd seeds) only.  Block A is untouched, because nothing here recalibrates.

Usage:  python run_coverage_finegrid.py [--workers N]
"""

from __future__ import annotations

import argparse
import itertools
import json
import os

import numpy as np

from chi2 import chi2_sf, het_chi2
from run_three_observable_campaign import (ARMS, COLOURS, OBSERVABLES, REPS,
                                           SEVERITIES, ALPHA, config,
                                           evaluation_job, rules, summary,
                                           wilson, ztest)
from runner import HERE, pmap, seed_for

STEP = 0.0625
C_GRID = tuple(round(STEP * i, 4) for i in range(1, 17))
OLD_GRID = (0.25, 0.5, 0.75, 1.0)
BASELINE_AGG = "three_observable_block_b_results_20260903.json"
OUT_AGG = "coverage_finegrid_block_b_20260904.json"
DETECTORS = ("latency",) + OBSERVABLES


def jobs_for(scenarios):
    """Block B jobs (odd seeds) over the fine grid.  Same parity filter as the
    2026-09-03 campaign, so the shared cells select the same replicates."""
    return [(c, arm, colour, sev, scen, rep)
            for c, arm, colour, sev, scen in itertools.product(
                C_GRID, ARMS, COLOURS, SEVERITIES, scenarios)
            for rep in range(REPS)
            if seed_for(f"ea|{c}|{arm}", rep) % 2 == 1]


def counts_by_cell(controls):
    out = {}
    for c in C_GRID:
        cell = [r for r in controls if r["c"] == c]
        out[c] = {d: (sum(r[d]["lead_fired"] for r in cell), len(cell))
                  for d in DETECTORS}
    return out


def spearman(xs, ys):
    """Rank correlation of rate against coverage.  Reported because 15 df can show a
    MONOTONE dependence where 3 df could only show inhomogeneity, and monotone is the
    claim the paper would actually like to make."""
    good = [(x, y) for x, y in zip(xs, ys) if np.isfinite(y)]
    if len(good) < 4:
        return {"rho": float("nan"), "n": len(good)}
    x = np.asarray([g[0] for g in good], float)
    y = np.asarray([g[1] for g in good], float)

    def rank(v):
        order = np.argsort(v, kind="stable")
        r = np.empty(v.size, float)
        r[order] = np.arange(1, v.size + 1, dtype=float)
        # average ranks within ties
        for value in np.unique(v):
            m = v == value
            if np.count_nonzero(m) > 1:
                r[m] = r[m].mean()
        return r

    rx, ry = rank(x), rank(y)
    if np.std(rx) < 1e-12 or np.std(ry) < 1e-12:
        return {"rho": float("nan"), "n": len(good)}
    return {"rho": float(np.corrcoef(rx, ry)[0, 1]), "n": len(good)}


def reproduction_check(cells):
    """The four shared cells must match the committed 2026-09-03 counts exactly."""
    path = os.path.join(HERE, "aggregates", BASELINE_AGG)
    with open(path, encoding="utf-8") as f:
        base = json.load(f)["false_alarms_by_coverage"]
    rows, mismatches = [], []
    for c in OLD_GRID:
        for d in DETECTORS:
            k, n = cells[c][d]
            ref = base[str(c)][d]
            bk, bn = int(ref["false_alarms"]), int(ref["n_controls"])
            ok = (k, n) == (bk, bn)
            rows.append({"c": c, "detector": d, "now": [k, n],
                         "committed_20260903": [bk, bn], "match": ok})
            if not ok:
                mismatches.append(f"c={c} {d}: {k}/{n} now vs {bk}/{bn} committed")
    return {"n_comparisons": len(rows), "n_matched": sum(r["match"] for r in rows),
            "rows": rows, "mismatches": mismatches}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=None)
    args = ap.parse_args()

    jobs = jobs_for(("degrade", "control"))
    print(f"{len(C_GRID)} coverage cells, {len(C_GRID) - 1} degrees of freedom, "
          f"{len(jobs)} block B runs")
    rows = pmap(evaluation_job, jobs, args.workers, "fine-grid block B")
    controls = [r for r in rows if r["scen"] == "control"]
    cells = counts_by_cell(controls)

    repro = reproduction_check(cells)
    print(f"reproduction of the four shared cells: "
          f"{repro['n_matched']}/{repro['n_comparisons']} field comparisons matched")
    for line in repro["mismatches"]:
        print("  MISMATCH " + line)

    het, pooled, trend = {}, {}, {}
    for d in DETECTORS:
        counts = [cells[c][d] for c in C_GRID]
        het[d] = het_chi2(counts)
        k = sum(a for a, _ in counts)
        n = sum(b for _, b in counts)
        pooled[d] = {"false_alarms": k, "n_controls": n, "rate": k / n,
                     "wilson95": wilson(k, n)}
        trend[d] = spearman(list(C_GRID), [a / b if b else float("nan")
                                           for a, b in counts])
    lk, ln = pooled["latency"]["false_alarms"], pooled["latency"]["n_controls"]
    per_cell_vs_latency = {}
    for d in OBSERVABLES:
        pooled[d]["versus_latency_ztest"] = ztest(
            pooled[d]["false_alarms"], pooled[d]["n_controls"], lk, ln)
        cellwise = []
        for c in C_GRID:
            k, n = cells[c][d]
            bk, bn = cells[c]["latency"]
            t = ztest(k, n, bk, bn)
            cellwise.append({"c": c, "detector": [k, n], "latency": [bk, bn],
                             "p_two_sided": t["p_two_sided"],
                             "fails_alpha_0_05": t["p_two_sided"] < ALPHA})
        per_cell_vs_latency[d] = {
            "cells": cellwise,
            "n_cells_failing": sum(x["fails_alpha_0_05"] for x in cellwise),
            "n_cells": len(cellwise)}

    payload = {
        "_experiment": ("coverage-heterogeneity of three observed-subgraph detectors on "
                        "a 16-cell grid, frozen 2026-09-03 rules, block B only"),
        "design": {
            "coverage_grid": list(C_GRID), "step": STEP,
            "degrees_of_freedom": len(C_GRID) - 1,
            "partition": "odd seeds only (block B)",
            "rules_recalibrated": False,
            "rule_source": {x: rules()[x]["leading"]["k"] for x in OBSERVABLES},
            "latency_rule": "detection_rule_recal.json",
            "n_block_b_runs": len(jobs), "n_controls": len(controls),
            "scenarios": ["control", "degrade"], "reps": REPS,
            "arms": list(ARMS), "colours": list(COLOURS), "severities": list(SEVERITIES),
            "supersedes_nothing": ("extends the four-cell grid of "
                                   "three_observable_block_b_results_20260903.json; "
                                   "that campaign's numbers stand unchanged"),
        },
        "reproduction_of_shared_cells": repro,
        "false_alarms_by_coverage": {str(c): {d: {"false_alarms": cells[c][d][0],
                                                  "n_controls": cells[c][d][1],
                                                  "rate": (cells[c][d][0] /
                                                           cells[c][d][1])}
                                              for d in DETECTORS}
                                     for c in C_GRID},
        "coverage_heterogeneity": het,
        "coverage_trend_spearman": trend,
        "pooled_false_alarms": pooled,
        "per_cell_versus_latency": per_cell_vs_latency,
        "per_run_rows": rows,
    }
    path = os.path.join(HERE, "aggregates", OUT_AGG)
    with open(path, "x", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=1)
        f.write("\n")

    print()
    for d in DETECTORS:
        h = het[d]
        print(f"{d:22s} pooled {pooled[d]['false_alarms']:4d}/{pooled[d]['n_controls']:5d} "
              f"= {100 * pooled[d]['rate']:5.2f}%   "
              f"chi-square {h['chi2']:8.2f} on {h['df']} df, p = {h['p_value']:.3e}   "
              f"rank corr with coverage {trend[d]['rho']:+.3f}")
    for d in OBSERVABLES:
        x = per_cell_vs_latency[d]
        print(f"{d:22s} fails latency match in {x['n_cells_failing']} of {x['n_cells']} cells")
    print(f"\nwritten: aggregates/{OUT_AGG}")
    if repro["mismatches"]:
        raise SystemExit("shared-cell reproduction FAILED; see mismatches above")


if __name__ == "__main__":
    main()
