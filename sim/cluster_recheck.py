"""Recompute the coverage-heterogeneity tests at the correct unit of analysis.

THE DEFECT.  The 2026-09-03 campaign and the 2026-09-04 fine-grid rerun both compute the
coverage-heterogeneity chi-square over RUNS, treating every control run in a coverage cell as
an independent Bernoulli trial.  They are not independent.  Each (coverage, arm, replicate)
triple is simulated once per (noise colour, severity) pair -- nine of them -- and the three
observed-subgraph detectors are INVARIANT to that pair: measured on the fine grid, all three
give an identical fired/not-fired verdict across all nine noise settings in 494 of 494 groups.
The noise instrument perturbs reported scalar streams; it does not perturb graph topology, so
a graph detector sees the same graph nine times and votes the same way nine times.

A chi-square built on nine copies of one observation is inflated by a design effect of
about nine.  This script recomputes each test with one observation per group.

For the graph detectors the collapse is exact: the group value is 0 or 1 with no information
lost.  For the scalar comparator it is not -- it varies within 70 of 494 groups -- so its
group value is the fraction of its nine runs that fired, and its corrected statistic is
reported as an approximation with that caveat attached rather than as an exact test.

Run:  python cluster_recheck.py
"""

from __future__ import annotations

import collections
import json
import os

from chi2 import chi2_sf
from runner import HERE

FINE = "coverage_finegrid_block_b_20260904.json"
COARSE = "three_observable_block_b_results_20260903.json"
DETECTORS = ("latency", "obs_giant_frac", "obs_mean_degree", "obs_susceptibility")


def load(name):
    with open(os.path.join(HERE, "aggregates", name), encoding="utf-8") as f:
        return json.load(f)


def het(cells):
    """Heterogeneity statistic over (k, n) pairs; k may be fractional after collapsing."""
    tk = sum(k for k, _ in cells)
    tn = sum(n for _, n in cells)
    p = tk / tn if tn else float("nan")
    df = len(cells) - 1
    if not (0.0 < p < 1.0):
        return {"chi2": float("nan"), "df": df, "p": float("nan"), "pooled": p}
    chi = sum((k - n * p) ** 2 / (n * p * (1 - p)) for k, n in cells)
    return {"chi2": chi, "df": df, "p": chi2_sf(chi, df), "pooled": p}


def report(label, agg_name):
    d = load(agg_name)
    rows = [r for r in d["per_run_rows"] if r["scen"] == "control"]
    groups = collections.defaultdict(list)
    for r in rows:
        groups[(r["c"], r["arm"], r["rep"])].append(r)
    sizes = sorted({len(v) for v in groups.values()})
    coverages = sorted({r["c"] for r in rows})

    print(f"\n=== {label} ===")
    print(f"{len(rows)} control runs, {len(groups)} (coverage, arm, replicate) groups, "
          f"group size(s) {sizes}, {len(coverages)} coverage cells")

    print(f"\n{'detector':22s} {'invariant':>12s}  "
          f"{'run-level chi2':>16s} {'p':>11s}   {'group-level chi2':>17s} {'p':>11s}  {'ratio':>6s}")
    out = {}
    for det in DETECTORS:
        invariant = sum(1 for v in groups.values()
                        if len({x[det]["lead_fired"] for x in v}) == 1)

        run_cells, grp_cells = [], []
        for c in coverages:
            rr = [r for r in rows if r["c"] == c]
            run_cells.append((sum(r[det]["lead_fired"] for r in rr), len(rr)))
            gg = [v for k, v in groups.items() if k[0] == c]
            grp_cells.append((sum(sum(x[det]["lead_fired"] for x in v) / len(v)
                                  for v in gg), len(gg)))
        a, b = het(run_cells), het(grp_cells)
        ratio = a["chi2"] / b["chi2"] if b["chi2"] else float("nan")
        exact = invariant == len(groups)
        print(f"{det:22s} {invariant:5d}/{len(groups):<6d}  "
              f"{a['chi2']:16.2f} {a['p']:11.2e}   {b['chi2']:17.2f} {b['p']:11.2e}  "
              f"{ratio:6.2f}{'' if exact else '  (approx: varies within groups)'}")
        out[det] = {"n_groups": len(groups), "groups_invariant": invariant,
                    "collapse_exact": exact,
                    "run_level": a, "group_level": b, "design_effect_ratio": ratio}
    return out


if __name__ == "__main__":
    result = {"_analysis": ("coverage-heterogeneity tests recomputed with one observation per "
                            "(coverage, arm, replicate) group instead of one per run"),
              "coarse_4_cell_20260903": report("4-cell grid, 2026-09-03 campaign", COARSE),
              "fine_16_cell_20260904": report("16-cell grid, 2026-09-04 rerun", FINE)}
    path = os.path.join(HERE, "aggregates", "cluster_recheck_20260904.json")
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        json.dump(result, f, indent=1)
        f.write("\n")
    print(f"\nwritten: aggregates/cluster_recheck_20260904.json")
