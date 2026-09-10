"""Recompute the graph-versus-latency LEAD estimates at the correct unit of analysis.

WHY THIS EXISTS.  Section 8.10 of RESULTS.md suppressed every graph-versus-latency lead
estimate, coverage-dependent lead curve and win/loss claim from the observed-graph detector
campaign.  The stated ground was that the detectors' false-alarm rates were unmatched: at the
RUN unit the giant-component interval lay wholly above the comparator's and the mean-degree
interval wholly below it.  The 2026-09-04 correction showed that statement is computed on
2,160 control rows that are really 240 groups of 9, design effect 9.04 for both graph
detectors and 3.79 for the comparator.  At the group unit all three intervals overlap, so the
rule's premise fails and the rule does not fire.

A RULE THAT DOES NOT FIRE REMOVES THE REASON GIVEN FOR SILENCE.  IT SUPPLIES NO REASON TO
SPEAK.  Section 8.12 says what is owed before any lead number may be quoted: "the censoring
rule, the bootstrap and the paired intervals all need the same clustering treatment before any
lead number is quoted, and none of that has been done."  This script is that treatment.  It
reports numbers; it decides nothing about the manuscript.

THE COLLAPSE.  Each (coverage, arm, replicate) triple is simulated once per (noise colour,
severity) pair, nine of them.  The noise instrument perturbs reported scalar streams, not graph
topology, so a graph detector sees the same graph nine times.  This script MEASURES that
invariance rather than assuming it, and refuses to collapse a detector that is not invariant
without saying so.  The comparator is not invariant, so its group value is the mean of its nine
runs and every interval on it is a cluster bootstrap over groups, never over runs.

THE CENSORING.  Unchanged from the campaign: a non-firing detector is assigned cycle 320.  The
censored lead is therefore horizon-dependent and is reported as such, alongside a paired
estimate restricted to groups in which both detectors fired.

THE CONTROL TWIN.  Section 9.2's 2026-09-04 correction found the paired lead contaminated by
runs whose matched control also alarmed.  Each degradation group has a control twin at the same
(coverage, arm, replicate).  This script reports the paired estimate both with and without the
contaminated groups, because the contaminated one is the quantity the earlier record printed.

Run:  python cluster_lead.py
"""

from __future__ import annotations

import collections
import json
import os
import random

from runner import HERE

HEADLINE = "graph_detector_headline.json"
OUT = "cluster_lead_20260906.json"

COMPARATOR = "latency"
GRAPH = ("obs_giant_frac", "obs_mean_degree")
DETECTORS = (COMPARATOR,) + GRAPH

CENSOR_AT = 320
BOOTSTRAP_DRAWS = 2000
BOOTSTRAP_SEED = 20260906


def load(name):
    with open(os.path.join(HERE, "aggregates", name), encoding="utf-8") as f:
        return json.load(f)


def alarm_cycle(row, det):
    """Censored alarm cycle for one detector on one run: t_lead, or the horizon if it never fired."""
    block = row[det]
    t = block.get("t_lead")
    if block.get("lead_fired") and t is not None:
        return float(t), True
    return float(CENSOR_AT), False


def group_key(row):
    return (row["c"], row["arm"], row["rep"])


def invariance(groups, det):
    """How many groups return an identical censored alarm cycle across all nine noise settings."""
    same = 0
    for rows in groups.values():
        vals = {alarm_cycle(r, det)[0] for r in rows}
        if len(vals) == 1:
            same += 1
    return same


def collapse(groups, det):
    """One value per group: the mean censored alarm cycle, and whether the detector fired in all nine."""
    out = {}
    for key, rows in groups.items():
        cycles = [alarm_cycle(r, det) for r in rows]
        mean_cycle = sum(c for c, _ in cycles) / len(cycles)
        fired_all = all(f for _, f in cycles)
        fired_any = any(f for _, f in cycles)
        out[key] = {"cycle": mean_cycle, "fired_all": fired_all, "fired_any": fired_any}
    return out


def boot(values, draws=BOOTSTRAP_DRAWS, seed=BOOTSTRAP_SEED):
    """Cluster bootstrap: resample GROUPS with replacement, never runs."""
    if not values:
        return {"mean": float("nan"), "se": float("nan"), "ci95": [float("nan"), float("nan")], "n": 0}
    rng = random.Random(seed)
    n = len(values)
    mean = sum(values) / n
    means = []
    for _ in range(draws):
        draw = [values[rng.randrange(n)] for _ in range(n)]
        means.append(sum(draw) / n)
    means.sort()
    lo = means[int(0.025 * draws)]
    hi = means[int(0.975 * draws) - 1]
    centred = sum((m - mean) ** 2 for m in means) / (draws - 1)
    return {"mean": mean, "se": centred ** 0.5, "ci95": [lo, hi], "n": n}


def main():
    d = load(HEADLINE)
    rows = d["per_run_rows"]

    degrade = [r for r in rows if r["scen"] == "degrade"]
    control = [r for r in rows if r["scen"] == "control"]

    dgroups = collections.defaultdict(list)
    for r in degrade:
        dgroups[group_key(r)].append(r)
    cgroups = collections.defaultdict(list)
    for r in control:
        cgroups[group_key(r)].append(r)

    sizes = sorted({len(v) for v in dgroups.values()})
    coverages = sorted({r["c"] for r in degrade})

    print(f"\n=== graph-versus-latency lead at the group unit ===")
    print(f"{len(rows)} runs in the headline aggregate: {len(degrade)} degradation, {len(control)} control")
    print(f"{len(dgroups)} (coverage, arm, replicate) degradation groups, group size(s) {sizes}, "
          f"{len(coverages)} coverage levels")

    # --- invariance, measured rather than assumed -------------------------------------
    print(f"\n-- invariance of the censored alarm cycle across the nine noise settings --")
    inv = {}
    for det in DETECTORS:
        same = invariance(dgroups, det)
        inv[det] = {"invariant_groups": same, "n_groups": len(dgroups)}
        note = "invariant" if same == len(dgroups) else "NOT invariant -- group value is a mean"
        print(f"  {det:18s} {same:3d} of {len(dgroups)} groups identical across all 9   {note}")

    # --- collapse ---------------------------------------------------------------------
    coll = {det: collapse(dgroups, det) for det in DETECTORS}
    ccoll = {det: collapse(cgroups, det) for det in DETECTORS}

    results = {
        "_experiment": "graph-versus-latency lead estimates recomputed at the (coverage, arm, replicate) group unit",
        "_why": "RESULTS.md 8.10 suppressed these on a run-unit matching rule that does not fire at the group unit; "
                "8.12 records that the censoring, bootstrap and paired intervals were never given the same treatment.",
        "_reports_only": "This file reports measurements. It does not authorise reporting any of them in the manuscript.",
        "design": {
            "n_runs": len(rows),
            "n_degradation_runs": len(degrade),
            "n_control_runs": len(control),
            "n_groups": len(dgroups),
            "group_sizes": sizes,
            "group_key": "(coverage, arm, replicate)",
            "censoring": f"non-firing detector assigned cycle {CENSOR_AT}",
            "bootstrap": f"{BOOTSTRAP_DRAWS} draws, cluster bootstrap resampling GROUPS not runs, seed {BOOTSTRAP_SEED}",
            "sign": "positive means the GRAPH detector alarms EARLIER than the comparator",
        },
        "invariance": inv,
        "pooled": {},
        "by_coverage": {},
    }

    # --- pooled lead, three estimands --------------------------------------------------
    print(f"\n-- pooled lead over {len(dgroups)} groups, positive means the graph detector alarms EARLIER --")
    header = f"{'detector':18s} {'estimand':26s} {'mean':>9s} {'se':>8s} {'95% CI':>22s} {'n':>5s}"
    print(header)
    for det in GRAPH:
        censored, paired, clean = [], [], []
        for key in dgroups:
            g, c = coll[det][key], coll[COMPARATOR][key]
            censored.append(c["cycle"] - g["cycle"])
            if g["fired_all"] and c["fired_all"]:
                paired.append(c["cycle"] - g["cycle"])
                twin = ccoll[det].get(key)
                twin_c = ccoll[COMPARATOR].get(key)
                twin_alarmed = (twin and twin["fired_any"]) or (twin_c and twin_c["fired_any"])
                if not twin_alarmed:
                    clean.append(c["cycle"] - g["cycle"])
        block = {}
        for label, vals in (("censored (tau=320)", censored),
                            ("paired, both fired", paired),
                            ("paired, clean twin", clean)):
            b = boot(vals)
            block[label] = b
            ci = f"[{b['ci95'][0]:+.2f}, {b['ci95'][1]:+.2f}]" if b["n"] else "[--, --]"
            print(f"{det:18s} {label:26s} {b['mean']:+9.2f} {b['se']:8.2f} {ci:>22s} {b['n']:5d}")
        results["pooled"][det] = block

    # --- coverage-dependent lead curve -------------------------------------------------
    print(f"\n-- coverage-dependent lead curve, censored (tau=320), cluster bootstrap per coverage --")
    print(f"{'coverage':>9s}  {'groups':>6s}  " + "  ".join(f"{d:>28s}" for d in GRAPH))
    for cov in coverages:
        keys = [k for k in dgroups if k[0] == cov]
        row_out = {}
        cells = []
        for det in GRAPH:
            vals = [coll[COMPARATOR][k]["cycle"] - coll[det][k]["cycle"] for k in keys]
            b = boot(vals)
            row_out[det] = b
            cells.append(f"{b['mean']:+7.2f} [{b['ci95'][0]:+6.2f},{b['ci95'][1]:+6.2f}]")
        results["by_coverage"][str(cov)] = {"n_groups": len(keys), "detectors": row_out}
        print(f"{cov:9.4f}  {len(keys):6d}  " + "  ".join(f"{c:>28s}" for c in cells))

    # --- win/loss at the group unit ----------------------------------------------------
    print(f"\n-- win/loss: groups in which the graph detector alarms strictly earlier --")
    for det in GRAPH:
        both = [k for k in dgroups if coll[det][k]["fired_all"] and coll[COMPARATOR][k]["fired_all"]]
        wins = sum(1 for k in both if coll[det][k]["cycle"] < coll[COMPARATOR][k]["cycle"])
        fired = sum(1 for k in dgroups if coll[det][k]["fired_any"])
        print(f"  {det:18s} fired in {fired:3d} of {len(dgroups)} groups; "
              f"both fired in {len(both):3d}; graph earlier in {wins:3d} of {len(both) if both else 0}")
        results["pooled"][det]["win_loss"] = {
            "groups_graph_fired": fired,
            "groups_both_fired": len(both),
            "groups_graph_earlier": wins,
            "n_groups": len(dgroups),
        }

    path = os.path.join(HERE, "aggregates", OUT)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=1)
    print(f"\nwrote aggregates/{OUT}")


if __name__ == "__main__":
    main()
