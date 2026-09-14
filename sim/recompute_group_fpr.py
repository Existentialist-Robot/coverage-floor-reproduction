"""Group-unit false-alarm rates for the leading and lagging detectors, with cluster intervals.

WHY THIS FILE EXISTS. The manuscript previously printed the leading detector's group-unit
false-alarm rate as "4.48% (32.22 of 720 groups; Wilson 95% CI 3.19-6.24%)". Two things were
wrong with that and one was right.

  RIGHT: the 4.48% point estimate. It is the mean per-group alarm fraction and it reproduces
  exactly from the committed rows.

  WRONG 1: "32.22 of 720 groups" reads as a count, and a count of groups cannot be fractional.
  The leading detector's outcome is invariant within 640 of the 720 groups and VARIES within the
  other 80, so the figure sums per-group fractions rather than counting alarmed groups whole.

  WRONG 2: a Wilson interval is a binomial-proportion interval and is not defined for a
  fractional success count. The interval printed beside that figure was computed as though 32.22
  were a count of successes.

WHY NOT JUST COLLAPSE TO AN INTEGER. Because a group is not a set of replicates. It is one seed
evaluated across a 3x3 grid of noise colour (red, systematic, white) and severity (high, mid,
low), so "at least one run in the group alarmed" is a union over nine different counterfactual
conditions rather than a false-alarm rate. That rule is measurably contaminated by how many
conditions were bundled: at group size nine it reads 14.17% and at group size three 10.00%, while
the underlying per-run rates are 4.26% and 4.58%. The any-rule number moves four points on
packaging alone. A majority rule and an all rule are equally arbitrary and are reported here only
so the choice is visible rather than asserted.

WHAT THIS PRINTS AND WHAT THE PAPER USES. The paper uses the CLUSTER MEAN with a CLUSTER
BOOTSTRAP interval, resampling whole groups, for BOTH detectors, so the two are computed the same
way. The lagging detector's 35 of 720 is a genuine integer count because it is invariant within
every group; its cluster interval is printed anyway so the pair is comparable.

Run:  python sim/recompute_group_fpr.py
"""

import collections
import json
import pathlib
import random

HERE = pathlib.Path(__file__).resolve().parent
SOURCE = HERE / "aggregates" / "ea_headline_recal.json"

GROUP_KEYS = ("c", "arm", "rep")
DETECTORS = (("lead_fired", "leading"), ("lag_fired", "lagging"))

# Declared before the numbers were looked at, and printed in the manuscript beside the interval.
BOOTSTRAP_DRAWS = 20000
BOOTSTRAP_SEED = 20260910


def load_control_groups():
    rows = json.loads(SOURCE.read_text(encoding="utf-8"))["per_run_rows"]
    control = [r for r in rows if r.get("scen") == "control"]
    groups = collections.defaultdict(list)
    for r in control:
        groups[tuple(r[k] for k in GROUP_KEYS)].append(r)
    return control, list(groups.values())


def cluster_bootstrap(fractions, draws=BOOTSTRAP_DRAWS, seed=BOOTSTRAP_SEED):
    """Percentile interval, resampling whole GROUPS with replacement.

    The unit of resampling is the group because the group is the independent unit: runs inside one
    group share a seed. Resampling runs would treat correlated copies as independent and return an
    interval that is too narrow, which is the error class this whole recomputation exists to fix.
    """
    rnd = random.Random(seed)
    n = len(fractions)
    point = sum(fractions) / n
    reps = []
    for _ in range(draws):
        reps.append(sum(fractions[rnd.randrange(n)] for _ in range(n)) / n)
    reps.sort()
    return point, reps[int(0.025 * draws)], reps[int(0.975 * draws)]


def main():
    control, groups = load_control_groups()
    sizes = collections.Counter(len(g) for g in groups)
    print(f"control rows: {len(control)}")
    print(f"groups: {len(groups)}  sizes: {dict(sizes)}")

    for key, label in DETECTORS:
        fracs = [sum(bool(r[key]) for r in g) / len(g) for g in groups]
        point, lo, hi = cluster_bootstrap(fracs)
        disagree = sum(1 for g in groups if 0 < sum(bool(r[key]) for r in g) < len(g))
        any_rule = sum(1 for g in groups if any(bool(r[key]) for r in g))
        maj_rule = sum(1 for g in groups if sum(bool(r[key]) for r in g) * 2 > len(g))
        all_rule = sum(1 for g in groups if all(bool(r[key]) for r in g))
        n = len(groups)
        print()
        print(f"=== {label} detector ===")
        print(f"  PAPER FIGURE  mean per-group alarm fraction: {point * 100:.2f}%")
        print(f"                95% cluster bootstrap: [{lo * 100:.2f}%, {hi * 100:.2f}%]"
              f"  ({BOOTSTRAP_DRAWS} draws, seed {BOOTSTRAP_SEED})")
        print(f"  groups whose runs disagree: {disagree} of {n}")
        print(f"  reported for visibility, NOT used by the paper --"
              f" any {any_rule}/{n}, majority {maj_rule}/{n}, all {all_rule}/{n}")

    print()
    print("bundle-size contamination of the any-rule, which is why no integer collapse is used:")
    for size in sorted(sizes, reverse=True):
        sub = [g for g in groups if len(g) == size]
        runs = sum(len(g) for g in sub)
        fired = sum(sum(bool(r["lead_fired"]) for r in g) for g in sub)
        anyr = sum(1 for g in sub if any(bool(r["lead_fired"]) for r in g))
        print(f"  size {size}: per-run {fired}/{runs} = {fired / runs * 100:.2f}%"
              f"   any-rule {anyr}/{len(sub)} = {anyr / len(sub) * 100:.2f}%")


if __name__ == "__main__":
    main()
