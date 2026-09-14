"""
verify_reported_numbers.py -- the reader contract (frozen spec section 7).

Every number reported in RESULTS.md (and therefore every number the paper's section 5
may use) must be derivable from the committed JSON aggregates.  This checker enforces
that, and enforces the discipline the aggregates are supposed to embody:

  recomputation  RECOMPUTATION -- every cell summary in every aggregate is recomputed from the
      per-run rows stored in the same file and must match.
  seed hygiene  SEED HYGIENE -- no scored run uses a calibration seed (>= 900000); and the
      follow-up campaign's calibration null is entirely >= 900000, internally unique,
      and DISJOINT from campaign 1's calibration block and from every scored seed.
  matched pairing  MATCHED PAIRING -- in coverage-and-noise campaign, each cell's degradation and control runs use the same
      seed set; in coverage-and-coupling campaign, coupled and severed runs use the same seed set.
  rule identity  RULE IDENTITY -- the detection rule embedded in each aggregate is byte-identical
      to the rule file that aggregate declares (configs/detection_rule.json for the
      first campaign, configs/detection_rule_recal.json for the follow-up), so no
      aggregate was produced under a rule it does not name.
  claim binding  CLAIM BINDING -- every entry in claims.json (the numbers quoted in RESULTS.md)
      resolves to the stated path in the stated aggregate and matches within tolerance.
  determinism  DETERMINISM -- a sample of scored runs is re-executed from its recorded seed and
      must reproduce its recorded per-run scalars exactly.
  aggregate verification  CALIBRATION (follow-up campaign) -- the large-null FPR curves are recomputed from
      the committed per-run sufficient statistics; the recalibrated thresholds are
      exactly what the pre-declared selection rule picks from those curves; and the
      recalibrated rule differs from the frozen one in the THRESHOLDS ONLY.
  iso-line  ISO-LINE -- every reported Delta t = 0 crossing and its +-1 s.e. band is
      re-derived from the cell means by the same linear interpolation.
  e-b decision rule  coverage-and-coupling campaign DECISION RULE -- the pooled paired coupling-loss CIs, the power analysis and
      the verdict are recomputed from the per-run rows, so the verdict is the
      pre-declared rule applied to the data and not a reading of the table.
  cross-campaign analysis CROSS-CAMPAIGN ANALYSIS -- recal_analysis.json is a pure function of two committed
      aggregates; it is recomputed from them, the Delta t decomposition identity must
      close, and the two campaigns must have scored identical seeds where they overlap.
  aggregate verification REMEDIATION SWEEP (2026-08-05) -- remediation_runs.json: seed hygiene and matched
      degradation/control pairing per (variant, c); both declared rules byte-identical to
      their config files; every variant's parameter overrides ACTUALLY present in its
      rows and every non-overridden parameter at its frozen default; the base variant
      reproducing ea_headline_recal.json run-for-run (which is what proves the
      diagnostic tap perturbs nothing); every cell summary, both scorings' Delta t
      identity, every iso-line and the comparator-gain table recomputed from the rows.
  remediation analysis REMEDIATION ANALYSIS -- remediation_analysis.json is a pure function of committed
      aggregates; it is recomputed and must reproduce the committed file, the
      restricted-mean-time-to-alarm identity must close against the censored mean, the
      censoring decomposition must sum to the censored mean, and every reported
      probability must be its own count over its own denominator.
  both-class calibration BOTH-CLASS CALIBRATION -- calibration_null_bothclass.json: the FPR curve is
      recomputed from the committed per-run rows, the chosen threshold is what the
      pre-declared selection gives, the rule on disk matches the aggregate, the
      both-class rule differs from the recalibrated rule in the CLASS SET and the
      THRESHOLD only, and it was calibrated on the identical 600-run null.

Exit code 0 if every check passes, 1 otherwise.

Run:  python verify_reported_numbers.py [--sample 6] [--skip-determinism]
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
AGG = os.path.join(HERE, "aggregates")
CFG = os.path.join(HERE, "configs")
SCORED_SEED_CEILING = 900_000

# scored aggregates, campaign 1 (frozen rule) then the 2026-08-05 follow-up (recalibrated
# rule).  Each aggregate names its own rule file in design.rule_file; campaign 1's
# aggregates predate that field and default to the frozen rule.
EA_FILES = ("ea_headline.json", "ea_headline_recal.json")
AUX_FILES = ("aux_arms.json", "aux_arms_recal.json")
EB_FILES = ("eb_joint_sweep.json", "eb_joint_sweep_n500.json")
SCORED_FILES = EA_FILES + AUX_FILES + EB_FILES
CAL_NULL = "calibration_null_large.json"
# remediation campaign 2026-08-05 (review findings)
REM_RUNS = "remediation_runs.json"
REM_ANALYSIS = "remediation_analysis.json"
CAL_BOTHCLASS = "calibration_null_bothclass.json"

FAILURES: list[str] = []
CHECKS = [0]


def ok(cond: bool, msg: str) -> bool:
    CHECKS[0] += 1
    if not cond:
        FAILURES.append(msg)
    return bool(cond)


def close(a, b, tol=1e-6) -> bool:
    if a is None and b is None:
        return True
    if a is None or b is None:
        return False
    if isinstance(a, bool) or isinstance(b, bool):
        return bool(a) == bool(b)
    try:
        fa, fb = float(a), float(b)
    except (TypeError, ValueError):
        return a == b
    if math.isnan(fa) and math.isnan(fb):
        return True
    return abs(fa - fb) <= tol + tol * max(abs(fa), abs(fb))


def load(name):
    p = os.path.join(AGG, name)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def _field_eq(got, want: str) -> bool:
    try:
        return abs(float(got) - float(want)) < 1e-9
    except (TypeError, ValueError):
        return str(got) == want


def resolve(obj, path: str):
    """Resolve a slash path with bracket selectors:

      cells[3]/dt_censored_mean                      positional index
      cells[c=0.5,arm=mnar,colour=red]/fpr_leading   FIELD SELECTOR -- names the cell by
                                                     its own design fields instead of by
                                                     position, so adding a grid point
                                                     cannot silently re-point a claim at
                                                     a different cell.
    """
    cur = obj
    for part in path.split("/"):
        if not part:
            continue
        head = part.split("[", 1)[0]
        brackets = re.findall(r"\[([^\]]*)\]", part)
        if head:
            cur = cur[head]
        for b in brackets:
            if "=" in b:
                sel = dict(kv.split("=", 1) for kv in b.split(","))
                hits = [row for row in cur
                        if all(_field_eq(row.get(k), v) for k, v in sel.items())]
                if len(hits) != 1:
                    raise KeyError(f"selector [{b}] matched {len(hits)} rows, not 1")
                cur = hits[0]
            elif b.startswith("'"):
                cur = cur[b.strip("'")]
            else:
                cur = cur[int(b)]
    return cur


# ------------------------------------------------------------------- seed hygiene / matched pairing / rule identity
def load_rule_file(name: str):
    p = os.path.join(CFG, name)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def scored_seeds(d) -> list[int]:
    rows = d.get("per_run_rows", []) + d.get("ensign_rows", [])
    return [r["seed"] for r in rows if "seed" in r]


def check_seeds_and_rule() -> None:
    if not ok(load_rule_file("detection_rule.json") is not None,
              "configs/detection_rule.json missing"):
        return
    all_scored: set[int] = set()
    for name in SCORED_FILES:
        d = load(name)
        if d is None:
            continue
        seeds = scored_seeds(d)
        all_scored |= set(seeds)
        ok(all(s < SCORED_SEED_CEILING for s in seeds),
           f"seed hygiene {name}: a scored run used a calibration seed (>= {SCORED_SEED_CEILING})")
        ok(len(seeds) > 0, f"seed hygiene {name}: no seeds recorded")
        rf = (d.get("design") or {}).get("rule_file", "detection_rule.json")
        rule = load_rule_file(rf)
        if not ok(rule is not None, f"rule identity {name}: declared rule file configs/{rf} missing"):
            continue
        if "detection_rule" in d:
            ok(json.dumps(d["detection_rule"], sort_keys=True)
               == json.dumps(rule, sort_keys=True),
               f"rule identity {name}: embedded detection rule differs from configs/{rf}")

    for name in EA_FILES:
        d = load(name)
        if not d:
            continue
        groups: dict[tuple, dict[str, set]] = {}
        for r in d["per_run_rows"]:
            key = (r["c"], r["arm"], r["colour"], r["sev"])
            groups.setdefault(key, {}).setdefault(r["scen"], set()).add(r["seed"])
        bad = [k for k, v in groups.items()
               if set(v) == {"degrade", "control"} and v["degrade"] != v["control"]]
        ok(not bad, f"matched pairing {name}: {len(bad)} cells have unmatched seed sets between the "
                    f"degradation and control arms, e.g. {bad[:2]}")

    for name in EB_FILES:
        d = load(name)
        if not d:
            continue
        groups = {}
        for r in d["per_run_rows"]:
            groups.setdefault((r["c"], r["g_name"]), {}).setdefault(
                r["mode"], set()).add(r["seed"])
        bad = [k for k, v in groups.items()
               if set(v) == {"coupled", "severed"} and v["coupled"] != v["severed"]]
        ok(not bad, f"matched pairing {name}: {len(bad)} cells have unmatched seed sets between the "
                    f"coupled and severed arms, e.g. {bad[:2]}")

    # --- follow-up campaign: the calibration null's own seed hygiene ---------
    cal = load(CAL_NULL)
    if cal is None:
        return
    des = cal["design"]
    nulls = [r["seed"] for r in cal["per_run_rows"] if r["scen"] == "control"]
    degs = [r["seed"] for r in cal["per_run_rows"] if r["scen"] == "degrade"]
    ok(len(nulls) >= 600,
       f"seed hygiene {CAL_NULL}: calibration null has {len(nulls)} runs, fewer than the "
       f"declared 600")
    ok(len(nulls) == des["n_null_control_runs"],
       f"seed hygiene {CAL_NULL}: declared n_null_control_runs != rows present")
    ok(len(set(nulls)) == len(nulls), f"seed hygiene {CAL_NULL}: duplicate null seeds")
    ok(all(s >= SCORED_SEED_CEILING for s in nulls + degs),
       f"seed hygiene {CAL_NULL}: a calibration run used a seed below the scored ceiling")
    ok(not (set(nulls) & all_scored),
       f"seed hygiene {CAL_NULL}: calibration null seeds intersect scored evaluation seeds")
    c1_block = set(range(des["calibration_degrade_seed_base"],
                         des["calibration_degrade_seed_base"]
                         + des["n_calibration_degrade_runs"]))
    ok(not (set(nulls) & c1_block),
       f"seed hygiene {CAL_NULL}: the fresh null overlaps campaign 1's calibration block")
    ok(min(nulls) == des["null_seed_base"]
       and max(nulls) == des["null_seed_base"] + len(nulls) - 1,
       f"seed hygiene {CAL_NULL}: null seed block is not the contiguous declared range")


# ------------------------------------------------------------------------- recomputation
def check_ea_recompute(name: str = "ea_headline.json") -> None:
    d = load(name)
    if d is None:
        return
    rows = d["per_run_rows"]
    n_cycles = d["design"]["n_cycles"]
    for cell in d["cells"]:
        sel = [r for r in rows
               if r["c"] == cell["c"] and r["arm"] == cell["arm"]
               and r["colour"] == cell["colour"] and r["sev"] == cell["severity"]]
        deg = [r for r in sel if r["scen"] == "degrade"]
        ctl = [r for r in sel if r["scen"] == "control"]
        tag = (f"recomputation {name} cell c={cell['c']} {cell['arm']} {cell['colour']}"
               f"/{cell['severity']}")
        ok(len(deg) == cell["n_reps"], f"{tag}: n_reps mismatch")
        ok(close(np.mean([r["dt_censored"] for r in deg]),
                 cell["dt_censored_mean"], 1e-9),
           f"{tag}: dt_censored_mean not reproducible from per-run rows")
        ok(close(np.mean([r["lead_fired"] for r in deg]),
                 cell["lead_detect_rate"], 1e-9), f"{tag}: lead_detect_rate")
        ok(close(np.mean([r["lead_fired"] for r in ctl]),
                 cell["fpr_leading"], 1e-9), f"{tag}: fpr_leading")
        ok(close(np.mean([r["achieved_edge_frac"] for r in deg]),
                 cell["achieved_edge_frac_mean"], 1e-9),
           f"{tag}: achieved_edge_frac_mean")
        # dt_censored must equal (lagging or censor) - (leading or censor)
        bad = [r for r in sel
               if r["dt_censored"] != ((n_cycles if r["t_lag"] is None else r["t_lag"])
                                       - (n_cycles if r["t_lead"] is None
                                          else r["t_lead"]))]
        ok(not bad, f"{tag}: dt_censored inconsistent with the alarm cycles "
                    f"in {len(bad)} runs")

    # ---- pooled false-alarm accounting (follow-up campaign aggregates) -------
    p = d.get("pooled_false_alarms")
    if p:
        ctl_all = [r for r in rows if r["scen"] == "control"]
        mid = [r for r in ctl_all if r["sev"] == d["design"]["headline_severity"]]
        ok(len(ctl_all) == p["n_control_runs"],
           f"recomputation {name} pooled: n_control_runs mismatch")
        ok(close(np.mean([r["lead_fired"] for r in ctl_all]),
                 p["fpr_leading_pooled"], 1e-9),
           f"recomputation {name} pooled: fpr_leading_pooled not reproducible from the rows")
        ok(close(np.mean([r["lag_fired"] for r in ctl_all]),
                 p["fpr_lagging_pooled"], 1e-9),
           f"recomputation {name} pooled: fpr_lagging_pooled")
        for col, v in p["fpr_leading_by_colour_mid_severity"].items():
            ok(close(np.mean([r["lead_fired"] for r in mid if r["colour"] == col]),
                     v, 1e-9), f"recomputation {name} pooled: FPR by colour {col}")
        for cs, v in p["fpr_leading_by_c_mid_severity"].items():
            ok(close(np.mean([r["lead_fired"] for r in mid
                              if abs(r["c"] - float(cs)) < 1e-12]), v, 1e-9),
               f"recomputation {name} pooled: FPR by c {cs}")
        ok(p["meets_external_bar"] == (p["fpr_leading_pooled"] <= p["external_bar"]),
           f"recomputation {name} pooled: meets_external_bar flag inconsistent with the number")
        ok(p["meets_declared_target"] == (p["fpr_leading_pooled"]
                                          <= p["declared_target"]),
           f"recomputation {name} pooled: meets_declared_target flag inconsistent")


# ------------------------------------------------------------------------- iso-line
def _interp_cross(cs: list[float], ys: list[float]):
    """The reported iso-line rule: first sign change of Delta t along ascending c,
    linearly interpolated.  Identical to the runner's implementation."""
    for i in range(len(cs) - 1):
        if (ys[i] <= 0.0) != (ys[i + 1] <= 0.0):
            return cs[i] + (cs[i + 1] - cs[i]) * (0.0 - ys[i]) / (ys[i + 1] - ys[i])
    return None


def check_iso_line(name: str = "ea_headline.json") -> None:
    d = load(name)
    if d is None:
        return
    sev = d["design"]["headline_severity"]
    for key, entry in d["iso_line_dt_zero"].items():
        arm, col = key.split("|")
        cs = sorted([x for x in d["cells"]
                     if x["arm"] == arm and x["colour"] == col
                     and x["severity"] == sev], key=lambda x: x["c"])
        tag = f"iso-line {name} iso {key}"
        ok([x["c"] for x in cs] == entry["c"], f"{tag}: c list")
        ok(all(close(a["dt_censored_mean"], b, 1e-9)
               for a, b in zip(cs, entry["dt"])), f"{tag}: dt list")
        ok(all(close(a["n_reps"], b, 1e-9)
               for a, b in zip(cs, entry.get("n_reps", [x["n_reps"] for x in cs]))),
           f"{tag}: n_reps list")
        ok(close(_interp_cross([x["c"] for x in cs],
                               [x["dt_censored_mean"] for x in cs]),
                 entry["c_at_dt_zero"], 1e-9), f"{tag}: c_at_dt_zero")
        band = [_interp_cross([x["c"] for x in cs],
                              [x["dt_censored_mean"] + s * x["dt_censored_se"]
                               for x in cs]) for s in (-1.0, 1.0)]
        ok(all(close(a, b, 1e-9) for a, b in zip(band, entry["c_at_dt_zero_band"])),
           f"{tag}: +-1 s.e. band")


def check_gate1_recompute() -> None:
    d = load("gate1_bifurcation.json")
    if d is None:
        return
    rows = d["per_seed_rows"]
    for cell in d["branch_table"]:
        sel = [r for r in rows if r["n"] == cell["n"] and r["a"] == cell["a"]
               and r["branch"] == cell["branch"]]
        tag = f"recomputation Gate1 N={cell['n']} A={cell['a']} {cell['branch']}"
        ok(len(sel) == cell["n_seeds"], f"{tag}: n_seeds mismatch")
        ok(close(np.mean([r["S"] for r in sel]), cell["S_mean"], 1e-9),
           f"{tag}: S_mean")
        ok(close(np.mean([r["z"] for r in sel]), cell["z_mean"], 1e-9),
           f"{tag}: z_mean")
    graded = {k: v for k, v in d["tests"].items() if v.get("pass") is not None}
    ok((d["verdict"] == "PASS") == all(v["pass"] for v in graded.values()),
       "recomputation Gate1: verdict inconsistent with its own graded tests")
    lr = d.get("long_run_table") or []
    lrows = d.get("long_run_rows") or []
    for cell in lr:
        sel = [r for r in lrows if r["n"] == cell["n"] and r["a"] == cell["a"]]
        ok(close(np.mean([r["z_end"] for r in sel]), cell["z_end_mean"], 1e-9),
           f"recomputation Gate1 long-run N={cell['n']} A={cell['a']}: z_end_mean")


def check_eb_recompute(name: str = "eb_joint_sweep.json") -> None:
    d = load(name)
    if d is None:
        return
    rows = d["per_run_rows"]
    for cell in d["cells"]:
        sel = [r for r in rows if r["c"] == cell["c"]
               and r["g_name"] == cell["g_name"] and r["mode"] == cell["mode"]]
        tag = f"recomputation {name} c={cell['c']} g={cell['g_name']} {cell['mode']}"
        ok(len(sel) == cell["n_reps"], f"{tag}: n_reps mismatch")
        vals = [r["corr_indicator_latent"] for r in sel
                if r["corr_indicator_latent"] is not None
                and np.isfinite(r["corr_indicator_latent"])]
        ok(close(np.mean(vals), cell["corr_mean"], 1e-9), f"{tag}: corr_mean")
        ok(close(np.mean([r["attrition_frac"] for r in sel]),
                 cell["attrition_frac_mean"], 1e-9), f"{tag}: attrition_frac_mean")
    for g, v in d["interaction"].items():
        for c, pc in v["per_c"].items():
            cf = float(c)
            cp = {r["rep"]: r for r in rows if r["c"] == cf
                  and r["g_name"] == g and r["mode"] == "coupled"}
            sv = {r["rep"]: r for r in rows if r["c"] == cf
                  and r["g_name"] == g and r["mode"] == "severed"}
            shared = sorted(set(cp) & set(sv))
            loss = [abs(sv[i]["corr_indicator_latent"])
                    - abs(cp[i]["corr_indicator_latent"]) for i in shared]
            ok(close(np.mean(loss), pc["tracking_loss_mean"], 1e-9),
               f"recomputation {name} interaction g={g} c={c}: tracking_loss_mean")
            ok(len(shared) == pc["n_pairs"],
               f"recomputation {name} interaction g={g} c={c}: n_pairs")

    # ------------------------------------------------------------------- e-b decision rule
    dec = d.get("decision")
    if not dec:
        return
    z80 = 1.959963985 + 0.841621234
    for g, v in dec["per_g"].items():
        pooled = []
        for c in sorted({r["c"] for r in rows}):
            cp = {r["rep"]: r for r in rows if r["c"] == c
                  and r["g_name"] == g and r["mode"] == "coupled"}
            sv = {r["rep"]: r for r in rows if r["c"] == c
                  and r["g_name"] == g and r["mode"] == "severed"}
            for i in sorted(set(cp) & set(sv)):
                pooled.append(abs(sv[i]["corr_indicator_latent"])
                              - abs(cp[i]["corr_indicator_latent"]))
        a = np.asarray([x for x in pooled if np.isfinite(x)], float)
        tag = f"e-b decision rule {name} decision g={g}"
        ok(a.size == v["n_pairs_pooled"], f"{tag}: n_pairs_pooled")
        ok(a.size >= 500 * 3 * 0.9,
           f"{tag}: pooled pairs {a.size} far short of 500 per cell x 3 coverage levels")
        ok(close(a.mean(), v["tracking_loss_mean"], 1e-9), f"{tag}: tracking_loss_mean")
        ok(close(a.std(ddof=1), v["paired_sd"], 1e-9), f"{tag}: paired_sd")
        ok(close(z80 * a.std(ddof=1) / math.sqrt(a.size), v["mde_80pct_power"], 1e-9),
           f"{tag}: minimum detectable effect at 80% power")
        ok(v["excludes_zero"] == (v["lo95"] > 0 or v["hi95"] < 0),
           f"{tag}: excludes_zero flag inconsistent with its own CI")
        ok(v["sign"] == int(np.sign(v["tracking_loss_mean"])), f"{tag}: sign")
    both = all(v["excludes_zero"] for v in dec["per_g"].values())
    signs = {v["sign"] for v in dec["per_g"].values()}
    want = ("RESOLVED" if both and len(signs) == 1 else
            "COUPLING-STRENGTH-DEPENDENT REGIME CHANGE" if both else "UNRESOLVED")
    ok(dec["verdict"] == want,
       f"e-b decision rule {name}: verdict {dec['verdict']!r} is not what the pre-declared rule gives "
       f"({want!r})")


def check_aux_recompute(name: str = "aux_arms.json") -> None:
    d = load(name)
    if d is None:
        return
    rows = d["per_run_rows"]
    for cell in d["revocation_arm"]["cells"]:
        deg = [r for r in rows if r["condition"] == cell["condition"]
               and r["scen"] == "degrade"]
        ctl = [r for r in rows if r["condition"] == cell["condition"]
               and r["scen"] == "control"]
        tag = f"recomputation {name} {cell['condition']}"
        ok(len(deg) == cell["n_reps"], f"{tag}: n_reps mismatch")
        ok(close(np.mean([r["dt_censored"] for r in deg]),
                 cell["dt_censored_mean"], 1e-9), f"{tag}: dt_censored_mean")
        ok(close(np.mean([r["lead_fired"] for r in ctl]),
                 cell["fpr_leading"], 1e-9), f"{tag}: fpr_leading")
    irows = d["ensign_rows"]
    for cond in ("static_full", "revocation"):
        sel = [r for r in irows if r["condition"] == cond]
        v = d["ensign_arm"][cond]
        ok(close(np.nanmean([r["bias_reduction"] for r in sel]),
                 v["bias_reduction_mean"], 1e-9),
           f"recomputation {name} Ensign {cond}: bias_reduction_mean")


# ------------------------------------------------------------------------- aggregate verification
def check_calibration_null() -> None:
    d = load(CAL_NULL)
    if d is None:
        return
    des = d["design"]
    rows = d["per_run_rows"]
    null = [r for r in rows if r["scen"] == "control"]
    deg = [r for r in rows if r["scen"] == "degrade"]
    target = des["target_fpr"]
    ok(target <= des["external_bar"],
       "aggregate verification calibration: the declared target is not at or below the external bar")

    # leading FPR curve, recomputed from the per-run sufficient statistic
    stats = np.asarray([r["lead_sustained_score"] for r in null], float)
    for x in d["leading_fpr_curve"]:
        ok(close(float(np.mean(stats >= x["k"])), x["fpr"], 1e-12),
           f"aggregate verification calibration: leading FPR at k={x['k']} not reproducible from the rows")
        ok(int(np.sum(stats >= x["k"])) == x["n_false_alarms"],
           f"aggregate verification calibration: n_false_alarms at k={x['k']}")
    adm = [x["k"] for x in d["leading_fpr_curve"] if x["fpr"] <= target]
    ok(adm and close(min(adm), d["chosen"]["leading_k"], 1e-12),
       f"aggregate verification calibration: leading k {d['chosen']['leading_k']} is not the smallest k "
       f"with null FPR <= {target}")
    ok(d["chosen"]["leading_calibration_fpr"] <= target,
       "aggregate verification calibration: chosen leading k does not meet its own target")

    # lagging grid, recomputed: FPR from the null rows, detection from the cal degrade
    # rows' committed per-candidate-rule alarm cycles
    if not ok(all("lag_alarm_cycle_by_rule" in r for r in deg),
              "aggregate verification calibration: the calibration degradation rows do not carry the "
              "per-candidate-rule lagging alarm cycles, so the lagging selection is "
              "not recomputable from the committed aggregate"):
        deg = []
    for x in d["lagging_grid"] if deg else []:
        key = f"{x['k_sd']}|{x['sustain']}"
        fpr = float(np.mean([r[f"lag_sustained_drop_s{x['sustain']}"] >= x["k_sd"]
                             for r in null]))
        ok(close(fpr, x["fpr"], 1e-12),
           f"aggregate verification calibration: lagging FPR at {key} not reproducible")
        hits = [r["lag_alarm_cycle_by_rule"][key] for r in deg]
        fired = [h for h in hits if h is not None]
        ok(close(len(fired) / len(deg), x["detect_rate"], 1e-12),
           f"aggregate verification calibration: lagging detect_rate at {key}")
        ok(close(float(np.median(fired)) if fired else None, x["median_t"], 1e-12),
           f"aggregate verification calibration: lagging median_t at {key}")
    lag_adm = [x for x in d["lagging_grid"] if x["fpr"] <= target
               and x["detect_rate"] >= 0.5 and x["median_t"] is not None]
    ok(bool(lag_adm), "aggregate verification calibration: no admissible lagging comparator on the grid")
    if lag_adm:
        best = min(x["median_t"] for x in lag_adm)
        chosen = [x for x in lag_adm
                  if close(x["k_sd"], d["chosen"]["lagging_k_sd"], 1e-12)
                  and x["sustain"] == d["chosen"]["lagging_sustain"]]
        ok(len(chosen) == 1, "aggregate verification calibration: chosen lagging rule is not on the grid "
                             "or is not admissible")
        if chosen:
            ok(close(chosen[0]["median_t"], best, 1e-12),
               "aggregate verification calibration: the chosen lagging rule is not the earliest-detecting "
               "admissible one -- the comparator lost an advantage it is owed")

    # the recalibrated rule must differ from the frozen one in THRESHOLDS ONLY
    fr = d["frozen_rule_campaign1"]
    rc = d["recalibrated_rule"]
    disk = load_rule_file("detection_rule_recal.json")
    ok(disk is not None and json.dumps(disk, sort_keys=True)
       == json.dumps(rc, sort_keys=True),
       "aggregate verification calibration: configs/detection_rule_recal.json differs from the rule "
       "recorded in the aggregate")
    for k, v in fr["leading"].items():
        if k == "k":
            continue
        ok(close(rc["leading"].get(k), v, 1e-12),
           f"aggregate verification calibration: leading rule FORM changed at '{k}' -- only the threshold "
           f"may move")
    for k, v in fr["lagging"].items():
        if k in ("k_sd", "sustain"):
            continue
        ok(close(rc["lagging"].get(k), v, 1e-12),
           f"aggregate verification calibration: lagging rule FORM changed at '{k}'")
    ok(close(rc["leading"]["k"], d["chosen"]["leading_k"], 1e-12),
       "aggregate verification calibration: rule file threshold != chosen threshold")
    ok(des["n_null_control_runs"] == rc["n_calibration_seeds"],
       "aggregate verification calibration: the rule file's declared null size != the null actually run")


# ------------------------------------------------------------------------- cross-campaign analysis
def _deep_diff(a, b, path: str = "") -> list[str]:
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                out.append(f"{path}/{k}: present in only one side")
            else:
                out += _deep_diff(a[k], b[k], f"{path}/{k}")
        return out
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [f"{path}: length {len(a)} vs {len(b)}"]
        out = []
        for i, (x, y) in enumerate(zip(a, b)):
            out += _deep_diff(x, y, f"{path}[{i}]")
        return out
    return [] if close(a, b, 1e-9) else [f"{path}: {a!r} vs {b!r}"]


def check_recal_analysis() -> None:
    """cross-campaign analysis -- the cross-campaign analysis is a pure function of two committed
    aggregates, so it is re-derived here and must reproduce the committed file."""
    d = load("recal_analysis.json")
    if d is None:
        return
    try:
        import recal_analysis
        fresh = recal_analysis.compute()
    except Exception as e:                                     # noqa: BLE001
        ok(False, f"cross-campaign analysis recal_analysis: recomputation raised {e}")
        return
    committed = {k: v for k, v in d.items() if k != "_generated_by"}
    diffs = _deep_diff(fresh, committed)
    ok(not diffs, f"cross-campaign analysis recal_analysis.json is not reproducible from its sources: "
                  f"{len(diffs)} differences, e.g. {diffs[:3]}")
    ok(d["n_seed_mismatches"] == 0,
       "cross-campaign analysis recal_analysis: the two campaigns did not score identical seeds where they "
       "overlap, so the paired decomposition is invalid")
    ok(d["delta_t_decomposition_pooled"]["max_abs_identity_residual"] < 1e-9,
       "cross-campaign analysis recal_analysis: the Delta t decomposition does not close "
       "(change != comparator effect + leading-threshold effect)")
    for b in d["pooled_false_alarm_rates"]:
        n, k = b["n_control_runs"], b["n_false_alarms"]
        ok(close(k / n, b["fpr_leading"], 1e-12),
           f"cross-campaign analysis recal_analysis: FPR for {b['campaign']} not k/n")
        ok(b["above_external_bar_ci_lower"]
           == (b["fpr_leading_wilson95"][0] > b["external_bar"]),
           f"cross-campaign analysis recal_analysis: {b['campaign']} bar flag inconsistent with its CI")


# ------------------------------------------------------------------------ aggregate verification
PARAM_FIELDS = ("output_knee", "output_sigma", "formation_driver",
                "partial_obs_sigma", "n_agents", "comp_ref")


def check_remediation_runs() -> None:
    d = load(REM_RUNS)
    if d is None:
        return
    from model import Config
    default = Config()
    rows = d["per_run_rows"]
    des = d["design"]
    tag0 = f"aggregate verification {REM_RUNS}"

    # --- the two declared rules must be the files they name ---------------------
    for key, fname in (("detection_rule", des["rule_file"]),
                       ("detection_rule_bothclass", des["rule_file_bothclass"])):
        rule = load_rule_file(fname)
        if ok(rule is not None, f"{tag0}: declared rule file configs/{fname} missing"):
            ok(json.dumps(d[key], sort_keys=True) == json.dumps(rule, sort_keys=True),
               f"{tag0}: embedded {key} differs from configs/{fname}")

    # --- seed hygiene and matched pairing --------------------------------------
    ok(all(r["seed"] < SCORED_SEED_CEILING for r in rows),
       f"{tag0}: a scored run used a calibration seed (>= {SCORED_SEED_CEILING})")
    ok(len(rows) == des["n_runs"], f"{tag0}: design.n_runs != rows present")
    groups: dict[tuple, dict[str, set]] = {}
    for r in rows:
        groups.setdefault((r["variant"], r["c"]), {}).setdefault(
            r["scen"], set()).add(r["seed"])
    bad = [k for k, v in groups.items()
           if set(v) == {"degrade", "control"} and v["degrade"] != v["control"]]
    ok(not bad, f"{tag0}: {len(bad)} (variant, c) groups have unmatched seed sets "
                f"between the degradation and control arms, e.g. {bad[:2]}")

    # --- every variant's overrides really are in its rows, and nothing else moved
    for v in des["variants"]:
        vr = [r for r in rows if r["variant"] == v["variant"]]
        tag = f"{tag0} variant {v['variant']}"
        if not ok(bool(vr), f"{tag}: declared in the design but has no rows"):
            continue
        for field in PARAM_FIELDS:
            want = v["overrides"].get(field, getattr(default, field))
            ok(all(_field_eq(r.get(field), want) for r in vr),
               f"{tag}: parameter '{field}' is not {want} in every row -- the sweep did "
               f"not sweep what the design says it swept")
        ok(all(k in PARAM_FIELDS or k == "programs_per_cycle"
               for k in v["overrides"]),
           f"{tag}: an override is not among the recorded parameter fields, so the "
           f"aggregate cannot evidence that it was applied")

    # --- the base variant must reproduce the committed headline runs ------------
    ea = load("ea_headline_recal.json")
    if ea is not None:
        idx = {(r["c"], r["scen"], r["rep"]): r for r in ea["per_run_rows"]
               if r["arm"] == des["arm"] and r["colour"] == des["colour"]
               and r["sev"] == des["severity"]}
        shared = 0
        mism = []
        for r in rows:
            if r["variant"] != des["base_variant"]:
                continue
            o = idx.get((r["c"], r["scen"], r["rep"]))
            if o is None:
                continue
            shared += 1
            if (r["seed"] != o["seed"] or r["t_lead"] != o["t_lead"]
                    or r["t_lag"] != o["t_lag"]
                    or not close(r["achieved_edge_frac"], o["achieved_edge_frac"],
                                 1e-12)):
                mism.append((r["c"], r["scen"], r["rep"]))
        ok(shared > 0, f"{tag0}: the base variant shares no run with "
                       f"ea_headline_recal.json, so it cannot be checked against it")
        ok(not mism, f"{tag0}: the base variant does not reproduce "
                     f"ea_headline_recal.json in {len(mism)} of {shared} shared runs "
                     f"(e.g. {mism[:3]}) -- the diagnostic tap perturbed the substrate")

    # --- cell summaries, both scorings ----------------------------------------
    n_cycles = rows[0]["n_cycles"]
    for cellrow in d["cells"]:
        deg = [r for r in rows if r["variant"] == cellrow["variant"]
               and r["c"] == cellrow["c"] and r["scen"] == "degrade"]
        ctl = [r for r in rows if r["variant"] == cellrow["variant"]
               and r["c"] == cellrow["c"] and r["scen"] == "control"]
        tag = f"{tag0} cell {cellrow['variant']} c={cellrow['c']}"
        ok(len(deg) == cellrow["n_reps"], f"{tag}: n_reps mismatch")
        for field, key, src in (
                ("dt_censored_mean", "dt_censored", deg),
                ("bc_dt_censored_mean", "bc_dt_censored", deg),
                ("lead_detect_rate", "lead_fired", deg),
                ("lag_detect_rate", "lag_fired", deg),
                ("bc_lead_detect_rate", "bc_lead_fired", deg),
                ("fpr_leading", "lead_fired", ctl),
                ("fpr_lagging", "lag_fired", ctl),
                ("bc_fpr_leading", "bc_lead_fired", ctl),
                ("achieved_edge_frac_mean", "achieved_edge_frac", deg),
                ("lag_sigma_drop_sustained_mean", "lag_sigma_drop_sustained", deg),
                ("lead_max_score_mean", "lead_max_score", deg)):
            ok(close(np.mean([r[key] for r in src]), cellrow[field], 1e-9),
               f"{tag}: {field} not reproducible from the per-run rows")
    for pre in ("", "bc_"):
        bad = [r for r in rows
               if r[pre + "dt_censored"]
               != ((n_cycles if r[pre + "t_lag"] is None else r[pre + "t_lag"])
                   - (n_cycles if r[pre + "t_lead"] is None else r[pre + "t_lead"]))]
        ok(not bad, f"{tag0}: {pre}dt_censored inconsistent with the alarm cycles in "
                    f"{len(bad)} runs")

    # --- every iso-line, per variant ------------------------------------------
    for name, entry in d["iso_line_dt_zero"].items():
        cs = sorted([x for x in d["cells"] if x["variant"] == name],
                    key=lambda x: x["c"])
        tag = f"aggregate verification {REM_RUNS} iso {name}"
        ok([x["c"] for x in cs] == entry["c"], f"{tag}: c list")
        ok(all(close(a["dt_censored_mean"], b, 1e-9)
               for a, b in zip(cs, entry["dt"])), f"{tag}: dt list")
        ok(close(_interp_cross([x["c"] for x in cs],
                               [x["dt_censored_mean"] for x in cs]),
                 entry["c_at_dt_zero"], 1e-9), f"{tag}: c_at_dt_zero")
        band = [_interp_cross([x["c"] for x in cs],
                              [x["dt_censored_mean"] + s * x["dt_censored_se"]
                               for x in cs]) for s in (-1.0, 1.0)]
        ok(all(close(a, b, 1e-9) for a, b in zip(band, entry["c_at_dt_zero_band"])),
           f"{tag}: +-1 s.e. band")

    # --- the comparator-gain table --------------------------------------------
    k_sd = float(d["detection_rule"]["lagging"]["k_sd"])
    for g in d["comparator_gain_table"]:
        deg = [r for r in rows if r["variant"] == g["variant"]
               and r["scen"] == "degrade"]
        tag = f"aggregate verification {REM_RUNS} gain {g['variant']}"
        ok(len(deg) == g["n_degrade_runs"], f"{tag}: n_degrade_runs")
        ok(close(np.mean([r["lag_sigma_drop_sustained"] for r in deg]),
                 g["comparator_gain_sigma_drop"], 1e-9),
           f"{tag}: comparator_gain_sigma_drop not reproducible")
        ok(close(np.mean([r["lag_fired"] for r in deg]),
                 g["lag_detect_rate_pooled"], 1e-9), f"{tag}: lag_detect_rate_pooled")
        ok(close(g["k_sd_required"], k_sd, 1e-12), f"{tag}: k_sd_required")
        ok(g["comparator_can_fire"]
           == (g["comparator_gain_sigma_drop"] >= k_sd),
           f"{tag}: comparator_can_fire flag inconsistent with its own numbers")


# ------------------------------------------------------------------------ remediation analysis
def check_remediation_analysis() -> None:
    d = load(REM_ANALYSIS)
    if d is None:
        return
    tag0 = f"remediation analysis {REM_ANALYSIS}"
    try:
        import remediation_analysis
        fresh = remediation_analysis.compute()
    except Exception as e:                                         # noqa: BLE001
        ok(False, f"{tag0}: recomputation raised {e}")
        return
    committed = {k: v for k, v in d.items() if k != "_generated_by"}
    diffs = _deep_diff(fresh, committed)
    ok(not diffs, f"{tag0} is not reproducible from its sources: {len(diffs)} "
                  f"differences, e.g. {diffs[:3]}")

    # the identity that makes the censored mean a legitimate estimand
    ok(d["rmta_identity"]["max_abs_discrepancy"] < 1e-9,
       f"{tag0}: mean(dt_censored) is not the restricted-mean-time-to-alarm difference "
       f"at the run-length horizon, so the de-censoring argument does not hold")
    for row in d["censoring_accounting"]:
        ok(close(row["contribution_sum"], row["dt_censored_mean"], 1e-9),
           f"{tag0}: the censoring decomposition at c={row['c']} {row['arm']}/"
           f"{row['colour']} does not sum to the censored mean")
        ok(close(row["frac_both_fired"] + row["frac_lead_only"]
                 + row["frac_lag_only"] + row["frac_neither_fired"], 1.0, 1e-9),
           f"{tag0}: firing-pattern fractions at c={row['c']} do not sum to 1")
    for row in d["detection_probability"]:
        ok(close(row["lead_detect"], row["lead_n_fired"] / row["n_reps"], 1e-12),
           f"{tag0}: lead_detect at c={row['c']} is not its own count over its own n")
        ok(close(row["lag_detect"], row["lag_n_fired"] / row["n_reps"], 1e-12),
           f"{tag0}: lag_detect at c={row['c']} is not its own count over its own n")
    for row in d["win_probability"]:
        if row["n_decided"]:
            ok(close(row["p_lead_first"], row["n_lead_first"] / row["n_decided"],
                     1e-12),
               f"{tag0}: p_lead_first at c={row['c']} is not n_lead_first/n_decided")
        ok(row["n_decided"] <= row["n_reps"],
           f"{tag0}: more decided runs than runs at c={row['c']}")
    mv = d["mnar_vs_mcar_verdict"]
    ok(mv["verdict"].startswith("DELETE")
       == (mv["gap_bootstrap"].get("straddles_zero", True)
           and not mv["achieved_edge_mechanism_supported"]),
       f"{tag0}: the MNAR verdict string is not what its own two tests give")
    ok(d["graph_object"]["is_cumulative"] is False
       and d["graph_object"]["tie_window_cycles"] > 0,
       f"{tag0}: the graph-object block does not record a finite tie window")


# ------------------------------------------------------------------------ both-class calibration
def check_bothclass_calibration() -> None:
    d = load(CAL_BOTHCLASS)
    if d is None:
        return
    tag = f"both-class calibration {CAL_BOTHCLASS}"
    des = d["design"]
    rows = d["per_run_rows"]
    stats = np.asarray([r["lead_sustained_score"] for r in rows], float)
    target = des["target_fpr"]
    ok(len(rows) == des["n_null_control_runs"],
       f"{tag}: declared n_null_control_runs != rows present")
    ok(len(rows) >= 600, f"{tag}: null smaller than the declared 600 runs")
    ok(all(r["seed"] >= SCORED_SEED_CEILING for r in rows),
       f"{tag}: a calibration run used a seed below the scored ceiling")
    big = load(CAL_NULL)
    if big is not None:
        want = {r["seed"] for r in big["per_run_rows"] if r["scen"] == "control"}
        ok({r["seed"] for r in rows} == want,
           f"{tag}: the both-class rule was not calibrated on the SAME null as large-null recalibration, so "
           f"the two rules are not comparable")
    for x in d["leading_fpr_curve"]:
        ok(close(float(np.mean(stats >= x["k"])), x["fpr"], 1e-12),
           f"{tag}: FPR at k={x['k']} not reproducible from the rows")
        ok(int(np.sum(stats >= x["k"])) == x["n_false_alarms"],
           f"{tag}: n_false_alarms at k={x['k']}")
    adm = [x["k"] for x in d["leading_fpr_curve"] if x["fpr"] <= target]
    ok(bool(adm) and close(min(adm), d["chosen"]["leading_k"], 1e-12),
       f"{tag}: chosen k {d['chosen']['leading_k']} is not the smallest k with null "
       f"FPR <= {target}")
    ok(d["chosen"]["leading_calibration_fpr"] <= target,
       f"{tag}: the chosen threshold does not meet its own target")
    disk = load_rule_file("detection_rule_bothclass.json")
    ok(disk is not None and json.dumps(disk, sort_keys=True)
       == json.dumps(d["bothclass_rule"], sort_keys=True),
       f"{tag}: configs/detection_rule_bothclass.json differs from the aggregate")
    fo = d["formation_only_rule"]
    bc = d["bothclass_rule"]
    for k, v in fo["leading"].items():
        if k in ("k", "classes"):
            continue
        ok(close(bc["leading"].get(k), v, 1e-12),
           f"{tag}: leading rule FORM changed at '{k}' -- only the class set and the "
           f"threshold may move")
    ok(json.dumps(bc["lagging"], sort_keys=True)
       == json.dumps(fo["lagging"], sort_keys=True),
       f"{tag}: the lagging comparator moved with the leading class set, so the "
       f"lead-time contrast is not a contrast")
    ok(list(bc["leading"]["classes"]) == list(des["classes"]),
       f"{tag}: the rule's class set is not the one the design declares")
    ok(close(bc["leading"]["k"], d["chosen"]["leading_k"], 1e-12),
       f"{tag}: rule file threshold != chosen threshold")


# ------------------------------------------------------------------------ aggregate verification
TR2 = "tranche2_runs.json"


def check_tranche2() -> None:
    """aggregate verification -- tranche-2 (detrending-window sensitivity ... non-monotone driver scenarios): seed hygiene, matched pairing, per-variant
    threshold selection reproducible from the null rows, and every summary table
    reproducible from its own per-run rows."""
    d = load(TR2)
    if d is None:
        return
    tag0 = f"aggregate verification {TR2}"
    des = d["design"]

    # --- seed hygiene ----------------------------------------------------------
    scored = ([r["seed"] for r in d["a20_detrend"]["per_run_rows"]]
              + [r["seed"] for r in d["a21_control_baseline"]["per_run_rows"]]
              + [r["seed"] for r in d["a22_bimodality"]["per_run_rows"]]
              + [r["seed"] for r in d["a24_per_run_rows"]]
              + [r["seed"] for r in d["a25_scenarios"]["per_run_rows"]])
    ok(all(s < SCORED_SEED_CEILING for s in scored),
       f"{tag0}: a scored run used a calibration seed (>= {SCORED_SEED_CEILING})")
    ok(all(r["seed"] >= SCORED_SEED_CEILING
           for r in d["a20_detrend"]["null_rows"]),
       f"{tag0}: an detrending-window sensitivity calibration run used a seed below the scored ceiling")

    # --- detrending-window sensitivity: thresholds are the frozen selection rule on the null rows --------
    import calibrate_recal as _cr
    for name, k in des["k_by_variant"].items():
        stats = np.asarray([r[name] for r in d["a20_detrend"]["null_rows"]], float)
        curve = _cr.fpr_curve(stats.tolist(), _cr.K_GRID)
        adm = [x["k"] for x in curve if x["fpr"] <= _cr.TARGET_FPR]
        ok(bool(adm) and close(min(adm), k, 1e-12),
           f"{tag0} detrending-window sensitivity {name}: threshold {k} is not the smallest admissible k on "
           f"the committed null rows")
    rows20 = d["a20_detrend"]["per_run_rows"]
    seeds_by = {}
    for r in rows20:
        seeds_by.setdefault((r["c"], r["scen"]), set()).add(r["seed"])
    bad = [c for c in des["c_grid"]
           if seeds_by.get((c, "degrade")) != seeds_by.get((c, "control"))]
    ok(not bad, f"{tag0} detrending-window sensitivity: degradation and control arms not matched seed-for-seed "
                f"at c={bad[:3]}")
    nc = rows20[0]["n_cycles"]
    for row in d["a20_detrend"]["table"]:
        deg = [r for r in rows20 if r["c"] == row["c"] and r["scen"] == "degrade"]
        ctl = [r for r in rows20 if r["c"] == row["c"] and r["scen"] == "control"]
        key = f"t_lead_{row['variant']}"
        dts = [(nc if r["t_lag"] is None else r["t_lag"])
               - (nc if r[key] is None else r[key]) for r in deg]
        tag = f"{tag0} detrending-window sensitivity {row['variant']} c={row['c']}"
        ok(close(float(np.mean(dts)), row["dt_censored_mean"], 1e-9),
           f"{tag}: dt_censored_mean not reproducible")
        ok(close(float(np.mean([r[key] is not None for r in deg])),
                 row["lead_detect"], 1e-12), f"{tag}: lead_detect")
        ok(close(float(np.mean([r[key] is not None for r in ctl])),
                 row["fpr_leading"], 1e-12), f"{tag}: fpr_leading")

    # --- control-baseline sensitivity table --------------------------------------------------------------
    rows21 = d["a21_control_baseline"]["per_run_rows"]
    for row in d["a21_control_baseline"]["table"]:
        rs = [r for r in rows21 if r["c"] == row["c"]]
        nc = rs[0]["n_cycles"]
        for field, key in (("dt_censored_ownbase", "t_lead_ownbase"),
                           ("dt_censored_ctlbase", "t_lead_ctlbase")):
            dts = [(nc if r["t_lag"] is None else r["t_lag"])
                   - (nc if r[key] is None else r[key]) for r in rs]
            ok(close(float(np.mean(dts)), row[field], 1e-9),
               f"{tag0} control-baseline sensitivity c={row['c']}: {field} not reproducible")

    # --- end-state mixture analysis occupancy + BIC ----------------------------------------------------
    import run_tranche2 as _t2
    g1 = load("gate1_bifurcation.json")
    if g1 is not None:
        for n in ("400", "800"):
            if n not in d["a22_bimodality"]["mixture_bic"]:
                continue
            pooled = [r["z_end"] for r in g1["long_run_rows"] if str(r["n"]) == n]
            fresh = _t2._mixture_bic(np.asarray(pooled))
            ok(close(fresh["delta_bic_1_minus_2"],
                     d["a22_bimodality"]["mixture_bic"][n]["delta_bic_1_minus_2"],
                     1e-6),
               f"{tag0} end-state mixture analysis: mixture BIC at N={n} not reproducible from the "
               f"committed Bifurcation-validation rows")

    # --- small-perturbation recovery table ---------------------------------------------------------------
    rows23 = d["a23_small_shock"]["per_run_rows"]
    for row in d["a23_small_shock"]["table"]:
        rs = [r for r in rows23 if r["a"] == row["a"] and r["frac"] == row["frac"]]
        ok(len(rs) == row["n_seeds"],
           f"{tag0} small-perturbation recovery a={row['a']} f={row['frac']}: n_seeds")
        ok(close(float(np.mean([r["returned"] for r in rs])),
                 row["returned_frac"], 1e-12),
           f"{tag0} small-perturbation recovery a={row['a']} f={row['frac']}: returned_frac")

    # --- record-revision decomposition: aggregates reproducible from rows (seeded bootstrap -> exact) ------
    from scoring import bootstrap_mean as _bm
    rows24 = d["a24_per_run_rows"]
    for leg in ("retroactive", "disclosed_equiv", "reconstruction_residual"):
        for stat in ("mean_abs", "level", "shape_sd", "alarm_flip_frac"):
            fresh = _bm([r[leg][stat] for r in rows24 if r[leg] is not None])
            ok(close(fresh["mean"],
                     d["a24_revocation_decomposition"][leg][stat]["mean"], 1e-9),
               f"{tag0} record-revision decomposition {leg}.{stat}: not reproducible from the rows")

    # --- non-monotone driver scenarios table -----------------------------------------------------------------
    rows25 = d["a25_scenarios"]["per_run_rows"]
    for row in d["a25_scenarios"]["table"]:
        rs = [r for r in rows25 if r["scen"] == row["scenario"]]
        tag = f"{tag0} non-monotone driver scenarios {row['scenario']}"
        ok(len(rs) == row["n_reps"], f"{tag}: n_reps")
        ok(close(float(np.mean([r["collapsed"] for r in rs])),
                 row["collapsed_frac"], 1e-12), f"{tag}: collapsed_frac")
        ok(close(float(np.mean([r["lead_fired_recal"] for r in rs])),
                 row["lead_alarm_rate_recal"], 1e-12),
           f"{tag}: lead_alarm_rate_recal")
        ok(close(float(np.mean([r["lag_fired_recal"] for r in rs])),
                 row["lag_alarm_rate"], 1e-12), f"{tag}: lag_alarm_rate")
    # dip_recover and oscillate must not have collapsed -- if they did, the "false
    # positive" framing of the scenario leg is wrong and the report must change
    for scen in ("dip_recover", "oscillate"):
        row = next((x for x in d["a25_scenarios"]["table"]
                    if x["scenario"] == scen), None)
        if row is not None:
            ok(row["collapsed_frac"] == 0.0,
               f"{tag0} non-monotone driver scenarios {scen}: trajectory collapsed in {row['collapsed_frac']:.0%}"
               f" of runs -- alarms there are not false positives and the scenario "
               f"does not test what non-monotone driver scenarios declares")


# ------------------------------------------------------------------------- claim binding
def check_claims() -> None:
    p = os.path.join(HERE, "claims.json")
    if not ok(os.path.exists(p), "claim binding claims.json missing -- RESULTS.md numbers are "
                                 "not bound to the aggregates"):
        return
    with open(p) as f:
        claims = json.load(f)
    rp = os.path.join(HERE, claims.get("results_file", "RESULTS.md"))
    if not ok(os.path.exists(rp), f"claim binding {claims.get('results_file')} missing"):
        return
    with open(rp, encoding="utf-8") as f:
        results_text = f.read()
    for cl in claims["claims"]:
        d = load(cl["file"])
        if not ok(d is not None, f"claim binding {cl['id']}: aggregate {cl['file']} missing"):
            continue
        try:
            got = resolve(d, cl["path"])
        except Exception as e:                            # noqa: BLE001
            ok(False, f"claim binding {cl['id']}: path {cl['path']} did not resolve ({e})")
            continue
        ok(close(got, cl["value"], cl.get("tol", 1e-6)),
           f"claim binding {cl['id']}: claim says {cl['value']}, aggregate says {got} "
           f"({cl['file']}:{cl['path']})")
        if cl.get("text"):
            ok(cl["text"] in results_text,
               f"claim-text binding {cl['id']}: bound text not present verbatim in "
               f"{claims.get('results_file', 'RESULTS.md')}: {cl['text']!r}")


# ------------------------------------------------------------------------- determinism
def check_determinism(sample: int) -> None:
    from model import Config, run_config
    from scoring import score_run
    from runner import load_rule
    for name in EA_FILES:
        d = load(name)
        if d is None:
            continue
        rule = load_rule((d.get("design") or {}).get("rule_file",
                                                    "detection_rule.json"))
        rows = d["per_run_rows"]
        rng = np.random.Generator(np.random.PCG64(4242))
        pick = rng.choice(len(rows), size=min(sample, len(rows)), replace=False)
        for i in pick:
            r = rows[int(i)]
            cfg = Config(consent_c=r["c"], consent_arm=r["arm"],
                         eps_colour=r["colour"], eps_severity=r["sev"],
                         scenario=r["scen"], seed=r["seed"])
            s = run_config(cfg).series()
            sc = score_run(s, rule, cfg.n_cycles)
            tag = (f"determinism {name} seed={r['seed']} c={r['c']} {r['arm']} "
                   f"{r['colour']}/{r['sev']} {r['scen']}")
            ok(sc["t_lead"] == r["t_lead"],
               f"{tag}: t_lead {sc['t_lead']} != {r['t_lead']}")
            ok(sc["t_lag"] == r["t_lag"], f"{tag}: t_lag {sc['t_lag']} != {r['t_lag']}")
            ok(close(float(s["achieved_edge_frac"][-1]), r["achieved_edge_frac"],
                     1e-12), f"{tag}: achieved_edge_frac")

    # the remediation sweep is deterministic too, and its runs go through a Simulation
    # SUBCLASS (the cohesion tap), so a sample is re-executed through the same code path
    # and must reproduce its committed scalars under BOTH scorings.
    d = load(REM_RUNS)
    if d is not None:
        import run_remediation as rr
        overrides = {v["variant"]: v["overrides"] for v in d["design"]["variants"]}
        rows = d["per_run_rows"]
        rng = np.random.Generator(np.random.PCG64(9191))
        pick = rng.choice(len(rows), size=min(max(2, sample // 2), len(rows)),
                          replace=False)
        for i in pick:
            r = rows[int(i)]
            got = rr._job((r["variant"], overrides[r["variant"]], r["c"], r["scen"],
                           r["rep"]))
            tag = (f"determinism {REM_RUNS} {r['variant']} c={r['c']} {r['scen']} "
                   f"rep={r['rep']} seed={r['seed']}")
            for key in ("seed", "t_lead", "t_lag", "bc_t_lead", "bc_t_lag"):
                ok(got[key] == r[key], f"{tag}: {key} {got[key]} != {r[key]}")
            for key in ("achieved_edge_frac", "lag_sigma_drop_sustained",
                        "lead_max_score", "coh_mean_a048"):
                ok(close(got[key], r[key], 1e-12), f"{tag}: {key}")

    # the calibration null is deterministic too: re-run one null and one calibration
    # degradation seed and reproduce the committed sufficient statistics exactly.
    d = load(CAL_NULL)
    if d is None:
        return
    import calibrate_recal as cr
    for scen in ("control", "degrade"):
        row = next((r for r in d["per_run_rows"] if r["scen"] == scen), None)
        if row is None:
            continue
        got = cr._job((scen, row["seed"]))
        tag = f"determinism {CAL_NULL} {scen} seed={row['seed']}"
        ok(close(got["lead_sustained_score"], row["lead_sustained_score"], 1e-12),
           f"{tag}: lead_sustained_score")
        for su in d["design"]["sustain_grid"]:
            ok(close(got[f"lag_sustained_drop_s{su}"],
                     row[f"lag_sustained_drop_s{su}"], 1e-12),
               f"{tag}: lag_sustained_drop_s{su}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, default=6)
    ap.add_argument("--skip-determinism", action="store_true")
    args = ap.parse_args()

    check_seeds_and_rule()
    check_gate1_recompute()
    check_calibration_null()
    for f in EA_FILES:
        check_ea_recompute(f)
        check_iso_line(f)
    for f in AUX_FILES:
        check_aux_recompute(f)
    for f in EB_FILES:
        check_eb_recompute(f)
    check_recal_analysis()
    check_bothclass_calibration()
    check_remediation_runs()
    check_remediation_analysis()
    check_tranche2()
    check_claims()
    if not args.skip_determinism:
        check_determinism(args.sample)

    print(f"verify_reported_numbers: {CHECKS[0]} checks, {len(FAILURES)} failures")
    for f in FAILURES[:40]:
        print("  FAIL", f)
    if len(FAILURES) > 40:
        print(f"  ... and {len(FAILURES) - 40} more")
    sys.exit(1 if FAILURES else 0)


if __name__ == "__main__":
    main()
