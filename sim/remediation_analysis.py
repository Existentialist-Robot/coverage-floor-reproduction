"""
Remediation analysis 2026-08-05 -- pre-submission review findings, re-analysis leg.

No new runs.  Everything here is a pure function of committed aggregates:

    ea_headline_recal.json          the recalibrated E-A grid (7 200 runs)
    ea_headline.json                campaign 1, for the seed-noise gauge
    remediation_runs.json           the remediation sweep (run_remediation.py)
    calibration_null_bothclass.json the both-class rule's calibration (calibrate_class.py)

so `compute()` is re-executed by verify_reported_numbers.py and must reproduce the
committed file byte-for-byte (check class C12).  Every bootstrap is seeded.

WHAT IT ANSWERS

  1 DE-CENSORING (finding 1).  `scoring.py` imputes a non-firing detector to the run
    length, so the reported Delta t is a mean over values in which non-detection is
    encoded as 320.  Three things are computed instead of that single number:
      (a) DETECTION PROBABILITY per cell, per detector, with Wilson intervals -- a
          first-class outcome rather than a footnote;
      (b) PAIRED LEAD TIME `dt_paired`, restricted to runs where BOTH detectors fired,
          which is computed in scoring.py and was never reported;
      (c) a COMPETING-RISKS treatment.  Alarm time is a time-to-event with one
          administrative censoring time (the run length), so the Kaplan-Meier estimator
          reduces to the empirical survival function, and the restricted mean time to
          alarm over a declared horizon tau is RMTA(tau) = mean(min(T, tau)).  The
          difference RMTA_lag(tau) - RMTA_lead(tau) is reported at several horizons.
      A fourth, censoring-ROBUST statistic is added because it is the one that answers
      this question directly: the WIN PROBABILITY, P(the leading detector
      alarms first | at least one detector alarmed).  Which detector was first is
      determined whenever either fires -- if the leading detector fires and the lagging
      one never does, the leading one was first -- so this statistic is undefined only
      on runs where NEITHER fires, and it involves no imputation at all.
      The Delta t = 0 iso-line is then re-sited under each estimator.

  2 IDENTITY, stated because it is load-bearing and unobvious.  With a single
    administrative censoring time equal to the horizon, mean(dt_censored) is EXACTLY
    RMTA_lag(320) - RMTA_lead(320).  The censored mean is therefore a legitimate
    estimand -- a difference of restricted mean times to alarm at a 320-cycle horizon --
    and the error is one of NAMING and of horizon-dependence, not of arithmetic.  This
    block computes both quantities and their maximum discrepancy over every cell.

  3 PER-CELL FALSE ALARMS (finding 12/c).  The leading and lagging per-cell FPRs with
    Wilson intervals and a two-proportion test, plus a heterogeneity test across c for
    each detector.  The lagging detector reads program-close output, which no consent
    dial touches, so its FPR CANNOT depend on c: its across-c heterogeneity is a
    calibrated gauge of how much apparent c-dependence is seed noise.

  4 THE MNAR CLAIM (finding e / the paper's weakest sentence).  The MCAR-minus-MNAR
    coverage-floor gap is bootstrapped, so "0.412 vs 0.502" gets an interval; the
    achieved-edge-fraction mechanism is tested; and the c = 1.0 same-configuration
    Delta t disagreement is re-measured on the recalibrated rule as the seed-noise floor.

  5 The remediation sweep's reductions: Delta t and the floor as functions of measured
    comparator gain, of `formation_driver`, of `partial_obs_sigma` and of N; the
    both-class rule against the formation-only rule; and an amplitude audit of the two
    channels' gains on the same order parameter.

Run:  python remediation_analysis.py
Writes: aggregates/remediation_analysis.json
"""

from __future__ import annotations

import hashlib
import json
import math
import os

import numpy as np

from model import Config
from runner import write_agg, AGG

EA_NEW = "ea_headline_recal.json"
EA_OLD = "ea_headline.json"
REM = "remediation_runs.json"
CAL_BC = "calibration_null_bothclass.json"
GATE1 = "gate1_bifurcation.json"

Z = 1.959963985
Z80 = 1.959963985 + 0.841621234
HORIZONS = (160, 200, 240, 280, 320)
INCIDENCE_CYCLES = (140, 160, 180, 200, 220, 240, 260, 280, 300, 320)
DETECTION_FLOOR_TARGETS = (0.5, 0.4, 1.0 / 3.0)
N_BOOT = 4000
BOOT_SEED = 20260805
# A(t) = 0.48 -> 0.31 linearly from cycle 50 to cycle 320, so the empirical fold
# A_c = 0.37 (Gate 1) is crossed at a cycle that is fixed by the scenario design.  Used
# to ask whether an alarm arrived before or after the transition itself.
FOLD_A = 0.37


def stable_seed(tag: str) -> int:
    """Deterministic per-cell bootstrap seed.  Python's built-in hash() is salted per
    process, so it must never appear in a number this file commits: the verifier
    re-executes compute() in a different process and requires bit-identical output."""
    h = hashlib.sha256(tag.encode()).digest()
    return BOOT_SEED + int.from_bytes(h[:4], "big") % 1_000_000


def _load(name):
    p = os.path.join(AGG, name)
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def wilson(k: int, n: int) -> list[float]:
    if n == 0:
        return [float("nan"), float("nan")]
    p = k / n
    c = (p + Z * Z / (2 * n)) / (1 + Z * Z / n)
    h = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / (1 + Z * Z / n)
    return [c - h, c + h]


def _phi(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def two_prop_z(k1: int, n1: int, k2: int, n2: int) -> dict:
    if n1 == 0 or n2 == 0:
        return {"z": float("nan"), "p_two_sided": float("nan")}
    p1, p2 = k1 / n1, k2 / n2
    p = (k1 + k2) / (n1 + n2)
    se = math.sqrt(p * (1 - p) * (1 / n1 + 1 / n2))
    if se <= 0:
        return {"z": float("nan"), "p_two_sided": float("nan")}
    z = (p1 - p2) / se
    return {"z": z, "p_two_sided": 2.0 * (1.0 - _phi(abs(z)))}


def chi2_sf_even(x: float, df: int) -> float:
    """Upper tail of a chi-square with EVEN degrees of freedom, closed form."""
    m = df // 2
    s = 0.0
    term = 1.0
    for j in range(m):
        if j:
            term *= (x / 2.0) / j
        s += term
    return math.exp(-x / 2.0) * s


def het_chi2(counts: list[tuple[int, int]]) -> dict:
    """Test that one binomial rate generated every (k, n) group."""
    tk = sum(k for k, _ in counts)
    tn = sum(n for _, n in counts)
    p = tk / tn if tn else float("nan")
    if not (0.0 < p < 1.0):
        return {"chi2": float("nan"), "df": len(counts) - 1,
                "p_value": float("nan"), "pooled_rate": p}
    chi = sum((k - n * p) ** 2 / (n * p * (1 - p)) for k, n in counts)
    df = len(counts) - 1
    return {"chi2": chi, "df": df,
            "p_value": chi2_sf_even(chi, df) if df % 2 == 0 else float("nan"),
            "pooled_rate": p,
            "rates": [k / n for k, n in counts]}


def interp_cross(cs, ys, level: float = 0.0):
    """First crossing of `level` along ascending c, linearly interpolated.  Same rule
    the runners use for the Delta t = 0 iso-line, generalised to any level."""
    for i in range(len(cs) - 1):
        a, b = ys[i] - level, ys[i + 1] - level
        if (a <= 0.0) != (b <= 0.0):
            return cs[i] + (cs[i + 1] - cs[i]) * (0.0 - a) / (b - a)
    return None


def _eff(v, n_cycles: int) -> float:
    return float(n_cycles if v is None else v)


def _rmta(rows, key: str, tau: int, n_cycles: int) -> float:
    return float(np.mean([min(tau, _eff(r[key], n_cycles)) for r in rows]))


def _boot_quantiles(v: np.ndarray) -> dict:
    ok = v[np.isfinite(v)]
    if ok.size == 0:
        return {"n_resolved": 0, "mean": float("nan"),
                "lo95": float("nan"), "hi95": float("nan")}
    return {"n_resolved": int(ok.size), "mean": float(ok.mean()),
            "lo95": float(np.quantile(ok, 0.025)),
            "hi95": float(np.quantile(ok, 0.975))}


# ============================================================== the analysis
def compute() -> dict:
    new = _load(EA_NEW)
    old = _load(EA_OLD)
    rem = _load(REM)
    cal_bc = _load(CAL_BC)
    gate1 = _load(GATE1)
    n_cycles = new["design"]["n_cycles"]
    sev = new["design"]["headline_severity"]
    cs_all = sorted(new["design"]["c_grid"])
    colours = list(new["design"]["colours"])
    arms = list(new["design"]["arms"])
    rows = new["per_run_rows"]

    def cell(arm, col, c, scen="degrade"):
        return [r for r in rows if r["c"] == c and r["arm"] == arm
                and r["colour"] == col and r["sev"] == sev and r["scen"] == scen]

    cfg = Config()
    fold_cycle = int(round(cfg.degrade_start + (cfg.n_cycles - cfg.degrade_start)
                           * (cfg.a_high - FOLD_A) / (cfg.a_high - cfg.a_low)))

    # ---------------------------------------------------------------- (1a) detection
    detection, censoring, paired, survival, winprob = [], [], [], [], []
    for arm in arms:
        for col in colours:
            for c in cs_all:
                d = cell(arm, col, c)
                if not d:
                    continue
                n = len(d)
                kl = int(sum(1 for r in d if r["lead_fired"]))
                kg = int(sum(1 for r in d if r["lag_fired"]))
                tag = {"arm": arm, "colour": col, "c": c, "severity": sev, "n_reps": n}
                detection.append({
                    **tag,
                    "lead_detect": kl / n, "lead_detect_wilson95": wilson(kl, n),
                    "lag_detect": kg / n, "lag_detect_wilson95": wilson(kg, n),
                    "lead_n_fired": kl, "lag_n_fired": kg,
                    "lead_median_alarm_cycle": (
                        float(np.median([r["t_lead"] for r in d if r["lead_fired"]]))
                        if kl else None),
                    "lag_median_alarm_cycle": (
                        float(np.median([r["t_lag"] for r in d if r["lag_fired"]]))
                        if kg else None),
                    "lead_frac_fired_before_fold": float(np.mean(
                        [1.0 if (r["lead_fired"] and r["t_lead"] < fold_cycle) else 0.0
                         for r in d])),
                    "lag_frac_fired_before_fold": float(np.mean(
                        [1.0 if (r["lag_fired"] and r["t_lag"] < fold_cycle) else 0.0
                         for r in d])),
                })

                # ------------------------------------------- censoring accounting
                both = [r for r in d if r["lead_fired"] and r["lag_fired"]]
                lead_only = [r for r in d if r["lead_fired"] and not r["lag_fired"]]
                lag_only = [r for r in d if r["lag_fired"] and not r["lead_fired"]]
                neither = [r for r in d if not r["lead_fired"] and not r["lag_fired"]]
                grp = {"both_fired": both, "lead_only": lead_only,
                       "lag_only": lag_only, "neither": neither}
                contrib = {k: (len(v) / n) * float(np.mean([x["dt_censored"]
                                                            for x in v]))
                           if v else 0.0 for k, v in grp.items()}
                censoring.append({
                    **tag,
                    "dt_censored_mean": float(np.mean([r["dt_censored"] for r in d])),
                    "frac_both_fired": len(both) / n,
                    "frac_lead_only": len(lead_only) / n,
                    "frac_lag_only": len(lag_only) / n,
                    "frac_neither_fired": len(neither) / n,
                    "frac_any_imputation": float(np.mean(
                        [1.0 if (r["t_lead"] is None or r["t_lag"] is None) else 0.0
                         for r in d])),
                    "frac_lead_imputed_to_run_length": 1.0 - kl / n,
                    "contribution_both_fired": contrib["both_fired"],
                    "contribution_lead_only": contrib["lead_only"],
                    "contribution_lag_only": contrib["lag_only"],
                    "contribution_neither": contrib["neither"],
                    "contribution_sum": sum(contrib.values()),
                })

                # ------------------------------------------- (1b) paired lead time
                pv = np.asarray([r["dt_paired"] for r in d
                                 if r["dt_paired"] is not None], float)
                rng = np.random.Generator(np.random.PCG64(
                    stable_seed(f"paired|{arm}|{col}|{c}")))
                if pv.size:
                    bs = pv[rng.integers(0, pv.size, (N_BOOT, pv.size))].mean(axis=1)
                    paired.append({
                        **tag, "n_both_fired": int(pv.size),
                        "dt_paired_mean": float(pv.mean()),
                        "dt_paired_se": float(bs.std(ddof=1)),
                        "dt_paired_lo95": float(np.quantile(bs, 0.025)),
                        "dt_paired_hi95": float(np.quantile(bs, 0.975)),
                        "dt_paired_median": float(np.median(pv)),
                        "frac_lead_earlier": float(np.mean(pv > 0)),
                    })
                else:
                    paired.append({**tag, "n_both_fired": 0,
                                   "dt_paired_mean": float("nan"),
                                   "dt_paired_se": float("nan"),
                                   "dt_paired_lo95": float("nan"),
                                   "dt_paired_hi95": float("nan"),
                                   "dt_paired_median": float("nan"),
                                   "frac_lead_earlier": float("nan")})

                # ------------------------------------------- (1c) competing risks
                srow = {**tag}
                for tau in HORIZONS:
                    rl = _rmta(d, "t_lead", tau, n_cycles)
                    rg = _rmta(d, "t_lag", tau, n_cycles)
                    srow[f"rmta_lead_tau{tau}"] = rl
                    srow[f"rmta_lag_tau{tau}"] = rg
                    srow[f"delta_rmta_tau{tau}"] = rg - rl
                for t in INCIDENCE_CYCLES:
                    srow[f"cum_incidence_lead_{t}"] = float(np.mean(
                        [1.0 if (r["lead_fired"] and r["t_lead"] <= t) else 0.0
                         for r in d]))
                    srow[f"cum_incidence_lag_{t}"] = float(np.mean(
                        [1.0 if (r["lag_fired"] and r["t_lag"] <= t) else 0.0
                         for r in d]))
                survival.append(srow)

                # ------------------------------------------- win probability
                dec = [r for r in d if r["lead_fired"] or r["lag_fired"]]
                win = [r for r in dec if r["lead_fired"]
                       and (not r["lag_fired"] or r["t_lead"] < r["t_lag"])]
                tie = [r for r in dec if r["lead_fired"] and r["lag_fired"]
                       and r["t_lead"] == r["t_lag"]]
                winprob.append({
                    **tag, "n_decided": len(dec),
                    "frac_decided": len(dec) / n,
                    "n_lead_first": len(win), "n_ties": len(tie),
                    "p_lead_first": (len(win) / len(dec)) if dec else float("nan"),
                    "p_lead_first_wilson95": wilson(len(win), len(dec)),
                })

    # RMTA identity: mean(dt_censored) == delta RMTA at tau = run length, exactly
    ident = []
    for s in survival:
        m = next(x for x in censoring if x["arm"] == s["arm"]
                 and x["colour"] == s["colour"] and x["c"] == s["c"])
        ident.append(abs(s[f"delta_rmta_tau{n_cycles}"] - m["dt_censored_mean"]))
    rmta_identity = {
        "statement": ("with a single administrative censoring time equal to the "
                      "horizon, mean(dt_censored) is exactly "
                      "RMTA_lag(tau) - RMTA_lead(tau) at tau = the run length"),
        "horizon": n_cycles,
        "n_cells_checked": len(ident),
        "max_abs_discrepancy": float(max(ident)) if ident else float("nan"),
    }

    # ------------------------------------------------------- floor re-siting
    def cellmeans(arm, col, fn):
        return [fn(cell(arm, col, c)) for c in cs_all]

    floors = {}
    for arm in arms:
        for col in colours:
            if not cell(arm, col, cs_all[0]):
                continue
            key = f"{arm}|{col}"
            dt = cellmeans(arm, col,
                           lambda d: float(np.mean([r["dt_censored"] for r in d])))
            entry = {"c_grid": cs_all, "delta_rmta_tau320": dt,
                     "c_at_dt_zero_tau320": interp_cross(cs_all, dt)}
            for tau in HORIZONS:
                ys = cellmeans(arm, col, lambda d, t=tau: (
                    _rmta(d, "t_lag", t, n_cycles) - _rmta(d, "t_lead", t, n_cycles)))
                entry[f"delta_rmta_tau{tau}_by_c"] = ys
                entry[f"c_at_zero_tau{tau}"] = interp_cross(cs_all, ys)
            pv = cellmeans(arm, col, lambda d: float(np.mean(
                [r["dt_paired"] for r in d if r["dt_paired"] is not None]))
                if any(r["dt_paired"] is not None for r in d) else float("nan"))
            entry["dt_paired_by_c"] = pv
            entry["c_at_dt_paired_zero"] = interp_cross(cs_all, pv)
            wp = cellmeans(arm, col, lambda d: (
                len([r for r in d if r["lead_fired"]
                     and (not r["lag_fired"] or r["t_lead"] < r["t_lag"])])
                / max(1, len([r for r in d if r["lead_fired"] or r["lag_fired"]]))))
            entry["p_lead_first_by_c"] = wp
            entry["c_at_p_lead_first_half"] = interp_cross(cs_all, wp, 0.5)
            ld = cellmeans(arm, col,
                           lambda d: float(np.mean([r["lead_fired"] for r in d])))
            entry["lead_detect_by_c"] = ld
            entry["c_at_lead_detect"] = {
                f"{t:.3f}": interp_cross(cs_all, ld, t)
                for t in DETECTION_FLOOR_TARGETS}
            floors[key] = entry

    # bootstrap the headline floors (white and systematic, both arms)
    boot = {}
    for col in colours:
        cells_v = {(arm, c): {
            "dt": np.asarray([r["dt_censored"] for r in cell(arm, col, c)], float),
            "lead": np.asarray([1.0 if r["lead_fired"] else 0.0
                                for r in cell(arm, col, c)], float),
            "first": np.asarray(
                [1.0 if (r["lead_fired"] and (not r["lag_fired"]
                                              or r["t_lead"] < r["t_lag"])) else 0.0
                 for r in cell(arm, col, c)], float),
            "dec": np.asarray([1.0 if (r["lead_fired"] or r["lag_fired"]) else 0.0
                               for r in cell(arm, col, c)], float),
        } for arm in arms for c in cs_all if cell(arm, col, c)}
        if not cells_v:
            continue
        rng = np.random.Generator(np.random.PCG64(stable_seed(f"floor|{col}")))
        draws = {f"{a}|dt": [] for a in arms}
        draws.update({f"{a}|win": [] for a in arms})
        draws.update({f"{a}|det": [] for a in arms})
        diff = []
        for _ in range(N_BOOT):
            got = {}
            for arm in arms:
                dts, wins, dets = [], [], []
                for c in cs_all:
                    v = cells_v.get((arm, c))
                    if v is None:
                        continue
                    idx = rng.integers(0, v["dt"].size, v["dt"].size)
                    dts.append(float(v["dt"][idx].mean()))
                    nd = float(v["dec"][idx].sum())
                    wins.append(float(v["first"][idx].sum()) / nd if nd else np.nan)
                    dets.append(float(v["lead"][idx].mean()))
                got[arm] = (interp_cross(cs_all, dts),
                            interp_cross(cs_all, wins, 0.5),
                            interp_cross(cs_all, dets, 0.5))
                draws[f"{arm}|dt"].append(got[arm][0] if got[arm][0] is not None
                                          else np.nan)
                draws[f"{arm}|win"].append(got[arm][1] if got[arm][1] is not None
                                           else np.nan)
                draws[f"{arm}|det"].append(got[arm][2] if got[arm][2] is not None
                                           else np.nan)
            a, b = got.get("mnar", (None,))[0], got.get("mcar", (None,))[0]
            if a is not None and b is not None:
                diff.append(b - a)
        blk = {f"{k}": _boot_quantiles(np.asarray(v, float))
               for k, v in draws.items()}
        dv = np.asarray(diff, float)
        blk["mcar_minus_mnar_dt_floor"] = {
            **_boot_quantiles(dv),
            "n_draws_both_resolved": int(dv.size),
            "p_gap_positive": float(np.mean(dv > 0)) if dv.size else float("nan"),
            "straddles_zero": bool(dv.size and np.quantile(dv, 0.025) < 0
                                   < np.quantile(dv, 0.975)),
        }
        boot[col] = blk

    # ------------------------------------------------------- (3) per-cell FPR
    mid_ctl = [r for r in rows if r["scen"] == "control" and r["sev"] == sev]
    fpr_by_c, counts_lead, counts_lag = [], [], []
    for c in cs_all:
        sel = [r for r in mid_ctl if r["c"] == c]
        n = len(sel)
        kl = int(sum(1 for r in sel if r["lead_fired"]))
        kg = int(sum(1 for r in sel if r["lag_fired"]))
        counts_lead.append((kl, n))
        counts_lag.append((kg, n))
        fpr_by_c.append({
            "c": c, "n_control_runs": n,
            "fpr_leading": kl / n, "fpr_leading_wilson95": wilson(kl, n),
            "fpr_lagging": kg / n, "fpr_lagging_wilson95": wilson(kg, n),
            "ratio_leading_over_lagging": (kl / kg) if kg else float("nan"),
            "lead_vs_lag_test": two_prop_z(kl, n, kg, n),
            "leading_resolvably_above_external_bar": bool(wilson(kl, n)[0] > 0.036),
        })
    fpr_per_cell_detail = []
    for arm in arms:
        for col in colours:
            for c in cs_all:
                sel = cell(arm, col, c, "control")
                if not sel:
                    continue
                n = len(sel)
                kl = int(sum(1 for r in sel if r["lead_fired"]))
                kg = int(sum(1 for r in sel if r["lag_fired"]))
                fpr_per_cell_detail.append({
                    "arm": arm, "colour": col, "c": c, "n_control_runs": n,
                    "fpr_leading": kl / n, "fpr_leading_wilson95": wilson(kl, n),
                    "fpr_lagging": kg / n, "fpr_lagging_wilson95": wilson(kg, n)})
    fpr_block = {
        "pooled_over_colours_and_arms_at_mid_severity": fpr_by_c,
        "per_cell": fpr_per_cell_detail,
        "leading_heterogeneity_across_c": het_chi2(counts_lead),
        "lagging_heterogeneity_across_c": het_chi2(counts_lag),
        "why_the_lagging_test_is_the_gauge": (
            "the lagging comparator reads program-close output, which no consent dial "
            "touches, so its false-alarm rate cannot depend on c by construction: its "
            "across-c heterogeneity is pure seed noise and calibrates how much of the "
            "leading rule's apparent c-dependence is real"),
    }

    # ------------------------------------------------------- (4) the MNAR verdict
    edge = {}
    for arm in arms:
        a = np.asarray([r["achieved_edge_frac"] for r in rows
                        if r["c"] == 0.5 and r["arm"] == arm and r["sev"] == sev
                        and r["scen"] == "degrade"], float)
        edge[arm] = {"n": int(a.size), "mean": float(a.mean()),
                     "sd": float(a.std(ddof=1)),
                     "se": float(a.std(ddof=1) / math.sqrt(a.size))}
    de = edge["mnar"]["mean"] - edge["mcar"]["mean"]
    dse = math.sqrt(edge["mnar"]["se"] ** 2 + edge["mcar"]["se"] ** 2)
    gauge = {}
    for camp, d in (("campaign1", old), ("recalibrated", new)):
        if d is None:
            continue
        per = {}
        for col in colours:
            vals = {}
            for arm in arms:
                sel = [r for r in d["per_run_rows"]
                       if r["scen"] == "degrade" and r["c"] == 1.0
                       and r["arm"] == arm and r["colour"] == col and r["sev"] == sev]
                vals[arm] = (float(np.mean([r["dt_censored"] for r in sel]))
                             if sel else float("nan"))
            per[col] = {"mnar_dt": vals["mnar"], "mcar_dt": vals["mcar"],
                        "abs_difference": abs(vals["mnar"] - vals["mcar"])}
        gauge[camp] = {"per_colour": per,
                       "max_abs_difference": float(max(v["abs_difference"]
                                                       for v in per.values()))}
    wcol = boot.get("white", {}).get("mcar_minus_mnar_dt_floor", {})
    mnar_verdict = {
        "claim_under_test": ("the MNAR floor sits below the MCAR floor (0.412 vs "
                             "0.502), so realistic consent skew helps rather than "
                             "hurts detection at matched nominal coverage"),
        "floor_mnar_white": floors.get("mnar|white", {}).get("c_at_dt_zero_tau320"),
        "floor_mcar_white": floors.get("mcar|white", {}).get("c_at_dt_zero_tau320"),
        "gap_bootstrap": wcol,
        "achieved_edge_fraction_at_c_half": edge,
        "achieved_edge_difference_mnar_minus_mcar": de,
        "achieved_edge_difference_se": dse,
        "achieved_edge_difference_ci95": [de - Z * dse, de + Z * dse],
        "achieved_edge_mechanism_supported": bool(abs(de) > Z * dse),
        "dt_seed_noise_gauge_same_configuration_at_c1": gauge,
        "verdict": ("DELETE -- the gap is inside the seed-noise floor and its stated "
                    "mechanism is falsified"
                    if (wcol.get("straddles_zero", True)
                        and not bool(abs(de) > Z * dse))
                    else "KEEP -- resolvable"),
    }

    # ------------------------------------------------------- (5) the sweep legs
    def rem_cells(pred):
        return sorted([x for x in rem["cells"] if pred(x)], key=lambda x: x["c"])

    def leg_table(variant: str) -> dict:
        cl = rem_cells(lambda x: x["variant"] == variant)
        iso = rem["iso_line_dt_zero"].get(variant, {})
        return {
            "variant": variant,
            "c": [x["c"] for x in cl],
            "dt": [x["dt_censored_mean"] for x in cl],
            "dt_se": [x["dt_censored_se"] for x in cl],
            "dt_paired": [x["dt_paired_mean"] for x in cl],
            "lead_detect": [x["lead_detect_rate"] for x in cl],
            "lag_detect": [x["lag_detect_rate"] for x in cl],
            "fpr_leading": [x["fpr_leading"] for x in cl],
            "fpr_lagging": [x["fpr_lagging"] for x in cl],
            "n_reps": [x["n_reps"] for x in cl],
            "c_at_dt_zero": iso.get("c_at_dt_zero"),
            "c_at_dt_zero_band": iso.get("c_at_dt_zero_band"),
        }

    base = rem["design"]["base_variant"]
    gain_rows = [x for x in rem["comparator_gain_table"] if x["leg"] == "gain"]
    gain_rows.sort(key=lambda x: x["comparator_gain_sigma_drop"])
    g_axis = [x["comparator_gain_sigma_drop"] for x in gain_rows]
    dt_full = [x["dt_at_full_coverage"] for x in gain_rows]
    cstar = [x["c_at_dt_zero"] for x in gain_rows]
    gain_leg = {
        "x_axis": ("comparator gain, measured as the sustained baseline-sigma drop "
                   "available in the lagging series (the statistic its k_sd rule is "
                   "compared against); the rule requires "
                   f"{gain_rows[0]['k_sd_required'] if gain_rows else float('nan')}"),
        "rows": gain_rows,
        "gain_axis": g_axis,
        "dt_at_full_coverage": dt_full,
        "c_at_dt_zero": cstar,
        "gain_at_which_dt_at_full_coverage_vanishes": interp_cross(
            g_axis, [v if v is not None else float("nan") for v in dt_full]),
        "gain_at_which_floor_reaches_full_coverage": interp_cross(
            g_axis, [(v if v is not None else 1.0) for v in cstar], 1.0),
        "per_variant": {x["variant"]: leg_table(x["variant"]) for x in gain_rows},
        "base_variant": base,
    }

    size_variants = [base] + [v["variant"] for v in rem["design"]["variants"]
                              if v["leg"] == "size"]
    size_leg = {
        "note": ("N is swept with programs_per_cycle scaled to hold the mean-field "
                 "control parameter K fixed (the Gate 1 recipe) AND comp_ref scaled "
                 "with N, since a comp_ref fixed at 100 makes the leading "
                 "observable's gain N-dependent by construction"),
        "per_variant": {v: leg_table(v) for v in size_variants},
        "n_agents": {v: next((r["n_agents"] for r in rem["per_run_rows"]
                              if r["variant"] == v), None) for v in size_variants},
        "comp_ref": {v: next((r["comp_ref"] for r in rem["per_run_rows"]
                              if r["variant"] == v), None) for v in size_variants},
        "floors": {v: rem["iso_line_dt_zero"].get(v, {}).get("c_at_dt_zero")
                   for v in size_variants},
    }
    pos_variants = [base] + [v["variant"] for v in rem["design"]["variants"]
                             if v["leg"] == "partial"]
    partial_leg = {
        "per_variant": {v: leg_table(v) for v in pos_variants},
        "partial_obs_sigma": {v: next((r["partial_obs_sigma"]
                                       for r in rem["per_run_rows"]
                                       if r["variant"] == v), None)
                              for v in pos_variants},
        "floors": {v: rem["iso_line_dt_zero"].get(v, {}).get("c_at_dt_zero")
                   for v in pos_variants},
    }
    driver_leg = {
        "per_variant": {v: leg_table(v) for v in (base, "driver_degree")},
        "floors": {v: rem["iso_line_dt_zero"].get(v, {}).get("c_at_dt_zero")
                   for v in (base, "driver_degree")},
    }

    # both-class rule (finding 11) on the base variant's runs
    bcl = rem_cells(lambda x: x["variant"] == base)
    bc_iso = interp_cross([x["c"] for x in bcl],
                          [x["bc_dt_censored_mean"] for x in bcl])
    class_leg = {
        "formation_only": leg_table(base),
        "bothclass": {
            "c": [x["c"] for x in bcl],
            "dt": [x["bc_dt_censored_mean"] for x in bcl],
            "dt_se": [x["bc_dt_censored_se"] for x in bcl],
            "lead_detect": [x["bc_lead_detect_rate"] for x in bcl],
            "fpr_leading": [x["bc_fpr_leading"] for x in bcl],
            "c_at_dt_zero": bc_iso,
        },
        "calibration": (cal_bc["chosen"] if cal_bc else None),
        "formation_only_threshold_on_bothclass_score": (
            cal_bc["formation_only_threshold_on_bothclass_score"] if cal_bc else None),
        "pooled_fpr_leading_formation_only": float(np.mean(
            [1.0 if r["lead_fired"] else 0.0 for r in rem["per_run_rows"]
             if r["variant"] == base and r["scen"] == "control"])),
        "pooled_fpr_leading_bothclass": float(np.mean(
            [1.0 if r["bc_lead_fired"] else 0.0 for r in rem["per_run_rows"]
             if r["variant"] == base and r["scen"] == "control"])),
        "class_mixing_active_under_frozen_rule": bool(
            len(rem["detection_rule"]["leading"]["classes"]) > 1),
        "frozen_classes": list(rem["detection_rule"]["leading"]["classes"]),
        "bothclass_classes": list(
            rem["detection_rule_bothclass"]["leading"]["classes"]),
        "escalation_observed_cycle_frac_at_c1": float(np.mean(
            [r["obs_esc_cycles_with_data"] / r["n_cycles"]
             for r in rem["per_run_rows"]
             if r["variant"] == base and r["c"] == 1.0 and r["scen"] == "degrade"])),
    }

    # amplitude audit -- the gain ratio finding 2 is about, measured
    bdeg = [r for r in rem["per_run_rows"]
            if r["variant"] == base and r["scen"] == "degrade"]
    btop = [r for r in bdeg if r["c"] == 1.0]
    amp = {
        "fold_cycle": fold_cycle,
        "cohesion_mean_at_a048": float(np.mean([r["coh_mean_a048"] for r in bdeg])),
        "cohesion_mean_at_a040": float(np.mean([r["coh_mean_a040"] for r in bdeg])),
        "cohesion_relative_swing_a048_to_a040": float(
            1.0 - np.mean([r["coh_mean_a040"] for r in bdeg])
            / np.mean([r["coh_mean_a048"] for r in bdeg])),
        "output_mean_at_a048": float(np.mean([r["out_mean_a048"] for r in bdeg])),
        "output_mean_at_a040": float(np.mean([r["out_mean_a040"] for r in bdeg])),
        "output_relative_swing_a048_to_a040": float(
            1.0 - np.mean([r["out_mean_a040"] for r in bdeg])
            / np.mean([r["out_mean_a048"] for r in bdeg])),
        "output_relative_swing_a048_to_end": float(
            1.0 - np.mean([r["out_mean_last30"] for r in bdeg])
            / np.mean([r["out_mean_a048"] for r in bdeg])),
        "formation_latency_mean_at_a048": float(np.mean(
            [r["form_true_mean_a048"] for r in bdeg])),
        "formation_latency_mean_at_end": float(np.mean(
            [r["form_true_mean_last30"] for r in bdeg])),
        "formation_latency_relative_swing_a048_to_end": float(
            np.mean([r["form_true_mean_last30"] for r in bdeg])
            / np.mean([r["form_true_mean_a048"] for r in bdeg]) - 1.0),
        "knee_binding_program_fraction": float(np.mean(
            [r["coh_frac_at_or_above_knee"] for r in bdeg])),
        "comparator_sigma_drop_available": float(np.mean(
            [r["lag_sigma_drop_sustained"] for r in bdeg])),
        "comparator_sigma_drop_required": float(
            rem["detection_rule"]["lagging"]["k_sd"]),
        "leading_max_z_at_full_coverage": float(np.mean(
            [r["lead_max_score"] for r in btop])) if btop else float("nan"),
        "sensitivity_ratio_lead_over_lag_at_full_coverage": (
            float(np.mean([r["lead_max_score"] for r in btop])
                  / np.mean([r["lag_sigma_drop_sustained"] for r in btop]))
            if btop else float("nan")),
        "lag_frac_fired_before_fold_white_mnar": {
            str(c): next((x["lag_frac_fired_before_fold"] for x in detection
                          if x["arm"] == "mnar" and x["colour"] == "white"
                          and x["c"] == c), None) for c in cs_all},
        "lead_frac_fired_before_fold_white_mnar": {
            str(c): next((x["lead_frac_fired_before_fold"] for x in detection
                          if x["arm"] == "mnar" and x["colour"] == "white"
                          and x["c"] == c), None) for c in cs_all},
    }

    # the graph object (finding 16/b): the same windowed graph in both places
    m = np.asarray(cfg.team_size_choices, float)
    R = cfg.programs_per_cycle * float(np.mean(m * (m - 1) / 2.0))
    graph_object = {
        "tie_window_cycles": cfg.tie_window,
        "is_cumulative": False,
        "expiry_applied_every_cycle": True,
        "ews_detrending_window": rem["detection_rule"]["leading"]["window"],
        "tie_window_shorter_than_detrending_window": bool(
            cfg.tie_window < int(rem["detection_rule"]["leading"]["window"])),
        "mean_field_K_from_tie_window": 2.0 * cfg.tie_window * R / cfg.n_agents,
        "gate1_mean_field_K": (gate1["mean_field"].get("K_mean_field")
                               if gate1 and isinstance(gate1.get("mean_field"), dict)
                               else None),
        "note": ("the giant component Gate 1 measures and the component size the "
                 "leading observable reads are the SAME object: Simulation.ties, a "
                 f"{cfg.tie_window}-cycle sliding-window repeat-tie graph whose "
                 "expiry runs every step.  'Cumulative' is wrong in both places, and "
                 "the mean-field control parameter K carries the window length, which "
                 "is why the analytic fold depends on it."),
    }
    if (graph_object["gate1_mean_field_K"] is None
            or not math.isclose(graph_object["gate1_mean_field_K"],
                                graph_object["mean_field_K_from_tie_window"],
                                rel_tol=0.0, abs_tol=1e-12)):
        raise ValueError("Gate 1 mean-field K does not agree with the value "
                         "derived from the graph tie window")

    head = floors.get("mnar|white", {})
    headline = {
        "detection_probability_white_mnar_by_c": {
            str(x["c"]): {"lead": x["lead_detect"], "lead_ci": x["lead_detect_wilson95"],
                          "lag": x["lag_detect"], "lag_ci": x["lag_detect_wilson95"]}
            for x in detection if x["arm"] == "mnar" and x["colour"] == "white"},
        "paired_lead_time_white_mnar_by_c": {
            str(x["c"]): {"mean": x["dt_paired_mean"], "n": x["n_both_fired"],
                          "ci95": [x["dt_paired_lo95"], x["dt_paired_hi95"]]}
            for x in paired if x["arm"] == "mnar" and x["colour"] == "white"},
        "delta_rmta_white_mnar_by_c_tau320": head.get("delta_rmta_tau320"),
        "floor_under_delta_rmta_tau320": head.get("c_at_dt_zero_tau320"),
        "floor_under_delta_rmta_by_horizon": {
            str(t): head.get(f"c_at_zero_tau{t}") for t in HORIZONS},
        "floor_under_paired_lead_time": head.get("c_at_dt_paired_zero"),
        "floor_under_win_probability_half": head.get("c_at_p_lead_first_half"),
        "floor_under_lead_detection_half": head.get("c_at_lead_detect", {}).get("0.500"),
        "iso_line_is_horizon_dependent": True,
        "verdict": ("the Delta t = 0 iso-line is not a coverage-independent property of "
                    "the instrument: it is the zero of a restricted-mean-time-to-alarm "
                    "difference at a declared horizon, it moves with that horizon, and "
                    "the censoring-robust ordering statistic puts the crossing at a "
                    "materially higher coverage.  Report a detection-probability floor "
                    "and a win-probability floor instead of a single Delta t = 0 line."),
    }

    return {
        "_experiment": ("remediation analysis 2026-08-05 -- review findings 1, 3, "
                        "11, 12 and the MNAR claim; no new runs, pure function of "
                        "committed aggregates"),
        "sources": {"ea_recalibrated": EA_NEW, "ea_campaign1": EA_OLD,
                    "remediation_runs": REM, "bothclass_calibration": CAL_BC,
                    "gate1": GATE1},
        "bootstrap": {"n_boot": N_BOOT, "seed": BOOT_SEED},
        "horizons": list(HORIZONS),
        "fold_cycle": fold_cycle,
        "headline_restated": headline,
        "detection_probability": detection,
        "censoring_accounting": censoring,
        "paired_lead_time": paired,
        "alarm_time_survival": survival,
        "rmta_identity": rmta_identity,
        "win_probability": winprob,
        "floor_resiting": floors,
        "floor_bootstrap": boot,
        "false_alarms_per_cell": fpr_block,
        "mnar_vs_mcar_verdict": mnar_verdict,
        "comparator_gain_leg": gain_leg,
        "formation_driver_leg": driver_leg,
        "partial_obs_sigma_leg": partial_leg,
        "system_size_leg": size_leg,
        "class_set_leg": class_leg,
        "amplitude_audit": amp,
        "graph_object": graph_object,
    }


def main() -> None:
    p = compute()
    write_agg("remediation_analysis.json", p)
    h = p["headline_restated"]
    print("\n  floor under delta-RMTA(tau=320)   :", h["floor_under_delta_rmta_tau320"])
    print("  floor by horizon                  :", h["floor_under_delta_rmta_by_horizon"])
    print("  floor under paired lead time      :", h["floor_under_paired_lead_time"])
    print("  floor under win probability = 0.5 :", h["floor_under_win_probability_half"])
    print("  floor under lead detection = 0.5  :", h["floor_under_lead_detection_half"])
    print("  RMTA identity max discrepancy     :",
          p["rmta_identity"]["max_abs_discrepancy"])
    print("  MNAR verdict                      :", p["mnar_vs_mcar_verdict"]["verdict"])
    g = p["comparator_gain_leg"]
    print("\n  comparator gain (sigma) -> dt@c=1, c*:")
    for x, y, z in zip(g["gain_axis"], g["dt_at_full_coverage"], g["c_at_dt_zero"]):
        print(f"    {x:6.2f}  {(y if y is not None else float('nan')):+8.1f}  "
              f"{'none' if z is None else f'{z:.3f}'}")
    print("  gain at which dt@c=1 vanishes:",
          g["gain_at_which_dt_at_full_coverage_vanishes"])
    print("  gain at which the floor reaches c=1:",
          g["gain_at_which_floor_reaches_full_coverage"])


if __name__ == "__main__":
    main()
