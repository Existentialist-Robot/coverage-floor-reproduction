"""
Tranche-2 remediation campaign 2026-08-05, addressing the residue from the first
remediation tranche.

Six legs, one aggregate.  Everything follows the campaign discipline: no claim,
kill criterion or frozen operating point changes; every scored seed is < 900000 and
derived by seed_for(tag, rep); the 950000+ block is calibration-only; variants exist
to REPORT the headline's dependence on a choice, never to re-tune it.

  Detrending/window sensitivity: rolling var/AC computed
                 under {in-window linear, Gaussian-residual} x window {20, 30, 45}.
                 Per variant the leading threshold is RE-PICKED by the frozen matched false-alarm calibration/large-null recalibration
                 procedure (smallest k on the 0.25 grid with FPR <= 3.5%) on the SAME
                 600-run null (seed block 950000+, re-simulated bit-for-bit), because
                 changing the statistic changes its null distribution.  The lagging
                 comparator is inherited verbatim.  coverage-and-noise campaign white/mid/MNAR, coverage-and-noise campaign seeds.
  Control-scenario baseline: the reference-period mean/sd
                 of each rolling statistic is taken from a MATCHED CONTROL RUN (same
                 seed, scenario="control") instead of from the degradation run's own
                 reference period, 40 of whose 61 cycles are post-onset.  Threshold
                 unchanged (recal k): on a control run the "twin" is the run itself,
                 so the calibration-null FPR is unchanged by construction.
  Bimodality check: (a) 1- vs 2-component Gaussian-mixture BIC on the
                 committed Bifurcation-validation long-run end states (pooled over the fold window,
                 per N); (b) NEW long runs at 2000 cycles (vs 700) around the crossing
                 A, with the middle-band occupancy recorded at four horizons -- if the
                 rising-with-N occupancy of RESULTS section 2 is incomplete relaxation,
                 it must FALL with horizon; if it is genuine, it persists.
  Small-shock relaxation: perturb-and-relax at shock_frac 0.05 (inside the linear
                 regime) alongside the frozen 0.35, with the log-linear fit's R^2 and
                 the fraction of runs that failed to return reported per (A, frac) --
                 the 0.35 recovery times are only relaxation times if the shock stayed
                 in the basin.
  Revocation decomposition: the 0.993-sigma retroactive shift decomposed into LEVEL
                 (mean displacement) and SHAPE (sd of the residual displacement), and
                 re-measured after the re-derived stream is RE-BASELINED against its
                 own reference period -- the operation available to an operator to
                 whom the coverage change was DISCLOSED.  What disclosure cannot
                 absorb is the honest residual of the semantics claim.
  Non-monotone A(t) scenario: a dip that recovers without crossing
                 the fold (dip_recover), a bounded oscillation (oscillate), and a ramp
                 that crosses, holds, and partially recovers (ramp_hold_recover).
                 Scored under the frozen recal AND both-class rules, thresholds
                 UNCHANGED: alarms on non-transitioning trajectories are false
                 positives the paper must own.

Run:   python run_tranche2.py [--workers N] [--quick]
Writes: aggregates/tranche2_runs.json
"""

from __future__ import annotations

import argparse
import math

import numpy as np

from calibrate_recal import (K_GRID, NULL_SEED_BASE, TARGET_FPR, _cfg as null_cfg,
                             fpr_curve, max_sustained_min)
from indicators import fill_gaps, rolling_var_ac, lagging_alarm, alarm_from_z, baseline_z
from model import Config, run_config
from scoring import score_run, bootstrap_mean
from runner import pmap, seed_for, write_agg, load_rule

C_GRID = (0.25, 0.333, 0.5, 0.75, 1.0)
ARM, COLOUR, SEVERITY = "mnar", "white", "mid"

# ---------------------------------------------------------------- detrending-window sensitivity variants
DETREND_WINDOWS = (20, 30, 45)
DETREND_MODES = ("linear", "gauss")
GAUSS_BANDWIDTH = 10.0          # cycles; declared, not tuned (Dakos et al. 2012 use
                                # kernel-residual detrending with a fixed bandwidth)
BASE_VARIANT = "linear_w30"     # the frozen statistic


def _gauss_residual(x: np.ndarray, bw: float) -> np.ndarray:
    """Residual after Gaussian-kernel smoothing of the WHOLE series (the Dakos
    construction), computed on the gap-filled series; nan stays nan."""
    y = np.asarray(x, float)
    n = y.size
    out = np.full(n, np.nan)
    idx = np.arange(n, dtype=float)
    ok = np.isfinite(y)
    if ok.sum() < 5:
        return out
    xi, yi = idx[ok], y[ok]
    for i in range(n):
        if not ok[i]:
            continue
        w = np.exp(-0.5 * ((xi - i) / bw) ** 2)
        out[i] = y[i] - float((w @ yi) / w.sum())
    return out


def _rolling_var_ac_gauss(x: np.ndarray, window: int, bw: float,
                          min_frac: float, max_gap: int) -> tuple[np.ndarray, np.ndarray]:
    """Rolling var / lag-1 AC of the KERNEL-RESIDUAL series, no in-window detrend
    (the trend was removed globally)."""
    y = fill_gaps(x, max_gap)
    r = _gauss_residual(y, bw)
    n = r.size
    var = np.full(n, np.nan)
    ac = np.full(n, np.nan)
    need = int(math.ceil(min_frac * window))
    for end in range(window - 1, n):
        seg = r[end - window + 1:end + 1]
        m = ~np.isnan(seg)
        if m.sum() < max(4, need):
            continue
        d = seg - np.nanmean(seg)
        dm = d[m]
        var[end] = float(np.var(dm, ddof=1))
        pair = m[:-1] & m[1:]
        if pair.sum() < 4:
            continue
        a, b = d[:-1][pair], d[1:][pair]
        sa, sb = a.std(), b.std()
        if sa < 1e-12 or sb < 1e-12:
            continue
        ac[end] = float(np.mean((a - a.mean()) * (b - b.mean())) / (sa * sb))
    return var, ac


def _variant_score(series: dict, mode: str, window: int, L: dict) -> np.ndarray:
    """The composite leading score (combine='mean', frozen) under a stat variant,
    standardised against the run's own reference period exactly as baseline_z does."""
    zs = []
    for cls in L["classes"]:
        if mode == "linear":
            var, ac = rolling_var_ac(series[cls], window,
                                     float(L["min_frac"]), int(L["max_gap"]))
        else:
            var, ac = _rolling_var_ac_gauss(series[cls], window, GAUSS_BANDWIDTH,
                                            float(L["min_frac"]), int(L["max_gap"]))
        for x in (var, ac):
            base = x[window - 1:int(L["start"])]
            base = base[~np.isnan(base)]
            if base.size < 8:
                zs.append(np.full_like(x, np.nan))
                continue
            mu, sd = float(base.mean()), float(base.std(ddof=1))
            zs.append((x - mu) / (sd if sd > 1e-9 else 1e-9))
    return np.nanmean(np.vstack(zs), axis=0) if zs else np.asarray([])


def _variants():
    return [f"{m}_w{w}" for m in DETREND_MODES for w in DETREND_WINDOWS]


def _detrend_null_job(job):
    """One 950000+ null run reduced to the sustained-score statistic per variant."""
    seed = job
    cfg = null_cfg("control", seed)
    s = run_config(cfg).series()
    L = load_rule("detection_rule_recal.json")["leading"]
    out = {"seed": int(seed)}
    for name in _variants():
        mode, w = name.split("_w")
        score = _variant_score(s, mode, int(w), L)
        out[name] = max_sustained_min(score, int(L["sustain"]), int(L["start"]))
    return out


def _detrend_deg_job(job):
    """One coverage-and-noise campaign degradation/control run scored under every variant at its re-picked
    threshold (passed in), plus the verbatim lagging comparator."""
    c, rep, scen, k_by_variant = job
    rule = load_rule("detection_rule_recal.json")
    L = rule["leading"]
    cfg = Config(consent_c=c, consent_arm=ARM, eps_colour=COLOUR,
                 eps_severity=SEVERITY, scenario=scen,
                 seed=seed_for(f"ea|{c}|{ARM}", rep))
    s = run_config(cfg).series()
    lag = lagging_alarm(s["output_close"], rule["lagging"])
    t_lag = lag["alarm_cycle"]
    row = {"c": c, "rep": rep, "scen": scen, "seed": cfg.seed,
           "t_lag": t_lag, "n_cycles": cfg.n_cycles}
    for name in _variants():
        mode, w = name.split("_w")
        score = _variant_score(s, mode, int(w), L)
        al = alarm_from_z(score, score, float(k_by_variant[name]),
                          int(L["sustain"]), "var_only", int(L["start"]))
        row[f"t_lead_{name}"] = al["alarm_cycle"]
    return row


# ---------------------------------------------------------------- control-baseline sensitivity baseline
def _cbase_job(job):
    """Score one coverage-and-noise campaign degradation run with reference stats taken from its matched
    control twin (same seed, scenario='control').  Threshold = frozen recal k."""
    c, rep = job
    rule = load_rule("detection_rule_recal.json")
    L = rule["leading"]
    seed = seed_for(f"ea|{c}|{ARM}", rep)
    base = dict(consent_c=c, consent_arm=ARM, eps_colour=COLOUR,
                eps_severity=SEVERITY, seed=seed)
    s_deg = run_config(Config(scenario="degrade", **base)).series()
    s_ctl = run_config(Config(scenario="control", **base)).series()
    w, start = int(L["window"]), int(L["start"])
    zs = []
    for cls in L["classes"]:
        var_d, ac_d = rolling_var_ac(s_deg[cls], w, float(L["min_frac"]),
                                     int(L["max_gap"]))
        var_c, ac_c = rolling_var_ac(s_ctl[cls], w, float(L["min_frac"]),
                                     int(L["max_gap"]))
        for xd, xc in ((var_d, var_c), (ac_d, ac_c)):
            ref = xc[w - 1:start]
            ref = ref[~np.isnan(ref)]
            if ref.size < 8:
                zs.append(np.full_like(xd, np.nan))
                continue
            mu, sd = float(ref.mean()), float(ref.std(ddof=1))
            zs.append((xd - mu) / (sd if sd > 1e-9 else 1e-9))
    score = np.nanmean(np.vstack(zs), axis=0)
    al = alarm_from_z(score, score, float(L["k"]), int(L["sustain"]),
                      "var_only", start)
    own = score_run(s_deg, rule, Config().n_cycles)     # own-baseline reference
    lag = lagging_alarm(s_deg["output_close"], rule["lagging"])
    return {"c": c, "rep": rep, "seed": seed,
            "t_lead_ctlbase": al["alarm_cycle"],
            "t_lead_ownbase": own["t_lead"],
            "t_lag": lag["alarm_cycle"], "n_cycles": Config().n_cycles}


# ---------------------------------------------------------------- end-state mixture analysis bimodality
LONG2_A = (0.38, 0.37, 0.36)
LONG2_SIZES = ((400, 12), (800, 24))
LONG2_CYCLES = 2000
LONG2_SEEDS = 12
LONG2_CHECKPOINTS = (600, 1000, 1500, 2000)


def _long2_job(job):
    n, p, a, rep = job
    cfg = Config(n_agents=n, programs_per_cycle=p, scenario="control", a_control=a,
                 burn_in_a=0.90, n_cycles=LONG2_CYCLES,
                 seed=seed_for(f"t2long|{n}|{a}", rep))
    s = run_config(cfg).series()
    z = s["mean_degree"]
    row = {"n": n, "a": a, "rep": rep, "seed": cfg.seed}
    for cp in LONG2_CHECKPOINTS:
        row[f"z_tail_{cp}"] = float(z[cp - 60:cp].mean())
    return row


def _mixture_bic(x: np.ndarray) -> dict:
    """1- vs 2-component Gaussian mixture by EM (deterministic init: split at the
    median), BIC for each.  Positive delta_bic = the 2-component fit is preferred."""
    x = np.sort(np.asarray(x, float))
    n = x.size
    mu, sd = float(x.mean()), float(x.std(ddof=1))
    ll1 = float(np.sum(-0.5 * math.log(2 * math.pi) - math.log(max(sd, 1e-9))
                       - 0.5 * ((x - mu) / max(sd, 1e-9)) ** 2))
    bic1 = -2 * ll1 + 2 * math.log(n)
    lo, hi = x[:n // 2], x[n // 2:]
    pi = 0.5
    m1, s1 = float(lo.mean()), max(float(lo.std()), 1e-3)
    m2, s2 = float(hi.mean()), max(float(hi.std()), 1e-3)
    for _ in range(500):
        p1 = pi * np.exp(-0.5 * ((x - m1) / s1) ** 2) / s1
        p2 = (1 - pi) * np.exp(-0.5 * ((x - m2) / s2) ** 2) / s2
        tot = np.maximum(p1 + p2, 1e-300)
        g = p1 / tot
        pi = float(g.mean())
        if g.sum() < 1e-9 or (1 - g).sum() < 1e-9:
            break
        m1 = float((g @ x) / g.sum())
        m2 = float(((1 - g) @ x) / (1 - g).sum())
        s1 = max(1e-3, math.sqrt(float((g @ (x - m1) ** 2) / g.sum())))
        s2 = max(1e-3, math.sqrt(float(((1 - g) @ (x - m2) ** 2) / (1 - g).sum())))
    ll2 = float(np.sum(np.log(np.maximum(
        pi * np.exp(-0.5 * ((x - m1) / s1) ** 2) / (s1 * math.sqrt(2 * math.pi))
        + (1 - pi) * np.exp(-0.5 * ((x - m2) / s2) ** 2) / (s2 * math.sqrt(2 * math.pi)),
        1e-300))))
    bic2 = -2 * ll2 + 5 * math.log(n)
    return {"n": n, "bic_1comp": bic1, "bic_2comp": bic2,
            "delta_bic_1_minus_2": bic1 - bic2,
            "two_component_preferred": bool(bic1 - bic2 > 10.0),
            "comp_means": [m1, m2], "comp_sds": [s1, s2], "comp_weight_1": pi}


# ---------------------------------------------------------------- small-perturbation recovery small shock
SHOCK2_A = (0.48, 0.44, 0.42, 0.40, 0.39, 0.38, 0.37)
SHOCK2_FRACS = (0.05, 0.35)
SHOCK2_SEEDS = 12
SHOCK2_AT = 40
SHOCK2_CYCLES = 140


def _shock2_job(job):
    a, frac, rep = job
    cfg = Config(scenario="control", a_control=a, burn_in_a=0.90,
                 n_cycles=SHOCK2_CYCLES, shock_cycle=SHOCK2_AT, shock_frac=frac,
                 seed=seed_for(f"gate1shock|{a}", rep))   # gate1's seeds: matched
    s = run_config(cfg).series()
    z = s["mean_degree"]
    z_star = float(z[max(0, SHOCK2_AT - 25):SHOCK2_AT].mean())
    z_end = float(z[-30:].mean())
    seg = z[SHOCK2_AT:SHOCK2_AT + 60]
    dev = np.abs(z_star - seg)
    ok = dev > 1e-3
    lam, r2 = float("nan"), float("nan")
    if ok.sum() >= 8:
        x = np.arange(seg.size, dtype=float)[ok]
        y = np.log(dev[ok])
        xm, ym = x.mean(), y.mean()
        sxx = float((x - xm) @ (x - xm))
        if sxx > 1e-9:
            slope = float((x - xm) @ (y - ym)) / sxx
            lam = -slope
            resid = y - (ym + slope * (x - xm))
            sst = float((y - ym) @ (y - ym))
            r2 = 1.0 - float(resid @ resid) / sst if sst > 1e-12 else float("nan")
    return {"a": a, "frac": frac, "rep": rep, "z_star": z_star, "z_end": z_end,
            "recovery_rate": lam, "fit_r2": r2,
            "returned": bool(z_end > 0.6 * z_star)}


# ---------------------------------------------------------------- record-revision decomposition revocation
REV_REPS = 60


def _rev_job(job):
    """Level/shape decomposition of the retroactive shift, and its comparison
    against the DISCLOSED coverage change of the same magnitude.

    Both streams are already self-baselined by baseline_z (the 0.993-sigma of the
    aux arm is a displacement between two self-standardised series), so the
    decomposition is on that same object.  The disclosed comparator is the
    static_equiv run on the SAME seed -- the same agents non-consenting from cycle
    0, i.e. the record an operator holds when a coverage reduction of identical
    magnitude was disclosed up front.  The reconstruction residual |re-derived -
    static_equiv| measures what re-derivation cannot reconstruct of the honest
    low-coverage record: the pre-revocation history was captured at HIGH coverage
    (different noise realisations, different partial-observation error), and
    deleting agents from it is not the same as never having observed them."""
    rep = job
    rule = load_rule("detection_rule_recal.json")
    L = rule["leading"]
    base = dict(consent_c=0.5, consent_arm="mnar", eps_colour="red",
                eps_severity="mid", scenario="degrade",
                seed=seed_for(f"aux|{rep}", rep))
    cfg = Config(revocation=True, revocation_retroactive=True, **base)
    s = run_config(cfg).series()
    sq = run_config(Config(revocation=True, revocation_retroactive=True,
                           revocation_cycle=0, **base)).series()
    start, t_rev = int(L["start"]), cfg.revocation_cycle
    args = (int(L["window"]), tuple(L["classes"]), start,
            float(L["min_frac"]), int(L["max_gap"]))
    za = baseline_z(s, *args)
    srev = dict(s)
    for cls in L["classes"]:
        srev[cls] = s[cls + "_rev"]
    zb = baseline_z(srev, *args)
    zq = baseline_z(sq, *args)
    ia = 0.5 * (za["z_var"] + za["z_ac"])       # as published (pre-revocation record)
    ib = 0.5 * (zb["z_var"] + zb["z_ac"])       # retroactively re-derived
    iq = 0.5 * (zq["z_var"] + zq["z_ac"])       # honest low-coverage-from-day-one
    seg = slice(start, max(start + 1, t_rev))
    da, db, dq = ia[seg], ib[seg], iq[seg]
    k = float(L["k"])

    def _cmp(x, y):
        m = np.isfinite(x) & np.isfinite(y)
        if not m.any():
            return None
        diff = y[m] - x[m]
        return {"mean_abs": float(np.mean(np.abs(diff))),
                "level": float(diff.mean()),
                "shape_sd": (float(diff.std(ddof=1)) if diff.size > 2
                             else float("nan")),
                "alarm_flip_frac": float(np.mean((x[m] >= k) != (y[m] >= k)))}

    retro = _cmp(da, db)          # the 0.993-sigma object, decomposed
    disclosed = _cmp(da, dq)      # published vs disclosed-from-day-one record
    recon = _cmp(dq, db)          # re-derivation vs the honest record
    if retro is None:
        return {"rep": rep, "seed": cfg.seed, "valid": False}
    return {"rep": rep, "seed": cfg.seed, "valid": True,
            "retroactive": retro, "disclosed_equiv": disclosed,
            "reconstruction_residual": recon}


# ---------------------------------------------------------------- non-monotone driver scenarios scenarios
class ConfigT2(Config):
    """Config with three declared non-monotone A(t) shapes.  The base scenarios
    ("degrade", "control") are untouched, so every existing run reproduces
    bit-for-bit; model.py is not edited."""

    def a_of_t(self, t: int) -> float:
        s = self.scenario
        if s == "dip_recover":
            # down to 0.40 (ABOVE the empirical fold A_c ~ 0.37) and back: no
            # transition occurs, so any alarm is a false positive
            if t <= 100:
                return self.a_high
            if t <= 140:
                return self.a_high + (0.40 - self.a_high) * (t - 100) / 40.0
            if t <= 200:
                return 0.40 + (self.a_high - 0.40) * (t - 140) / 60.0
            return self.a_high
        if s == "oscillate":
            # bounded oscillation a_high .. 0.40, period 80, from cycle 50 on
            if t <= 50:
                return self.a_high
            amp = 0.5 * (self.a_high - 0.40)
            return self.a_high - amp * (1.0 - math.cos(2.0 * math.pi
                                                       * (t - 50) / 80.0))
        if s == "ramp_hold_recover":
            # crosses the fold, holds, then recovers -- the transition happens,
            # then the driver (not the system) comes back
            if t <= 50:
                return self.a_high
            if t <= 180:
                return self.a_high + (self.a_low - self.a_high) * (t - 50) / 130.0
            if t <= 230:
                return self.a_low
            if t <= 300:
                return self.a_low + (self.a_high - self.a_low) * (t - 230) / 70.0
            return self.a_high
        return super().a_of_t(t)


SCEN2 = ("dip_recover", "oscillate", "ramp_hold_recover")
SCEN2_REPS = 40
SCEN2_C = 0.5


def _scen_job(job):
    scen, rep = job
    recal = load_rule("detection_rule_recal.json")
    both = load_rule("detection_rule_bothclass.json")
    cfg = ConfigT2(consent_c=SCEN2_C, consent_arm=ARM, eps_colour=COLOUR,
                   eps_severity=SEVERITY, scenario=scen,
                   seed=seed_for(f"t2scen|{scen}", rep))
    s = run_config(cfg).series()
    a = s["a_t"]
    z = s["mean_degree"]
    r1 = score_run(s, recal, cfg.n_cycles)
    r2 = score_run(s, both, cfg.n_cycles)
    return {"scen": scen, "rep": rep, "seed": cfg.seed,
            "a_min": float(a.min()), "z_min": float(z.min()),
            "z_end": float(z[-30:].mean()),
            "collapsed": bool(z[-30:].mean() < 0.5 * z[:30].mean()),
            "t_lead_recal": r1["t_lead"], "t_lag_recal": r1["t_lag"],
            "lead_fired_recal": r1["lead_fired"], "lag_fired_recal": r1["lag_fired"],
            "t_lead_both": r2["t_lead"], "lead_fired_both": r2["lead_fired"],
            "dt_censored_recal": r1["dt_censored"]}


# ---------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()
    W = args.workers
    q = args.quick
    rule = load_rule("detection_rule_recal.json")
    L = rule["leading"]

    # ---- detrending-window sensitivity: null sweep -> per-variant thresholds -> coverage-and-noise campaign sweep -------------
    n_null = 60 if q else 600
    null_rows = pmap(_detrend_null_job,
                     [NULL_SEED_BASE + r for r in range(n_null)],
                     W, "detrending-window sensitivity variant calibration null")
    k_by_variant, curves = {}, {}
    for name in _variants():
        curve = fpr_curve([r[name] for r in null_rows], K_GRID)
        adm = [x for x in curve if x["fpr"] <= TARGET_FPR]
        k_by_variant[name] = min(x["k"] for x in adm) if adm else float("inf")
        curves[name] = next((x for x in curve if x["k"] == k_by_variant[name]), None)
    reps20 = 8 if q else 40
    djobs = [(c, r, scen, k_by_variant) for c in C_GRID for r in range(reps20)
             for scen in ("degrade", "control")]
    drows = pmap(_detrend_deg_job, djobs, W, "detrending-window sensitivity detrend/window coverage-and-noise campaign sweep")

    detrend_tab = []
    for name in _variants():
        for c in C_GRID:
            d = [r for r in drows if r["c"] == c and r["scen"] == "degrade"]
            k_ = [r for r in drows if r["c"] == c and r["scen"] == "control"]
            nc = d[0]["n_cycles"]
            dts = [(nc if r["t_lag"] is None else r["t_lag"])
                   - (nc if r[f"t_lead_{name}"] is None else r[f"t_lead_{name}"])
                   for r in d]
            detrend_tab.append({
                "variant": name, "k": k_by_variant[name], "c": c,
                "n_reps": len(d),
                "dt_censored_mean": float(np.mean(dts)),
                "lead_detect": float(np.mean(
                    [r[f"t_lead_{name}"] is not None for r in d])),
                "fpr_leading": float(np.mean(
                    [r[f"t_lead_{name}"] is not None for r in k_])),
            })

    # ---- control-baseline sensitivity: control-scenario baseline -------------------------------------
    reps21 = 8 if q else 40
    crows = pmap(_cbase_job, [(c, r) for c in C_GRID for r in range(reps21)],
                 W, "control-baseline sensitivity control-baseline standardisation")
    cbase_tab = []
    for c in C_GRID:
        rs = [r for r in crows if r["c"] == c]
        nc = rs[0]["n_cycles"]
        def _dt(key):
            return [(nc if r["t_lag"] is None else r["t_lag"])
                    - (nc if r[key] is None else r[key]) for r in rs]
        cbase_tab.append({
            "c": c, "n_reps": len(rs),
            "dt_censored_ownbase": float(np.mean(_dt("t_lead_ownbase"))),
            "dt_censored_ctlbase": float(np.mean(_dt("t_lead_ctlbase"))),
            "lead_detect_ownbase": float(np.mean(
                [r["t_lead_ownbase"] is not None for r in rs])),
            "lead_detect_ctlbase": float(np.mean(
                [r["t_lead_ctlbase"] is not None for r in rs])),
            "median_advance_cycles": float(np.median(
                [r["t_lead_ownbase"] - r["t_lead_ctlbase"] for r in rs
                 if r["t_lead_ownbase"] is not None
                 and r["t_lead_ctlbase"] is not None] or [float("nan")])),
        })

    # ---- end-state mixture analysis: bimodality -----------------------------------------------------
    import json as _json
    import os as _os
    from runner import AGG as _AGG
    with open(_os.path.join(_AGG, "gate1_bifurcation.json")) as f:
        gate1 = _json.load(f)
    bic_committed = {}
    for n in ("400", "800"):
        pooled = [r["z_end"] for r in gate1["long_run_rows"] if str(r["n"]) == n]
        if len(pooled) >= 20:
            bic_committed[n] = _mixture_bic(np.asarray(pooled))
    sizes22 = LONG2_SIZES[:1] if q else LONG2_SIZES
    seeds22 = 3 if q else LONG2_SEEDS
    ljobs = [(n, p, a, r) for (n, p) in sizes22 for a in LONG2_A
             for r in range(seeds22)]
    lrows = pmap(_long2_job, ljobs, W, "end-state mixture analysis long-horizon runs (2000 cycles)")
    occ_tab = []
    for (n, p) in sizes22:
        rs = [r for r in lrows if r["n"] == n]
        z_hi = float(np.mean([r[f"z_tail_{LONG2_CHECKPOINTS[-1]}"] for r in rs
                              if r["a"] == LONG2_A[0]]))
        z_lo = float(np.mean([r[f"z_tail_{LONG2_CHECKPOINTS[-1]}"] for r in rs
                              if r["a"] == LONG2_A[-1]]))
        rngz = max(1e-9, z_hi - z_lo)
        mid = 0.5 * (z_hi + z_lo)
        for cp in LONG2_CHECKPOINTS:
            zs = np.asarray([r[f"z_tail_{cp}"] for r in rs if r["a"] == 0.37])
            occ_tab.append({"n": n, "horizon": cp, "a": 0.37,
                            "n_seeds": int(zs.size),
                            "mid_band_occupancy": float(np.mean(
                                np.abs(zs - mid) < 0.20 * rngz))})
        pooled = np.asarray([r[f"z_tail_{LONG2_CHECKPOINTS[-1]}"] for r in rs])
        bic_committed[f"{n}_horizon2000"] = _mixture_bic(pooled)

    # ---- small-perturbation recovery: small shock ----------------------------------------------------
    seeds23 = 4 if q else SHOCK2_SEEDS
    sjobs = [(a, f, r) for a in SHOCK2_A for f in SHOCK2_FRACS
             for r in range(seeds23)]
    srows = pmap(_shock2_job, sjobs, W, "small-perturbation recovery perturb-and-relax, two shock sizes")
    shock_tab = []
    for f in SHOCK2_FRACS:
        for a in SHOCK2_A:
            rs = [r for r in srows if r["a"] == a and r["frac"] == f]
            lam = [r["recovery_rate"] for r in rs if np.isfinite(r["recovery_rate"])]
            r2s = [r["fit_r2"] for r in rs if np.isfinite(r["fit_r2"])]
            shock_tab.append({
                "frac": f, "a": a, "n_seeds": len(rs),
                "recovery_rate_mean": float(np.mean(lam)) if lam else None,
                "recovery_time_cycles": (float(1.0 / np.mean(lam))
                                         if lam and np.mean(lam) > 1e-6 else None),
                "fit_r2_mean": float(np.mean(r2s)) if r2s else None,
                "returned_frac": float(np.mean([r["returned"] for r in rs])),
            })

    # ---- record-revision decomposition: revocation decomposition ---------------------------------------
    reps24 = 10 if q else REV_REPS
    rrows = [r for r in pmap(_rev_job, list(range(reps24)),
                             W, "record-revision decomposition revocation decomposition") if r["valid"]]
    rev = {"n_reps": len(rrows)}
    for leg in ("retroactive", "disclosed_equiv", "reconstruction_residual"):
        rev[leg] = {
            stat: bootstrap_mean([r[leg][stat] for r in rrows
                                  if r[leg] is not None])
            for stat in ("mean_abs", "level", "shape_sd", "alarm_flip_frac")}

    # ---- non-monotone driver scenarios: non-monotone scenarios -----------------------------------------
    reps25 = 8 if q else SCEN2_REPS
    scrows = pmap(_scen_job, [(s_, r) for s_ in SCEN2 for r in range(reps25)],
                  W, "non-monotone driver scenarios non-monotone A(t) scenarios")
    scen_tab = []
    for s_ in SCEN2:
        rs = [r for r in scrows if r["scen"] == s_]
        scen_tab.append({
            "scenario": s_, "n_reps": len(rs),
            "collapsed_frac": float(np.mean([r["collapsed"] for r in rs])),
            "a_min": float(np.mean([r["a_min"] for r in rs])),
            "lead_alarm_rate_recal": float(np.mean(
                [r["lead_fired_recal"] for r in rs])),
            "lead_alarm_rate_bothclass": float(np.mean(
                [r["lead_fired_both"] for r in rs])),
            "lag_alarm_rate": float(np.mean([r["lag_fired_recal"] for r in rs])),
            "dt_censored_recal_mean": float(np.mean(
                [r["dt_censored_recal"] for r in rs])),
            "median_t_lead_recal": (float(np.median(
                [r["t_lead_recal"] for r in rs if r["lead_fired_recal"]]))
                if any(r["lead_fired_recal"] for r in rs) else None),
        })

    write_agg("tranche2_runs.json", {
        "_experiment": ("tranche-2 remediation 2026-08-05: "
                        "detrending/window sensitivity, control-scenario baseline, "
                        "bimodality (mixture BIC + 2000-cycle horizon probe), "
                        "small-shock relaxation, revocation level/shape "
                        "decomposition, non-monotone A(t) scenarios"),
        "design": {
            "detrend_variants": _variants(), "gauss_bandwidth": GAUSS_BANDWIDTH,
            "base_variant": BASE_VARIANT, "k_by_variant": k_by_variant,
            "n_null": n_null, "c_grid": list(C_GRID),
            "arm": ARM, "colour": COLOUR, "severity": SEVERITY,
            "long2": {"a": list(LONG2_A), "sizes": [list(x) for x in sizes22],
                      "cycles": LONG2_CYCLES, "seeds": seeds22,
                      "checkpoints": list(LONG2_CHECKPOINTS)},
            "shock2": {"a": list(SHOCK2_A), "fracs": list(SHOCK2_FRACS),
                       "seeds": seeds23, "shock_cycle": SHOCK2_AT,
                       "seed_tag": "gate1shock|{a} -- gate1's own seeds, matched"},
            "rev": {"reps": reps24, "cell": "c=0.5, red/mid, MNAR (the fixed auxiliary-condition scope bound)"},
            "scenarios": {"names": list(SCEN2), "reps": reps25, "c": SCEN2_C,
                          "thresholds": "frozen recal + bothclass, UNCHANGED"},
            "seed_hygiene": ("all scored seeds seed_for(tag, rep) < 900000; the "
                             "950000+ block is calibration-only (detrending-window sensitivity nulls)"),
        },
        "a20_detrend": {"per_variant_thresholds": {n: curves[n] for n in _variants()},
                        "table": detrend_tab, "null_rows": null_rows,
                        "per_run_rows": drows},
        "a21_control_baseline": {"table": cbase_tab, "per_run_rows": crows},
        "a22_bimodality": {"mixture_bic": bic_committed, "occupancy_by_horizon": occ_tab,
                           "per_run_rows": lrows},
        "a23_small_shock": {"table": shock_tab, "per_run_rows": srows},
        "a24_revocation_decomposition": rev,
        "a24_per_run_rows": rrows,
        "a25_scenarios": {"table": scen_tab, "per_run_rows": scrows},
    })

    print("\n  detrending-window sensitivity dt@c=1 by variant:")
    for name in _variants():
        row = next(r for r in detrend_tab if r["variant"] == name and r["c"] == 1.0)
        print(f"    {name:12s} k={k_by_variant[name]:5.2f}  "
              f"dt={row['dt_censored_mean']:+7.1f}  DR={row['lead_detect']:.2f}  "
              f"FPR={row['fpr_leading']:.3f}")
    print("  control-baseline sensitivity own-vs-control baseline dt@c=0.5:",
          next((f"{r['dt_censored_ownbase']:+.1f} -> {r['dt_censored_ctlbase']:+.1f}"
                for r in cbase_tab if r["c"] == 0.5), "n/a"))
    print("  end-state mixture analysis BIC:", {k: round(v["delta_bic_1_minus_2"], 1)
                          for k, v in bic_committed.items()})
    print("  record-revision decomposition retroactive:", round(rev["retroactive"]["mean_abs"]["mean"], 3),
          "(level", round(rev["retroactive"]["level"]["mean"], 3),
          "shape", round(rev["retroactive"]["shape_sd"]["mean"], 3), ")",
          "disclosed:", round(rev["disclosed_equiv"]["mean_abs"]["mean"], 3),
          "recon-resid:", round(rev["reconstruction_residual"]["mean_abs"]["mean"], 3))
    for r in scen_tab:
        print(f"  non-monotone driver scenarios {r['scenario']:18s} collapsed={r['collapsed_frac']:.2f} "
              f"leadFP={r['lead_alarm_rate_recal']:.2f} "
              f"lagFP={r['lag_alarm_rate']:.2f}")


if __name__ == "__main__":
    main()
