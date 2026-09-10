# Simulation results — COMPLEX NETWORKS 2026

*Report of 2026-08-04. Implements the frozen simulation specification `sim-spec-FROZEN-2026-08-04`, which is not distributed with this bundle;
seven implementation-forced amendments are logged in that spec's §9. Every number below is
bound to a committed aggregate through [`claims.json`](claims.json) and checked by
[`verify_reported_numbers.py`](verify_reported_numbers.py). **§7 is written for the
kill/fallback judgment and should be read before §2–6.***

> **CORRECTED 2026-09-04.** "Bound through `claims.json`" is not one guarantee. Of the 128
> claims, **61 of 128** are independently recomputed from per-run rows, **27 of 128** are
> reproduced by re-running the same analysis module rather than checked independently, and
> **40 of 128** are not recomputed at all. §6 has the full breakdown, the aggregate-file
> coverage, and what changes about the verifier's assurance.

> **Read §8 before quoting anything from §3–§5.** The follow-up campaign of 2026-08-05
> (§8, spec §9 amendments A-8 … A-10) recalibrated the detection threshold on a 600-run
> null, firmed the iso-line, and raised E-B to ≥ 500 matched pairs. It **supersedes §3's
> lead times, §3.1's iso-line, §3.2's false-alarm rates, §4.1's timing contrast and §5's
> E-B numbers.** §§2–6 are kept verbatim as the campaign-1 record — both operating points
> stay independently verifiable — but the numbers the paper reports are §8's.

---

## 1. What ran

| stage | script | scale | wall time |
|---|---|---|---|
| rule calibration | [`calibrate_rule.py`](calibrate_rule.py) | 80 runs on calibration seeds ≥ 900000 | ~5 min |
| **Gate 1** — bifurcation obligation | [`run_gate1_bifurcation.py`](run_gate1_bifurcation.py) | 832 equilibrium + 84 shock + 144 long runs | ~11 min |
| **E-A** headline sweep | [`run_ea_headline.py`](run_ea_headline.py) | 4 320 runs (4 c × 3 colours × 3 severities × 2 arms × 2 scenarios × 30 reps) | 1 199 s |
| **auxiliary arms** | [`run_aux_arms.py`](run_aux_arms.py) | 360 + 120 runs at (c = 0.5, mid red noise) | 145 s |
| **E-B** joint sweep | [`run_eb.py`](run_eb.py) | 2 160 runs (3 c × 3 g × 2 modes × 120 reps) | 684 s |
| figures | [`figures.py`](figures.py) | 6 figures, PNG + SVG | seconds |

Substrate: N = 400 agents, 320 recorded cycles after a 170-cycle unrecorded burn-in,
12 programs per cycle, teams of 3–5 on the Guimerà newcomer(p = 0.28)/repeat-tie(q = 0.55)
rule with live turnover, repeat-tie window 20 cycles. Total compute well under two hours
on one desktop CPU. All configs are in the `design` block of each aggregate.

**The frozen detection rule** (`configs/detection_rule.json`, set on calibration seeds only
and applied unchanged in every cell): rolling variance and rolling lag-1 autocorrelation of
the in-window linearly detrended observed formation-latency series, window 30, each
standardised against its own distribution over the run's reference period (cycles 29–90),
combined as the mean of the two z-scores, alarm when the composite ≥ **5.0** for **5**
consecutive cycles. **The rule *form* itself — window 30, formation-only classes, sustain 5,
mean combine — is the argmax of a 36-candidate search** (window ∈ {30, 40, 50} × classes ∈
{both, formation only} × sustain ∈ {3, 5} × combine ∈ {and, mean, var_only}) on median lead
time at full coverage over the 40 calibration degradation runs, ties broken on AUC
([`calibrate_rule.py`](calibrate_rule.py) lines 100–128); what follows describes only how the
*threshold* was calibrated within the chosen form. Lagging comparator: program-close output on a 10-cycle trailing mean,
alarm below baseline − **3.5** σ for 5 cycles. Both thresholds were set to the smallest
value with calibration control-arm FPR ≤ 5%, so lead time is compared at matched
false-alarm rate; among admissible lagging rules the earliest-detecting was taken, i.e. the
comparator gets every advantage.

---

## 2. GATE 1 — bifurcation obligation: **PASS**

No early-warning signal may be reported before this gate (the bifurcation gate defined in the specification).
Four independent lines of evidence, all on the cumulative repeat-tie graph:

1. **Mean-field saddle-node (analytic).** The substrate's own update rule is
   `z = K·A·g(z)` with `K = 7.6`; it has a tangency — one fixed point below the fold,
   three above — at **A = 0.29925**, `z_fold = 1.109`, with a bistable window
   **0.0735** wide in A. This is a fold in the dynamics, not a percolation threshold in a
   parameter, which is what Boettiger & Hastings (2013) requires before any EWS is claimed.
2. **End-state bimodality (700-cycle runs).** At A = 0.38 every seed ends on the upper
   branch; at A = 0.36 every seed ends on the lower branch; at **A = 0.37 the seeds split**
   (z_end s.d. 0.41 at N = 400, 0.50 at N = 800). Middle-band occupancy at the crossing is
   **0.167** (N = 400) and 0.25 (N = 800), and the largest gap in the pooled end-state
   distribution is 28% (N = 400) / 50% (N = 800) of the branch range. Individual runs sit on
   one branch or the other; the ensemble mean interpolates two modes.
3. **Hysteresis at four system sizes.** With programs-per-cycle scaled with N to hold K
   fixed, the bistable A-window is **0.14–0.15** wide at every N ∈ {200, 400, 800, 1600},
   and the upper-vs-lower branch gap in z reaches **1.73** at N = 1600. Neither a continuous
   transition nor noise alone can produce this.
4. **Critical slowing down on the latent structure.** Post-shock (35% of repeat ties
   deleted) recovery time of the mean degree diverges from **31.5 cycles at A = 0.48 to
   333.6 cycles at A = 0.37** — a factor of **10.6** — measured directly on the substrate,
   independent of any observed indicator.

**Mixture-model cells.** Deterministic expectation-maximization fits of one- and
two-component Gaussian mixtures to pooled long-run end states give
ΔBIC(1−2) = **204.5 at N = 400** and **258.9 at N = 800** on 72 end states per size
from the committed 700-cycle Gate-1 rows. The corresponding 2,000-cycle cells give
**119.9 at N = 400** and **134.0 at N = 800** on 36 end states per size. Provenance:
`a22_bimodality.mixture_bic` in
[`aggregates/tranche2_runs.json`](aggregates/tranche2_runs.json), with the 700-cycle
rows inherited from [`aggregates/gate1_bifurcation.json`](aggregates/gate1_bifurcation.json).

> **CORRECTED 2026-09-05.** This used to say the rounded range **ΔBIC 120–259** is what the
> manuscript prints. The manuscript no longer prints any ΔBIC value at all — the bimodality
> claim built on it was cut entirely, and the built PDF contains zero occurrences of "BIC".
> Cause: the pooled fit above mixes **6 of 6** driver values at 700 cycles and **3 of 3** at
> 2,000 cycles on both sides of the fold, and pooling unimodal distributions with different
> means makes two components preferable by BIC regardless of whether the substrate is
> bistable — so the pooled ΔBIC is a property of the pooling, not evidence of bimodality.
> Refit **within** fixed driver values, this repository's own mixture family and its own
> ΔBIC > 10 threshold, 12 of 12 seeds per cell: two components win in **3 of 12 cells at 700
> cycles** (A = 0.37 at both N = 400 and N = 800, A = 0.36 at N = 400) and **2 of 6 cells at
> 2,000 cycles** (A = 0.36 at N = 400, A = 0.37 at N = 800). Neither the within-driver counts
> nor a per-driver-value breakdown is stored in any aggregate in this sim directory — no file
> here can bind them — so treat the two counts above as reported, not independently
> reproduced by this record. The **hysteresis** clause below is untouched by any of this: it
> rests on branch-mean separation above 0.30 over 832 committed equilibrium rows at 4 of 4
> sizes (§2 item 3) and never used the mixture fit.

**All-returned recovery cells.** Under a 35% repeat-tie shock, all 12 of 12 seeds
return at both A = 0.48 and A = 0.39. Mean-degree recovery time rises from
**31.504 cycles** in the `frac=0.35,A=0.48` cell to **115.642 cycles** in the
`frac=0.35,A=0.39` cell, a **3.67-fold** increase, printed as 31.5 to 115.6 cycles
and 3.7-fold. Provenance: `a23_small_shock.table` in
[`aggregates/tranche2_runs.json`](aggregates/tranche2_runs.json), 12 seeds per cell.

**Recovery-scenario cell.** In the `ramp_hold_recover` cell at N = 400 and c = 0.5,
the driver crosses the fold, holds, and returns; the repeat-tie graph remains collapsed
in **40 of 40 runs**. Provenance: `a25_scenarios.table` in
[`aggregates/tranche2_runs.json`](aggregates/tranche2_runs.json), 40 runs in the cell.

Empirical fold **A_c ≈ 0.37**, above the duplicate-free mean-field estimate of 0.299, as
expected (re-tying an already-tied pair adds no edge). Figure:
[`figures/figA1_gate1.png`](figures/figA1_gate1.png). Aggregate:
[`aggregates/gate1_bifurcation.json`](aggregates/gate1_bifurcation.json).

The gate's second obligation — that the indicator does not fire in a non-transitioning
control — is discharged by E-A's paired control arm (§3, FPR column), not here.

---

## 3. E-A — lead time under consent-bounded observation

Primary arm (centrality-correlated consent, MNAR), mid proxy-error severity. Δt is the
lead-time advantage in cycles over the equal-FPR lagging comparator, censoring a
non-firing detector at the run length; s.e. from a 2 000-sample bootstrap over the 30
matched-seed replicates.

| c | white Δt | red Δt | systematic Δt | lead detect | lagging detect | achieved edge fraction |
|---|---|---|---|---|---|---|
| 1.00 | **+88.8** ± 12.4 | **+16.6** ± 13.4 | **+94.9** ± 11.7 | 0.83 / 0.37 / 0.87 | 0.73 | 1.000 |
| 0.75 | +65.7 ± 12.7 | +13.4 ± 16.0 | +53.6 ± 14.0 | 0.87 / 0.37 / 0.77 | 0.63 | 0.560 |
| 0.50 | **+48.8** ± 11.6 | +6.6 ± 11.5 | +59.2 ± 10.7 | 0.60 / 0.20 / 0.70 | 0.57 | 0.256 |
| 0.25 | +17.3 ± 13.1 | **−3.3** ± 11.0 | −2.7 ± 8.6 | 0.30 / 0.17 / 0.17 | 0.70 | 0.062 |

**Is there a lead-time advantage? Yes, and it degrades gracefully rather than falling off a
cliff.** Under white noise the leading indicator beats the lagging output metric by ~89
cycles at full coverage and still by **+48.8 cycles at half coverage** — roughly 15% of the
run length, and about 60 cycles before the fold. Discrimination (AUC of the leading rule's
score, transition vs matched control) is 0.89 at c = 1.0, **0.93 at c = 0.5**, and falls to
0.72 at c = 0.25.

**Noise colour dominates amplitude.** Red (AR(1), φ = 0.75) noise costs far more than white
noise at the same amplitude: Δt collapses from +88.8 to **+16.6** at full coverage, and
detection rate from 0.83 to 0.37. This reproduces Dakos et al. (2012)'s finding on the
substrate rather than assuming it, and it is the arm that sets the binding coverage floor.

### 3.1 The Δt = 0 iso-line — the coverage floor

The differentiating deliverable. Linear interpolation of the cell means in c:

| noise colour | MNAR (primary) | MCAR (reference) |
|---|---|---|
| white | no crossing — Δt > 0 at every swept c, so the floor lies **below 0.25** | no crossing |
| red | **c = 0.333** | c = 0.395 |
| systematic | **c = 0.261** | no crossing |

So the coverage floor is **c ≈ 0.26–0.33 under the primary MNAR arm**, set by noise colour,
and below 0.25 for white noise. The ±1 s.e. seed-variation band on the iso-line is wide —
for MNAR/red it runs from 0.864 to off-scale — because Δt near the floor is small relative
to seed variation. Figure: [`figures/fig1_ea_headline_mnar.png`](figures/fig1_ea_headline_mnar.png)
(panels = noise colour; severity on the y-axis; solid Δt = 0 iso-line,
dashed ±1 s.e. band). Confusion time-series with cost axis:
[`figures/fig2_confusion_mnar_white.png`](figures/fig2_confusion_mnar_white.png).

### 3.2 False-alarm rate against the 3.6% external bar

Pooled over all **2 160** control runs, the leading rule's FPR is **8.98%** — above its own
5% nominal and **2.5× Falmagne, Stephenson & Levin (2026)'s 3.6% bar**. By noise colour at
mid severity: white 10.8%, red 4.6%, systematic 9.6%. Per-cell FPR ranges 0.00–0.233
(median 0.067) on 30 control runs each. The lagging comparator's pooled FPR is **0.83%**.

This is an honest failure to hit the external bar, and its cause is identifiable: the
threshold was chosen as the *smallest* k with FPR ≤ 5% on 40 calibration runs per arm, which
overfits the calibration sample. Raising k would trade detection rate for FPR along the
curve already committed in `aggregates/calibration_rule_grids.json`; the paper should either
re-calibrate on a larger null sample or report the operating point honestly with its
achieved FPR, not the nominal one.

**Bounded-oscillation control.** In the `oscillate` cell at N = 400 and c = 0.5,
the formation-only rule alarms in 11 of 40 non-transitioning runs, the both-class rule
alarms in 17 of 40, and the lagging comparator alarms in 3 of 40: **27.5%**, **42.5%**,
and **7.5%**, respectively. The manuscript prints the latter two as 43% and 7.5%.
Thresholds are unchanged from the calibrated rules. Provenance: the
`lead_alarm_rate_recal`, `lead_alarm_rate_bothclass`, and `lag_alarm_rate` fields in
`a25_scenarios.table[scenario=oscillate]` in
[`aggregates/tranche2_runs.json`](aggregates/tranche2_runs.json), 40 runs in the cell.

### 3.3 Nominal c versus achieved observed-edge fraction

Reported separately per the frozen specification. A formation event is recorded only through a dyad whose
**both** endpoints consent, so the achieved fraction is close to c²: nominal 0.25 → achieved
**0.062**, nominal 0.50 → **0.256**, nominal 0.75 → 0.560, nominal 1.00 → 1.000. A consent
regime described as "half the participants" delivers **one quarter of the edges**. That gap
is the reportable number.

**Negative result on the MNAR/MCAR contrast.** The two consent arms give nearly identical
achieved coverage (at c = 0.5: **0.2564** MNAR vs **0.2533** MCAR) and their Δt values differ
by no more than seed noise. The scale of that noise is directly measurable here: at c = 1.0
the two arms are *the same physical configuration* (consent is universal, so the assignment
rule is irrelevant) run on different seed sets, and they differ by 26 cycles (white), 11
(red) and 22 (systematic) — i.e. differences of that size are noise, and every MNAR/MCAR
difference in the table is smaller. **E-A does not resolve a non-random-consent penalty.**

The recalibrated white-noise floor comparison uses **4,000 paired bootstrap
resamples** of runs. The MCAR-minus-MNAR floor gap is **0.080**, with a 95% interval
of **[−0.133, +0.223]**, printed as [−0.13, +0.22]. At c = 0.5, the achieved-edge
comparison pools 300 run rows per consent arm across the three noise colours; the
MNAR-minus-MCAR difference is **0.0012 [−0.0060, +0.0084]**. Provenance:
`mnar_vs_mcar_verdict.gap_bootstrap` and
`mnar_vs_mcar_verdict.achieved_edge_difference_*` in
[`aggregates/remediation_analysis.json`](aggregates/remediation_analysis.json),
derived from [`aggregates/ea_headline_recal.json`](aggregates/ea_headline_recal.json)
at N = 400.

> **CORRECTED 2026-09-05.** This used to print the difference's interval as **[−0.0030,
> +0.0053]**, taking the 300 rows per arm (100 replicates × 3 proxy-noise colours) as 300
> independent draws. They are not: achieved edge fraction depends only on the consent
> assignment, which never reads the noise stream, so all three colours return the identical
> value within a run pair — verified directly in `mnar_vs_mcar_verdict.achieved_edge_fraction_at_c_half`,
> where the per-arm mean and sd already sit at their replicate-level values on n = 300, and
> confirmed at the (arm, replicate) level in **200 of 200** groups. The effective sample size
> is 100 per arm, not 300, a design effect of 3; the interval widens by √3, from
> `mnar_vs_mcar_verdict.achieved_edge_difference_se` **0.002115** to **0.003663**, giving the
> [−0.0060, +0.0084] interval above. The point difference (+0.00117) and the "does not
> support the mechanism" verdict are unchanged — it still straddles zero either way.

---

## 4. Auxiliary arms at (c = 0.5, mid red noise)

### 4.1 Revocation vs equivalent static coverage

Three conditions on matched seeds; `static_equiv` removes *the same agents* from cycle 0, so
the final consent set is identical and only the timing differs. Revocation removes a
centrality-weighted 40% of consenters (mean **80.8 agents**) at cycle 140 and the instrument
re-derives its whole history from the surviving set.

| condition | achieved edge fraction | Δt | lead detect | FPR | aggregate invalidation |
|---|---|---|---|---|---|
| static_full (c = 0.5 throughout) | 0.256 | +9.1 ± 8.7 | 0.35 | 0.067 | 0 by construction |
| revocation (retroactive at cycle 140) | 0.164 | +23.2 ± 9.6 | 0.47 | **0.250** | **0.993** |
| static_equiv (same agents, from cycle 0) | 0.092 | −5.5 ± 7.0 | 0.25 | 0.117 | 0 by construction |

**Semantics leg — holds, but on the invalidation measure, not on timing.** Re-deriving the
already-published pre-revocation indicator series shifts it by **0.993 baseline σ** on
average, and flips its alarm status in **0.93%** of pre-revocation cycles. A static coverage
reduction is exactly zero on this measure by construction: there is nothing to re-derive.
Revocation also **more than doubles the false-alarm rate** (0.250 vs 0.117), because the
retroactive rewrite injects a discontinuity into the observer's own reference statistics.
Both effects are unavailable to a static reduction, so revocation is distinguishable from an
equivalent static reduction — spec §6 kill criterion 1 is **not** triggered.

**Caveat that must travel with the timing number.** The paired Δt difference is +28.8
cycles, 95% CI [+8.6, +47.9], nominally "distinguishable" — but revocation enjoyed *full*
c = 0.5 coverage for the first 140 of 320 cycles, so its cumulative achieved coverage (0.164)
is nearly double static_equiv's (0.092). The timing advantage is largely a coverage
confound and **should not be reported as evidence for the semantics claim.** The
invalidation measure is the clean evidence.

**Record-displacement cells.** Each contrast below is the mean over 60 matched seeds
at c = 0.5 under mid-severity red noise. For each seed, the level is the mean
difference between the two standardized indicator series, shape is the standard
deviation of their de-meaned difference, and mean absolute displacement is the mean
absolute pointwise difference over the published pre-revocation series.

| contrast | mean absolute displacement | level shift | shape s.d. |
|---|---:|---:|---:|
| published vs re-derived | **0.993** | **0.024** | 0.840 |
| published vs from-day-one static-equivalent | **1.453** | 0.222 | 1.254 |
| re-derived vs from-day-one static-equivalent | **1.462** | −0.198 | 1.322 |

The units are baseline standard deviations. Provenance:
`a24_revocation_decomposition` and `a24_per_run_rows` in
[`aggregates/tranche2_runs.json`](aggregates/tranche2_runs.json), 60 matched seeds.

### 4.2 The Ensign fix (inverse-propensity correction)

Propensity estimated from observed consent status alone: `logistic(consent ~ 1 + degree in
the observed subgraph)`. The observer sees an observed program's consented members but not
its non-consenting ones, so it infers the team size from the estimated propensity
(`m̂ = k/p̂`), undoes the partial-observation shrinkage there, and inverse-weights by the
probability a team of `m̂` would have been observed at all.

| condition | \|bias\| raw | \|bias\| after IPW | bias reduction (95% CI) | mean agents mislabelled |
|---|---|---|---|---|
| static consent | 0.986 | 0.285 | **+0.702** [+0.673, +0.730] | 0 |
| under revocation | 1.166 | 0.865 | **+0.302** [+0.274, +0.331] | **80.8** |

**Where it fails, quantified.** Under static consent the correction removes **71%** of the
MNAR level bias. Under revocation it removes only **26%** — the reduction falls to
**43% of its static value**, and the raw bias is *larger* to begin with. **Two estimators appear
in that sentence and both are named here, because the paper quotes the percentages and this record
quotes the fractions.** The 71% and 26% are one minus the ratio of pooled absolute biases over 60
runs per arm, 71.12% and 25.86%, with 95% bootstrap intervals of 69–73 and 24–28 over 2,000
resamples of runs. The table's +0.702 and +0.302, and the 43% ratio above, are the per-run mean of
that same fraction. Under the pooled estimator the revocation-over-static ratio is 36% rather than
43%. The mechanism is
exactly the one predicted: revocation erases the consent record, so 80.8 agents per run
enter the fit labelled as never-consenting; the estimated consent rate drops from
**p̂ = 0.505 to p̂ = 0.303**, and that post-revocation rate is then applied to
pre-revocation records generated at the higher rate. The estimator has no way to detect
this from the record it is given. Aggregate:
[`aggregates/aux_arms.json`](aggregates/aux_arms.json).

---

## 5. E-B — the joint (c × g) sweep

Statistic: Spearman rank correlation between the published (20-cycle trailing mean) observed
indicator and the latent adaptive capacity A(t), over the pre-fold stretch (cycles 90 to the
A(t) = 0.37 crossing). More negative = better tracking. 120 matched-seed replicates per
cell; coupling-severed control on identical seeds with the response drawn from the same RNG
stream and discarded.

| g | c = 0.25 | c = 0.50 | c = 1.00 |
|---|---|---|---|
| 0 (reference) | **−0.449** [−0.515, −0.376] | −0.561 | **−0.602** [−0.664, −0.536] |
| mid, coupled | −0.416 | −0.548 | −0.627 |
| mid, severed | −0.395 | −0.555 | −0.595 |
| high, coupled | −0.408 | −0.554 | −0.581 |
| high, severed | −0.395 | −0.555 | −0.595 |

**The coverage main effect is resolved:** tracking quality falls monotonically with
coverage, and the c = 1.0 and c = 0.25 confidence intervals do not overlap. Detection rate
falls with it, 0.81 → 0.73 → 0.36.

**The coupling effect and the (c × g) interaction are NOT resolved.** Every one of the six
paired coupling-loss estimates has a 95% CI containing zero (largest |loss| = 0.053, CI
half-widths ≈ 0.07), and the two coupling levels point in **opposite directions**: at
g = mid the loss is larger at full coverage (reads as "compound"), at g = high it is larger
at low coverage (reads as "protect"). Neither reading is supported. Figure:
[`figures/fig3_eb_joint.png`](figures/fig3_eb_joint.png).

Mechanism note on why coupling has so little purchase: agents can game *which ties get
recorded*, but the indicator is a **fluctuation** statistic of the structure they sit inside,
and the proxy-response channel does not move it much. The Bol-style attrition channel
(mean **483** withdrawals at g = mid, **933** at g = high, all replaced so the ecosystem does
not empty) pushes the other way — the churn it injects destroys repeat ties and slightly
*raises* the indicator's alarm rate. The two channels partly cancel. Whether that is a
finding about variance-family indicators or an under-powered coupling implementation is not
decidable from these runs.

---

## 6. Verification

> **CORRECTED 2026-09-04.** This section previously said the verifier "enforces six classes of
> check". The source defines **fourteen**, C1 through C14. And "every number bound" (§1's
> phrasing) conflated three different guarantees, broken out below with denominators.

[`verify_reported_numbers.py`](verify_reported_numbers.py) enforces fourteen classes of check
(C1–C14): a cell summary recomputed from the per-run rows in the same aggregate (C1, C8, C11,
C12, C14); a count or set property of seeds, e.g. no scored run on a calibration seed (C2, C3);
byte identity of a rule file, e.g. the detection rule embedded in each aggregate against
`configs/detection_rule.json` (C4, C7, C11, C13); a stored flag compared to the stored numbers
it summarises (C1 `meets_external_bar`, C9 `excludes_zero`, C10 `above_external_bar_ci_lower`,
C12 the MNAR verdict string); every number in this file bound to an aggregate path through
[`claims.json`](claims.json) *and* present verbatim in this file's prose; and a random sample of
scored runs re-executed from its recorded seed reproducing its per-run scalars exactly.

**No check performs an independent Wilson interval, bootstrap standard error, p-value,
chi-square or design-effect check.** This is not the same as saying the verifier computes no
interval or p-value anywhere — it calls `recal_analysis.compute()` (C10) and
`remediation_analysis.compute()` (C12), and both of those modules compute Wilson intervals,
bootstrap quantiles and a heterogeneity p-value. What C10 and C12 do with that output is diff it
against the committed file, which proves the analysis module is deterministic and tests nothing
about whether its unit of analysis is independent. Re-running the same implementation cannot
catch a wrong independence assumption baked into that implementation, which is why the ninefold
design-effect defect in §8.10–§8.12 passed a green verifier.

The claim map contains **128 claims** whose `results_file` is `RESULTS.md`. Of those: **61 of
128** are independently recomputed from per-run rows; **27 of 128** bind fields of
`recal_analysis.json` that C10 reproduces by re-running that same module (self-reproduction, not
independent recomputation); **40 of 128** bind a stored field no check recomputes at all,
including the headline's bootstrap standard error, both `auc_discrimination` claims, and all
seven Gate 1 test-statistic claims. The verifier opens **13 of 104** aggregate files in
`sim/aggregates/`; the other **91 of 104** — including every `graph_detector_*`,
`three_observable_*` and fine-grid aggregate behind §8.10–§8.12 — are never opened, and none of
the 13 postdates 2026-08-25, so no campaign run in September is covered by this verifier at all.
Provenance: the `claims` array and `results_file` fields in [`claims.json`](claims.json),
counted over 128 entries; the 61/27/40 partition and the 13/104 file-open count from the
round-2 adversarial-review recomputation of 2026-09-04 (tracked outside this repository).

**The manuscript's "3,873 automated verifier checks" is this script's own printed total, not
a count of the 128 claims above.** Running `verify_reported_numbers.py` end to end prints
`verify_reported_numbers: 3873 checks, 0 failures` as its last line — one increment per
assertion evaluated across the fourteen classes C1–C14, so the same run can and does contribute
many checks to a single bound claim. Provenance: the script's own stdout, reproduced 2026-09-05.

---

## 7. What looks weak — for the kill/fallback judgment

Ordered by how much they threaten the paper, not by how easy they are to fix.

1. **C-III has no resolved evidence.** E-B's coupling main effect and its (c × g)
   interaction are both null at 120 paired replicates, and the two g levels disagree in
   sign. Conflict 2 forbids a pointer-only E-B while §1 carries C-III, so the paper as
   frozen currently asserts a conjunct its own experiment does not support. Three live
   options: (a) strengthen the coupling channel so it actually bites a variance-family
   indicator and re-run; (b) report the null as the result and re-word C-III as a bounded
   negative ("restricted observation neither protects nor compounds measurably at this
   coupling strength"); (c) drop C-III and take the pointer-only E-B that conflict 2 then
   permits. **This is the call I cannot make.**
2. **FPR misses the external bar by 2.5×** (8.98% pooled vs 3.6%). Fixable — re-calibrate
   the threshold on a much larger null sample rather than 40 runs — but it must be fixed
   or reported at its achieved value, because the frozen specification names that bar in the paper.
3. **The lagging comparator's insensitivity is a modelling assumption, and it is the single
   most attackable choice in the design.** Program output depends on team cohesion through
   a saturating function plus a 6-cycle close delay and a 10-cycle reporting period, and
   cohesion is nearly flat on the upper branch (0.125 at A = 0.48 vs 0.118 at A = 0.40),
   only collapsing at the fold. That is §3's premise — coordination strain shows in *how*
   teams form long before it shows in *what* programs ship — implemented, not tested. One
   could reasonably say the lead time measures the assumption. Mitigations already
   in place: the comparator gets the same statistical strictness, the earliest-detecting
   admissible rule, and a smoothed (lower-variance) series. The paper must state the
   assumption explicitly rather than let it be discovered.
4. **The MNAR arm buys nothing measurable.** Achieved coverage is indistinguishable from
   MCAR (0.2564 vs 0.2533 at c = 0.5) and no Δt difference exceeds seed noise. The
   non-randomness story is currently carried by the *mechanism* and the literature (de Man
   2023), not by a simulation result. Either the consent–centrality coupling needs to be
   made structurally consequential (e.g. consent correlated with something the observable
   actually depends on) or §4.3 should stop implying the sim demonstrates an MNAR penalty.
5. **Red noise nearly kills the instrument even at full coverage** (Δt +16.6, detection
   0.37). Honest, and arguably the most useful engineering number in the paper, but it
   means the headline "+88.8 cycles" is a white-noise number and must never be quoted
   without its colour.
6. **The revocation timing contrast is coverage-confounded** (§4.1). The semantics claim
   rests entirely on the aggregate-invalidation measure and the doubled FPR. That is
   defensible but narrower than the frozen spec's phrasing implies.
7. **Seed variation is large relative to several reported effects.** Δt s.e. is 9–16 cycles
   on 30 replicates; the c = 1.0 duplicate-configuration check shows 11–26 cycles of pure
   seed noise. Cells at c = 0.25 and the whole red-noise column are individually
   under-powered. Raising E-A to ~100 replicates costs about an hour and would firm up the
   iso-line, which is the paper's differentiating deliverable and currently has a very wide
   band.
8. **A(t) is a linear ramp and degradation is the only scenario shape.** No shock,
   oscillatory or partially-recovering scenario was run, so nothing here speaks to
   Looker et al. (2026)'s non-normal false-positive generator beyond the flat control arm.
9. **Single substrate size for all reported experiments.** Gate 1 sweeps N to 1600, but
   E-A, E-B and the auxiliaries are all at N = 400. The coverage floor's dependence on
   system size is unmeasured, and it is the number most likely to be asked about.

**Status of these nine after the follow-up campaign (§8):** item 1 → **§8.8** (addressed:
E-B run to ≥ 500 matched pairs, verdict by a pre-declared rule); item 2 → **§8.2–8.3**
(addressed and partly reversed: the threshold was worse than reported, the recalibrated
rate is 4.39% and still above the 3.6% bar); item 5 → **§8.4** (worse: red noise now has
*no* lead-time advantage at any coverage); item 6 → **§8.7** (the revocation timing
contrast is no longer even nominally distinguishable — the semantics leg now rests
entirely on aggregate invalidation); item 7 → **§8.5–8.6** (addressed: 100 seeds/cell on
the iso-line neighbourhood, band now 0.13 wide in c). Items 3, 4, 8 and 9 are untouched
and stand as written.

---

## 8. Follow-up campaign — 2026-08-05

*Three computation jobs, no design changes. Logged as spec §9 amendments **A-8**
(threshold recalibrated on a dedicated 600-run null), **A-9** (one added E-A grid point,
c = 0.333, at the headline severity, for iso-line localisation; replication raised to 100
seeds/cell on the iso-line neighbourhood) and **A-10** (E-B raised to ≥ 500 matched pairs
with a decision rule declared before the runs). The substrate, the grid semantics, the
statistics, the claim set and the kill criteria are untouched. Campaign 1's rule file and
aggregates are retained, so the two operating points can be compared directly — which
§8.5 does.*

### 8.1 What ran

| stage | script | scale | wall time |
|---|---|---|---|
| large-null recalibration | [`calibrate_recal.py`](calibrate_recal.py) | 600 control runs on fresh seeds 950000+, plus campaign 1's 40 calibration degradation runs | 188 s |
| **E-A** re-evaluation + iso-line replication | [`run_ea_headline.py`](run_ea_headline.py) `--rule detection_rule_recal.json` | **7 200 runs** (5 c × 3 colours × 3 severities × 2 arms × 2 scenarios; 100 reps on c ∈ {0.25, 0.333, 0.5} at mid severity, 30 elsewhere) | 2 020 s |
| **auxiliary arms** re-evaluation | [`run_aux_arms.py`](run_aux_arms.py) `--rule …` | 360 + 120 runs, unchanged design | 625 s |
| **E-B** at ≥ 500 matched pairs | [`run_eb.py`](run_eb.py) `--rule … --reps 500` | **9 000 runs** (3 c × 3 g × 2 modes × 500 reps) | 3 083 s |
| cross-campaign analysis | [`recal_analysis.py`](recal_analysis.py) | no runs — derived from two committed aggregates | seconds |
| figures | [`figures.py`](figures.py) | fig1–fig3 re-emitted; fig4, fig5, fig6 new | seconds |

Every E-A and auxiliary run in this campaign is the *same run* campaign 1 scored: a seed
depends on (c, arm, replicate) only, so re-running reproduces campaign 1's substrate
bit-for-bit and only the scoring rule differs. The verifier checks that directly — 4 320
runs overlap between the two E-A aggregates and **0** of them disagree on the seed.

### 8.2 The recalibration — the campaign-1 threshold's real false-alarm rate was 11.0%

The rule **form** is untouched (statistics, window, signal class, mixing, combination,
sustain length, reference period, and the lagging metric's definition all as frozen); only
the two thresholds were re-picked, on **600 control runs at seed block 950000+** — fresh,
disjoint from campaign 1's 900000–900039 calibration block and from every scored seed.
Calibration stays at c = 1.0, white/mid, MNAR: the instrument may not re-tune itself per
consent regime.

Measured on that null, campaign 1's threshold **k = 5.0 has a false-alarm rate of
11.0%** (66 of 600), not the 5.0% it was selected at on 40 runs. That is the overfit,
quantified: the diagnosis in §7 item 2 was right about the cause and understated the size.

Re-picking per the same pre-declared procedure with the target tightened to ≤ 3.5% (just
under the external bar) gives:

| | campaign 1 | recalibrated | selection |
|---|---|---|---|
| leading threshold k | 5.0 | **6.5** | smallest k with null FPR ≤ 3.5% → **3.50%** (21/600) |
| lagging k_sd / sustain | 3.5 σ / 5 | **3.5 σ / 3** | earliest-detecting rule with null FPR ≤ 3.5% → 2.67% (16/600), median calibration detection 276 → **262** cycles |

The comparator therefore got *more* generous, not less: the same σ threshold sustained for
3 cycles instead of 5 fires 14 cycles earlier at the median. That is the frozen A-2
discipline (both detectors at one target, comparator takes every admissible advantage)
applied to a null large enough to resolve 3.5%. Figure:
[`figures/fig5_calibration_null.png`](figures/fig5_calibration_null.png). Aggregate:
[`aggregates/calibration_null_large.json`](aggregates/calibration_null_large.json), whose
committed per-run rows carry the sufficient statistic for every point on both FPR curves,
so the verifier re-derives both selections rather than trusting them.

### 8.3 False alarms against both bars — halved, still above

Pooled over **3 600** evaluation control runs, the recalibrated leading rule's FPR is
**4.39%**, Wilson 95% CI **[3.77%, 5.11%]** at the run unit, corrected to
**[3.45%, 5.56%]** at the group unit below — down from campaign 1's **8.98%**
[7.85%, 10.26%] on 2 160 controls. At the run unit the CI lower bound sat **above** both
the 3.5% declared target and the **3.6% external bar** (Falmagne, Stephenson & Levin
2026); **at the corrected group unit it does not**, and the correction below carries the
arithmetic. The point estimate misses both bars either way, and the honest headline is:
**recalibration halves the excess and does not eliminate it.** What does NOT survive the
correction is the stronger reading this sentence used to carry, that the miss is real
rather than seed noise: that was an exclusion argument, and the corrected interval
contains both bars.

> **CORRECTED 2026-09-05.** This used to print the Wilson interval as **[3.77%, 5.11%]**,
> computed over the 3,600 control runs as independent trials. They are not: those 3,600 are
> **400 independent control groups** (one per (coverage, arm, replicate) triple) each scored
> under nine correlated noise-colour/severity copies, the same nine-fold pseudo-replication
> already named for the graph detectors in §8.10–§8.12. Recomputed at the group unit over
> the 400 groups, the interval is **[3.45%, 5.56%]**. The point estimate, **4.39%**, is
> unchanged — it is a mean, not a variance statistic — and it remains above both bars.
>
> **THE LOWER-BOUND READING IN THE SENTENCE ABOVE DOES NOT SURVIVE THIS CORRECTION, and
> saying so is the point of writing the correction at all.** That sentence's claim is
> specifically that the CI *lower bound* clears both bars. At the group unit the lower
> bound is **3.45%**, which sits BELOW the 3.5% declared target and below the 3.6%
> external bar. So the honest reading is narrower than the one this record has been
> carrying: the point estimate misses both bars, and the interval no longer excludes
> the declared target. "The miss is real and not seed noise" holds for the estimate
> and not for the bound.
>
> An earlier draft of this very correction block said the verdict "stands with a wider,
> more honest interval around it", which contradicts the interval it had just printed.
> Recorded rather than quietly replaced: a correction that widens an interval has to be
> re-read against every inference the old width supported, and this one was not.

Where the residual sits, at mid severity: by noise colour, white **6.25%**, systematic
3.61%, red 3.19%; by coverage, c = 0.25 → **6.33%**, c = 0.333 → 4.50%, c = 0.50 → 2.50%,
c = 0.75 → 2.22%, c = 1.00 → 5.56%. Two readings, and the second is a result rather than a
defect:

1. Per-cell FPRs on 30 control runs cannot resolve a 3.5% rate. At c = 1.0 the MNAR and
   MCAR arms are *the same physical configuration* on disjoint seed sets, and their FPRs
   differ by **0.100** (white: 0.133 vs 0.033; systematic: 0.133 vs 0.033) — pure seed
   noise of 10 percentage points. Only the pooled number is quotable.
2. The **only coverage level whose excess is individually resolvable** on these sample
   sizes is the lowest: c = 0.25 at **6.33%** on 600 controls, Wilson 95%
   **[4.65%, 8.57%]**. c = 0.333 (4.50%, [3.11%, 6.47%]) and c = 1.00 (5.56% on 180
   controls, [3.05%, 9.92%]) both include the bar.

   > **CORRECTED 2026-09-04.** This used to read the series as a directional cost of
   > restricted coverage — "a rule calibrated at full coverage runs hotter than its nominal
   > rate when coverage is restricted," a *second*, independent cost of consent-bounded
   > observation. The full series is **6.33%, 4.50%, 2.50%, 2.22%, 5.56%** at c = 0.25,
   > 0.333, 0.5, 0.75, 1.00 — non-monotone, with two of five levels below the 3.5% target and
   > the calibration level itself (c = 1.0) above it. At the (arm, replicate) group level,
   > heterogeneity across the five levels is **p = 0.010 over runs and p = 0.15 over groups**
   > (chi-square 13.2 versus 6.8 on 4 df). The correct statement: **the false-alarm rate
   > varies with coverage in a way these run counts do not resolve** — not a directional cost
   > the operator cannot correct.

The lagging comparator now runs at a pooled **4.92%**, so the two detectors are genuinely
matched (4.39% vs 4.92%). **In campaign 1 they were not:** 8.98% against 0.83%, a
comparator operating at a tenth of the leading rule's achieved false-alarm rate. Campaign
1's lead times were measured against a much stricter — therefore much later — comparator,
which is a large part of why they were large. §8.5 quantifies how much.

### 8.4 E-A lead time under the recalibrated rule

Primary arm (MNAR), mid severity. Δt is the lead-time advantage in cycles over the
equal-FPR lagging comparator, censoring a non-firing detector at the run length; ± is a
bootstrap s.e. over matched-seed replicates.

| c | n/cell | white Δt | red Δt | systematic Δt | lead detect (w/r/s) | lag detect | AUC (w/r/s) |
|---|---|---|---|---|---|---|---|
| 1.00 | 30 | **+58.1** ± 14.3 | **−9.8** ± 12.3 | **+60.4** ± 12.3 | 0.70 / 0.23 / 0.83 | 0.73 | 0.89 / 0.84 / 0.89 |
| 0.75 | 30 | +12.3 ± 12.8 | −22.9 ± 13.3 | +22.9 ± 15.1 | 0.57 / 0.20 / 0.63 | 0.83 | 0.90 / 0.79 / 0.93 |
| 0.50 | 100 | **+9.2** ± 6.8 | −25.5 ± 5.6 | +7.6 ± 6.9 | 0.45 / 0.12 / 0.42 | 0.77 | 0.91 / 0.81 / 0.91 |
| 0.333 | 100 | −8.3 ± 7.1 | −16.8 ± 6.9 | −12.3 ± 7.4 | 0.31 / 0.23 / 0.27 | 0.83 | 0.83 / 0.72 / 0.79 |
| 0.25 | 100 | −22.0 ± 6.6 | −32.7 ± 4.7 | −25.9 ± 5.8 | 0.17 / 0.07 / 0.13 | 0.73 | 0.72 / 0.67 / 0.73 |

> **CORRECTED 2026-09-05.** The manuscript's headline lead time is no longer this table's
> c = 1.0/white cell at n = 30 (**+58.1 ± 14.3**). It now reads **"53 ± 11 cycles (bootstrap
> s.e., n = 40)"**, taken from a 40-replicate extension of the same (c = 1.0, MNAR, white,
> mid) cell — these 30 seeds are nested inside that 40, nothing was swapped. The n = 40
> figure is **+53.125, bootstrap s.e. 11.3960** in `comparator_gain_leg.per_variant["knee0.30_sig0.100"]`
> in [`aggregates/remediation_analysis.json`](aggregates/remediation_analysis.json) (the
> frozen comparator-gain reference point, `comp_ref.knee0.30_sig0.100 = 100.0`) — the same
> row already quoted as "formation-only, c = 1.00: +53.1" in §9.3's table and as
> "linear_w30 (frozen), Δt @ c=1: +53.1" in §10.1's. The reported sample s.d. of 71.6076
> for that cell is not itself stored in this aggregate (only the bootstrap s.e. is); it is
> carried here as the manuscript's own figure, not independently re-derived. **This table's
> n = 30 row is kept as the campaign-2 record and as the input to §8.5's decomposition
> (which is defined over the 72 cells the two campaigns share at 30 replicates each) — it is
> not the number to quote as the paper's headline any more.**

**The lead-time advantage survives, much smaller, and only above roughly half coverage.**
Under white noise the leading indicator beats the lagging metric by **+58.1 cycles at full
coverage** (campaign 1: +88.8) and by **+9.2 cycles at half coverage** (campaign 1:
+48.8). The systematic-bias column behaves the same way. **Red noise now has no advantage
at any swept coverage at the headline mid severity** — Δt is negative in all five cells,
from −9.8 at c = 1.0 to −32.7 at c = 0.25 — which sharpens §7 item 5 from "nearly kills the
instrument" to "kills it at this amplitude".

**Severity is not incidental, and the headline severity is not the worst case.** At **low**
proxy-error severity red noise recovers to **+44.8** cycles at full coverage and **+14.1**
at c = 0.5; at **high** severity even white noise is negative below full coverage
(**+3.2** at c = 1.0, **−14.6** at c = 0.5). The full (c × colour × severity) surface is in
fig1. The one-line summary is that the instrument needs *either* high coverage *or*
low-amplitude uncorrelated proxy error, and buys nothing when both are unfavourable — which
is a sharper engineering statement than campaign 1 could make, and it is the sentence a
practitioner will actually use.

**Discrimination is nearly untouched:** AUC of the leading score (transition vs matched
control) is **0.89 at c = 1.0** and **0.91 at c = 0.5** under white noise, against 0.89 and
0.93 in campaign 1. That separation matters for the paper's story: what the honest
threshold costs is *timing at a matched false-alarm rate*, not the indicator's ability to
tell a transitioning system from a stable one. The signal is there; converting it into
early warning at an operating point a regulator would accept is what consent-bounded
coverage makes expensive.

Figures: [`figures/fig1_ea_headline_mnar_recal.png`](figures/fig1_ea_headline_mnar_recal.png)
(Ferretti form, panels = noise colour), MCAR reference in
[`figures/fig1_ea_headline_mcar_recal.png`](figures/fig1_ea_headline_mcar_recal.png),
confusion time-series in
[`figures/fig2_confusion_mnar_white_recal.png`](figures/fig2_confusion_mnar_white_recal.png)
and [`figures/fig2_confusion_mnar_red_recal.png`](figures/fig2_confusion_mnar_red_recal.png).
Aggregate: [`aggregates/ea_headline_recal.json`](aggregates/ea_headline_recal.json).

### 8.5 Where the lead time went — decomposition, not hand-waving

Two thresholds moved at once and both push Δt down, so the drop is attributed exactly. For
every run, Δt = (lagging alarm cycle) − (leading alarm cycle), each censored at the run
length, so on matched seeds

&nbsp;&nbsp;&nbsp;&nbsp;Δt(recalibrated) − Δt(campaign 1) = **comparator effect** + **leading-threshold effect**

with the comparator effect = mean change in the lagging alarm cycle, and the
leading-threshold effect = −(mean change in the leading alarm cycle). Over the **72** cells
the two campaigns share (30 replicates each, the overlap), the mean change is
**−34.15 cycles**, of which **−12.95** is the comparator getting earlier and **−21.20** is
the stricter alarm threshold. The comparator accounts for **37.9%** of the loss, and the
identity closes to 7.1 × 10⁻¹⁵. Aggregate:
[`aggregates/recal_analysis.json`](aggregates/recal_analysis.json).

So roughly **two-fifths of the shrinkage is the price of an honestly matched comparator**
and three-fifths the price of an honestly calibrated threshold. Neither is a change of
design; both are what campaign 1's numbers would have been had its null been large enough.

A second effect in the same aggregate justifies A-9's added replication: on the 30 shared
replicates the white/MNAR c = 0.5 cell reads **+23.3** cycles, and on all 100 replicates it
reads **+9.2**. Thirty seeds cannot site a cell relative to zero, which is exactly the §7
item 7 complaint.

### 8.6 The iso-line, firmed

With 100 seeds/cell on c ∈ {0.25, 0.333, 0.5} and a measured point at campaign 1's
interpolated crossing, the Δt = 0 coverage floor is:

| noise colour | MNAR (primary) | ±1 s.e. band | MCAR (reference) | campaign 1 (MNAR) |
|---|---|---|---|---|
| white | **c = 0.412** | 0.345 – 0.477 | c = 0.502 | no crossing (floor below 0.25) |
| systematic | **c = 0.436** | 0.375 – 0.494 | c = 0.473 | c = 0.261 |
| red | **no crossing — Δt ≤ 0 at every swept c** | +1 s.e. crossing at 0.948 | no crossing | c = 0.333 |

**The coverage floor is c ≈ 0.41–0.44 under the primary MNAR arm, and its band is now
0.13 wide in c instead of running off-scale.** That is the differentiating deliverable
at decision grade, and it is a much less comfortable number than campaign 1's:
the floor is not a quarter of participants, it is **more than four in ten** — and under
autocorrelated proxy error no coverage in the swept range buys any warning at all.
Nominal-versus-achieved (§3.3) makes it starker: c = 0.412 nominal is an achieved
observed-edge fraction of roughly **0.17**, because a formation event is recorded only when
both endpoints consent. Figure:
[`figures/fig4_isoline_recal.png`](figures/fig4_isoline_recal.png).

**Estimator and horizon cells.** These floors use the N = 400, white-noise MNAR arm.
At τ = 240, ΔRMTA is −3.83 cycles at c = 0.25 and +6.15 cycles at c = 0.333;
linear interpolation gives c* = **0.2819**, printed as the lower endpoint 0.28.
Those two cells contain 100 runs each. Across τ ∈ {160, 200, 240, 280}, the same
interpolation gives c* = 0.303, 0.286, 0.282, and 0.296. At τ = 320, interpolation
of lead-detection probability gives **c ≈ 0.61**, and interpolation of the probability
that the leading detector fires first conditional on either detector firing gives
**c ≈ 0.69**. The τ = 320 ΔRMTA floor uses 100 runs per cell at c ∈ {0.25, 0.333,
0.5} and 30 runs per cell at c ∈ {0.75, 1.0}; 4,000 run-resampling bootstrap draws
give **0.417 [0.315, 0.595]**. Provenance: `floor_resiting.mnar|white` and
`floor_bootstrap.white.mnar|dt` in
[`aggregates/remediation_analysis.json`](aggregates/remediation_analysis.json), derived
from [`aggregates/ea_headline_recal.json`](aggregates/ea_headline_recal.json).

**Comparator-gain cells.** Nine declared combinations of `output_knee` and
`output_sigma`, each evaluated on 40 runs at each of five coverage levels at N = 400,
give τ = 320 ΔRMTA crossings from **0.328 to 0.413**, printed as a comparator-gain
band of **0.33–0.41**. Provenance: `comparator_gain_leg.rows` and
`comparator_gain_leg.by_variant` in
[`aggregates/remediation_analysis.json`](aggregates/remediation_analysis.json), derived
from [`aggregates/remediation_runs.json`](aggregates/remediation_runs.json).

### 8.7 Auxiliary arms under the recalibrated rule

Same design, same seeds, same (c = 0.5, mid red noise) scope bound (D16); only the scoring
rule changed. Red noise at c = 0.5 has no lead-time advantage at all now (§8.4), so every
Δt in this arm is negative, and the timing contrast is the only thing that moves:

| condition | achieved edge fraction | Δt | lead detect | FPR | aggregate invalidation |
|---|---|---|---|---|---|
| static_full (c = 0.5 throughout) | 0.256 | −37.6 ± 7.5 | 0.13 | 0.033 | 0 by construction |
| revocation (retroactive at cycle 140) | 0.164 | −23.8 ± 9.3 | 0.22 | 0.083 | **0.993** |
| static_equiv (same agents, from cycle 0) | 0.092 | −34.2 ± 8.0 | 0.17 | 0.067 | 0 by construction |

Three changes from §4.1, all in the direction of *less* claim:

1. **The timing leg is no longer distinguishable.** The paired Δt difference
   (revocation − static_equiv) is **+10.5 cycles, 95% CI [−6.6, +27.6]** — it straddles
   zero, where campaign 1 had +28.8 [+8.6, +47.9]. §4.1 already said this number should not
   be reported as evidence because it is coverage-confounded; it is now not evidence at
   all, which removes the temptation.
2. **The doubled false-alarm rate does not survive.** Revocation 0.083 vs static_equiv
   0.067 (campaign 1: 0.250 vs 0.117). That effect was mostly the overfitted threshold
   sitting where the retroactive rewrite's discontinuity could trip it.
3. **The semantics leg holds, unchanged, on the invalidation measure.** Re-deriving the
   already-published pre-revocation series still shifts it by **0.993 baseline σ** and
   flips its alarm status in 0.67% of pre-revocation cycles, against exactly zero for a
   static reduction by construction. That quantity does not depend on the threshold, and it
   is now the *only* evidence for the leg — spec §6 kill criterion 1 is still not
   triggered, but C-II's semantics sub-claim rests on a single measure and the paper must
   say so.

**The Ensign-fix arm is threshold-independent and unchanged:** the correction removes
**+0.702** [+0.673, +0.730] of the MNAR level bias under static consent and only
**+0.302** [+0.274, +0.331] under revocation — **43% of its static value** — with **80.8**
agents per run entering the propensity fit mislabelled as never-consenting and p̂ falling
from 0.505 to 0.303. Aggregate:
[`aggregates/aux_arms_recal.json`](aggregates/aux_arms_recal.json).

### 8.8 E-B at 500 matched pairs per cell — the C-III evidence gap

Design unchanged in every respect (3 × 3 grid, N = 400, the same Spearman statistic on the
published 20-cycle trailing mean against latent A(t) over the pre-fold stretch, the same
proxy-response and Bol-style attrition channels, the same coupling-severed control on
identical seeds with the response drawn from the same RNG stream and discarded). Only the
replicate count rose: **500 matched pairs per cell, 9 000 runs, 3 083 s**.

**The coverage main effect is resolved and now tight.** At g = 0 the tracking correlation
is **−0.406** [−0.443, −0.371] at c = 0.25, −0.565 at c = 0.50 and **−0.605**
[−0.633, −0.577] at c = 1.00 — monotone, with non-overlapping intervals at the ends and
s.e. ≈ 0.015 instead of 0.035. Leading-detection rate over the same runs rises
**0.194 → 0.462 → 0.666**.

**The pre-declared decision rule returns UNRESOLVED**, and the aggregate applies it
mechanically rather than leaving it to prose. The paired coupling-loss
(|corr|<sub>severed</sub> − |corr|<sub>coupled</sub>, pooled over the three coverage levels
at fixed g, 1 500 pairs each):

| g | paired coupling loss | 95% CI | excludes 0? | MDE at 80% power |
|---|---|---|---|---|
| mid | **−0.0083** | [−0.0264, +0.0110] | no | 0.0265 |
| high | **−0.0223** | [−0.0414, −0.0030] | **yes** | 0.0268 |

Since one of the two CIs still straddles zero, the rule's third branch applies: **report
unresolved at this substrate with the power analysis.** Three things must travel with that
verdict, and all three are stronger than campaign 1's bare null:

1. **The one interval that does exclude zero points the *opposite* way to the "gaming
   degrades the indicator" reading.** At g = high the loss is **negative** — coupling makes
   the published indicator track the latent state *better*, by 0.022 in |corr| — and
   leading-detection rate rises with coupling too (0.666 → 0.698 at c = 1.0,
   0.194 → 0.236 at c = 0.25). The mechanism is the one §5 already suspected, now with the
   power to see it: the attrition channel withdraws **933** participants per run at
   g = high, all replaced, and that churn destroys accumulated repeat ties and *raises* the
   fluctuation statistic the indicator is built from. The proxy-response channel — agents
   preferentially recruiting partners already in the record — has no purchase on a variance
   statistic. Whatever else E-B shows, it does not show restricted observation protecting an
   indicator from gaming, because at this substrate gaming does not damage this indicator.
2. **The interaction is bounded, not merely unmeasured.** The c-dependence of the coupling
   loss — the quantity C-III is about — is **−0.0047** at g = mid (95% CI
   [−0.051, +0.041]) and **−0.0148** at g = high ([−0.062, +0.032]), against a minimum
   detectable range of **0.066** and **0.068** at 80% power. So the honest statement is:
   *any* protect-or-compound interaction at these coupling strengths is smaller than
   **0.07** tracking-correlation units, on a scale where the whole coverage main effect is
   0.20. That is a quantified bounded negative, and it is reportable.
3. **What it would take to resolve it.** The observed pooled g = mid effect would need
   **15 257** matched pairs to reach 80% power — 34× the runs, about 29 hours of the same
   compute, for an effect of 0.008. Resolving the *interaction* is worse still. The
   conclusion is not "run more seeds"; it is that a variance-family indicator on this
   substrate is close to insensitive to this coupling channel, and that a coupling
   implementation which actually bites a fluctuation statistic would be a **design change**
   — out of scope for this campaign and a decision for the paper, not the executor.

**One more number the extra replication bought, worth reporting on its own.** The cell
statistic is the mean *signed* correlation; the paired statistic is a difference of |corr|
per run. They differ because a fraction of runs **anti-track** — the published indicator
moves the *wrong way* against the latent state — and that fraction is a direct cost of
restricted coverage: **18.2%** of runs at c = 0.25, 8.0% at c = 0.50, **6.4%** at
c = 1.00 (g = 0, 500 runs per cell). At quarter coverage nearly one published indicator in
five would have told its reader the ecosystem was improving while adaptive capacity fell.

Figures: [`figures/fig3_eb_joint_n500.png`](figures/fig3_eb_joint_n500.png) (the (c × g)
plane with the severed control as reference contour) and
[`figures/fig6_eb_coupling_forest.png`](figures/fig6_eb_coupling_forest.png) (every paired
coupling-loss CI, per cell and pooled). Aggregates:
[`aggregates/eb_joint_sweep_n500.json`](aggregates/eb_joint_sweep_n500.json) (verdict in its
`decision` block) and [`aggregates/recal_analysis.json`](aggregates/recal_analysis.json)
(the interaction bound and the anti-tracking fractions).

### 8.9 What the follow-up campaign changes for the paper

Facts, in the order they matter for the Aug 20 kill/fallback decision:

1. **The instrument still beats its lagging comparator, at a far higher coverage price.**
   +58.1 cycles at full coverage and +9.2 at half coverage under white noise, at a matched
   false-alarm rate — against +88.8 and +48.8 when the comparator was running ten times
   stricter. Kill criterion 2 ("no lead-time advantage under any regime") is **not**
   triggered.
2. **The coverage floor is the paper's number and it moved a long way: c ≈ 0.41–0.44 under
   MNAR**, band 0.13 wide, versus 0.26–0.33 before. Anything in §1/§5 that implies a
   quarter of participants suffices must be rewritten.
3. **Red noise at mid or high severity removes the advantage entirely.** The white-noise
   headline may never be quoted without its colour *and* its severity (§8.4).
4. **The false-alarm rate still misses the external bar** — 4.39% [3.45, 5.56] against
   3.6% (§8.3's group-unit interval, corrected 2026-09-05) — and the paper must report
   the achieved rate, not the nominal one.

   > **CORRECTED 2026-09-04.** This used to read the residual as "concentrated at low
   > coverage, which is itself a reportable second cost of consent-bounded observation." The
   > by-coverage series is non-monotone (6.33%, 4.50%, 2.50%, 2.22%, 5.56% at c = 0.25, 0.333,
   > 0.5, 0.75, 1.0), and heterogeneity across coverage is p = 0.010 over runs but p = 0.15 at
   > the (arm, replicate) group level (§8.3). **The false-alarm rate varies with coverage in a
   > way these run counts do not resolve** — it is not a directional cost.
5. **C-II's semantics leg now rests on one measure** (aggregate invalidation, 0.993 σ). The
   timing contrast and the doubled-FPR effect are both gone (§8.7). Kill criterion 1 is
   still not triggered, but the sub-claim is narrower than the frozen spec's phrasing.
6. **C-III cannot be asserted as protect-or-compound.** At 500 matched pairs the
   interaction is bounded below ≈ 0.07 correlation units and the only resolved coupling
   effect runs the *unexpected* way. Of §7 item 1's three options, (a) is a design change,
   and the evidence now supports (b) — re-word C-III as a bounded negative with the bound
   quantified — over (c). **That re-wording is a decision for the paper, not the executor's.**
7. **What did not change:** Gate 1 (§2) is untouched — the bifurcation obligation is
   discharged on the substrate, independent of any detection rule. The Ensign-fix arm
   (§4.2) is threshold-independent and stands as reported. §7 items 3, 4, 8 and 9 (the
   comparator-insensitivity assumption, the null MNAR/MCAR contrast, the single scenario
   shape, the single substrate size) are untouched by this campaign and remain the
   standing weaknesses.

### 8.10 Observed-graph detectors under the matched-calibration protocol

Two observed-graph series were scored with the existing leading-detector form. The form
retained the 30-cycle window, in-window linear detrending, within-run reference-period
standardisation, mean of variance and lag-1-autocorrelation z scores, five-cycle sustain,
cycle-90 alarm start, and the existing quarter-step threshold grid. Each graph observable
was calibrated independently on the same fresh dedicated null of 600 control runs at full
coverage, white/mid error, and MNAR consent. The seeds 960000 through 960599 are disjoint
from every scored seed and earlier calibration block. No result-dependent choice was made
between the two observables.

The smallest admissible thresholds were **k = 10.50** for observed giant-component
fraction and **k = 10.75** for observed mean repeat-tie degree. Each produced 21 false
alarms among 600 calibration controls, or **3.50%** (Wilson 95% CI **2.30% to 5.29%**).

The frozen rules were then scored on **4,320 of 4,320 headline runs**: four coverage
levels, three error colours, three severities, two consent arms, two scenarios, and 30
replicates in every cell. The 2,160-control pooled denominator is the same convention used
for the latency detector on this grid. Evaluated false-alarm rates were:

| detector | false alarms / controls | evaluated rate | Wilson 95% CI |
|---|---:|---:|---:|
| observed giant-component fraction | 153 / 2,160 | **7.08%** | **6.08% to 8.24%** |
| formation-latency detector | 92 / 2,160 | **4.26%** | **3.49% to 5.20%** |
| observed mean repeat-tie degree | 36 / 2,160 | **1.67%** | **1.21% to 2.30%** |

**THE FOUR SENTENCES THAT FOLLOWED HERE ARE RETIRED, NOT QUALIFIED, AND THE CORRECTION
BELOW SAYS WHY.** They read: neither graph detector achieved a false-alarm rate matched to
the latency detector at the run unit; the giant-component interval lay wholly above the
latency interval and the mean-degree interval wholly below it at the run unit; therefore no
graph-versus-latency lead estimate, coverage-dependent lead curve or win/loss claim was
reportable from this campaign; and the stop rule fired before timing was interpreted.

At the group unit all three intervals overlap, so none of the four holds. **The consequence
is larger than a wording fix: the stated reason for suppressing every timing result from
this campaign was that the detectors were unmatched, and at the corrected unit they are not
distinguishable from the comparator in either direction.** Whether anything from this
campaign is therefore reportable is a separate question, and this record does not settle
it — a rule that does not fire is not a licence to report, it removes the reason given for
not reporting.

> **CORRECTED 2026-09-04.** This was called "the predeclared stop rule." No prior declaration
> of the matching criterion was found — in the runner, the frozen spec, or the 2026-09-03 spike
> and results notes, four files searched — so call it **post hoc** unless one turns up. And
> "wholly above / wholly below" is a run-unit statement over 2,160 control rows that are 240
> groups of 9 (one seed per group in 240 of 240 groups): design effect 9.04 for both graph
> detectors, 3.79 for the comparator. At the group unit (n_eff): giant-component fraction 7.08%
> [4.46%, 11.06%], mean repeat-tie degree 1.67% [0.65%, 4.21%], comparator 4.26% [2.88%, 6.25%]
> — **both graph intervals overlap the comparator's**, so the rule's premise does not hold at
> the unit the detectors actually vary at, and the rule does not fire there.

Both raw graph series had 0 NaN cycles among 345,600 cycles at each of four coverage
levels. At coverage 0.25, the derived giant-component detector score was NaN in 693 of
248,400 scored cycles because sparse windows sometimes had undefined lag-1
autocorrelation; the alarm code treated each such cycle as not firing. It had 0 of 248,400
score NaN cycles at each higher coverage. Mean-degree had 0 of 248,400 score NaN cycles
at each of four coverage levels. For comparison, formation latency had 9,693 raw NaN
cycles among 345,600 cycles at coverage 0.25, but gap handling left 0 of 248,400 scored
cycles undefined.

The low-coverage stability diagnostic used 120 of 120 white/mid runs and 38,400 of 38,400
cycles across both consent arms and both scenarios. Consented-node counts ranged from 86
to 120 nodes among 400 nodes per run. Observed edges averaged 25.961 per cycle, ranged
from 0 to 72, and were zero in 1 of 38,400 cycles. Observed giant-component fraction was
defined in 38,400 of 38,400 cycles, but its median within-run standard deviation was
0.02437, its run means ranged from 0.04412 to 0.08300 across 120 runs, and its median
within-run correlation with the full-graph giant-component fraction was only 0.313 across
120 runs. The quarter-coverage estimator is finite but not stable as a proxy for the
full-graph quantity.

Censoring and uncertainty were implemented exactly as for formation latency: every
non-firing detector was assigned cycle 320, and cell means used 2,000 deterministic
run-resampling bootstrap draws. The bootstrap standard error is the standard deviation of
the bootstrap means; the stored 95% interval is the 2.5th to 97.5th percentile interval.
Those timing summaries remain in the aggregate for audit but are not interpreted because
the false-alarm matching condition failed.

Aggregates:
[graph detector calibration](aggregates/graph_detector_calibration.json),
[graph detector headline](aggregates/graph_detector_headline.json), and
[low-coverage stability](aggregates/graph_detector_low_coverage_stability.json).
Runner: [graph detector campaign](run_graph_detector_campaign.py).

### 8.11 Three observed-subgraph observables under a disjoint split

> **CORRECTED 2026-09-04. EVERY SIGNIFICANCE TEST IN THIS SUBSECTION IS COMPUTED OVER RUNS
> AND IS INFLATED BY A DESIGN EFFECT OF EXACTLY NINE.** Each (coverage, arm, replicate)
> triple is simulated once per noise colour and severity, nine settings in all, and the three
> observed-subgraph detectors are invariant to that grid: they return an identical
> fired/not-fired verdict across all nine settings in 494 of 494 groups on the finer grid, and
> in 127 of 127 on this one. The noise instrument perturbs reported scalar streams, not graph
> topology, so a graph detector sees the same graph nine times and votes the same way nine
> times. Treating those nine as independent trials understates every standard error by a
> factor of three and multiplies every chi-square by nine.
>
> **THE CONCLUSIONS BELOW DO NOT SURVIVE THE CORRECTION.** At the correct unit of analysis,
> on this four-cell grid, no detector's false-alarm rate is significantly inhomogeneous across
> coverage and all three observables are MATCHED to the comparator: the "4 of 4 coverage
> cells fail" and "all 3 of 3 observables are unmatched" findings are artifacts. The numbers
> below are left standing rather than edited so the correction is auditable. **Read
> §8.12 for the corrected analysis**, which also reports what survives on a sixteen-cell grid.
>
> None of the figures in this subsection was ever bound by the verifier. Its 128 claims cover
> none of them, so its 3,873 checks passing said nothing about any number here.

Even seeds set thresholds on block A controls; odd seeds supplied every block B rate and
lead. The blocks shared 0 of 240 unique seeds. The fixed 2.00--16.00 quarter-step grid
selected **13.00** for observed giant-component fraction, **7.50** for observed mean
repeat-tie degree, and **14.75** for observed susceptibility.

> **STATED 2026-09-04 — this criterion was not previously written down here.** The selection
> rule is **minimum absolute difference between the observable's false-alarm count and the
> comparator's count on the same block A controls (52 of 1,017), larger k on ties**
> ([`run_three_observable_campaign.py`](run_three_observable_campaign.py) lines 128–135). Each
> observable landed at 54 of 1,017 rows — the closest a nine-step curve can get to 52 — which
> is **6 of 113** independent (coverage, arm, replicate) groups, not 54 independent trials.
> This corrects an earlier characterisation of these as "the smallest admissible value at
> 3.5%": they were not selected against a 3.5% target at all, but by count-matching to the
> comparator.

Each graph rule fired on
54 of 1,017 block A controls; latency fired on 52 of 1,017 controls. Susceptibility is the
mean finite-cluster size after excluding every largest component on the consented-node,
consented-dyad repeat-tie subgraph.

On block B, latency fired on 40 of 1,143 controls, **3.50%**. Giant-component fraction
fired on 63 of 1,143 controls, **5.51%** (two-sided two-proportion z-test p = 0.0204), and
mean degree fired on 81 of 1,143 controls, **7.09%** (p = 0.000128). Both are unmatched,
so their lead times are not reported. Susceptibility fired on 54 of 1,143 controls,
**4.72%**, Wilson 95% CI **3.64%--6.11%**, versus latency's **2.58%--4.73%**; p = 0.1403
POOLED. **THE POOLED TEST IS THE WRONG TEST HERE AND SUSCEPTIBILITY IS NOT MATCHED.** Its
rate is not constant across coverage: 4 of 4 coverage cells fail a two-sided
two-proportion z-test against latency on the same controls, at p = 0.00140 (c = 0.25),
p = 0.00027 (c = 0.50), p = 0.01240 (c = 0.75) and p = 0.00014 (c = 1.00). The pooled
figure passes only because two cells in which it never fires cancel two in which it fires
roughly three times too often. A homogeneity test across the four cells gives chi-square
**67.53** on 3 degrees of freedom for susceptibility against **2.68** on 3 degrees of
freedom for latency, whose rate is constant. Its full-coverage paired censored lead
difference against latency of **+13.20 cycles**, percentile 95% CI **-19.61 to +48.87**,
across 15 of 15 block B replicates, therefore **is not reportable**: it sits on a
false-alarm rate of 36 of 288 against latency's 11 of 288 in that same cell, a factor of
3.3, and a detector firing 3.3 times as often fires earlier. **All 3 of 3 observables are
unmatched and no lead comparison from this campaign is reportable.**

Coverage-specific false alarms for latency, giant fraction, mean degree, and
susceptibility were respectively 10/252, 27/252, 18/252, and 0/252 at c = 0.25;
13/306, 27/306, 45/306, and 0/306 at c = 0.50; 6/297, 9/297, 0/297, and 18/297 at
c = 0.75; and 11/288, 0/288, 18/288, and 36/288 at c = 1.00. At c = 0.25,
undefined detector-score cycles were 333/115,920 for giant fraction, 0/115,920 for mean
degree, and 270/115,920 for susceptibility. Their median within-run Pearson correlations
with the matching full-graph series were 0.316, 0.604, and -0.108 across 504 of 504 runs.
Mean degree is usable at quarter coverage; the other 2 of 3 graph observables are not.
At c = 1.00 the observed subgraph is the full graph, so every observable correlates 1.000
with its full-graph counterpart by construction; that row is definitional and is not
evidence of tracking.

Every non-firing detector was censored at cycle 320, the same rule used for latency.
Uncertainty used 2,000 deterministic run-resampling draws; reported intervals are
percentile 95% confidence intervals and standard errors are bootstrap standard
deviations. Aggregate: [block A calibration](aggregates/three_observable_block_a_calibration_20260903.json)
and [block B results](aggregates/three_observable_block_b_results_20260903.json). Runner:
[three-observable campaign](run_three_observable_campaign.py).

### 8.12 The same detectors on a sixteen-cell coverage grid, clustering corrected

Two changes from §8.11, and the second is why the first was worth running.

**The grid.** Coverage runs on 0.0625 steps from 0.0625 to 1.00, sixteen cells, fifteen
degrees of freedom against the earlier three. The frozen 2026-09-03 thresholds are reused
unchanged — 13.00, 7.50 and 14.75 — because recalibrating on a new grid would produce a
pre-declared rule written after the first campaign's outcome was known, which is the error
that campaign already made once with its pooled criterion. The grid is a strict superset of
the old one, so the four shared cells draw the same seeds under the same rules and must
reproduce their committed counts. They do, at **16 of 16 field comparisons**. 8,892 block B
runs, 609 s on 26 workers.

**The unit of analysis.** Every test is computed once per (coverage, arm, replicate) group
rather than once per run. The measured design effect is **9.00** in every cell of every
detector on both grids, which is what the invariance predicts analytically.

Coverage heterogeneity, comparator against the three observables:

| grid | detector | chi-square over runs | corrected chi-square | df | corrected p |
|---|---|---|---|---|---|
| 4-cell | comparator | 2.68 | **0.30** | 3 | 0.960 |
| 4-cell | giant fraction | 39.85 | **4.43** | 3 | 0.219 |
| 4-cell | mean degree | 49.94 | **5.55** | 3 | 0.136 |
| 4-cell | susceptibility | 67.53 | **7.50** | 3 | 0.058 |
| 16-cell | comparator | 164.73 | **18.30** | 15 | 0.247 |
| 16-cell | giant fraction | 162.99 | **18.11** | 15 | 0.257 |
| 16-cell | mean degree | 145.68 | **16.19** | 15 | 0.370 |
| 16-cell | susceptibility | 383.52 | **42.61** | 15 | 0.00018 |

**One finding survives and the extra degrees of freedom are what recovered it.** Observed
susceptibility's false-alarm rate depends on observation coverage, chi-square 42.61 on 15
degrees of freedom, p = 0.00018, while the scalar comparator's does not, 18.30 on 15 degrees
of freedom, p = 0.247. On the four-cell grid the same corrected test gives p = 0.058 and
decides nothing. The four-cell campaign could not have found this: the effect is real and
three degrees of freedom had no power against it.

**The direction is opposite for two of the three observables, which three degrees of freedom
could not have shown either.** Susceptibility's false-alarm rate rises with coverage, rank
correlation +0.484 over the sixteen cells; giant-component fraction's falls, −0.602. On the
four-cell grid these read +0.949 and −1.000 on four points each.

Matching against the comparator, pooled, at the group level:

| grid | detector | pooled rate | p over runs | verdict over runs | corrected p | corrected verdict | cells failing |
|---|---|---|---|---|---|---|---|
| 4-cell | giant fraction | 5.51% | 2.0e-02 | unmatched | 0.4395 | **matched** | 0 of 4, was 3 |
| 4-cell | mean degree | 7.09% | 1.3e-04 | unmatched | 0.2017 | **matched** | 0 of 4, was 2 |
| 4-cell | susceptibility | 4.72% | 1.4e-01 | matched | 0.6230 | **matched** | 0 of 4, was 4 |
| 16-cell | giant fraction | 4.45% | 4.6e-01 | matched | 0.8072 | **matched** | 0 of 16, was 9 |
| 16-cell | mean degree | 8.30% | 4.4e-16 | unmatched | 0.0068 | unmatched | 1 of 16, was 9 |
| 16-cell | susceptibility | 7.89% | 9.5e-14 | unmatched | 0.0130 | unmatched | 2 of 16, was 11 |

> **FOOTNOTE, 2026-09-04 — coverage-mixture effect, not a separate finding.** The 16-cell
> "unmatched" verdicts for mean degree (8.30%) and susceptibility (7.89%) are driven entirely by
> cells the thresholds (13.00, 7.50, 14.75) were never matched on. On the **four shared cells**
> the fine grid reproduces the earlier block B counts exactly (16 of 16 field comparisons,
> above), matched by construction. Over the **12 new cells only**: latency 4.36%, giant fraction
> 4.09%, mean degree **8.72%**, susceptibility **8.99%**. The thresholds were set to match the
> comparator's count on block A's four-cell coverage mixture (§8.11); a pooled-count match on
> one mixture does not hold on another. **"Unmatched on the sixteen-cell grid" is not carried
> forward as a finding separate from the heterogeneity result above** — susceptibility's
> false-alarm rate depending on coverage while the comparator's does not already says this.

**§8.11's headline conclusion is withdrawn.** It reads "all 3 of 3 observables are unmatched
and no lead comparison from this campaign is reportable." On its own data at the correct unit
of analysis, all three are matched on that grid, and the per-cell rejections that drove the
conclusion — three, two and four cells — are zero. Whether the suppressed lead comparisons
should now be reported is a separate question and is not answered here: the censoring rule,
the bootstrap and the paired intervals all need the same clustering treatment before any lead
number is quoted, and none of that has been done.

The comparator's rate is homogeneous across coverage on both grids once corrected, so the
contrast §8.11 wanted to draw does hold at sixteen cells. It does not hold at four.

Aggregates: [fine grid](aggregates/coverage_finegrid_block_b_20260904.json),
[heterogeneity recheck](aggregates/cluster_recheck_20260904.json),
[matching recheck](aggregates/cluster_matching_20260904.json). Runners:
[fine-grid campaign](run_coverage_finegrid.py), [heterogeneity](cluster_recheck.py),
[matching](cluster_matching.py). The chi-square upper tail is
[chi2.py](chi2.py), validated by [test_chi2.py](test_chi2.py) at 22 checks against the
even-degree closed form already in this tree, published critical values for 1 through 20
degrees of freedom, and a Monte Carlo check on odd degrees of freedom. **Before it existed no
p-value in §8.11 could be computed at all**: the only chi-square tail function here is closed
form for even degrees of freedom and returns NaN for odd, and these tests run on
`cells - 1`, which was 3.

### 8.13 The suppressed lead estimates, recomputed at the group unit — nothing survives

§8.10 suppressed every graph-versus-latency lead estimate, coverage-dependent lead curve and
win/loss claim from this campaign, on a matching rule computed at the run unit. The 2026-09-04
correction retired that rule: at the group unit all three false-alarm intervals overlap and the
rule does not fire. §8.12 recorded what was still owed before any lead number could be quoted —
"the censoring rule, the bootstrap and the paired intervals all need the same clustering
treatment" — and that none of it had been done. This is that treatment. Runner:
[cluster lead](cluster_lead.py). Aggregate:
[group-unit lead](aggregates/cluster_lead_20260906.json).

**The reproduction gate first.** Recomputing detection rates from `per_run_rows` reproduces the
committed cell aggregates exactly: formation latency 0.3704 both ways, observed giant-component
fraction 0.1542 both ways, observed mean repeat-tie degree 0.0000 both ways, over 2,160 of 2,160
degradation runs and 72 of 72 cells. No new number below is reported against an aggregate that
did not first reproduce.

**Invariance, measured rather than assumed.** Over 240 of 240 groups, the censored alarm cycle is
identical across all nine noise settings for observed giant-component fraction (240 of 240) and
for observed mean repeat-tie degree (240 of 240). The comparator is identical in **73 of 240**
groups and is therefore collapsed by mean, never by assumption.

**Observed mean repeat-tie degree fires on 0 of 2,160 degradation runs and on 36 of 2,160
controls.** It has no detections anywhere on this design, so it has no lead to estimate at any
unit. Its false-alarm rate matches the comparator's at the group unit because it is inert, and a
matched inert detector is not a comparison the paper can make.

**Observed giant-component fraction has detections and still yields no reportable lead.** It fires
in 37 of 240 groups against the comparator's noise-dependent firing. The paired estimate's sign
and its significance both move with the pairing rule, which is the finding:

| pairing rule at the group unit | lead | 95% CI | n of 240 |
|---|---:|---:|---:|
| comparator fired in at least 1 of 9 | +12.11 | [−1.83, +24.94] | 35 |
| comparator fired in at least 5 of 9 | +2.12 | [−10.83, +14.41] | 27 |
| comparator fired in 9 of 9 | +11.74 | [+0.56, +34.00] | **3** |
| comparator mean-censored at τ = 320 | +47.50 | [+30.86, +65.02] | 37 |

Positive means the graph detector alarms earlier. Intervals are cluster bootstraps over groups,
2,000 draws. **Two rules exclude zero and neither is usable.** The 9-of-9 rule rests on 3 groups.
The mean-censored rule credits the graph detector for the comparator's non-firing runs being
imputed to cycle 320, which is the censoring artifact §9.2 already identified in a different
column, not a lead. The two rules that neither rest on three groups nor impute the comparator's
misses both cross zero.

**The coverage-dependent lead curve is censoring, not signal.** Censored at τ = 320 the pooled
lead is **−23.61** [−29.97, −17.04] for giant-component fraction and **−41.50** [−46.42, −36.43]
for mean repeat-tie degree over 240 of 240 groups, and it grows more negative with coverage for
mean degree (−17.64 at c = 0.25 to −62.22 at c = 1.00) because the comparator's detection rate
rises with coverage while the graph detectors' does not. A curve driven by which detector fails to
fire is not a lead curve.

**Win/loss at the group unit is unreportable on its face.** Both detectors fire in 3 of 240 groups
for giant-component fraction and in 0 of 240 for mean repeat-tie degree.

**The conclusion, and it is narrower than §8.10's.** §8.10 said these comparisons were not
reportable because the detectors were unmatched. That reason is gone. They remain not reportable
for a different and stronger reason: on this design the two graph detectors produce three usable
paired observations between them, and every estimate that reaches significance does so through
either a three-group sample or the censoring imputation. **No graph-versus-latency lead estimate,
coverage-dependent lead curve or win/loss claim is reportable from this campaign at the correct
unit of analysis.** The manuscript's existing sentence, that lead-time comparisons are withheld
pending this correction, is correct as written and now rests on the measurement rather than on the
retired rule. Nothing in the manuscript changes, and this costs the page budget nothing.

## 9. Remediation campaign — 2026-08-05

*Nine targeted analyses, answering pre-submission review findings. Five are new runs; four
are estimator/reporting changes with
**no new runs** — pure functions of committed aggregates. Nothing here re-tunes anything:
the frozen operating point, every campaign-1/2 aggregate, the claim set and the kill
criteria are retained unchanged, and every variant run is matched **seed-for-seed** to the
headline run it is compared against (`seed_for("ea|{c}|mnar", rep)`), so the base variant
reproduces `ea_headline_recal.json` bit-for-bit where they overlap. The verifier checks
that directly and now runs **3,696 checks, 0 failures** (classes C11–C12 added; C12
re-executes `remediation_analysis.compute()` and requires byte-identical output).*

### 9.1 What ran

| stage | script | scale | wall time |
|---|---|---|---|
| both-class threshold calibration (A-11) | [`calibrate_class.py`](calibrate_class.py) | same 600-run null as A-8, same ≤ 3.5% target; only `leading.classes` and `leading.k` differ from the frozen rule | 412 s |
| remediation sweep (A-11 … A-15) | [`run_remediation.py`](run_remediation.py) | **5,100 runs** — 14 variants (9 gain × 40 reps, driver/partial/size × 30 reps) × 5 c, degradation + control, every run scored under BOTH rules | 1,677 s |
| re-analysis (A-16 … A-18) | [`remediation_analysis.py`](remediation_analysis.py) | no runs — pure function of `ea_headline_recal.json`, `ea_headline.json`, `remediation_runs.json`, `calibration_null_bothclass.json`, `gate1_bifurcation.json` | seconds |

Evidence: [`aggregates/remediation_runs.json`](aggregates/remediation_runs.json) +
[`aggregates/remediation_analysis.json`](aggregates/remediation_analysis.json) +
[`aggregates/calibration_null_bothclass.json`](aggregates/calibration_null_bothclass.json).

### 9.2 De-censoring the lead time — the primary objection, answered

The objection: `scoring.py` imputes a non-firing detector to the run length (320),
so at c = 0.5 the "+9.2 cycles" is 55% imputation, and the Δt = 0 iso-line sits exactly
where the imputation dominates. All four replacement quantities are now computed
(white/MNAR/mid, 100 reps on the iso-line neighbourhood):

| c | lead detect [Wilson 95] | lag detect | Δt censored (τ=320) | **Δt paired** (n both fired) | P(lead first \| any) |
|---|---|---|---|---|---|
| 0.25 | **0.17** [0.11, 0.26] | 0.73 | −22.0 | **+66.1** (11) [43.0, 89.3] | 0.20 |
| 0.333 | 0.31 [0.23, 0.41] | 0.83 | −8.3 | **+55.3** (27) [32.4, 78.7] | — |
| 0.50 | 0.45 [0.36, 0.55] | 0.77 | +9.2 | **+49.9** (38) [31.9, 67.1] | — |
| 0.75 | 0.57 [0.39, 0.73] | 0.83 | +12.3 | **+39.5** (15) [10.5, 65.8] | — |
| 1.00 | 0.70 [0.52, 0.83] | 0.73 | +58.1 | **+83.9** (17) [56.7, 110.3] | 0.5 crossing at c ≈ 0.69 |

> **CORRECTED 2026-09-04 — the paired lead is contaminated by control twins, and the
> contamination is at FULL coverage rather than low.** The paired column counts runs where
> both detectors fired, without asking whether that run's matched control also alarmed. Joining
> each degradation run to its control twin on (coverage, arm, colour, severity, replicate):
> at c = 1.00 the cell's control false-alarm rate is 13.3% (4 of 30) and **4 of 17** both-fired
> runs have a twin that also alarmed; excluding them the paired lead is **+70.0 (n = 13)**
> against the +83.9 (n = 17) printed above, and the contaminated runs alarm at median cycle 151
> against 195. At c = 0.25 the contamination is 2 of 11 and the clean mean is **72.4** against
> 66.1; at c = 0.50 it is 2 of 38 and **47.5** against 49.9. This was checked expecting the
> LOW-coverage paired lead to be inflated by early false alarms. It is not — the full-coverage
> one is, and finding 1 below, that the low-coverage failure mode is not firing rather than
> firing late, survives the check unchanged and is if anything strengthened by it.

Three findings, in decreasing order of comfort:

1. **The paired lead time is positive at every coverage** — its lower 95% bound never
   crosses zero, and at c = 0.25, where the censored mean reads −22, the leading detector
   *when it fires* still leads by 66 cycles (and fires first in 10 of its 11 both-fired
   runs). The instrument's low-coverage failure mode is **not firing**, not firing late.
2. **The RMTA identity closes to 5 × 10⁻¹⁴ over all 30 cells**: mean(dt_censored) is
   *exactly* RMTA_lag(τ) − RMTA_lead(τ) at τ = 320. The censored mean was a legitimate
   estimand mis-named "lead time"; the error was naming and horizon-dependence, not
   arithmetic.
3. **The Δt = 0 iso-line is not a coverage-independent property of the instrument.** Its
   location under ΔRMTA moves with the horizon — c\* = 0.28–0.30 at τ ∈ {160 … 280},
   0.412 at τ = 320 — and the censoring-robust estimators put the crossing materially
   higher: **lead-detection-probability = 0.5 at c ≈ 0.61**, **win-probability = 0.5 at
   c ≈ 0.69**. The paper must retire the single iso-line and report a detection floor and
   a win-probability floor alongside the ΔRMTA zero at a declared horizon. Bootstrap on
   the τ = 320 floor: MNAR 0.417 [0.315, 0.595].

### 9.3 The two-channel rule the manuscript described but never ran (A-11)

The frozen and recalibrated rules score `["obs_formation"]` only; §5's claim of "formation
**and** escalation latency" was false as written (spec §9's note (b) is corrected in
place). A-11 calibrates the genuine both-class rule on the **same 600-run null** at the
same ≤ 3.5% target (k = 5.0 vs 6.5; `calibration_null_bothclass.json`) and scores every
remediation run under both rules. Result — the misstatement was hiding a better rule:

| c | Δt censored, formation-only | Δt censored, **both-class** | lead detect f-only → both |
|---|---|---|---|
| 0.25 | −15.0 | **+10.9** | 0.20 → **0.55** |
| 0.333 | −5.6 | **+14.6** | 0.38 → 0.68 |
| 0.50 | +17.1 | **+17.2** | 0.45 → 0.65 |
| 0.75 | +22.1 | **+36.0** | 0.58 → 0.78 |
| 1.00 | +53.1 | **+45.2** | 0.70 → 0.83 |

**Under the both-class rule the censored Δt is positive at every swept coverage — there is
no Δt = 0 crossing in range** — at a pooled control FPR of 6.5% vs the formation-only
rule's 7.0% on the same control arms (small-n; the headline FPR remains §8.3's 4.39% on
3,600 controls). Adding the escalation channel chiefly rescues **detection probability at
low coverage** (0.17 → 0.55 at c = 0.25). The escalation series is observed at c = 1 in
100% of scored cycles, so this is signal the instrument always had and the frozen rule
discarded. §5 should report the both-class rule as the primary operating point going
forward — it is the rule the prose always described — with the formation-only numbers
retained as the frozen-rule record.

### 9.4 Comparator gain (A-12) — the crux survives the sweep

`output_knee` ∈ {0.135 … 0.60} × `output_sigma` ∈ {0.0125 … 0.10}, 9 variants, x-axis the
**measured** sustained σ-drop available in the lagging series (3.94–4.67σ across the
sweep; the 3.5σ rule is within reach everywhere — the "arithmetically incapable
of firing" objection was wrong at the frozen point, see §9.8). Result: Δt at full coverage stays in
**+50.2 … +56.8** across the entire grid — there is **no swept gain at which it vanishes**
— and the τ = 320 floor wanders **0.328–0.413 with no monotone trend in gain**. The
conditional-deliverable worst case (floor tracks gain linearly) does not
materialise; the honest statement is a floor *band* of roughly 0.33–0.41 under comparator
gain uncertainty, quoted alongside the estimator dependence of §9.2.

### 9.5 Formation driver (A-13) — the sensitivity check the code carried and never ran

`formation_driver = "degree"` wires the leading observable to mean repeat-tie degree (a
fast local quantity) instead of mean component size (the slow order parameter). Result:
**the degree-wired detector detects almost nothing** — lead detection 0.03 / 0.03 / 0.00 /
0.00 / 0.00 across c, censored Δt negative everywhere (−30 … −49), no floor crossing
anywhere. This variant was read as the check that would expose the lead time as
manufactured; what it actually shows is sharper and must be stated in §5 as a
conditionality: **the early-warning signal lives in slow, component-scale structure, and
an instrument that reads only fast local degree sees nothing**. That is a design
requirement on the capture layer (it must resolve cross-team component structure, which is
exactly what consent-bounded aggregation makes hard), not a free choice that happened to
flatter the headline.

### 9.6 Partial-observation amplitude (A-14) — the floor is one-sided in the assumed constant

`partial_obs_sigma` ∈ {0.55, **1.1**, 2.2}: τ = 320 floor at **0.389 / 0.374 / 0.631**.
Halving the assumed error amplitude barely moves the floor; doubling it pushes the floor
from ≈ 0.37 to ≈ 0.63. The floor is therefore robust on the optimistic side and sensitive
on the pessimistic side, and c\* must be reported as **c\*(partial_obs_sigma)** with the
frozen value named as an assumption — the concern stands in one direction.

### 9.7 System size (A-15) — magnitude moves, the floor band does not

N ∈ {200, **400**, 800} with `programs_per_cycle` scaled to hold K fixed and `comp_ref`
scaled with N. Censored Δt at full coverage: **+11.7 / +53.1 / +35.3** — the lead-time
*magnitude* is strongly N-dependent and the headline number must be stated as an N = 400
result. The τ = 320 floor: **0.317 / 0.374 / 0.327** — a ≈ 0.32–0.37 band with no
monotone N-trend at these run counts. Neither confirmation size destroys the effect; the
smallest system nearly does (11.7 censored, though its paired lead time remains positive).

> **CORRECTED 2026-09-05.** This used to call the headline "58 ± 14" and ask that it be
> *captioned* as an N = 400 result. The manuscript's headline is no longer 58 ± 14 at all —
> it now reads **"53 ± 11 cycles (bootstrap s.e., n = 40)"** (§8.4's correction), and this
> row's own printed **+53.1** at N = 400 is that same figure, not a caveat sitting beside a
> still-unmoved 58 ± 14.

### 9.8 The amplitude audit — re-scoring the original arithmetic

The originally proposed model of the lead time (a ~40× gain mismatch; a comparator that
"cannot fire before the fold"; a knee that "never binds") is now measured rather than
argued (`amplitude_audit`):

- **Sensitivity ratio lead : lag at full coverage = 2.28**, not ~40. The original estimate compared
  the leading observable's *end-state* swing (+82% formation latency) against the lagging
  observable's *pre-fold* swing (3.5% at A = 0.40); like-for-like, the comparator's
  available signal is 4.44σ of its reporting noise against a 3.5σ threshold.
- **The comparator does fire** — detection 0.73–0.83 across the comparator cells — and fires *before the
  fold* (cycle 225) in 10–13% of runs; the leading rule manages 53% at c = 1. The lead
  time is a race the comparator loses on speed, not a walkover it never enters.
- **The knee binds in 18.9% of program draws**, not never.

The original directional instinct (§9.5's driver result; §9.6's one-sided floor) was
right where it was right, and the specific arithmetic was wrong where it was wrong. Both
are recorded here for the manuscript.

### 9.9 The MNAR sentence (A-18) — DELETE, now on an interval rather than a judgement

The claim under test: "the MNAR floor sits below the MCAR floor (0.412 vs 0.502) …
realistic consent skew helps rather than hurts detection." Bootstrap of the MCAR-minus-
MNAR floor gap over 4,000 paired resamples: **mean +0.080, 95% CI [−0.133, +0.223] —
straddles zero** (P(gap > 0) = 0.87, far short of resolution). The stated mechanism is
not supported by the data: the achieved observed-edge difference (MNAR − MCAR at c = 0.5) is
**0.0012 [−0.0060, +0.0084]** (§3.3's correction of 2026-09-05: the interval widens by √3
once the 300 pooled rows are recognised as 100 replicates × 3 noise-invariant colours, a
design effect of 3, not 300 independent rows) — no advantage for a consenting central
participant is resolved in observed edges. And the seed-noise gauge (MNAR vs MCAR at c = 1, the same physical
configuration on disjoint seeds) shows |Δt| disagreements up to **39.2 cycles** under the
recalibrated rule — larger than campaign 1's 26.2 — so no 0.09 difference in c can be read
off these grids. §5's sentence is deleted and the null reported, per RESULTS §7 item 4,
which said this all along.

### 9.10 Terminology and the graph object (A-19)

"Cumulative repeat-tie network/graph" is corrected to **"20-cycle sliding-window
repeat-tie graph"** everywhere including the abstract; `tie_window = 20` and the 30-cycle
detrending window are stated in §5 (`graph_object`: the gate's giant component and the
leading observable read the *same* windowed object, whose expiry runs every cycle, and the
mean-field K carries the window length). A genuinely cumulative graph densifies
monotonically and cannot collapse; the old descriptor contradicted the mechanism that
produces the headline phenomenon.

### 9.11 What the remediation campaign changes for the paper

1. **Retire the unconditional Δt = 0 iso-line.** Report detection probability per cell
   (first-class), Δt|both-fired, ΔRMTA at declared horizons, and the win-probability
   ordering; quote the floor as estimator- and horizon-dependent (τ = 320 ΔRMTA: ≈ 0.41;
   detection-0.5: ≈ 0.61; win-0.5: ≈ 0.69) with the §9.4 gain band and §9.6's
   c\*(partial_obs_sigma) attached (§9.2, §9.4, §9.6).
2. **Promote the both-class rule to the primary operating point** — it is the rule §5
   already described, it is calibrated on the same null at the same target, and it removes
   the in-range Δt = 0 crossing while nearly tripling low-coverage detection (§9.3).
   Frozen-rule numbers stay as the audit trail.
3. **State the driver conditionality as a design requirement:** the signal exists only for
   an instrument that resolves slow component-scale structure (§9.5).
4. ~~**"58 ± 14" becomes "58 ± 14 at N = 400 under the frozen rule"** with the N-sweep
   quoted (§9.7).~~ **CORRECTED 2026-09-05: this recommendation was overtaken by a bigger
   change.** The manuscript did not caption 58 ± 14 as N = 400; it replaced the headline
   itself with the nested 40-replicate figure, **"53 ± 11 cycles (bootstrap s.e., n = 40)"**
   (§8.4). The N-sweep of §9.7 is still quoted alongside it.
5. **The MNAR-helps sentence is deleted and the null reported** (§9.9).
6. **Abstract/§5 terminology corrected** per A-19; the verifier count in §5 updates from
   2,262 to **3,696** (§9.10, this campaign).
7. **What did not change:** Gate 1, the red-noise kill of the advantage, the achieved-FPR
   obligation (§8.3), the revocation semantics leg, and E-B's bounded-null reading are all
   untouched. §7 items 8 (single scenario shape) and the §8.9 residue not covered by
   the first remediation tranche (detrending/window sensitivity, control-scenario baseline, bimodality test,
   small-shock relaxation, revocation decomposition) remain the standing weaknesses —
   triaged as the
   candidate second tranche.

## 10. Tranche-2 remediation — 2026-08-05

*The residue from the first remediation tranche, addressed in full. Same discipline as §9: sensitivity
and robustness reporting only, matched seeds throughout, scored seeds < 900000, the
950000+ block calibration-only. 1,620 runs, ~14 min wall at 8 workers. Evidence:
[`aggregates/tranche2_runs.json`](aggregates/tranche2_runs.json), runner
[`run_tranche2.py`](run_tranche2.py). Verifier class **C14** covers the aggregate;
the full suite now runs **3,861 checks, 0 failures**.*

### 10.1 Detrending and window — the lead time survives the filter

Six variants — {in-window linear, Gaussian-residual (bandwidth 10)} × window
{20, 30, 45} — each with its threshold re-picked by the frozen procedure on the same
600-run null (re-simulated bit-for-bit; the frozen `linear_w30` variant re-selects
k = 6.50, the recal rule's own threshold, which is the sanity anchor):

| variant | k | Δt @ c=1 | Δt @ c=0.5 |
|---|---|---|---|
| linear_w20 | 5.00 | +44.2 | — |
| **linear_w30 (frozen)** | **6.50** | **+53.1** | **+17.1** |
| linear_w45 | 10.25 | +45.1 | +12.9 |
| gauss_w20 | 5.50 | +32.1 | +3.1 |
| gauss_w30 | 6.50 | +46.6 | — |
| gauss_w45 | 9.50 | +41.9 | — |

Dakos-style filter sensitivity is real but bounded: the worst variant (Gaussian
residual at window 20, i.e. the window equal to the structural correlation time)
attenuates the full-coverage lead time by ~40%; **no variant kills or reverses it**,
and the requested window bracketing (20 vs 45, both sides of the
`tie_window` = 20 coincidence) moves it by ≈ ±9 cycles. §5 gains one sentence and this
table's range; Dakos et al. (2012) gets engaged in text as required.

### 10.2 Control-scenario baseline — no directional bias found

Re-standardising every rolling statistic against a matched control twin (same seed)
instead of the run's own partially-post-onset reference period:

| c | Δt own-baseline | Δt control-baseline | lead DR own → ctl |
|---|---|---|---|
| 0.25 | −14.9 | −14.8 | 0.20 → 0.20 |
| 0.333 | −5.6 | −12.8 | 0.38 → 0.33 |
| 0.50 | +17.1 | +5.0 | 0.45 → 0.35 |
| 0.75 | +22.1 | +32.5 | 0.58 → 0.65 |
| 1.00 | +53.1 | +50.9 | 0.70 → 0.78 |

The differences carry **no consistent sign across coverage** (down at mid-c, up at
0.75, flat at the extremes; per-run median alarm-time change is 0 cycles at every c),
so the contamination concern does not manifest as a directional bias at
these run counts — baseline inflation neither systematically manufactures nor
systematically suppresses the lead. Reported as a robustness note; the frozen
own-baseline rule (the one an operator can actually execute, which needs no
counterfactual control run) stays primary, now on evidence rather than assertion.

**Disclosure added 2026-09-05, and this note predates it.** The control twins used
above as a standardisation baseline were themselves measured alarming at **13.3 per
cent (4 of 30) at full coverage** by section 9's correction of 2026-09-04, filed
thirty days after this section was written and against which this section was never
re-read. The control-baseline column at c = 1.00 is therefore computed against a
population in which 13.3 per cent of controls alarm. Whether that moves the sign
pattern this section's conclusion rests on is **unestablished and is not asserted
either way**; the recomputation with alarming twins excluded is available from the
existing per-cell flags and has not been run.

### 10.3 Bimodality — confirmed, and the occupancy anomaly resolved

**The mixture test the paper lacked:** 1- vs 2-component Gaussian-mixture BIC
(deterministic EM) on the pooled long-run end states. ΔBIC(1−2) = **+204.5 (N=400)**
and **+258.9 (N=800)** on the committed 700-cycle Gate-1 rows, and **+119.9 / +134.0**
on new 2000-cycle runs — two components decisively preferred everywhere (threshold 10).

> **CORRECTED 2026-09-05 — see §2's fuller correction, which this repeats.** "Decisively
> preferred everywhere" describes the **pooled** fit only, and the pool mixes 6 of 6 driver
> values at 700 cycles and 3 of 3 at 2,000 on both sides of the fold — pooling unimodal
> distributions with different means prefers two components by BIC regardless of
> bistability, so this is not evidence the mixture test the paper lacked now confirms
> anything. The manuscript no longer prints a ΔBIC value at all. Refit within fixed driver
> values (§2): two components win in 3 of 12 cells at 700 cycles and 2 of 6 at 2,000. The
> **hysteresis** evidence (branch-mean separation, §2 item 3) does not depend on this test
> and is unaffected.

**The anomaly** — middle-band occupancy *rising* 0.167 → 0.25 from N=400 to N=800 on
700-cycle runs, the opposite of the finite-size expectation — **is incomplete
relaxation, exactly the alternative proposed.** Occupancy by horizon on
the new runs, at the crossing A = 0.37: N=400: 0.167 (h=600) → **0.000** (h ≥ 1000);
N=800: 0.167 (h=600) → **0.083** (h ≥ 1000, one seed of twelve). The larger system
relaxes more slowly, so at a fixed 700-cycle horizon it *looks* more mid-band; given
time, the middle band empties at both sizes. Gate 1's bimodality conclusion stands;
its §2 occupancy line must be re-written to say this.

### 10.4 Small-shock relaxation — the caveat is real and must be printed

- **At shock_frac = 0.05 the relaxation rate is unmeasurable** — fitted rates
  straddle zero with R² ≈ 0.3–0.5: a 5% tie deletion is inside the equilibrium
  fluctuation band, so the requested small-perturbation protocol cannot
  produce a rate at these system sizes (and even it fails to return in 25% of runs at
  A = 0.37 — the basin is genuinely shallow at the fold).
- **At the frozen 0.35 shock, the basin-escape audit changes the reportable claim:**
  returned-fraction is 1.00 for A ≥ 0.39, 0.75 at 0.38, and **0.58 at 0.37** — so the
  headline "10.6× divergence (31.5 → 333.6 cycles)" mixes true relaxation with runs
  that crossed the separatrix and never came back. The clean statement is the
  divergence over the all-returned range: **31.5 → 115.6 cycles (3.7×) across
  A = 0.48 → 0.39**, with the 333.6 at A = 0.37 quoted as an escape-contaminated
  upper anchor (42% non-returning). Critical slowing down survives; the multiplier
  the abstract can honestly carry is the 3.7× clean figure, with the fold-point value
  caveated.

### 10.5 Revocation decomposed — the original guess was wrong in an informative way

The 0.993σ retroactive displacement (reproduced exactly on the aux seeds), decomposed
against the static-equivalent run on the same seed (the record an operator holds when
the same-magnitude reduction is **disclosed** up front):

| contrast | mean abs | level | shape sd | alarm flips |
|---|---|---|---|---|
| retroactive (published vs re-derived) | **0.993** | **0.024** | 0.84 | 0.7% |
| disclosed (published vs from-day-one record) | 1.453 | 0.222 | 1.254 | 0.1% |
| reconstruction (re-derived vs from-day-one) | 1.462 | −0.198 | 1.322 | 0.8% |

Three sentences for the paper. The retroactive shift is **almost entirely shape, not
level** (level 0.024 of 0.993) — the original "roughly what deleting 40% of a sample
does to a re-derived level" is quantitatively wrong. It is **smaller than the
difference a disclosed same-magnitude reduction produces** (0.99 vs 1.45), so
revocation is not an outsized version of coverage loss. And the reconstruction
residual (1.46) shows **re-derivation does not recover the record an honestly
low-coverage instrument would have produced** — the pre-revocation history was
captured at high coverage with different error realisations, and deleting agents from
it is not the same as never having observed them. That last number is the sharpest
statement of "revocation changes the epistemic status of the record" the simulation
has produced.

### 10.6 Non-monotone A(t) — a real false-positive channel, and a caveat on §9's rule promotion

Three declared shapes at c = 0.5, thresholds unchanged, 40 reps each:

| scenario | collapsed | lead alarms (recal) | lead alarms (both-class) | lag alarms |
|---|---|---|---|---|
| dip_recover (to A=0.40, above the fold, and back) | 0.00 | 0.10 | 0.23 | 0.20 |
| oscillate (bounded, period 80) | 0.00 | **0.28** | **0.43** | 0.08 |
| ramp_hold_recover (crosses, holds, recovers) | **1.00** | 0.33 | 0.80 | 1.00 |

1. **Bounded oscillation is a genuine false-positive generator for variance-family
   rules**: 28% of non-transitioning oscillating runs alarm (vs 7.5% for the lagging
   comparator, and vs the ~4–6% flat-control FPR). A deployed instrument watching a
   seasonally-driven ecosystem needs this failure mode in §5's limitations.
2. **The both-class rule is markedly more FP-prone on non-monotone drivers** (43%
   oscillate, 23% dip) — §9.11's promotion of the both-class rule must therefore be
   conditioned: it buys low-coverage detection on monotone degradation at the price
   of oscillation robustness, and the paper should present the two rules as a
   declared trade rather than a strict upgrade. (This is what running one more leg
   buys: §9.3 alone read as a free lunch.)
3. **ramp_hold_recover is a hysteresis demonstration the paper gets for free**: the
   driver recovers fully and the repeat-tie graph stays collapsed in 40/40 runs —
   driver recovery is not system recovery, which is both a Gate-1 corollary and a
   policy-relevant sentence for §7.

### 10.7 What tranche 2 changes for the paper

1. §5 detrending sentence + Dakos engagement + the 32–53-cycle filter band (§10.1).
2. Baseline-contamination robustness note; own-baseline stays primary (§10.2).
3. ~~Gate-1 text: add the mixture-BIC line~~; re-write the occupancy sentence as a
   horizon effect (§10.3). **CORRECTED 2026-09-05: the mixture-BIC line was cut rather than
   added** — the pooled fit conflates driver values and does not survive a within-driver
   refit everywhere (§2, §10.3) — but the occupancy-as-horizon-effect rewrite still stands.
4. Critical-slowing-down multiplier: quote **3.7× clean (A 0.48 → 0.39)**, caveat the
   10.6×/333-cycle figure as escape-contaminated; state that 5% shocks are
   sub-noise (§10.4).
5. Revocation paragraph: shape-not-level, smaller-than-disclosed, and the
   reconstruction residual as the semantics-claim carrier (§10.5) — replaces the
   two by-construction-zero legs as the lead evidence, alongside the IPW result.
6. §5 limitations: oscillation false-positive channel; both-class promotion becomes
   a declared trade (§10.6) — **amends §9.11 item 2.**
7. Verifier count in §5: **3,861**.
8. Residue now empty except: one item (bistable-window width, analytic-vs-empirical —
   prose reconciliation), three mandatory-prose items, and the §C prose/architecture
   fixes from the review — all manuscript work, no further runs required.

## 11. Binding three headline numbers to their aggregate cells — 2026-09-02

Three internal checks flagged a gap on 2026-09-02, each having
searched this file for the values below and found nothing. All three are committed and
reproducible in aggregates this file already names; they were simply never quoted here.

| value as the manuscript prints it | aggregate | cell | raw figure |
|---|---|---|---|
| median alarm cycle 195 vs 286 at full coverage | [`aggregates/ea_headline_recal.json`](aggregates/ea_headline_recal.json) | `cells[c=1.0,arm=mnar,colour=white,severity=mid]` | `lead_median_t` **195.0**, `lag_median_t` **286.0** (n = 30 reps) |
| binds 124 claims in this file | [`claims.json`](claims.json) | `claims` array; `results_file` field | **124** entries; `results_file` = `"RESULTS.md"` |
| edge fraction 0.069 vs 0.062 at c = 0.25 | [`aggregates/eb_joint_sweep_n500.json`](aggregates/eb_joint_sweep_n500.json) | `cells[c=0.25,g_name=high,mode=coupled]` / `[…mode=severed]` | `achieved_edge_frac_mean` **0.06867274840750463** coupled, **0.06222346971260022** severed (n = 500 reps) |

Row 1 is §8.4's own primary operating cell — the cell that also supplies the **+58.1** ±
14.3 figure still bound in `claims.json` as `ea2.white_c1_dt`/`ea2.white_c1_se`
(`dt_censored_mean`/`dt_censored_se`) — read here on its median-time fields instead of its
mean. **This is no longer the manuscript's headline** (§8.4's correction of 2026-09-05: the
headline moved to the nested 40-replicate cell, 53 ± 11); `claims.json` itself was not
touched by that move and still binds the n = 30 figure. Row 3's cells sit at the same
(c, g, mode) coordinates as `withdrawn_total_mean` (bound as `eb2.withdrawn_high`) in the
same aggregate.

**Verifier coverage.** C5 (§6) binds every `claims.json` entry to its aggregate path and
checks the bound text appears verbatim in this file. Four entries were added on 2026-09-02 so
that binding now covers rows 1 and 3: `ea2.white_c1_lead_median` and `ea2.white_c1_lag_median`
on `lead_median_t`/`lag_median_t`, and `eb2.edgefrac_c025_high_coupled` and
`eb2.edgefrac_c025_high_severed` on `achieved_edge_frac_mean` at both modes. The claim count
moves from 124 to 128 with them, which is why row 2 above reads 124 and this paragraph reads
128: the table records the figure the manuscript prints, and this paragraph records the file.

Row 2 has no comparable binding and cannot get one. `load()` in `verify_reported_numbers.py`
reads only from `aggregates/`, and `claims.json` sits one directory above it, so no claim entry
can ever reference the claim file's own length. That number is checkable by opening the file,
which is what a referee would do.

## 12. Multiplicity across reported significance tests — 2026-09-05

The clustering corrections of 2026-09-04/05 (§8.12, §9.9, §3.3) fix a different defect than
the one this section names: each corrects a single test's own p-value for correlated
replicates within one comparison. None of them, and nothing else in this file, corrects for
running many such comparisons across coverage cells and detector families.

`paper/_qa/multiplicity-check.py` -- an internal quality-assurance script that lives beside the
manuscript and DOES NOT SHIP IN THIS BUNDLE, which carries `sim/` only -- extracts every
explicit printed p-value in this file: 35 raw occurrences, collapsing to **20 distinct tests** once recap
sentences (the same test restated in a later "what this changes for the paper" summary) and
the document's own naive-versus-corrected pairs (same comparison, kept once at its corrected
value) are removed. The method and every exclusion are stated in the script's docstring and
printed by a run of it. A further population is not included: at least three distinct
comparisons in this file are reported only as a confidence interval with no printed p (§3.3's
MNAR-vs-MCAR edge-fraction difference, and §8.8's two paired coupling-loss cells), and
converting a CI to a p-value would assume a sampling distribution the file does not state
uniformly across every interval. main.tex §5.4 additionally names "60 per-cell tests" for the
coverage-by-gaming analysis; no per-cell table in this file enumerates them, so their
individual p-values are not extractable from this record as written, and are not part of the
20 below.

Benjamini-Hochberg FDR at q = 0.05 and Holm's procedure at alpha = 0.05, applied to the 20
distinct p-values: **7 of 20 survive BH; 4 of 20 survive Holm.** In both procedures the
survivors are dominated by the four uncorrected per-cell z-tests of §8.11 (susceptibility
against the comparator at each of the four coverage cells, p = 0.00140/0.00027/0.01240/
0.00014) whose own significance §8.12 already withdraws for a different reason: those four
cells' own design effect, not multiplicity, stating "the per-cell rejections... are zero"
without printing a replacement per-cell figure. The two §8.12 16-cell "unmatched" matching
results (mean degree p = 0.0068, susceptibility p = 0.0130) survive BH but not Holm, and
§8.12's own footnote already attributes both to a coverage-mixture/threshold artefact rather
than treating them as independent findings.

**One test is this file's standing claim, and it survives both procedures.** Susceptibility's
coverage-dependent false-alarm rate on the sixteen-cell grid (chi-square 42.61, 15 df,
p = 0.00018, §8.12) is also the only p-value main.tex itself prints (Sec. 5, "Graph-valued
detectors"). Under Holm, the stricter of the two corrections, it is one of four survivors from
the family of 20 defined above; under BH it is one of seven. Correcting for multiplicity does
not overturn the paper's one quantitative graph-detector claim. It does mean that claim has
never, before this measurement, been checked against the number of looks this evidentiary
record actually took.
