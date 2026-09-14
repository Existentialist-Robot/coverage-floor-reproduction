# Results for the revised coverage-floor manuscript

Revised 2026-09-11 using retained synthetic observations; no new simulation campaign was run.
This is the current interpretation and reproduction route. The [historical ledger](history/RESULTS-before-review-2026-09-11.md)
is retained as provenance, including superseded interpretations and the numerical anchors used by
the legacy verifier. Its original relative links are relative to `sim/`. It is not the current
claim record. The public version 1.0.0 archive predates these corrections.

## Reproduction and scope

Run `python sim/review_reanalysis.py --check` from the repository root to regenerate the complete
revised analysis in memory and compare it exactly with
[`aggregates/review_reanalysis.json`](aggregates/review_reanalysis.json). That file contains all
input SHA-256 digests, retained sample sizes, confidence intervals, test statistics, the full
64-comparison detector family and its exclusions. Run without `--check` to write a new output.
`python sim/recompute_group_fpr.py` reproduces the paper's separate 20,000-draw cluster intervals.
`python sim/verify_reported_numbers.py` checks historical aggregates and historical claim bindings;
it does not certify every number in the manuscript or validate model assumptions.

## Model and observation channel

The scored substrate is one sliding-window repeat-tie graph, not the proposed full multilayer
architecture. Mean degree z controls re-engagement through A times
`0.05 + 0.95*z*z/(1.2*1.2 + z*z)`. The mean-field closure `z=7.6*A*g(z)` has its lower fold
at A=0.29925 and bistable width 0.0735. This approximation ignores repeated pairs and finite-graph
correlations. The empirical collapse near A=0.37 and empirical window 0.14–0.15 are distinct
finite-network measurements, not confirmations of the analytic location.

Formation latency uses the team's mean full-graph component size. Observed latency applies
the code's consenting-member timestamp-span shrinkage, missing-dyad noise and proxy error.
Only teams with at least two consenters yield an observation; graph edges require both endpoints.
The main warning score is the mean of standardized rolling variance and lag-1 autocorrelation,
not their maximum. Threshold 6.5, sustain 5, in-window linear detrending over 30 cycles;
the lagging program-output rule uses 3.5 baseline standard deviations and sustain 3.
Calibration uses 600 full-coverage white/mid null seeds with at most 3.5% alarms, taking the
smallest admissible leading threshold and earliest admissible lagging detection with at least
50% detection on calibration degradation runs. The auxiliary both-class rule uses threshold 5.0.

These equations and choices are synthetic. The cited early-warning, privacy, federation and
unlearning literature does not certify this model, its governance guarantees, or this estimator.

## Alarm timing, coverage and noise

The headline is a difference of restricted mean alarm times at horizon 320, assigning a non-firing
detector the horizon. It is not a conditional warning lead among issued alarms or time before
the fold. Full-coverage mid-white-noise gain is 53.125 cycles (bootstrap s.e. 11.396; 40 matched
seeds). Figure 1's full-coverage point is from the overlapping 30-seed threshold campaign and
is about 58 cycles; it was not changed to the 40-seed mean.

The MNAR white-noise zero crossing is 0.412. The bootstrap mean is 0.41667 with percentile
95% interval [0.31473,0.59549]. Four thousand draws independently resample seeds within each
coverage/arm cell, use a first linear crossing and omit out-of-range crossings; all 4,000
white-noise alarm-time draws cross in range. A zero crossing means zero estimated advantage,
not preservation of the full-coverage magnitude or reliable detection.

The leading-detection probability threshold at 0.5 is about 0.61; its conditional bootstrap
interval is [0.46730,0.85000], 3,987/4,000 in-range crossings. The conditional win probability
`P(leading first | either fires)` crosses 0.5 near 0.69; its interval is [0.46788,0.85365],
3,998/4,000 crossings. At horizon 240 the alarm-time point crossing is 0.28. The alternative
threshold intervals are exploratory, conditional on an in-range crossing, and not simultaneous.

Under mid-severity red noise the estimated advantage is absent across the tested coverages;
the low-severity full-coverage cell is about +45 cycles (30 runs, descriptive estimate).
This is conditional on the chosen correlation 0.75 and noise amplitudes. It does not establish
that noise control dominates recruitment in general.

## Flat-control rates and seed-paired comparator sensitivity

The original scalar warning and lagging-output rates are 4.48% and 4.86%, respectively. These
equally weight 720 independent (coverage, arm, replicate) groups: 240 with nine noise copies,
480 with three, 3,600 control runs total. Leading group fractions can be fractional; the lagging
output is invariant within each group (35/720). Their 20,000-draw percentile cluster intervals
are [3.47,5.56]% and [3.33,6.53]%, seed 20260910. These are mixture rates, not per-cell guarantees.
The 30-run full-coverage white-noise threshold cell has 4/30 leading control alarms; it is not
the 40-run headline arm. Bounded period-80 controls alarm in 17/40 both-class, 11/40
formation-only and 3/40 lagging-output runs. Those seasonal controls are not part of the flat-null claim.

Nine comparator variants retain mean full-coverage gains between 50.225 and 56.775 cycles.
The range is descriptive. New matched-seed differences use 40 pairs and 4,000 percentile
bootstrap draws, seed 20260911. Halving reporting noise gives a difference of −1.325 cycles,
95% interval [−2.625,−0.250], versus the base. The other seven non-base intervals include zero.
These unadjusted contrasts are sensitivity analyses, not equivalence tests. In particular,
overlap of nine point estimates with one base interval never established a null effect.

## Missingness and propensity correction

Independent node consent gives expected dyad retention c², but centrality-correlated consent
need not. The MCAR-minus-MNAR white-noise alarm-time threshold interval [−0.13264,+0.22300]
compares thresholds, not edge retention with c². Its resampling is independent within each
arm and coverage cell, not seed pairing between different arms. Likewise the achieved-edge
difference at c=0.5, 0.00117 with group-corrected Wald interval [−0.0060,+0.0084], compares arms,
not either arm with a square law. Replacing component size with mean degree in the latency
mapping yields detection 0.03/0.03/0/0/0 (30 runs per cell), versus 0.17–0.70 for the original
mapping; this is not a general result about all degree-based warning methods.

Generating consent depends on cumulative degree measured before observation, while the fitted
correction uses degree in the final consented graph. Those are not the same covariate; a shared
logistic form does not establish correct specification. The correction further averages clipped
propensities, estimates team size as observed members divided by that average, reverses the known
latency shrinkage and weights by estimated probability of observing at least two team members.
It is a model-specific combined adjustment, not a generic inverse-propensity guarantee.

Bias in each run is the mean recorded-minus-true formation latency over its finite-cycle
intersection; corrected bias uses its own analogous intersection. Reduction is one minus the
ratio of pooled absolute corrected to raw bias across 60 runs. Reanalysis gives static reduction
0.711218 [0.692076,0.731355] and revocation reduction 0.258607 [0.237552,0.279995], using 4,000
run-bootstrap draws at seed 20260911. Estimated mean consent falls from 0.505164 to 0.303084;
applying the final rate to earlier observations is the recorded failure mechanism. This average
does not establish a positivity failure. Per-agent fitted propensities were not retained, so
empirical positivity diagnostics cannot be reconstructed without rerunning simulation.

Record-revision mean absolute displacement is 0.993 baseline standard deviations; signed mean
shift is 0.024. They are separate summaries. A disclosed reduction to the same final consenting
set gives 1.453, and reconstructed versus observed-from-day-one displacement is 1.462. Observation
errors differ between these histories; the residual does not isolate deletion semantics.

## Recovery and coupling

The recovery calculation reproduces 31.50444 to 115.64234 cycles, ratio 3.670669, percentile
95% interval [1.870992,9.234169]. Twelve rates per driver cell are independently resampled,
4,000 draws, NumPy PCG64 seed 20260904, one index matrix per cell. All original rates and all
draws are retained. At the stressed cell one rate is negative and four fits have R² below 0.12.
This is an unstable finite increase, not proof of divergence. The 333.6-cycle fold-point
estimate mixes relaxation with escape (42% non-returning) and remains an upper anchor only.

The gaming statistic is correlation with a deterministic monotone driver and thus measures
trend strength, not latent-state prediction. Off-coupling signed means are −0.406
[−0.443,−0.371] at quarter coverage and −0.605 [−0.633,−0.577] at full coverage. Each coupling
level pools three coverages with 500 matched pairs each. Pooled coupling losses are −0.008324
[−0.026355,+0.010971] at mid coupling and −0.022336 [−0.041365,−0.002996] at high coupling.
Low-minus-full-coverage interaction intervals are [−0.050570,+0.041185] and
[−0.062142,+0.032448]. Their endpoints lie within ±0.07, but their widths exceed 0.07;
unresolved effects do not establish separability or an equivalence bound.

The existing normal-approximation power estimate for the observed pooled mid-coupling effect
is 15,257 pairs at 80% power and two-sided alpha 0.05: 10.2 times the existing 1,500 pooled
pairs, not thirty times that sample. This is effect-based planning, not measured future power.
At quarter coverage with coupling off, 91/500 runs (18.2%; Wilson interval 15.1–21.8%) have
the opposite-signed trend. Coupling has no matched flat-driver controls and changes achieved
edge coverage, so neither a causal mechanism nor robustness is established.

## Graph detectors and multiplicity

The graph experiments compare graph-valued detectors with the formation-latency warning rule;
their scalar reference is not the lagging program-output metric. The old prose conflated these.
All graph verdicts repeat exactly within 494 fine-grid seed groups; the latency reference is
invariant in only 424. The current script directly resamples group fractions rather than using
an approximate common design-effect divisor. The former 44.01 comparator statistic depended
on an inadequately documented estimated divisor; it is replaced by the explicit permutation test.

On the coarse grid, giant-component, mean-degree and susceptibility rates are 5.51%, 7.09% and
4.72%, versus formation-latency 3.50%. New paired group-mean tests give p=0.39086, 0.16516 and
0.54335. No difference is resolved; no equivalence margin or test was specified. All rate
intervals and per-coverage comparisons are emitted in the reanalysis JSON.

Coverage heterogeneity permutes whole seed groups within missingness arm, requiring
exchangeability under the null; 19,999 permutations, seed 20260911, with the plus-one correction.
On the fine grid susceptibility gives p=0.00015 (2 exceedances), formation-latency p=0.00030
(5 exceedances). Comparing their p-values does not compare effect strength. Frozen graph
thresholds remain 13.00, 7.50 and 14.75.

The corrected post-hoc detector family has 64 tests: ten coverage-heterogeneity tests (four
detectors on each of two graph grids, two detectors on the primary grid), six pooled graph
versus latency contrasts, and 48 distinct per-coverage graph contrasts. The twelve coarse-grid
per-cell contrasts recur identically on the fine grid and are counted once. Withdrawn run-level
tests are replaced, not retained. The fraction of bootstrap threshold gaps above zero is not a
null p-value and is excluded. Holm controls familywise error without independence; Benjamini–
Hochberg is a descriptive sensitivity under dependence. **Five survive Benjamini–Hochberg;
two survive Holm.** The Holm survivors are susceptibility and formation-latency fine-grid
heterogeneity, adjusted p=0.0096 and 0.0189. Neither the old seven/four nor eight/five count
describes this family.

This is not comprehensive multiplicity control across all past analyses. Coupling, recovery,
propensity correction, comparator sensitivity and threshold intervals remain exploratory,
unadjusted summaries. The former claim of 60 gaming tests had no enumerated test registry and
has been removed from the manuscript rather than represented as a verified family size.

The current script also recomputes group-bootstrap restricted-mean alarm-time gains versus
formation latency on the same four-cell graph grid, with 127 independent degradation groups.
Giant-component fraction gives −29.035 cycles [−36.493,−21.686], mean degree −44.405
[−51.385,−37.490], susceptibility −21.570 [−29.726,−13.127] (unadjusted 95% percentile
intervals; 4,000 group resamples, seed 20260911). All are negative: no graph-detector
alarm-time advantage is claimed. These exploratory effect intervals are outside the
64-test false-alarm family. The older two-detector lead campaign is separately preserved in
[cluster_lead_20260906.json](aggregates/cluster_lead_20260906.json); it is not the current
three-detector comparison and does not need to be described as unfinished.

## Finite-horizon recovery scenario and publication boundary

All 40/40 ramp–hold–recover runs fail to regain pre-shock mean degree within the run horizon
(Wilson interval for this proportion 91.2–100%). This does not establish permanent
irreversibility or failure of all possible interventions. The real-program capture audit
remains a workflow classification with no populated participant records; the simulation
does not validate real prediction, consent uptake or governance implementation.

The old public archive is retained unchanged. The revised supplementary bundle must carry
the current scripts, aggregates and file manifest before it can represent this numerical
record; a future archive identifier must not be invented or silently substituted.
