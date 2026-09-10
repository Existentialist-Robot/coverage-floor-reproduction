# `sim/` — the COMPLEX NETWORKS 2026 simulation

Implements the frozen simulation specification `sim-spec-FROZEN-2026-08-04`, which is not distributed with this bundle. Read
that spec first — it is the contract; every design decision lives there, and every
implementation-forced deviation is a numbered amendment in its §9.

Results and honest caveats: **[RESULTS.md](RESULTS.md)**.

## Dependencies

Pure Python **stdlib + numpy** for everything that produces a number. `matplotlib` is
used for figures if importable, and [`svgwrite_min.py`](svgwrite_min.py) is a
dependency-free SVG fallback if it is not. No LLM calls, no network access, no scipy
(Kendall tau, the AUC/Mann-Whitney statistic and the logistic propensity fit are
implemented in [`indicators.py`](indicators.py)).

## Reproduction, in order

```bash
# ---- campaign 1 (2026-08-04) -------------------------------------------------
python calibrate_rule.py          # freezes configs/detection_rule.json (~5 min)
python run_gate1_bifurcation.py   # GATE 1 -- must PASS before any EWS is reported
python run_ea_headline.py         # E-A headline sweep (~20 min, 18 workers)
python run_aux_arms.py            # revocation + Ensign arms
python run_eb.py                  # E-B joint (c x g) sweep
# ---- follow-up campaign (2026-08-05): the operating point the paper reports ---
python calibrate_recal.py --workers 18 --n-null 600      # -> detection_rule_recal.json
python run_ea_headline.py --workers 18 --reps 30 --rule detection_rule_recal.json \
       --out ea_headline_recal.json --c-extra 0.333 --iso-c 0.25,0.333,0.5 --iso-reps 100
python run_aux_arms.py --workers 18 --reps 60 --rule detection_rule_recal.json \
       --out aux_arms_recal.json
python run_eb.py --workers 18 --reps 500 --rule detection_rule_recal.json \
       --out eb_joint_sweep_n500.json
# ---- figures + the reader contract ------------------------------------------
python figures.py                 # figures/*.png + *.svg
python verify_reported_numbers.py # the reader contract; exits non-zero on any failure
```

The follow-up campaign re-scores the *same* runs under a threshold recalibrated on a
large dedicated null, adds replicates (and one c grid point) around the Δt = 0 iso-line,
and raises E-B to ≥ 500 matched pairs per cell — spec §9 amendments A-8 / A-9 / A-10.
Campaign 1's aggregates and rule file are kept, so both operating points stay
verifiable side by side.

Every script takes `--workers N`; the default is `cpu_count - 2`. `--quick` (where
offered) runs a reduced grid for smoke-testing only and must not be used for reported
numbers. Runtimes above are for an 22-core desktop CPU; the full set is well under two
hours.

## Determinism

Every run is keyed by an explicit integer seed derived from
`sha256(experiment-tag | replicate)` in [`runner.py`](runner.py), so any single cell can
be re-run in isolation and reproduces bit-for-bit. Paired arms are given the *same* seed
by construction — that is what "matched seeds" means here. Within a run, each subsystem
draws from its own stream spawned from one `SeedSequence`, so switching an instrument
dial off silences that dial without re-phasing the substrate.

**Scored seeds are always `< 900000`; calibration seeds are `>= 900000`.** No scored run
informs the detection rule. `verify_reported_numbers.py` enforces this, and for the
follow-up campaign it also enforces that the 600-run recalibration null (block
`950000+`) is disjoint from campaign 1's calibration block *and* from every scored seed.

## Files

| file | role |
|---|---|
| [`model.py`](model.py) | substrate (six-friction-point program lifecycle, Guimerà assembly, fold bifurcation) + instrument (the four frozen dials) |
| [`indicators.py`](indicators.py) | variance-family statistics, the level and trend rules, Kendall tau, AUC, logistic fit |
| [`scoring.py`](scoring.py) | applies the frozen rule and reduces a run to scalars — the single code path all three experiments score through |
| [`runner.py`](runner.py) | seed derivation, process-pool map, aggregate writer |
| [`calibrate_rule.py`](calibrate_rule.py) | freezes the detection rule on the disjoint calibration seed block |
| [`calibrate_recal.py`](calibrate_recal.py) | re-picks the two thresholds on a 600-run dedicated null (seed block 950000+); rule form untouched |
| [`run_gate1_bifurcation.py`](run_gate1_bifurcation.py) | Gate 1: mean-field saddle-node, end-state bimodality, hysteresis, relaxation-time divergence |
| [`run_ea_headline.py`](run_ea_headline.py) | E-A: (c × ε colour × severity × consent arm), paired controls |
| [`run_aux_arms.py`](run_aux_arms.py) | revocation-vs-static-equivalent and Ensign-fix arms at (c=0.5, mid red noise) |
| [`run_eb.py`](run_eb.py) | E-B: joint (c × g) sweep with the coupling-severed control |
| [`recal_analysis.py`](recal_analysis.py) | cross-campaign derivations (Δt decomposition, FPR intervals, seed-noise gauge, E-B interaction bound + anti-tracking) — no new runs |
| [`figures.py`](figures.py) | the three frozen figures + a Gate-1 appendix figure + fig4 (iso-line), fig5 (calibration null), fig6 (E-B coupling forest) |
| [`verify_reported_numbers.py`](verify_reported_numbers.py) | ten check classes: recomputation, seed hygiene, pairing, rule identity, claim binding, determinism, calibration selection, iso-line re-derivation, the E-B decision rule, cross-campaign analysis |
| `configs/detection_rule.json` | the campaign-1 frozen rule — **written by calibration, never hand-edited** |
| `configs/detection_rule_recal.json` | the recalibrated rule (same form, thresholds re-picked on the large null) — the operating point the paper reports |
| `aggregates/*.json` | committed aggregates; every reported number is derivable from these |
| `claims.json` | each number quoted in RESULTS.md bound to an aggregate path |

`_calcache.npz` is a local cache of the calibration runs (gitignored, regenerable —
delete it to force a fresh calibration).

## What is *not* here

No real-data validation, and no claim of any (frozen spec §7; that is future work).
The substrate is deliberately standard; the novelty budget is spent on the instrument
layer alone.
