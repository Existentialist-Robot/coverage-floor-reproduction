# Reproduction bundle — *A Coverage Floor for Early Warning on Consent-Bounded Temporal Networks*

This bundle holds the simulation code, the committed aggregates, the frozen operating point and the
verifier behind every number reported in the paper

> Eden Redman, Yashar Talebirad and Osmar R. Zaïane.
> *A Coverage Floor for Early Warning on Consent-Bounded Temporal Networks.*
> COMPLEX NETWORKS 2026, Granada.

The paper argues that a measurable floor of observation coverage is enough to raise early warnings
on a temporal network whose edges may only be aggregated under consent. Everything the paper prints
as a number is recomputable from what is here.

## What is in here

| Path | What it is |
|---|---|
| [`sim/model.py`](sim/model.py) | The agent model: team formation, repeat ties, attrition, and the coverage-limited observation channel. |
| [`sim/indicators.py`](sim/indicators.py) | The observation instrument. Variance-family early-warning statistics, Kendall tau, the AUC statistic and the logistic propensity fit, all implemented without scipy. |
| [`sim/runner.py`](sim/runner.py), [`sim/scoring.py`](sim/scoring.py) | The execution harness and the scoring rule. |
| `sim/run_*.py` | The experimental arms: the gate-1 bifurcation test, the headline coverage sweep, the auxiliary revocation and propensity-correction arms, and the joint coverage-by-coupling sweep. |
| `sim/calibrate_*.py` | Threshold calibration, including the dedicated-null recalibration that fixes the operating point the paper reports. |
| [`sim/verify_reported_numbers.py`](sim/verify_reported_numbers.py) | The verifier. It re-derives every reported number from the committed aggregates and exits non-zero on any mismatch. |
| [`sim/claims.json`](sim/claims.json) | The claim binding: each manuscript claim tied to the aggregate cell that produces it. |
| `sim/aggregates/` | The committed run aggregates. Every reported number is computed from these. |
| `sim/configs/` | The frozen operating point, which the paper's parameter table reports. |
| [`sim/RESULTS.md`](sim/RESULTS.md) | The results record, including the caveats and the negative results. |
| `sim/figures/` | The figures as published, plus their SVG sources. |

## Reproducing the reported numbers

The fastest check does not re-run the simulation. It re-derives every printed number from the
committed aggregates:

```bash
cd sim
python verify_reported_numbers.py
```

That exits zero and reports the number of checks performed. A full re-run of the campaigns, which
takes hours on a many-core machine, is documented in order in [`sim/README.md`](sim/README.md).

Dependencies are the Python standard library and numpy for everything that produces a number.
matplotlib is used for figures when it is importable, and [`sim/svgwrite_min.py`](sim/svgwrite_min.py)
is a dependency-free fallback when it is not. There are no network calls and no model API calls
anywhere in this bundle.

## What is deliberately not here

The manuscript itself is not distributed with the bundle, because the publisher takes copyright
transfer at camera-ready. Cite the paper and use these artifacts.

The frozen simulation specification and the project's internal planning and review material are also
excluded. Where the code refers to the frozen specification by name, that document is the design
contract the implementation was written against; its content is described in the paper's methods
section and in [`sim/RESULTS.md`](sim/RESULTS.md).

## Redaction

This bundle is published from a redacted copy of a private working tree. The redaction removed
internal work-item identifiers, two internal personal attributions, and references to a private
pre-submission review including reviewer identities. **No number, statistic, confidence interval,
table or figure was altered**, and no public scholarly citation was removed. The complete
before-and-after is published alongside this bundle as the redaction diff, so the change is
checkable without access to the private tree.

The verifier passes on this redacted copy at the same check count as on the source tree.

## Licence

Two licences, because this bundle holds two different kinds of thing.

- **Code** — the Python modules and scripts — is under the MIT licence. See [`LICENSE`](LICENSE).
- **Data** — the aggregates, the claim binding, the frozen configuration, the results record and
  the figures — is under Creative Commons Attribution 4.0 International. See
  [`LICENSE-DATA`](LICENSE-DATA).
