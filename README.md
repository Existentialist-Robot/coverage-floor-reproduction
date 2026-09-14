# Reproduction bundle — *A Coverage Floor for Early Warning on Consent-Bounded Temporal Networks*

This bundle holds the simulation code, the committed aggregates, the frozen operating point and the
verification apparatus behind the numbers reported in the paper

> Eden Redman, Yashar Talebirad and Osmar R. Zaïane.
> *A Coverage Floor for Early Warning on Consent-Bounded Temporal Networks.*
> COMPLEX NETWORKS 2026, Granada.

The paper argues that a measurable floor of observation coverage is enough to raise early warnings
on a temporal network whose edges may only be aggregated under consent. Everything the paper prints
as a number is recomputable from what is here. **All data in this bundle are synthetic simulation
output. There are no participant records of any kind.**

## This supersedes v1.0.0 / DOI 10.5281/zenodo.22699236

The archived release [10.5281/zenodo.22699236](https://doi.org/10.5281/zenodo.22699236) — this
repository's v1.0.0 — **predates the corrected analyses.** The manuscript was revised on 2026-09-11
after a pre-submission methods review, and the inferential record was recomputed. The Zenodo record
is deliberately left unchanged; a new archived version will be deposited separately. **Until then,
this repository is the current material and the DOI is not.**

### What changed since v1.0.0

**New analysis, not present in v1.0.0**

| File | What it adds |
|---|---|
| [`sim/review_reanalysis.py`](sim/review_reanalysis.py) | The revised inferential record. Enumerates the family of 64 post-hoc detector comparisons, applies Holm and Benjamini–Hochberg, and recomputes the recovery-ratio and bias-reduction intervals on independent seed groups with explicit bootstrap constructions. Writes `sim/aggregates/review_reanalysis.json`. |
| [`sim/test_review_reanalysis.py`](sim/test_review_reanalysis.py) | Unit tests for the multiplicity and bootstrap machinery above. |
| [`sim/recompute_group_fpr.py`](sim/recompute_group_fpr.py) | The flat-control false-alarm rates as the paper now reports them: cluster bootstrap over 720 independent (coverage, arm, replicate) groups rather than over correlated runs. |
| [`sim/emit_missingness_arm_interval.py`](sim/emit_missingness_arm_interval.py) | Persists the MNAR-minus-MCAR achieved-edge interval the paper prints as an aggregate cell, so it can be bound like every other reported number. Recomputes it from committed per-run rows; changes no reported value. |
| [`sim/aggregates/review_reanalysis.json`](sim/aggregates/review_reanalysis.json), [`sim/aggregates/missingness_arm_interval.json`](sim/aggregates/missingness_arm_interval.json) | The aggregates those two scripts produce. |
| [`sim/REPRODUCE.md`](sim/REPRODUCE.md) | The reproduction route for the revised analysis. |
| [`sim/history/RESULTS-before-review-2026-09-11.md`](sim/history/RESULTS-before-review-2026-09-11.md) | The pre-review results record, kept verbatim as provenance so the correction is inspectable rather than asserted. |
| [`MANIFEST.md`](MANIFEST.md), [`MANIFEST.json`](MANIFEST.json) | Every file with a purpose and a SHA-256. |

**Corrected in place**

`sim/RESULTS.md` now carries the revised claims and the corrections that produced them.
`sim/claims.json`, `sim/figures.py` and `fig4_isoline_recal` were updated with the reanalysis.
Several aggregates and run scripts changed only in wording (see *Redaction* below), not in any value.

**What the revision changed substantively.** The earlier record scored detector comparisons without
enumerating the post-hoc family or correcting for multiplicity, and computed several rates and
intervals over correlated runs rather than over independent groups. The revision enumerates the
family — 64 comparisons — applies Holm and Benjamini–Hochberg, and rebuilds the affected intervals
on independent seed groups. Two comparisons survive Holm. Conclusions that did not survive were
withdrawn rather than softened, and the pre-review record is retained above so both can be read side
by side.

**Provenance strings are not current conclusions.** `sim/verify_reported_numbers.py`,
`sim/claims.json` and the historical ledger prose bind the *legacy* record. They are retained
because they are the audit trail of the earlier release. They do not bind the revised manuscript.

## Reproducing the revised analysis

Requirements: Python 3.11 or compatible, NumPy, and Matplotlib for figures only. No network calls,
no model API calls, and no new simulation is needed for these checks.

From the repository root:

```bash
python sim/review_reanalysis.py --check          # the revised inferential record
python sim/recompute_group_fpr.py                # group-clustered false-alarm rates
python sim/emit_missingness_arm_interval.py --check
python sim/test_review_reanalysis.py             # unit tests
python sim/verify_reported_numbers.py            # legacy verifier: 3,873 checks, 128 ledger claims
```

Each exits zero. Run `review_reanalysis.py` and `emit_missingness_arm_interval.py` without `--check`
to regenerate their JSON. To rebuild the retained coverage figure only:

```bash
python -c "import sys; sys.path.insert(0,'sim'); import figures; figures.fig4_isoline()"
```

Read [`sim/RESULTS.md`](sim/RESULTS.md) for the current interpretation, sample units and exclusions,
and [`sim/REPRODUCE.md`](sim/REPRODUCE.md) for the same commands with their scope caveats. The full
simulation runners are included for provenance; they take hours and are unnecessary for the checks
above. Their order is documented in [`sim/README.md`](sim/README.md).

## What is in here

| Path | What it is |
|---|---|
| [`sim/model.py`](sim/model.py) | The agent model: team formation, repeat ties, attrition, and the coverage-limited observation channel. |
| [`sim/indicators.py`](sim/indicators.py) | The observation instrument. Variance-family early-warning statistics, Kendall tau, the AUC statistic and the logistic propensity fit, all implemented without scipy. |
| [`sim/runner.py`](sim/runner.py), [`sim/scoring.py`](sim/scoring.py) | The execution harness and the scoring rule. |
| `sim/run_*.py` | The experimental arms: the bifurcation validation, the headline coverage sweep, the auxiliary revocation and propensity-correction arms, the joint coverage-by-coupling sweep, the fine coverage grid and the graph-valued detector campaigns. |
| `sim/calibrate_*.py` | Threshold calibration, including the dedicated-null recalibration that fixes the operating point the paper reports. |
| [`sim/review_reanalysis.py`](sim/review_reanalysis.py) | The revised inferential record (see above). |
| [`sim/verify_reported_numbers.py`](sim/verify_reported_numbers.py) | The legacy verifier. Re-derives the earlier record's numbers from the committed aggregates and exits non-zero on any mismatch. |
| [`sim/claims.json`](sim/claims.json) | The legacy claim binding: 128 historical ledger claims tied to the aggregate cells that produce them. |
| `sim/aggregates/` | The committed run aggregates. Every reported number is computed from these. |
| `sim/configs/` | The frozen operating point, which the paper's parameter table reports. |
| [`sim/RESULTS.md`](sim/RESULTS.md) | The results record, including the caveats and the negative results. |
| `sim/figures/` | The figures as published, plus their SVG sources. |
| [`MANIFEST.md`](MANIFEST.md) | Every file, one line of purpose each, with SHA-256. |

## What is deliberately not here

The manuscript itself is not distributed, because the publisher takes copyright transfer at
camera-ready. Cite the paper and use these artifacts.

The frozen simulation specification, the private pre-submission review, and the project's internal
planning material are excluded. Where the code refers to the frozen specification by name, that
document is the design contract the implementation was written against; its content is described in
the paper's methods section and in [`sim/RESULTS.md`](sim/RESULTS.md).

`sim/audit_manuscript_numbers.py` — a read-only audit that binds every numeral printed in the
manuscript source to an aggregate cell — is also excluded, because it takes the LaTeX source as its
input and that source is not distributed. The one binding gap it found has been closed here, by
`sim/emit_missingness_arm_interval.py`.

## Redaction

This bundle is published from a redacted copy of a private working tree. The redaction replaced
internal work-item and amendment identifiers with the descriptive names used in the paper, and
removed references to the internal specification document and to the private pre-submission review.
**No number, statistic, confidence interval, table or figure was altered.** No public scholarly
citation was removed. The legacy verifier passes on this redacted copy at the same check count as on
the source tree, and `emit_missingness_arm_interval.py --check` reproduces its interval from the
redacted aggregates byte-for-byte.

## Licence

Two licences, because this bundle holds two different kinds of thing.

- **Code** — the Python modules and scripts — is under the MIT licence. See [`LICENSE`](LICENSE).
- **Data** — the aggregates, the claim binding, the frozen configuration, the results record and
  the figures — is under Creative Commons Attribution 4.0 International. See
  [`LICENSE-DATA`](LICENSE-DATA).
