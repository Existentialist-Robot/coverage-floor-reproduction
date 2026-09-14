# Reproduce the revised analysis

Requirements: Python 3.11 or compatible, NumPy and Matplotlib. No paid model calls, network
requests or new simulation are required for the retained-data checks below.

From the bundle or paper-repository root:

```text
python sim/review_reanalysis.py --check
python sim/recompute_group_fpr.py
python sim/test_review_reanalysis.py
python sim/verify_reported_numbers.py
```

The first command checks the current statistical record; the final command checks legacy
aggregates and historical textual bindings. Read sim/RESULTS.md for current interpretation,
sample units and exclusions. Run review_reanalysis.py without --check to regenerate its JSON.
To rebuild the retained coverage figure only, run `python -c "import sys; sys.path.insert(0,
'sim'); import figures; figures.fig4_isoline()"` as one line. Full simulation runners are
included for provenance but are unnecessary for these checks and will take substantially longer.

The original public release predates this analysis. The revised bundle is a supplement to
the revised manuscript; its manifest identifies the exact bytes. Historical ledger prose
and original derived-analysis verdict strings are provenance and may describe superseded
inferences. They are not the current conclusions.
