# Denial of Access as an Unmeasured Outcome in Nonprescription Statin Self-Selection Studies — analysis code

Analysis code for the manuscript by Kyle P. Molinari (in preparation for
Therapeutic Innovation & Regulatory Science).

## Contents

- `analysis.py` — consolidated certification script. Recomputes every
  derived quantity reported in the manuscript (CREST accuracy measures and
  Clopper-Pearson intervals, the bootstrap interval for the difference
  between positive predictive value and sensitivity, decision curve
  crossovers and net benefit over the prespecified 5% to 20% range, the
  TACTiC co-primary endpoints with and without mitigations, minimum
  mitigation counts, and screening-flow arithmetic) from numerical inputs
  transcribed from the cited sources, and asserts equality with the
  manuscript's reported values at reported precision. A clean run prints
  `ALL 43 ASSERTIONS PASSED`.
- `make_figure4.py` — generates Figure 4 (decision curves) as
  `figure4_decision_curves.png`.

## Requirements

Python 3.12.3, NumPy 2.4.4, SciPy 1.17.1, Matplotlib 3.10.8 (figure only).
Later versions are expected to work; the versions above are those under
which the reported values were certified.

## Run

    python analysis.py
    python make_figure4.py

`analysis.py` exits nonzero on the first failed assertion. One check is
marked SKIP pending one input (the correct-selection component count of the
TACTiC first co-primary endpoint); enter it as `TACTIC_CORRECT_SELECTORS`
when confirmed against the publication's outcomes table.

The bootstrap uses 200,000 multinomial resamples with a fixed seed
(20260826); a 5,000,000-resample reference run (seed 42) gives
(11.82, 47.37), confirming the reported interval at one decimal place.

## Data provenance

All inputs are counts transcribed from publicly available sources cited in
the manuscript's reference list (the CREST and TACTiC publications and
supplemental material, the TACTiC statistical analysis plan, and
ClinicalTrials.gov records). No source text is reproduced.

## License

To be selected by the author before public release.
