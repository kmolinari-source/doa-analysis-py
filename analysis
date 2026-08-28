#!/usr/bin/env python3
"""
Consolidated analysis for:

  Denial of Access as an Unmeasured Outcome in Nonprescription Statin
  Self-Selection Studies (Molinari, manuscript in preparation)

Purpose. Recomputes every derived quantity reported in the manuscript from
the numerical inputs transcribed from the cited sources, and asserts that
each recomputed value matches the value reported in the manuscript at the
reported precision. The script exits with an error if any assertion fails,
so a clean run certifies the reported numbers under the environment below.

Environment certified. Python 3.12.3, NumPy 2.4.4, SciPy 1.17.1.
Figure 4 is generated separately by make_figure4.py (Matplotlib 3.10.8).

Data provenance. All inputs are counts transcribed from public sources
cited in the manuscript's reference list:
  [7]  Nissen SE, et al. J Am Coll Cardiol. 2021;78(11):1114-1123 (CREST).
  [22] Nissen SE, et al. J Am Coll Cardiol. 2024;83(21):2080-2088 (TACTiC),
       including its supplemental material.
  [24] AstraZeneca. TACTiC statistical analysis plan, edition 5.0 (success
       criterion and interval method).
No text from any source is reproduced here; numerical facts only.
"""

import sys
import numpy as np
from scipy import stats

PASS = 0


def check(label, ok, detail=""):
    global PASS
    status = "PASS" if ok else "FAIL"
    print(f"[{status}] {label}" + (f"  ({detail})" if detail else ""))
    if not ok:
        sys.exit(f"ASSERTION FAILED: {label} {detail}")
    PASS += 1


def cp_interval(k, n, alpha=0.05):
    """Two-sided Clopper-Pearson exact interval for k successes in n."""
    lo = 0.0 if k == 0 else stats.beta.ppf(alpha / 2, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(1 - alpha / 2, k + 1, n - k)
    return lo, hi


def cp_check(label, k, n, pct, lo_r, hi_r, pct_dp=1, ci_dp=1):
    lo, hi = cp_interval(k, n)
    p = round(100 * k / n, pct_dp)
    l, h = round(100 * lo, ci_dp), round(100 * hi, ci_dp)
    check(label, (p, l, h) == (pct, lo_r, hi_r),
          f"{k}/{n} -> {p}% ({l}-{h}); manuscript {pct}% ({lo_r}-{hi_r})")


print("=" * 72)
print("Section 1. CREST cross-tabulation (Table 2, Panel A; ref 7)")
print("=" * 72)
# Rows: participant outcome (OK to use, Ask a doctor, Not right for you)
# Cols: clinician outcome  (OK to use, Ask a doctor, Not right for you)
# Cells reconstructed from the published cross-tabulation; they reconcile
# exactly with all published marginals and the published discordance
# decomposition (3 incorrect selections, 14 incorrect rejections,
# 2 incorrect "ask a doctor" outcomes), as verified below.
X = np.array([[23, 0, 3],
              [2, 0, 0],
              [14, 1, 457]])
N = int(X.sum())
check("Grand total = 500", N == 500)
check("Clinician marginals = 39 / 1 / 460", list(X.sum(axis=0)) == [39, 1, 460])
check("Participant marginals = 26 / 2 / 472", list(X.sum(axis=1)) == [26, 2, 472])
check("Unmitigated concordance = 480 (96.0%)",
      int(np.trace(X)) == 480 and round(100 * np.trace(X) / N, 1) == 96.0)
cp_check("Unmitigated concordance CI", 480, 500, 96.0, 93.9, 97.5)
# Prespecified mitigation reclassifies the single participant-NRFY /
# clinician-AAD cell as correct: 480 + 1 = 481.
check("Mitigated concordance = 481 (96.2%)",
      480 + int(X[2, 1]) == 481 and round(100 * 481 / N, 1) == 96.2)
cp_check("Mitigated concordance CI", 481, 500, 96.2, 94.1, 97.7)
check("Discordance decomposition 3 + 14 + 2 = 19",
      int(X[0, 2]) == 3 and int(X[2, 0]) == 14 and int(X[1, 0]) == 2
      and N - 481 == 19)

print()
print("=" * 72)
print("Section 2. CREST accuracy measures (Table 2, Panel B; Section 3.2)")
print("=" * 72)
# Handling rule (primary, Section 2.3): positive = participant "OK to use";
# clinician-eligible = clinician "OK to use" (n = 39).
elig = int(X[:, 0].sum())          # 39
check("Clinician-eligible n = 39", elig == 39)
cp_check("Prevalence", 39, 500, 7.8, 5.6, 10.5)
tp = int(X[0, 0])                  # 23
cp_check("Sensitivity (primary rule)", tp, 39, 59.0, 42.1, 74.4)
cp_check("Sensitivity (secondary rule, direct outcomes only)", tp, 37, 62.2, 44.8, 77.5)
# Specificity: clinician-ineligible = 460 + 1 = 461; app-positive among them = 3
cp_check("Specificity", 458, 461, 99.35, 98.11, 99.87, pct_dp=2, ci_dp=2)
cp_check("Positive predictive value", 23, 26, 88.5, 69.8, 97.6)
cp_check("Negative predictive value, before mitigation", 457, 472, 96.8, 94.8, 98.2)
cp_check("Negative predictive value, with mitigation", 458, 472, 97.0, 95.1, 98.4)
cp_check("Null classifier ('not right for you' to all)", 460, 500, 92.0, 89.3, 94.2)
check("Reported 96.2% exceeds null baseline by 4.2 pp",
      round(96.2 - 92.0, 1) == 4.2)
ppv = 100 * 23 / 26
sen = 100 * 23 / 39
check("PPV - sensitivity = 29.5 pp", round(ppv - sen, 1) == 29.5,
      f"{ppv:.2f} - {sen:.2f} = {ppv - sen:.2f}")

print()
print("=" * 72)
print("Section 3. Bootstrap CI for PPV - sensitivity (Section 2.3)")
print("=" * 72)
# 200,000 multinomial resamples of 500 participants from the observed cell
# proportions of the published cross-tabulation; percentile method.
B = 200_000
rng = np.random.default_rng(20260826)
p = (X / N).ravel()
draws = rng.multinomial(N, p, size=B).reshape(B, 3, 3)
tp_b = draws[:, 0, 0]
pos_b = draws[:, 0, :].sum(axis=1)
elig_b = draws[:, :, 0].sum(axis=1)
valid = (pos_b > 0) & (elig_b > 0)
diff = 100 * (tp_b[valid] / pos_b[valid] - tp_b[valid] / elig_b[valid])
lo_b, hi_b = np.percentile(diff, [2.5, 97.5])
print(f"  bootstrap percentile CI: ({lo_b:.2f}, {hi_b:.2f}); B = {B}, "
      f"valid resamples = {int(valid.sum())}")
# Large-B reference (5,000,000 resamples, seed 42): (11.824, 47.368).
check("Bootstrap CI matches reported (11.8 to 47.4)",
      round(lo_b, 1) == 11.8 and round(hi_b, 1) == 47.4,
      f"({lo_b:.1f}, {hi_b:.1f})")

print()
print("=" * 72)
print("Section 4. Decision curve analysis (Sections 2.5, 3.2; Figure 4)")
print("=" * 72)
# Application: TP = 23, FP = 3. Treat-all: TP = 39, FP = 461. N = 500.
def nb(tp_, fp_, pt):
    return tp_ / N - fp_ / N * (pt / (1 - pt))

cross_all = 16 / 474      # solves nb_app(pt) = nb_all(pt)
cross_none = 23 / 26      # solves nb_app(pt) = 0  (equals PPV)
print(f"  crossover vs unscreened access: {100*cross_all:.2f}%")
print(f"  crossover vs no access:         {100*cross_none:.2f}%")
check("Crossover vs unscreened access = 3.38% (text: 3.4%)",
      round(100 * cross_all, 2) == 3.38 and round(100 * cross_all, 1) == 3.4)
check("Crossover vs no access = 88.46% (text: 88.5%)",
      round(100 * cross_none, 2) == 88.46 and round(100 * cross_none, 1) == 88.5)
check("Crossover vs no access equals PPV", abs(cross_none - 23 / 26) < 1e-12)
nb5, nb20 = nb(23, 3, 0.05), nb(23, 3, 0.20)
check("Net benefit of application at 5% = 0.0457", round(nb5, 4) == 0.0457,
      f"{nb5:.6f}")
check("Net benefit of application at 20% = 0.0445", round(nb20, 4) == 0.0445,
      f"{nb20:.6f}")
grid = np.linspace(0.05, 0.20, 301)
nb_app = nb(23, 3, grid)
nb_all = nb(39, 461, grid)
check("Application dominates both comparators across 5-20% band",
      bool(np.all(nb_app > nb_all) and np.all(nb_app > 0)))
check("Exchange rates: 5% -> 1:19 and 20% -> 1:4",
      round((1 - 0.05) / 0.05) == 19 and round((1 - 0.20) / 0.20) == 4)

print()
print("=" * 72)
print("Section 5. TACTiC co-primary endpoints (Table 3; Section 3.4; ref 22, 24)")
print("=" * 72)
cp_check("First co-primary, with mitigations", 1085, 1196, 90.7, 88.9, 92.3)
cp_check("First co-primary, without mitigations", 1007, 1196, 84.2, 82.0, 86.2)
cp_check("Second co-primary, with mitigations", 1132, 1154, 98.1, 97.1, 98.8)
cp_check("Second co-primary, without mitigations", 934, 1154, 80.9, 78.5, 83.2)
check("Mitigation counts: 1085-1007 = 78 and 1132-934 = 198",
      1085 - 1007 == 78 and 1132 - 934 == 198)
check("Success criterion met with mitigations (lower bounds > 85%)",
      cp_interval(1085, 1196)[0] > 0.85 and cp_interval(1132, 1154)[0] > 0.85)
check("Success criterion NOT met without mitigations",
      cp_interval(1007, 1196)[0] < 0.85 and cp_interval(934, 1154)[0] < 0.85)

def k_min(base, n):
    for k in range(0, n - base + 1):
        if cp_interval(base + k, n)[0] > 0.85:
            return k
    return None

k1, k2 = k_min(1007, 1196), k_min(934, 1154)
check("Minimum mitigations for first endpoint to pass = 34", k1 == 34, f"k={k1}")
check("Minimum mitigations for second endpoint to pass = 71", k2 == 71, f"k={k2}")
# Composition of second-endpoint mitigations (ref 22, supplemental tables)
stopped_eligible, stopped_unconf, muscle, liver, missing_ldl = 60, 73, 54, 12, 90
check("Composition percentages of analysis population (n = 1,154)",
      round(100 * stopped_eligible / 1154, 1) == 5.2
      and round(100 * stopped_unconf / 1154, 1) == 6.3
      and round(100 * missing_ldl / 1154, 1) == 7.8)
check("Unconfirmed-symptom split 54 + 12 <= 73", muscle + liver <= stopped_unconf)
check("Discontinuation mitigations 60 + 73 = 133 of 198",
      stopped_eligible + stopped_unconf == 133 and 133 <= 198)

print()
print("=" * 72)
print("Section 6. TACTiC screening flow (Section 3.3; ref 22)")
print("=" * 72)
check("12,624 - 10,332 = 2,292 not rejected", 12624 - 10332 == 2292)
check("Consent proportion 1,196/2,292 = 52.2%",
      round(100 * 1196 / 2292, 1) == 52.2)
check("Reference standard coverage 1,196/12,624 = 9.5%",
      round(100 * 1196 / 12624, 1) == 9.5)

print()
print("=" * 72)
print("Section 7. Introduction and Discussion derived values")
print("=" * 72)
check("36.6% untreated = complement of 63.4% treated (ref 3)",
      round(100 - 63.4, 1) == 36.6)

print()
print("=" * 72)
print("PENDING INPUT (assertion skipped, not failed)")
print("=" * 72)
print("  [SKIP] Null classifier floor for TACTiC ('at least 80.3%'):")
print("         requires the correct-selection component count of the first")
print("         co-primary endpoint from the ref 22 outcomes table. Enter it")
print("         as TACTIC_CORRECT_SELECTORS below when confirmed; the floor")
print("         is then that count / 1,196, reported as a lower bound.")
TACTIC_CORRECT_SELECTORS = None
if TACTIC_CORRECT_SELECTORS is not None:
    floor = 100 * TACTIC_CORRECT_SELECTORS / 1196
    check("TACTiC null-classifier floor >= 80.3%", round(floor, 1) >= 80.3,
          f"{floor:.1f}%")

print()
print("=" * 72)
print(f"ALL {PASS} ASSERTIONS PASSED")
print("=" * 72)
