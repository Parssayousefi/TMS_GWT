# SUDMEX FPN adjustment diagnostics and NuMIT assessment

These are post-hoc diagnostics requested after inspection of the completed
primary analysis. They neither replace that analysis nor constitute a
preregistration. All 40 eligible SUDMEX T0/T1 pairs were retained: 22 active,
18 sham. No production processing settings, source data, frozen HCP workspace,
or earlier result tables were modified.

Results directory: `output/tms_diagnostics/20261125T152110_809161Z/`.
Scripts: `setup/analyse_sudmex_fpn_diagnostics.py` and
`setup/assess_tms_numit_feasibility.py`. The diagnostic plan is retained in
`PROVENANCE/SUDMEX_FPN_DIAGNOSTIC_PLAN_20261125.json`.

## Baseline imbalance accounts for most statistical attenuation

| Arm | Baseline FPN gateways | Two-week FPN gateways | Change |
| --- | ---: | ---: | ---: |
| Active | 13.9545 | 13.8182 | -0.1364 |
| Sham | 16.5556 | 13.2778 | -3.2778 |

The sham group began 2.6010 gateways higher; the standardized baseline
difference (active minus sham) was -0.7666. At follow-up the difference was
only +0.5404. The positive change-score contrast therefore reflects different
starting points and convergence, rather than an observed active-arm increase.
This is compatible with regression toward the mean but does not prove its
biological or measurement-error mechanism.

All rows below use identical participants and HC3 uncertainty with t intervals,
so the covariate effects can be compared without changing error estimators.

| Adjustment | Active-minus-sham effect | 95% CI | p |
| --- | ---: | --- | ---: |
| None, HC3 | 3.1414 | -0.0579 to 6.3407 | .05408 |
| Baseline FPN only | 1.2067 | -1.9904 to 4.4037 | .44928 |
| Change in mean FD only | 3.0682 | -0.3323 to 6.4686 | .07559 |
| Baseline FPN and change in mean FD | 1.0851 | -2.1301 to 4.3003 | .49806 |
| Baseline FPN, baseline FD and change in FD | 1.3094 | -1.8873 to 4.5061 | .41130 |

The original Welch test remains p = .04861, CI 0.0202 to 6.2626. The HC3
unadjusted p = .05408 is a different uncertainty calculation for the same
3.1414-point estimate; it is not a data change. Statistical significance at
.05 was already sensitive to a reasonable error-estimation choice.

For the joint model, the exact algebraic identity is:

`raw contrast = adjusted arm coefficient + baseline slope × baseline arm difference + FD-change slope × FD-change arm difference`.

The raw-to-adjusted reduction of 2.0563 gateways comprises 1.9430 from the
baseline term and 0.1133 from the FD-change term: approximately 94.5% and 5.5%
of the attenuation. These percentages describe this fitted regression; they
are not causal shares or proportions of brain activity explained. Baseline-only
adjustment already yields an inconclusive treatment contrast, so the finding
does not depend on adjusting for a post-treatment motion variable.

The joint FD-change slope is -3.5644 gateways/mm (CI -24.9524 to 17.8236).
Mean FD changes were -0.0182 mm in active and +0.0136 mm in sham. Including
baseline FD as well does not restore the FPN contrast. This addresses measured
mean FD; it does not rule out all residual motion or imaging artifacts.

The baseline/follow-up correlation was .168 overall and .274 in sham. These
are descriptive between-visit correlations, not a validated test-retest
reliability estimate. The negative baseline/change correlation cannot by
itself establish regression toward the mean because change subtracts baseline.
Baseline ANCOVA is appropriate for the baseline-adjustment question; change
regressed on baseline gives exactly the same arm coefficient as follow-up
regressed on baseline. See [Vickers and Altman](https://pmc.ncbi.nlm.nih.gov/articles/PMC1121605/)
and [Twisk et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC5898524/).

## Influence, resampling, and clinical outcomes

All 40 leave-one-participant-out joint-model estimates were positive but
inconclusive: +0.6233 to +1.7382 gateways, p = .2684 to .6938. None became
significant; the attenuation is not caused by one participant. The original
Welch test remained below .05 in only 13 of the 40 deletions. These are
influence diagnostics, not grounds to exclude anyone.

A 5,000-replicate bootstrap, resampling within the original treatment arms,
gave a percentile interval of -1.6884 to +3.9749 for the joint-adjusted
coefficient and -1.5807 to +4.1378 for baseline-only adjustment. Both include
zero. The unadjusted bootstrap interval was +0.1510 to +6.0455; this supports
describing the unadjusted result as fragile rather than claiming that every
reasonable calculation is nonsignificant.

As a separate exploratory check, all seven clinical outcomes were adjusted
for their own baseline value, with the same seven-outcome BH correction and
source-quality masks. Craving VAS remained lower under active treatment
(-1.5821 points, CI -2.8857 to -0.2786; q = .03276; n = 40). BIS total also
remained lower (-11.5640 points, CI -20.2602 to -2.8678; q = .02468; n = 39).
Cognitive and motor BIS subscales survived correction; CCQ-General, CCQ-Now
and nonplanning did not. This supports clinical differences within the imaging
cohort, without establishing their mediation by workspace metrics. The earlier
35 adjusted brain-clinical associations remain unsupported after correction.

## What the NuMIT paper supports

[Liardi et al. (2025)](https://doi.org/10.1371/journal.pcbi.1013629) compare
information atoms with null systems matched on total mutual information,
expressing each observed atom as a percentile in its corresponding null
distribution. Their demonstrations concern PID, including VAR models of MEG;
their discussion explicitly identifies PhiID as a future extension. The
supplement also describes empirical, data-informed nulls. A PhiID application
is therefore a methodological extension requiring validation, not an already
validated drop-in step for our fMRI gateway endpoint.

Normalization changes the question from information quantity to its position
relative to a chosen matched-information null. It does not correct baseline
arm imbalance, head motion, or a weak treatment contrast. A percentile is not
an amount of information in nats or a randomized treatment-effect p-value.

Participation coefficients are invariant to uniform positive rescaling of
all weights within a given matrix. A numerical check on a TMS scan confirmed
unchanged PCs and gateway/broadcaster counts after rescaling synergy by 100
and redundancy by .01. Thus dividing a whole scan by a single global number
cannot fix its FPN gateway count. Edge-specific nonlinear NuMIT transformations
can change relative weights and roles, defining a distinct sensitivity endpoint.

## TMS temporal-information audit and bounded pilot

The audit covered all 140 scans contributing to primary paired analyses:
80 SUDMEX scans and 60 Aachen scans. All five analytic Gaussian PhiID summaries
matched production matrices on 256 fixed sampled edges per scan, with maximum
absolute discrepancy 1.55e-14. The independent identity
`mean Phi = mean pairwise TDMI - 2 × mean single-region self-information`
held within 4e-16. No old psychedelic driver or private result analysis ran.

Mean TDMI was .8013 nats in SUDMEX and .7378 in Aachen; mean lag-one
autocorrelation was .7303 and .7117 respectively. Raw synergy and TDMI were
strongly correlated across scans within both datasets (r approximately .95),
a descriptive association with repeated people, not an independent-subject
significance test. This makes conditional-information robustness scientifically
relevant. Nevertheless, the SUDMEX active-minus-sham TDMI change was only
.00223 nats (CI -.02426 to .02873, p = .8653), providing no clear evidence of
a gross mean-TDMI treatment shift behind the unadjusted FPN result. Mean TDMI
cannot exclude heterogeneous edge-specific changes.

The method pilot used eight scans selected by subject order: the first SUDMEX
subject in each arm at both visits, and the first Aachen subject at all four
visits. It used 256 fixed sampled edges and 512 null draws per edge, with no
normalized group-effect tests. An empirical VAR coefficient pool was built
from all eligible baseline scans, independently of treatment assignment;
the same dataset-specific pool was used across arms and visits. A Gaussian
coefficient-pool comparison was also evaluated.

Under the generic Gaussian coefficient pool, 36.43% of sampled synergy
percentiles were exactly 0 or 1; the shared baseline pool reduced this to
1.12%. Phi_R saturation fell from 18.36% to .93%. This supports an empirical
null as a candidate for further validation. Rare draws missed the fixed TMI
tolerance of 1e-6 nats (at most about .0023% in an individual pilot batch);
the pilot excluded those draws from percentile denominators and recorded
their counts through valid fractions. A production implementation should
explicitly reject and replenish invalid draws.

Two independent seeds at 512 draws produced mean absolute per-edge quantile
differences of .0146 for synergy and .0176 for redundancy, with maxima .0859
and .1055. Full-graph gateway-count convergence is therefore an essential
remaining check. Population-level synthetic-null calibration was broadly
consistent with uniform quantiles (KS D .039-.127); individual nominal
calibration p-values were sometimes below .05 across the 30 correlated checks.
This is preliminary calibration, not proof of validity for finite fMRI records.
Finite-record and full-graph gateway stability were not tested by this pilot.

## Implementation cautions and recommended next step

The copied NuMIT helpers are explicitly marked experimental. The old
`run_numit.py` driver is unsuitable to run unchanged: it is configured for
private drug cohorts, uses only 100 draws, and passes each current edge chunk
as the empirical A pool instead of the full eligible pool it constructs.
Consequently its null reference depends on chunk boundaries and ROI ordering.
The new bounded TMS pilot avoids that dependency by passing a fixed shared
baseline pool. No edits were made to the old driver or mathematical helpers.

Recommended sequence before a full TMS NuMIT sensitivity analysis:

1. Freeze the TMS-only null construction and estimands. Keep one common null
   reference across compared visits/arms, one fixed lag of 2 s, the existing
   235-region workspace, and the existing statistical contrasts. Justify the
   shared baseline-pool adaptation; use alternative defensible nulls as
   declared method sensitivities, not select a null by treatment significance.
2. Validate finite-record estimation at the observed scan lengths and temporal
   dependence. Check null matching, invalid-draw handling, calibration, ROI
   ordering/chunk invariance, and Monte Carlo convergence of full PC/rank/count
   outputs. Common random draws across paired visits can reduce numerical noise.
3. Once validated, run the normalized branch consistently for all paired TMS
   scans. Normalize each composite (Phi/Phi_R) against its own null distribution;
   do not sum component percentiles or expect the raw additive identities to hold.
   Retain both raw and normalized branches, full baseline adjustment, the original
   comparisons and multiplicity families. Report normalized outcomes as a new
   post-hoc sensitivity analysis, not as a replacement of the Yousefi endpoint.

The manuscript should report the original nominal FPN result and its loss of
support under baseline adjustment prominently. A defensible current claim is
clinical improvement with uncertain workspace reconfiguration. It is premature
to claim a demonstrated workspace mechanism, a reproduction of the psychedelic
signature, a proven clinical-neural dissociation, or absence of a meaningful
neural effect. The adjusted intervals still allow a gateway increase of about
four regions. NuMIT would strengthen a robustness assessment; it cannot supply
missing evidence merely by changing the outcome definition.

## Verification and artifacts

- `fpn_model_comparison.csv`, `coefficient_attenuation_decomposition.csv`,
  `leave_one_out.csv`, `bootstrap_intervals.csv`, and
  `clinical_baseline_adjusted.csv` contain the diagnostic estimates.
- `adjustment_comparison.png/.pdf` and
  `baseline_motion_diagnostics.png/.pdf` were visually inspected.
- `numit_feasibility/` contains per-scan temporal-information metrics,
  pilot null quality/repeatability/calibration tables and input hashes.
- All five diagnostic regression estimates and intervals were checked with an
  independent NumPy HC3 covariance calculation. The original joint-adjusted
  result was reproduced exactly; ANCOVA equivalence and the algebraic
  attenuation identity passed. All source input hashes remained unchanged.
