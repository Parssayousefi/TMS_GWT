# Baseline versus the average of all post measurements

The requested analysis is complete for both datasets, using raw and normalized
metrics. It gives a simpler summary of the available follow-ups. In SUDMEX it
answers a different question from the acute randomized treatment comparison;
in Aachen it is exactly equivalent to the existing analysis.

## Averaging and participants

For every metric and participant, the response is:

`mean(all eligible post-visit values) - mean(all eligible pre-visit values)`.

Each post visit receives equal weight within a person. Each person then receives
equal weight in the group analysis, regardless of how many scans they have.
The averages are of already computed scan metrics, not time series or matrices.
No missing visits are imputed.

| Dataset | Participants | Pre scans | Post scans | Visits included |
|---|---:|---:|---:|---|
| SUDMEX | 41: 22 originally active, 19 originally sham | 41 | 88 | t0 versus every available t1, t14, t2, t3 and t4 |
| Aachen | 15 | 30 | 30 | Mean(pre90, pre120) versus mean(post90, post120) |

SUDMEX contributors include 40 t1 scans, 12 t14 scans, 18 t2 scans, 13 t3 scans
and 5 t4 scans. The other four usable follow-up scans in the full dataset lack
a usable baseline for this paired analysis and therefore cannot contribute.
Participant 046, originally sham, adds a valid t0/t14 pair despite an excluded
t1 scan. Repeating the analysis in the original 40-person acute sample isolates
the averaging change from this one-person expansion.

| Number of post scans averaged | Initially active | Initially sham |
|---:|---:|---:|
| 1 | 9 | 6 |
| 2 | 5 | 8 |
| 3 | 6 | 1 |
| 4 | 2 | 2 |
| 5 | 0 | 2 |

Aachen has complete data at both visits, so subtracting the two condition means
is algebraically identical to averaging the two matched pre/post differences.
It still provides 15 independent participant responses.

## Neural results

The existing scientific families and two-sided family-wise correction were
retained. Aachen uses all 32,768 joint sign flips. SUDMEX uses 199,999 sampled
sign vectors plus the observed configuration; within each participant, all
metrics receive the same sign. Max-T controls within-family error under the
sign-symmetry assumptions. Within-family Holm results are also retained.

Three SUDMEX raw pre/post effects survive in the full 41-person sample:

| Measure | Mean change | Marginal 95% CI | Family max-T p | Family Holm p |
|---|---:|---|---:|---:|
| Persistent synergy | +0.018150 nats | 0.008950 to 0.027349 | .001015 | .001381 |
| Total synergy | +0.008454 nats | 0.004046 to 0.012862 | .001535 | .001542 |
| Mean workspace synergy participation | −0.000282 | −0.000486 to −0.000077 | .015455 | .016309 |

None survives in the normalized full-sample analysis. Global Phi changes by
−0.000497 nats (−0.002847 to 0.001854; family max-T p=.9687), and Phi_R by
−0.000214 (−0.002792 to 0.002364; p=.9980). FPN gateway count changes by
−1.020 (−2.433 to 0.393; p=.1521). No gateway or broadcaster endpoint survives
its main family correction. These full-cohort temporal changes combine the
original arms and cannot be attributed specifically to active stimulation.

The same three raw findings survive when restricted to the original 40 people.
That sensitivity also yields a borderline decrease in normalized total synergy:
−0.012727 conditional-null quantile units, max-T p=.04644, Holm p=.06164.
It does not pass max-T in the full 41-person sample (p=.10936). This is
cohort-sensitive evidence, not a robust new finding or a reason to exclude
the additional eligible participant.

Aachen is unchanged: raw Phi increases by 0.001792 nats (family max-T p=.02393),
and normalized Phi by 0.007123 quantile units, or 0.7123 percentile points
(p=.04828). Neither Phi_R nor the role-count/participation outcomes passes
its main family correction. The previously documented lack of sham control
and sensitivity to the inferential model remain relevant.

## Original-arm comparison and clinical outcomes

In SUDMEX, the unadjusted difference between the original active and sham
groups in FPN changes is +3.534 gateways (p=.00940). Adjusting the all-post
outcome for its baseline reduces that to +1.636 (95% CI −0.711 to 3.983,
p=.1664); adding mean FD change gives +1.739 (−0.671 to 4.149, p=.1523).
The original 40-person sample gives the same qualitative conclusion.
No baseline-adjusted original-arm neural contrast survives its outcome-family
correction in either the raw or normalized analysis.

Clinical measures were averaged at the same MRI visits for each participant.
For an outcome to contribute, baseline and every selected follow-up score
had to be present and pass the existing source-quality rules. All seven
clinical scales decrease in the overall pre/post analysis after Holm correction
(37–40 participants depending on the scale). For example, craving VAS decreases
by 1.629 points (Holm p=.00651), and total BIS impulsivity by 8.775
(p=.00139). No baseline-adjusted contrast between the original arms survives
the seven-scale correction. No new brain–clinical association search was run.

## Why this is not a cleaner treatment test

The [SUDMEX protocol](https://openneuro.org/datasets/ds003037/versions/2.1.0)
has a two-week randomized phase, followed by open-label treatment and maintenance.
The t14 visit occurs after active treatment for participants who originally
received sham. Averaging all follow-ups therefore mixes treatment histories
and scheduled times from two weeks to twelve months. People with one follow-up
and people followed for a year do not represent the same post-treatment window.
Retention can also depend on clinical progress.

Thus this mean estimates a change over each person's observed follow-ups. It
does not estimate a common-duration treatment effect, and the label “originally
sham” does not describe a continuing sham exposure throughout follow-up.
The baseline-to-two-week active-versus-sham comparison remains the clearer
test of the randomized intervention. The all-post average is a useful
longitudinal sensitivity and descriptive summary alongside it.

The temporal-information audit also clarifies the raw synergy increases.
Average pairwise temporal mutual information rises by 0.022894 nats in SUDMEX;
twice mean single-region self-information rises by approximately 0.023390.
Their difference is the nearly unchanged Phi (−0.000497). Normalization removes
support for the main raw synergy/participation effects. This is consistent
with a change in the amount of predictable temporal information rather than
clear evidence for stronger relative integration; it does not establish the
biological or artifactual cause of that change.

## Files and verification

All outputs are isolated in
`output/tms_all_post_average/20261125T162521_698094Z/`:

- `scan_contributions.csv` records each scan and its within-person weight.
- `raw_participant_changes.csv` and `normalized_participant_changes.csv`
  contain all participant averages and differences.
- `neural_pre_post_results.csv` retains every endpoint, correction, and cohort.
- `original_arm_contrasts.csv` retains unadjusted and adjusted comparisons.
- `clinical_all_post_results.csv` and `clinical_matched_visit_coverage.csv`
  retain all clinical results and reasons for missing outcomes.
- `temporal_information_descriptive_summary.csv` contains the information
  balance audit; no additional hypothesis tests were added by that audit.
- `independent_verification.json` records checks of all 150 paired tests,
  300 original-arm contrasts, 14 clinical models, input hashes, visit weights,
  and an independent SciPy reproduction of the main max-T results.

Normalization was extended to 49 additional needed scans using the existing
reference distributions and unchanged baseline pools. None of the original
production outputs or earlier analyses was replaced. The executed exploratory
notebook includes this analysis as a separate final section.
