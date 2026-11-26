# TMS exploratory findings — 25 November 2026

The subsequently requested baseline-versus-average-of-all-follow-ups analysis
is documented separately in [TMS_ALL_POST_AVERAGE_20261125.md](TMS_ALL_POST_AVERAGE_20261125.md).
The results below describe the acute SUDMEX and paired Aachen comparisons.

The clearest new neural finding is an increase in whole-brain mean pairwise
whole-minus-sum Phi (Φ) in Aachen, averaged across the two iTBS intensity visits.
It passes family-wise max-T correction across the five global information
measures in both the raw and conditional-null-normalized analyses. This is a
candidate exploratory paper result, not evidence that TMS increases consciousness
or reproduces the psychedelic workspace signature.

## What was averaged

| Dataset | Main participants | Pre scans | Post scans | Analysis |
|---|---:|---:|---:|---|
| SUDMEX | 40: 22 active, 18 sham | 40 t0 | 40 t1 | One pre/post pair per person; active–sham comparison adjusted for that outcome's baseline |
| Aachen | 15 | 15 pre90 + 15 pre120 | 15 post90 + 15 post120 | Compute both matched changes, then average within each participant |

Aachen's main response is `((post90-pre90)+(post120-pre120))/2`. Each post scan
is paired with the baseline from its own visit. The sample size is 15 people,
not 30 independent observations. No fMRI time series, covariance matrices, or
PhiID matrices were averaged before participant metrics were computed. Averaging
the two metric changes equals subtracting the average of the two baseline
metrics from the average of the two post metrics, because all 15 people have
complete observations. The repeated-measures time effect gives the identical
parametric test: F(1,14)=9.0556, p=.009378.

For Phi, mean changes were +0.002043 nats at 90% (95% CI −0.000254 to 0.004341,
uncorrected p=.0771) and +0.001541 at 120% (−0.000624 to 0.003706, p=.1492).
The intensity difference was −0.000503 (−0.004164 to 0.003159, p=.7728).
The two visit changes correlated −.346 across participants; their average is
less variable than either visit alone. This supports reporting an average
across visits, not claiming independent replication at both intensities or
established equivalence of intensities. Eleven of 15 average changes were positive.

All usable completed scans were: SUDMEX t0=41, t1=43, t14=12, t2=19, t3=13,
t4=5 (133 total), and Aachen 15 at each of four sessions (60 total).
The SUDMEX main comparison requires both t0 and t1 and therefore uses 80 scans.
Later visits were not averaged into it. Eleven originally sham participants had
usable t1/t14 pairs and complete three-visit trajectories; delayed-active
comparisons produced no supported neural result. This open-label phase is
confounded with period, expectancy, and retention.

## Correction and surviving findings

All analyses are exploratory. Family-wise error is controlled within scientific
outcome families, separately for each dataset's main comparison: five global
information measures, two mean participation coefficients, the original FPN
gateway endpoint, seven other network gateway counts, eight broadcaster counts,
and two total-role counts. This controls error within each family, not across
the entire manuscript. The original FPN focus comes from the initial TMS
proposal. Secondary intensity discoveries are corrected across the 90%, 120%
and intensity-difference comparisons together within their outcome family.

Aachen uses two-sided, exact single-step max-|t| correction over all 32,768 joint
participant sign flips. The same sign applies to all of a person's outcomes,
preserving their dependence. Validity requires appropriate joint symmetry of
the paired errors; the sign flips do not supply a randomized sham comparison.
Holm correction of parametric p values is shown alongside max-T. SUDMEX uses
baseline-adjusted HC3 regression with within-family Holm correction. The seven
clinical scales form one family; the 35 brain–clinical associations form another.
We do not select the smaller correction separately for each endpoint.

| Aachen global measure | Average change (95% marginal CI) | Family max-T p | Family Holm p |
|---|---|---:|---:|
| Raw Φ | +0.001792 nats (0.000515 to 0.003070) | **.02393** | **.04689** |
| Raw ΦR | +0.002220 nats (0.000356 to 0.004085) | .06580 | .09169 |
| Normalized Φ | +0.7123 percentile points (0.1630 to 1.2617) | **.04828** | .07360 |
| Normalized ΦR | +0.7980 percentile points (0.1274 to 1.4687) | .07690 | .09206 |

Only Phi survives the global family max-T correction. No other main neural
family contains a surviving effect, including gateway counts and participation
coefficients. No extended network, subcortical, clinical-moderation, or multivariate
finding survives its family correction. With all 25 Aachen neural metrics in
one family, raw Phi max-T p=.1434 and normalized Phi p=.2246; family scope must
therefore be clear. The raw Phi paired Hedges' g_z is 0.735.

For baseline-adjusted SUDMEX clinical outcomes, cognitive impulsivity decreased
by 4.776 points (Holm p=.03560) and motor impulsivity by 4.931 (p=.04198).
Total BIS impulsivity decreased by 11.564 (p=.05290) and craving VAS by 1.582
(p=.07489); neither passes family-wise correction. These clinical data and
benefits were already reported by the original trial authors. Their reanalysis
provides clinical context, not the new paper's principal novelty. Brain–clinical
associations remain unsupported.

## Baseline, motion, and robustness

The SUDMEX FPN change contrast was +3.141 gateways before adjustment. Active
participants started with 13.955 FPN gateways versus 16.556 in sham. Adjustment
for baseline alone reduced the contrast to +1.207 (95% CI −1.990 to 4.404,
p=.449); motion alone left it at +3.068. Joint adjustment gave +1.085
(−2.130 to 4.300, p=.498). Baseline imbalance accounts algebraically for 94.5%
of the fitted attenuation, versus 5.5% for measured mean-FD change. Those are
regression contributions, not causal percentages. Regression to the mean is
compatible with this pattern but is not proved. Mean FD cannot exclude every
form of residual motion or physiological noise.

Aachen's raw Phi average stayed positive in every leave-one-person-out analysis
(+0.001392 to +0.002030 nats). Each uncorrected leave-one-out p was below .021;
the two-sided sign-flip p was .00452 and the Wilcoxon p was .00537. These are
dependent robustness checks, not additional independent replications.

Motion adjustment barely changed the raw Phi estimate (+0.001787), but its
HC3 uncertainty was larger: nominal p=.03329, five-global-measure Holm p=.16647.
Normalized Phi behaved similarly (estimate +0.7191 percentile points, nominal
p=.03417, Holm p=.17087). Thus there is little measured-motion attenuation of
the effect size, but corrected support is not retained by every model. No sham
condition in Aachen distinguishes stimulation effects from common pre/post
effects such as time, arousal, scanner drift, or repeated measurement.

## What normalization adds

The [NuMIT paper](https://doi.org/10.1371/journal.pcbi.1013629) compares information
components with a null distribution matched for total mutual information.
It demonstrates PID; our pairwise PhiID application is a methodological
extension. It belongs after PhiID estimation and before participation/role
calculation, not as a rerun of MRI preprocessing or spatial normalization.

We used a shared baseline-derived VAR coefficient pool within each dataset,
the same reference for both conditions, and two independent 32,768-draw null
references per TDMI grid point. All 140 main scans and all 25 endpoints were
processed. Interpolated CDFs differed from direct held-out CDFs by at most
0.000845. Independent scipy Lyapunov solutions matched, and 1,536 finite-record
simulations assessed calibration at retained lengths of 251 and 295 samples.
Finite-record bias remains: for example, average synergy percentile shifts
relative to population values were about −.020 to −.026. This validates useful
numerical properties, not the full generative adequacy of a VAR(1) null for
filtered, deconvolved fMRI. Null-pool uncertainty is not propagated.

Normalized Phi's max-T p values were .04852 and .04822 in the two independent
reference seeds, versus .04828 combined. Its paired-change Monte Carlo error
was small relative to between-person variation. Some discrete role counts
were less stable. Normalization therefore supplies a consistent, borderline
sensitivity for Phi, not an independently replicated discovery. Quantiles are
not nats, and separately normalized atoms/composites are not additive.

## Comparison and proposed manuscript

[Luppi et al.](https://elifesciences.org/articles/88173/figures) emphasized
redundancy-corrected Phi_R and its loss in workspace connections during
unconsciousness. The original whole-minus-sum Phi could instead increase during
unconsciousness. Our whole-brain Phi increase is consequently not evidence of
enhanced consciousness. The closer Phi_R comparison is positive but does not
pass the chosen family correction, and our spatial averaging also differs.

[Yousefi et al.](https://www.biorxiv.org/content/10.64898/2026.09.28.755172v1.full)
reported FPN/DAN gateway recruitment and changes in the distribution of
synergistic and redundant participation. We have not established that pattern
in these TMS samples. A difference in statistical significance does not prove
a difference between interventions; the data do not establish a neural–clinical
dissociation or a shared mechanism. We use published comparisons only, without
private psychedelic data or maps.

**Working title:** Temporal information balance following prefrontal stimulation:
an exploratory PhiID study.

**Central claim:** Across two iTBS visits in healthy volunteers, pre/post changes
in mean pairwise whole-minus-sum Phi provide exploratory evidence of a shift in
temporal information balance, including under conditional-null normalization.
Parallel analysis of a randomized clinical rTMS dataset does not establish
frontoparietal gateway recruitment. These results motivate separating global
information quantities from the regional organization of the synergistic
workspace when studying stimulation.

This supports a cautious exploratory report. It is not yet a strong causal
mechanism paper: the clearest signal is in the small uncontrolled dataset,
depends on averaging visits and the correction family, is weaker for Phi_R,
and loses corrected support in the HC3 motion sensitivity. A sham-controlled
replication of the stated global family and both intensities would be the
decisive next study. Further searching this dataset for a preferred significant
endpoint is not needed to make the present result interpretable.

## Review artifacts

- Original completed publication run: `output/tms_publication/20261125T145928_921789Z`.
- Baseline/motion diagnostics: `output/tms_diagnostics/20261125T152110_809161Z`.
- Normalized branch: `output/tms_numit/20261125T155252_384543Z`.
- Final family-wise tables and numerical verification: `output/tms_familywise/20261125T161051_102285Z`.
- Figures and scan accounting: `output/tms_exploratory_report/20261125`.
- Results notebook: `output/jupyter-notebook/TMS_Exploratory_Findings.ipynb`.

Independent checks verified 977 estimable neural contrasts, 100 delayed-active
contrasts and seven baseline-adjusted clinical models. One secondary amygdala
intensity model had undefined HC3 covariance due to unit leverage and is marked
unestimable in the final table. None of the main findings depends on that row.
Production processing, original outputs, frozen cohorts and D:/projects/GWT
were not modified.
