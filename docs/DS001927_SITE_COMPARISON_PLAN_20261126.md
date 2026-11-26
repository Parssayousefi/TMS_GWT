# Does stimulation site change brain information processing?

Analysis decision, 26 November 2026. The user selected this question before
ds001927 denoising, deconvolution or PID outcomes. This is a dated analysis
plan, not an externally registered protocol. Earlier Aachen and SUDMEX results
have already been inspected and remain exploratory context.

## Main comparison

For each person, estimate

`(FRO_post − FRO_pre) − (CTR_post − CTR_pre)`.

This asks whether the change after prefrontal stimulation differs from change
after stimulation of the temporoparietal control site. CTR receives active
stimulation; it is not sham or an untreated control. The original study used
three counterbalanced sessions and individually identified targets. FRO is a
dorsoanterior prefrontal target and is not assumed to equal the earlier cohorts'
DLPFC targets. See [Castrillon et al.](https://doi.org/10.1126/sciadv.aay2739).

The experimental unit is the participant, not a scan, region or region pair.
Retain a person for this contrast only when all four FRO/CTR scans pass the
declared scan QC and anatomical/signal review. An OCC failure alone does not
exclude a valid FRO/CTR comparison. Do not substitute a different scan or
impute missing outcomes. Preserve session and task identities explicitly.

## Two primary endpoints, one family

Use the existing Gaussian MMI PhiID estimator on CONN-denoised, rsHRF-deconvolved
regional signals at lag 1 sample (1.25 s). Average over the same unordered,
non-self regional pairs in every included scan:

1. **Redundancy-corrected integration, ΦR:**
   `I([X,Y]past;[X,Y]future) − I(Xpast;Xfuture) − I(Ypast;Yfuture) + rtr`.
   This is the pipeline's `phi_r`, distinct from uncorrected Φ.
2. **Six-mode information storage:** `rtr + xtx + yty + sts`.
   This is storage as a PhiID processing mode, not just a region's ordinary
   autocorrelation or the sum of its own-past mutual information.

These quantities are in nats before conditional-null normalization. Require
finite, nonconstant inputs and an identical complete participant set for both
primary endpoints. No directional hypothesis is imposed.

Use the participant's paired contrast for a two-sided studentized joint
sign-flip test with single-step maximum-|t| family-wise correction across
**these two endpoints** (alpha .05). Enumerate all sign patterns in chunks for
at most 23 participants. Apply the same sign to both endpoints of a person.
The sign-flip reference requires symmetric/exchangeable participant contrasts
under the null; it is not a claim of access to the original randomization
sequence. Report the mean contrast, marginal 95% t interval, unadjusted and
family-wise p-values, and all four condition means. Label marginal intervals
as such. Fewer than six valid complete participants is descriptive only.

No additional continuous FD predictor is included. Motion is handled by the
declared scan exclusions and the production nuisance regressors.

## Coverage is part of validity

The atlas begins with the existing 454 parcels. Before computing PID outcomes,
freeze one common parcel set using the intersection of parcels with at least
50% functional-mask overlap in **every required FRO/CTR scan of every retained
participant**. Use exactly that set in all four conditions and all people.
Report the surviving count, locations and pair count; describe the result as
integration/storage across the analyzed parcels, rather than implying complete
anatomical coverage. This is a study QC rule, not a universal coverage cutoff.

Inspect low-coverage raw and processed images to distinguish acquisition signal
loss from extraction or alignment mistakes. Complete denoised time-series
variance, finite-value and HRF checks before releasing PID. Nonzero background
noise does not establish valid signal. Any unexpected zero-variance or nonfinite
signal in the common set blocks release for investigation; do not silently
replace it with zero or silently change the set after outcome inspection.

Use a stricter 80% coverage intersection as a sensitivity with its parcel count
reported. Do not repeatedly choose coverage thresholds to recover significance.
An excessively depleted or anatomically implausible set requires a documented
methods review before outcomes, rather than automatic publication inference.
The original 454-parcel summaries may be retained as diagnostics only when
computable; known unsupported parcels do not qualify as a primary result.

## Bounded secondary questions

- **Information-processing profile:** the six disjoint modes (storage, transfer,
  copy, erasure, joint-to-individual, individual-to-joint) for FRO−CTR and
  OCC−CTR differences in pre/post change. Use the same all-six-run-complete
  participants and coverage set for this secondary analysis and one maximum-|t|
  family of 12 tests. The primary storage endpoint is repeated here for context;
  significance in this family does not substitute for the primary family.
- **Aachen-motivated follow-up:** total source synergy (`str+stx+sty+sts`)
  outside the frozen HCP workspace and the direct outside-minus-inside
  interaction, each for FRO−CTR change. One two-test family. Intersect the fixed
  workspace labels with the coverage-eligible parcel set; do not re-estimate a
  workspace from treatment effects. This is a new-cohort test of an existing
  exploratory candidate, not proof that the different protocols are identical.

OCC, workspace subdivisions and individual atoms cannot be promoted to primary
endpoints because their results are more attractive. Report the defined
families together. Clinical/subjective associations are outside this plan;
these data do not supply the earlier craving/impulsivity outcomes.

## Sensitivities and interpretation

For the two primary endpoints retain, in fixed order: the stricter mean-FD
gate (<=.20 mm), the 80% coverage set, lag 2 samples (2.5 s), single-participant
omissions, and the existing conditional-null normalization procedure rebuilt
from passing baseline records with the actual TR/retained lengths. Record null
calibration, acceptance and interpolation checks before using normalized
scores. Do not pool reference distributions across cohorts or select between
raw and normalized results by p-value. Normalization changes the estimand;
disagreement limits a claim about the allocation of information independent
of total predictability. Sensitivity results do not replace a failed primary
test or generate independent headline claims.

Lag units differ from the earlier 2-s cohorts. Their effects will be shown
separately, without pooling raw magnitudes or calling a shared direction a
replication. The information-decomposition atoms describe lagged statistical
dependence; even a positive experimental site contrast does not identify a
direct neuronal causal pathway. A null contrast will not be described as
workspace stability or equivalence.

If supported, the brief-report story is that **the site of stimulation matters
for integration or storage**, with a direct within-person active-site
comparison. That is the research question, not an observed ds001927 finding.
