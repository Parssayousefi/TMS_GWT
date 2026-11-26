# Current analysis status

The manuscript remains in development. The goal is one information-decomposition
treatment story supported across Aachen and SUDMEX. Haemodynamic associations
and repeatability are not the selected narrative. Nonsignificant change does
not establish an unchanged workspace or equivalence.

## DAN gateway audit

All 60 corrected Aachen scans were used to reconstruct the two matched visit
changes for each of 15 people, then average within person. The four-outcome
family is FPN gateways, DAN gateways, redundancy participation and synergy
participation. All counts, t intervals, Holm values and 32,768-pattern exact
joint sign-flip results reproduce the report.

| Representation | Mean DAN change | 95% marginal interval | t p | Exact sign-flip p | Holm4 p | maxT4 p |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| Raw | −0.400 | [−0.893, 0.093] | .1038 | .1367 | .3114 | .3888 |
| Normalized | −0.467 | [−0.904, −0.029] | .0382 | .0566 | .1527 | .1621 |

Normalized DAN gateway counts average 1.433 before and 0.967 after stimulation.
Nine participants decrease, five are unchanged and one increases. The apparent
significance of the plotted interval is nominal and specific to the t test.
The interval is not adjusted for the four-outcome family.

The SUDMEX randomized comparison includes 22 active and 18 sham participants.
The normalized active-minus-sham contrast adjusted for baseline and motion is
+0.400 gateways, 95% CI [−0.599, 1.399], p=.422, Holm4 p=1. This does not
reproduce the Aachen decrease. Counts represent classified parcels in a fixed
reference workspace, not loss of tissue or proof of reduced attention.

Audit code: `setup/audit_dan_gateway_result.py`. Aggregate results and an
explicitly labelled figure are in `PROVENANCE/dan_gateway_audit_20261126/`.

## Current publication scope

The corrected Aachen manuscript is a working draft. The two-dataset
haemodynamic draft is retained with a status notice and is not the selected
submission. The additional information-decomposition and spatial analyses below
are now complete. They retain matched participant contrasts and SUDMEX's
randomized comparison.
Later open-label SUDMEX measurements remain descriptive sensitivities.

The repeatability pilot is a separate exploratory analysis, not evidence of
an absence of stimulation effects and not the intended treatment paper.

## Deeper PID follow-up: completed

We tested six disjoint processing modes: storage, transfer, copying, erasure,
joint-to-individual conversion and individual-to-joint conversion. Together,
their raw values sum to time-delayed mutual information. The taxonomy follows
[Mediano et al. (2025)](https://doi.org/10.1073/pnas.2423297122). These are
descriptions of predictive information, not proof of neuronal causal flow.

Every mode was evaluated for seven fixed sets of region pairs: whole brain,
inside the reference workspace, outside it, across its boundary, FPN–striatum,
FPN–DAN and FPN–DMN. The boundary analysis addresses the requested workspace
question; the other scopes remain visible. Six modes form a family within a
scope, and the complete 42-outcome correction is also reported. These are
exploratory follow-ups; writing their specifications before extraction does
not make the overall project preregistered.

Aachen uses the mean of two matched visit changes for each of 15 participants,
with all 32,768 participant sign patterns. SUDMEX retains 22 active and 18 sham
participants in the two-week follow-up model, adjusting for baseline neural
values and motion change with HC3 uncertainty. Baseline-only models are shown
as sensitivities. We did not pool these different experimental contrasts.

| Follow-up | Result | Implication |
| --- | --- | --- |
| Raw modes, 42 endpoints per cohort | Smallest SUDMEX baseline+motion p=.3468; all its Holm6 p=1 | No shared effect is hidden merely by multiple-comparison correction |
| Workspace boundary in Aachen | Individual-to-joint increase +.007555 nats, exact p=.0243, maxT6 p=.0729; reverse conversion exact p=.0255, maxT6 p=.0779 | Nominal Aachen signals, without corrected support or SUDMEX replication |
| Separately normalized modes | Smallest SUDMEX baseline+motion p=.3504; all cross-dataset Holm6 and Holm42 p=1 | Conditional normalization does not establish a replicated treatment effect |
| Regional change patterns | Each cohort tested against a template learned only in the other; all six primary conjunction p=1 | No shared regional redistribution established |

The regional test uses participant-level projections onto independent-cohort
templates, not parcels as independent samples. All 454 regions and a cortical
400-region sensitivity are reported, together with every participant omission.
The descriptive map correlations are close to zero (approximately −.006 to
−.036 for the full atlas); they do not support a claim of opposite responses.

The normalized analysis forms each mode in nats **before** calculating its own
conditional-null percentile, then averages the pairwise percentiles. Percentile
atoms were never added. Each cohort's original baseline-only VAR pool was
reconstructed from the current time series and verified, including the corrected
Aachen routing. Two seeds each generated 32,768 draws at each .005-nat grid step.
The largest held-out interpolation error was .000701 on the 0–1 scale. Population
and finite-record simulations used independent Lyapunov solves and atom
reconstruction; the maximum mean deviations from .5 were .051 and .079,
respectively. These diagnostics do not make the conditional percentiles
distribution-free, and normalized modes no longer form an additive decomposition.
All seed-specific estimates remain available; they are numerical sensitivities.

The raw extraction reproduced every pair's five saved MATLAB composites within
1.85×10⁻¹³ and independently reconstructed sampled complete atom decompositions
within 2.04×10⁻¹⁴. Independent implementations of the HC3 models agreed. All input
hash checks passed. Numerical correctness does not imply a positive biological
result.

Code and local runs:

- `setup/explore_pid_processing_modes.py`:
  `output/pid_processing_modes/20261126T130022_906563Z`.
- `setup/explore_pid_change_patterns.py`:
  `output/pid_change_patterns/20261126T130710_947110Z`.
- `setup/normalize_pid_processing_modes.py`:
  `output/pid_modes_normalized/20261126T131608_596458Z`.

The dashboard includes all aggregate tests and normalization diagnostics. Their
scope-specific families do not erase earlier analyses or create project-wide
family-wise error control across every exploratory round. No positive shared
treatment narrative is established by the current results. A new dataset can
test a frozen question independently; it cannot retroactively supply replication
within Aachen and SUDMEX.

## Additional data

[The new dataset review](TMS_ADDITIONAL_DATASETS_20261126.md) recommends
OpenNeuro ds001927 as the most practical additional resting-state fMRI cohort.
Public metadata and sample imaging headers were downloaded and verified. The
full version 2.1.0 release is now downloading to `data/ds001927/` under the
user's instruction, with pinned content checksums. It has not been analysed.
