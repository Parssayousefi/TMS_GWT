# Target-based TMS information-flow analysis

The sham-controlled SUDMEX comparison does not support a stimulation-related change in net outward information flow, persistent synergy or persistent redundancy at the DLPFC target neighbourhood. All ten primary tests have family-wise p ≥ 0.406. This is insufficient evidence for a change, not proof of no effect. All original Aachen results are **suspended pending statistical reanalysis of the corrected inputs**; their present numbers must not be used as evidence of absence, of incorrect stimulation placement, or of a validated global effect.

## What was tested

SUDMEX contributed 40 complete acute-phase MRI pairs: 22 active and 18 sham, baseline versus the two-week visit only. The fixed left target is a 10-mm neighbourhood around the mean of all 21 published baseline MNI localizations: (−33.14, 43.71, 37.14) mm. It was defined before extracting targeted information outcomes. The mirrored right region is a spatial comparator, not an untreated control. Individual baseline coordinates were available for 19 usable participants (12 active, 7 sham); these form a separate localization sensitivity. The source specifies MNI but not its precise template variant, so localization is approximate. [Data descriptor](https://www.nature.com/articles/s41597-024-03242-y); [public localization workbook](https://zenodo.org/records/10409461).

The unchanged deconvolved 454-region time series were analysed at one TR (2 s). For target parcel X and remote parcel Y, outward transfer entropy is I(Xpast;Yfuture | Ypast), and inward transfer entropy reverses X and Y. Their difference describes net directional predictability. In the ΦID lattice:

- TEout = xtr + xty + str + sty; TEin = ytr + ytx + str + stx.
- Net unique-source contribution = (xtr + xty) − (ytr + ytx).
- Net synergistic contribution = sty − stx.
- Net TE equals those two net contributions summed, at the scan level.
- Persistent synergy is sts and persistent redundancy is rtr. These last two are symmetric pairwise measures; they do not encode a direction of flow.

The synergistic TE contribution requires the joint past of both regions; it cannot be interpreted as information originating solely in the stimulated region. This follows the [ΦID framework](https://doi.org/10.1073/pnas.2423297122). Seed parcels were weighted by anatomical overlap; all remote destinations received equal weight, excluding both seed sets. Both sides used the same destination pool. Effects are participant-level averages, not independent tests of thousands of edges. We did not search destination parcels or alter regions according to their effects.

## Primary randomized comparison

Post-treatment outcomes were regressed on baseline and active/sham assignment, with HC3 uncertainty. Ten tests combine five measures with two questions: change at the left target, and change at the left target relative to the right comparator. Shared treatment-label permutations (32,768 draws, preserving 22/18 counts) give single-step max-|t| family-wise p values. Label exchangeability among imaging completers is assumed; this is not a reconstruction of undisclosed trial randomization blocks. Marginal confidence intervals are not multiplicity-adjusted.

| Scope | Measure | Active−sham (nats) | 95% CI | Uncorrected permutation p | FWER p (10 tests) |
| --- | --- | --- | --- | --- | --- |
| Left target | Net outward transfer entropy | -0.000059 | [-0.001404, +0.001286] | 0.9263 | 1.0 |
| Left target | Net unique-source contribution | +0.000006 | [-0.000011, +0.000024] | 0.5148 | 0.993 |
| Left target | Net synergistic contribution | -0.000050 | [-0.001403, +0.001303] | 0.9388 | 1.0 |
| Left target | Persistent synergy | +0.013357 | [-0.016921, +0.043635] | 0.3589 | 0.9714 |
| Left target | Persistent redundancy | +0.001324 | [-0.001881, +0.004528] | 0.3957 | 0.9814 |
| Left minus right | Net outward transfer entropy | -0.000319 | [-0.002030, +0.001391] | 0.7039 | 1.0 |
| Left minus right | Net unique-source contribution | +0.000022 | [-0.000010, +0.000055] | 0.1591 | 0.6973 |
| Left minus right | Net synergistic contribution | -0.000358 | [-0.002062, +0.001347] | 0.6673 | 0.9999 |
| Left minus right | Persistent synergy | +0.012002 | [-0.012861, +0.036865] | 0.3197 | 0.9496 |
| Left minus right | Persistent redundancy | -0.002553 | [-0.005457, +0.000351] | 0.073 | 0.4057 |

Even before correction, none of the ten primary tests has p < .05. The smallest p concerns lower persistent redundancy in the left-minus-right comparison: p=0.073 uncorrected and 0.406 corrected. It is not a supported positive finding.

## Robustness checks

No test survives ten-test correction after adding change in mean framewise displacement, changing the lag to 4 s, matching pre/post record lengths within participant, or using individual baseline locations. Equal-length checks reuse the already deconvolved signals; they do not refit HRFs. Motion adjustment is a diagnostic sensitivity, because post-treatment motion could itself be influenced by treatment.

| definition | variant | model | n | minimum_maxT10_p |
| --- | --- | --- | --- | --- |
| common_baseline_target | lag1_equal_length | baseline | 40 | 0.4375 |
| common_baseline_target | lag1_full | baseline | 40 | 0.4057 |
| common_baseline_target | lag1_full | baseline_motion | 40 | 0.3433 |
| common_baseline_target | lag2_full | baseline | 40 | 0.6106 |
| individual_baseline_target | lag1_equal_length | baseline | 19 | 0.2587 |
| individual_baseline_target | lag1_full | baseline | 19 | 0.2588 |
| individual_baseline_target | lag1_full | baseline_motion | 19 | 0.3058 |
| individual_baseline_target | lag2_full | baseline | 19 | 0.4698 |

Leave-one-participant-out coefficients are saved separately as influence diagnostics; no omission was selected as the main analysis. All five measures, both scopes, both localization definitions and the three signal/lag variants remain available in the inference table. These bounded correction families do not correct the entire historical project search.

## Aachen alignment status

Public masks were found for eight usable participants. The independent audit established that the deposited Aachen T1s were already skull-stripped but underwent forced repeat extraction. The resulting masks removed substantial outer cortex and compromised anatomy-dependent processing. The deposited target masks align with the full native anatomy; their appearance outside the clipped standardized T1 did not establish incorrect stimulation placement. Numerical reproduction of a saved image using its own transform was insufficient to validate registration. All 15 participants and 60 runs have now passed independent anatomical review. The isolated corrected branch has completed CONN, parcellation, deconvolution and ΦID with output validation. Corrected native support retention is 99.63–100%; regional functional-coverage limitations remain documented. All old Aachen estimates, including whole-brain Φ, remain provisional until recomputed. See `PROVENANCE/AACHEN_ANATOMICAL_REPROCESSING_20261125.md` and `output/aachen_alignment_review/corrected_cohort/cohort_release_receipt.json`. The deposited seed is an intended analysis ROI; exact coil or peak-field location is not independently verified.

The separate SUDMEX audit covered 44 anatomical images and 142 functional runs. It found no comparable T1 over-extraction, so a blanket anatomical rerun is not indicated. Inferior functional coverage and some local alignment limitations remain; these were not resolved by the numerical checks on information estimates. See `output/sudmex_registration_audit/SUDMEX_FULL_PREPROCESSING_QC.md`.

## What can be called causal

Randomized active-versus-sham assignment provides a basis for estimating a treatment effect on these measured information summaries, subject to complete-case selection and other trial limitations. This analysis found no convincing such effect. Pairwise lagged fMRI predictability does not identify direct neuronal causal connections: unmeasured common inputs, indirect paths, haemodynamic modelling and slow sampling remain relevant. A net-flow average can also hide offsetting changes across destinations or simultaneous changes in both directions. The results therefore neither establish nor rule out every possible circuit effect.

## Verification and reproducibility

Original MATLAB PhiIDFull agreed for 36 fixed validation pairs across both datasets and both lags, all 16 atoms per pair: maximum error 1e-15 nats. Reconstructed five original matrix measures agreed within 2.9e-14 nats. Transfer entropy independently matched Gaussian conditional mutual information within 1.9e-15 nats. Synthetic unidirectional data and source swapping verified the direction convention. Independent full-design HC3/permutation calculations reproduced inference, and 676 input hashes were checked.

Run: `output/targeted_information_flow/20261125T193623_492441Z`. Reproduction scripts: `setup/analyse_targeted_information_flow.py`, `setup/verify_targeted_information_flow.py`, and `setup/verify_target_flow_matlab.m`. All calculations read frozen TMS signals only; HCP/psychedelic cohorts and production processing code were not rerun or changed. New directed measures were not normalized using the old symmetric-measure NuMIT references.

## Control workspace quality check

The workspace was derived independently from the control scans with unchanged GWT core functions. The HCP mask and gradient are byte-identical to the original GWT project reference. The matched sham sample has 18 participants at both visits; later open-label treatment is excluded. This is descriptive pattern QC, not proof of valid registration at every parcel or a treatment-effect result. Corrected Aachen baseline-workspace results are reported separately below.

| Condition | N | Workspace regions | Shared with HCP | Dice | Gradient rho |
| --- | --- | --- | --- | --- | --- |
| all_baseline | 41 | 227 | 179 | 0.775 | 0.676 |
| sham_t0 | 18 | 234 | 182 | 0.776 | 0.652 |
| sham_t1 | 18 | 232 | 173 | 0.741 | 0.615 |

## Corrected Aachen control workspace

Aachen has no sham arm. These workspaces use the two prestimulation visits, separately and averaged within each of the same 15 participants. The HCP mask remains fixed for treatment inference. Bootstraps, split halves and uniform leave-one-participant-out outputs describe sampling stability; no omission selects the primary map.

| Condition | N | Workspace regions | Shared with HCP | Dice | Gradient rho |
| --- | --- | --- | --- | --- | --- |
| baseline_mean | 15 | 239 | 176 | 0.743 | 0.596 |
| pre90 | 15 | 234 | 180 | 0.768 | 0.607 |
| pre120 | 15 | 229 | 158 | 0.681 | 0.512 |
