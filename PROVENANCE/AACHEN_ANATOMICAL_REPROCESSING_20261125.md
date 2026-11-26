# Aachen anatomical correction

The target-localization audit exposed a problem upstream of the target masks. The masks are within the native functional brain. The deposited T1-weighted anatomy is already skull-stripped, whereas the original fMRIPrep 23.2.1 run used `skull_strip_t1w = "force"`. Its resulting anatomical brain masks exclude substantial outer cortex and cerebellum. Independent audit files are in `output/aachen_alignment_review/`; parent cross-checks are in `output/aachen_target_alignment_audit/`.

This is not evidence that the original stimulation was outside the brain. Numerical reproduction of saved standardized BOLD images verified transform application but did not independently validate anatomical registration. Initial target-localized Aachen statistics in `output/targeted_information_flow/20261125T193623_492441Z` are suspended pending correction. Existing atlas-parcellated Aachen outputs require review and regeneration if the corrected anatomical preprocessing changes registrations or nuisance masks. The whole-brain Φ result must not currently be presented as a validated final result.

## Authorized correction and isolation

The user requested robust targeted analysis, challenged the invalid anatomical inference, and requested Aachen work in parallel. Cloud mode remains authorized. No HCP or psychedelic processing is performed.

- Only the permitted per-dataset setting in `pipeline/configs/datasets/aachen_itbs.yaml` changed: `fmriprep_extra_args: "--skull-strip-t1w skip"`.
- Production processing code is unchanged. The original derivatives and downstream results are retained.
- `setup/run_aachen_skullstrip_pilot.ps1` launches the same fMRIPrep 23.2.1 image, with no surface reconstruction or distortion correction, in a fresh work/output branch. Participant 06 is the initial pilot because it showed a clear anatomical mismatch, selected without reference to information-flow effects.
- New output: `data/aachen_itbs/derivatives/fmriprep_skullstrip_skip_20261125/`.
- Fresh work directory: `tmp/fmriprep_aachen_skullstrip_skip_20261125/sub-06/`.
- Pilot receipt: `PROVENANCE/aachen_skullstrip_pilot_06.json`.
- Before expanding or promoting outputs, verify native anatomical mask preservation, anatomical/functional alignment in native and standard space, cortical-atlas coverage and deposited target overlays. Do not move targets to force a desired result.
- After validation, all affected Aachen stages need rerunning using the unchanged production scripts: confound QC and preparation, CONN, parcellation, rsHRF and ΦID, followed by the same participant-level statistics. Corrections cannot be limited to moving the final target mask because tissue-derived nuisance regressors can also change.

The fMRIPrep 23.2.1 documentation explicitly supports skipping skull stripping for already masked anatomy: https://fmriprep.org/en/23.2.1/usage.html#specific-options-for-ants-registrations.

## SUDMEX independence

The completed independent SUDMEX review covers all 44 anatomical scans and all 142 functional runs. All raw T1s retain head tissue and the reviewed masks preserve peripheral cortex; no Aachen-like repeat-extraction failure was found. BOLD/T1 mask Dice is 0.883–0.956. A blanket T1 rerun is not supported. Inferior functional coverage and focal mismatch remain, however: some orbitofrontal/temporal-pole parcels have poor coverage despite passing motion QC. See `output/sudmex_registration_audit/SUDMEX_FULL_PREPROCESSING_QC.md` and its per-run/parcel evidence. The new SUDMEX information-flow analysis has independent numerical and original-MATLAB verification; no main ten-test comparison survives family-wise correction.

The independently derived SUDMEX baseline and sham workspaces resemble the frozen HCP reference, as documented in `docs/TMS_CONTROL_WORKSPACE_QC_20261125.md`. That descriptive pattern check does not supersede anatomical QC or validate the old Aachen preprocessing.

## Completed correction and control QC

The corrected 15-participant, 60-run cohort completed the full chain at 23:41 UTC on 25 November (26 November locally). All 15 participants passed independent anatomical review; native support retention is 99.63–100%, versus 44.2–52.9% originally. The release is `output/aachen_alignment_review/corrected_cohort/cohort_release_receipt.json`, SHA-256 `d02401c42be2ea20c9e6de008a30c9665d92e7151467bb0aa4c873a8dfd81895`. Its integrity review verified 762 linked files. Regional functional-coverage limitations and the resolved reference-enhancement investigations for subjects 10 and 14 remain in the cohort report.

The isolated downstream branch is `reruns/aachen_skullstrip_skip_20261125/pipeline/`. All 60 parcellated and deconvolved signals passed checks. All 300 PhiID matrices have the expected 102,831 finite unique pairs per scan, with checked symmetry, metadata, measure identities and hashes. No HRF or PhiID scan was skipped or failed. Numerical processing scripts were unchanged; permitted dataset paths/settings and logged runtime worker counts differ. See `logs/aachen_corrected_20261125/downstream/status.json` and its verification files.

The independently derived corrected Aachen baseline-mean workspace shares 176 regions with the frozen 235-region HCP/GWT mask (Dice 0.743, gradient rho 0.596). Separate pre90 and pre120 Dice values are 0.768 and 0.681; between-baseline Dice is 0.782. No post-stimulation scan enters this QC. Full results and participant-resampling diagnostics are in `output/control_workspace_qc/aachen_corrected_20261125T234156_788094Z/` and `docs/TMS_CONTROL_WORKSPACE_QC_20261125.md`.

Original derivatives, signals and effect tables are retained. The original Aachen treatment statistics, target-flow estimates and manuscript claims remain on hold until recomputed from the corrected branch. Old Aachen NuMIT references were learned from the original baseline inputs and cannot be treated as corrected normalization references. Control-workspace resemblance does not reinstate those historical treatment claims.
