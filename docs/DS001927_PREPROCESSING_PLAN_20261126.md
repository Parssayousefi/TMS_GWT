# ds001927: input review and preprocessing decisions

26 November 2026. This plan is for the downloaded v2.1.0 release. It records
choices before any new PID treatment results. It is not a registered protocol.

## What the deposited images contain

All 23 participants have three site sessions (FRO, OCC, CTR), each with a T1w
and pre/post resting runs: 69 anatomical images and 138 BOLD runs. BOLD has
600 volumes at TR 1.25 s (750 seconds, 12.5 minutes), with approximately
3 × 3 × 3.3 mm sampling. The paper's nominal slice thickness is 3 mm; header
spacing includes the gap. The three T1s per person are distinct files.

The [acquisition paper](https://doi.org/10.1126/sciadv.aay2739) describes an
individually localized dorsoanterior prefrontal target. FRO should not silently
be equated with the previous cohorts' DLPFC coordinates or network membership.

The [source README](https://raw.githubusercontent.com/OpenNeuroDatasets/ds001927/2.1.0/README)
documents pydeface. Visual inspection of central sagittal sections of all 69
T1s shows retained scalp/skull/neck: **defacing is not skull stripping**.
Normal brain extraction is needed, unlike the corrected Aachen branch.
No raw images have been modified.

One deposited anatomical scan, **sub-s022/ses-ctr**, has a flat superior cut
through brain tissue, confirmed in sagittal and coronal views. The upper
52 voxel planes are exactly zero in its RAS-oriented volume. Its other two
T1s appear complete. Exclude this anatomical scan from reference construction;
retain the participant's six BOLD runs for independent motion/alignment QC.
This is an image-quality decision made before treatment statistics, not a
participant exclusion. Do not infer whether acquisition or anonymization
caused the missing tissue.

All full T1 arrays and the first 24 frames of each BOLD contain finite values.
These checks and raw contact sheets do not establish acceptable motion,
registration or regional coverage. BIDS Validator 1.14.1 reports no errors;
its warning concerns missing slice timing. All 138 functional sidecars lack
SliceTiming, PhaseEncodingDirection, TotalReadoutTime and EffectiveEchoSpacing.
There are no deposited field maps.

Local audit: `output/ds001927_input_audit_20261126/`; reproduction:
`setup/audit_ds001927_inputs.py`. Download content verification is separately
recorded in `PROVENANCE/additional_datasets_20261126/download_complete.json`.

## First-stage processing and pilot

Use installed fMRIPrep **23.2.1**, matching the established cohorts. This is a
deliberate version-control choice, not a claim that it is the newest release.
The later [release notes](https://fmriprep.org/en/25.2.7/changes.html) were checked;
new versions should be evaluated as a separate version sensitivity if needed,
not silently substituted during a cohort. Retain container identity and commands.

1. Start with sub-s001, all three T1s and six BOLD runs. Review the completed
   pilot before cohort dispatch. Use sub-s022 as a second anatomy challenge
   case with its damaged CTR T1 explicitly filtered out.
2. Use `--skull-strip-t1w force` because the inputs retain the skull, plus
   bias-field correction and tissue segmentation. Construct one common,
   unbiased within-person anatomical reference (`--longitudinal`) from all
   usable T1s. All six runs register to that same reference. No choice depends
   on observed pre/post effects.
3. Perform motion correction and BOLD-to-T1 boundary-based registration with
   fMRIPrep's normal fallback checks, then T1-to-MNI normalization and combined
   BOLD resampling to MNI152NLin2009cAsym at 2 mm. No surface reconstruction
   or added spatial smoothing is needed for the volumetric parcel analysis.
4. Explicitly omit slice-timing correction: the acquisition times are unknown.
   Multiband factor 2 does not recover them. Omit susceptibility correction:
   no field maps or reliable PE direction/readout metadata are available.
   DICOM `COL` is not sufficient signed BIDS PE information for SyN correction.
   Do not invent metadata. Inspect orbitofrontal/temporal distortion and signal
   coverage and acknowledge this limitation. This is offline TMS, so no
   simultaneous pulse-interference removal is required.

The [fMRIPrep workflow documentation](https://fmriprep.org/en/23.2.1/workflows.html)
describes the extraction, registration and slice-timing behavior; the
[usage reference](https://fmriprep.org/en/23.2.1/usage.html) defines the options.
The pilot launcher is `setup/run_ds001927_pilot.ps1`, with isolated work and
derivative folders. It never replaces existing containers or treats an HTML
report alone as success. Raw data are mounted read-only.

## Release gates before denoising

Require all six expected images, masks and confound tables, 600 frames and
TR 1.25; verify site/task identity explicitly. Inspect extraction contours,
GM/WM/CSF segmentation, each BOLD-to-T1 registration and T1-to-MNI alignment.
Measure atlas and target coverage, including the inferior field of view.
Missing parcels are not filled with zero or treated as neural effects.
Pilot registration approval is distinct from scan-level motion inclusion.

Keep the established 3 mm translation / 3.15 degree rotation screen relative
to the first retained frame. The primary new-dataset gate additionally requires
mean FD ≤0.30 mm, no more than 20% of retained frames with FD >0.50 mm, and at
least 480 seconds of uninterrupted data. Report a stricter mean-FD ≤0.20 mm
sensitivity, without selecting the threshold on the resulting effects. DVARS
>1.5 is a diagnostic flag. These are study choices, not universal cutoffs.
They are frozen in `PROVENANCE/ds001927_preprocessing_20261126/qc_policy.json`
before pilot confound inspection. Add an explicit FD/DVARS audit;
the existing QC script does **not** enforce the YAML FD thresholds. In
particular, its salvage function describes an FD rule which is not actually
implemented. A ds001927 adapter must implement, verify and document any new
motion rule before results. Do not claim the old cohorts met an unimplemented
criterion. No final new-dataset motion exclusions have yet been assigned.
The [motion benchmark](https://doi.org/10.1016/j.neuroimage.2017.03.020)
shows that nuisance regression choices have different residual artifacts;
matching an existing method is not a substitute for QC.

For time-lag analysis, do not delete scattered frames and join the remaining
time points as if they were adjacent. Any censoring sensitivity must exclude
lagged pairs touching a bad frame, or use a declared contiguous clean window,
and report retained duration. Censoring/denoising behavior must also be
compatible with HRF estimation. Downstream release remains gated on this
decision and actual pilot confounds.

## Downstream signal and statistical plan

Retain the verified CONN definition for the comparable branch: 5 WM and 5 CSF
components, motion plus derivatives, linear detrending, 0.008–0.09 Hz filtering,
no GSR or added spatial smoothing. Remove at least eight initial frames
(10 seconds), extending to a longer leading nonsteady-state prefix when
flagged by fMRIPrep. Non-prefix flags require review. Five frames here would
mean only 6.25 s. Use true TR 1.25 in rsHRF. Keep the 454-region atlas and frozen HCP mask.
Document extra QC/sensitivity code separately; no frozen cohort is rerun.

An isolated dataset registry prevents accidental processing of other cohorts.
Pre/post is encoded in **task** and stimulation site in **session**. Preserve
both: join participant × site × time explicitly. The generic group loader
cannot safely represent this two-factor design; pooling all `restpre` scans
would count each person three times. Subject-level inference requires a
dataset-specific contrast adapter.

The proposed primary contrast is within-person
`(FRO_post − FRO_pre) − (CTR_post − CTR_pre)`, for mean pairwise ΦR and explicitly
defined six-mode storage, with two-sided family-wise correction across those
two endpoints. OCC and other modes need declared secondary families. CTR is
an active stimulation site, not sham. These definitions should be frozen
before inference; they do not guarantee a positive result.

One native sample means 1.25 s here and 2 s in Aachen/SUDMEX. Keep lag units
explicit and cohort estimates separate; never pool raw magnitudes as if the
lag and acquisition were identical. Specify lag sensitivities before analysis.
Unknown slice timing and a BOLD-derived estimator limit claims about directed
neuronal causality. The experimental contrast can support an effect of
stimulation site; a PhiID transfer atom alone cannot establish a causal pathway.

## Resources and current status

At input review, D: had approximately 95 GiB free. Do not launch all 23 people
at once. Measure pilot derivative/work size, then use bounded sequential
processing and verified cleanup of newly generated scratch only, or a suitable
documented storage location. Existing cohorts and source imaging stay intact.
Pilot receipts in `PROVENANCE/ds001927_preprocessing_20261126/` record actual
launch state; this document does not itself certify a completed pilot. The
sub-s001 six-run pilot finished successfully at 15:10 UTC on 26 November. All six
runs passed the prespecified motion gate and stricter mean-FD sensitivity,
retaining 592 volumes (740 s) each. Native extraction, tissue segmentation,
standard-space anatomy and all six BOLD reference overlays were visually
reviewed without gross clipping or misregistration. Inferior frontal/temporal
dropout remains: 121 parcel/run combinations have less than 50% functional-mask
coverage, including two with zero mask overlap. Nonzero reference support is
present in all parcels, which alone does not certify useful temporal signal.
The completion review is `pilot_s001_completion_review.json` in the provenance
directory. The sub-s022 pilot is now running with its damaged control-session
T1 excluded from reference construction. Full-cohort and outcome-analysis
release remain pending. MATLAB initially returned a
licence-server error; after VPN reconnection it started successfully and its
MATLAB, Signal Processing and Wavelet availability checks returned true.

The downstream copy is prepared in `reruns/ds001927_20261126/pipeline/`, with
104 processing scripts unchanged and three copies changed only for run paths,
dataset selection or TR lookup. The registry contains ds001927 only. No
denoising, deconvolution or PID analysis has been launched for this dataset.

Group-level outcome models will not include an additional continuous FD
covariate. Motion remains handled by the explicit scan-QC gate and the existing
denoising regressors.
