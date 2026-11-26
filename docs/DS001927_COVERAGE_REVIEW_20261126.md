# ds001927: tracing inferior signal loss

The two completed pilots do **not** show the repeat-extraction problem found in
Aachen. Weak inferior frontal and temporal BOLD signal is visible in the raw
acquisitions, before our preprocessing. The evidence is most consistent with
acquisition-related signal dropout. It does not establish the precise physical
cause of every affected voxel or validate every future participant.

## What was checked

Both pilots (s001 and s022) completed fMRIPrep successfully. Native brain
extraction, segmentation, standard-space anatomy and all twelve functional
reference overlays were inspected. Unlike Aachen's already skull-stripped T1s,
these inputs retain skull/scalp; normal extraction is appropriate. The damaged
s022 control-session T1 was excluded from reference construction, using its
other two anatomies. The resulting common anatomy has no obvious gross clipping
in the reviewed sections.

The audit read all 600 frames of each raw and processed pilot BOLD file, rather
than relying on the initial-frame input check. All values were finite. It then
applied the **existing** spatial transformations to raw temporal means and to
an all-ones image marking the acquired field of view. No registration was
re-estimated and no production image was changed. Applying the same transforms
to the native reference reproduced all twelve saved MNI references exactly
(maximum absolute voxel difference zero). This verifies the transform chain's
application, not anatomical registration ground truth.

Across 12 runs × 454 parcels:

| Check | Finding |
| --- | --- |
| Parcel/run cases with <50% functional-mask overlap | 220 of 5,448 (4.0%) |
| Those cases with >95% geometric acquisition coverage | All 220 |
| Median raw signal in the low-mask-coverage cases | 19.0% of the run's median signal inside the mask and atlas |
| Median MNI z coordinate of those cases | −21.3 mm |
| Acquired field of view covering the atlas, per run | 99.888%–100% |
| Parcels with zero temporal variance throughout | None; some individual voxels have zero variance |

Thus, “inferior brain coverage” was too imprecise: most of this is **weak signal
within the acquired field of view**, accompanied by restrictive functional
masks, rather than brain tissue absent because of an incorrect anatomical
extraction. Changing the extraction to `skip` would not restore signal that
is already weak in the original BOLD scans. A few geometric edge voxels also
lack coverage, but they do not explain the 220 low-mask-coverage cases.

The raw temporal mean is not motion corrected; its transformed image is a
diagnostic display, not a voxelwise before/after processing test. Small residual
alignment and distortion effects remain possible. No field maps or usable
phase-encoding/readout metadata were deposited. Finite, nonzero background
noise must not be counted as valid regional neural signal.

## Consequences for analysis

Keep the existing anatomical preprocessing configuration. Before PID inference,
freeze a common coverage-eligible parcel set across the required conditions and
participants, then validate the denoised time series. The new
[site-comparison plan](DS001927_SITE_COMPARISON_PLAN_20261126.md) specifies a
50% mask-overlap intersection and an 80% sensitivity. For s001 alone, its four
FRO/CTR runs retain 413/454 parcels at 50% and 361 at 80%; these are **pilot
counts**, not the eventual cohort mask. The same eligible pair set must be used
for all conditions to prevent a change in included regions from creating a
pre/post difference. Workspace labels remain frozen, intersected with that set.

All six s001 scans pass the primary and stricter motion gates. Five of six s022
scans pass the primary gate; CTR pre exceeds the 3-mm maximum-axis translation
limit (3.1135 mm). Therefore s022 does not enter the main frontal-versus-control
comparison under the existing rule. No threshold was relaxed to retain it.

## Processing status

The reviewed configuration is released for **fMRIPrep processing** of the
remaining 21 participants, one at a time. It does not release their anatomy,
denoised signals or PID outcomes automatically. The sequential queue stops on
a failed process or incomplete derivative inventory and generates motion tables
and anatomical review figures for each completed person.

D: had only about 47 GiB free after the pilots. All 651 pilot derivative files
(19.55 GB) were copied to
`E:/data/TMS_GWT/ds001927_20261126/derivatives/fmriprep_23_2_1`, with every checksum
verified. The D: originals remain intact. New work and derivatives use the
dedicated E: tree; no raw input, existing cohort or frozen HCP result was altered.
The queue started s002 on 26 November at 17:10 UTC. No ds001927 PID results exist
at this release.

Reproduction: `setup/audit_ds001927_coverage.py` (prepare, warp, summarize),
`setup/stage_ds001927_storage.py`, and `setup/run_ds001927_cohort.py`.
Local evidence is in `output/ds001927_coverage_trace_20261126/`.
Public decision/checksum receipts are in
`PROVENANCE/ds001927_preprocessing_20261126/`, particularly
`pilot_coverage_release.json`, `analysis_plan_v1.json`, `storage.json`, and
`cohort_queue.json`. The queue receipt is a point-in-time snapshot in Git;
its local copy records the current running state.
