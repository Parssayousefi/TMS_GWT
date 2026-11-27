# ds001927 preprocessing restart, 27 November 2026

The user explicitly resumed the five unprocessed participants after a pause
to free the workstation. The sequential queue restarted at 14:03 Zurich time
with s018; s019, s020, s021 and s023 follow automatically. One participant is
processed at a time, with the same pinned fMRIPrep image and acquisition/QC
settings. The local `cohort_queue.json` is the live status; Git stores a
point-in-time copy.

Before this restart, 18/23 participants had completed fMRIPrep: the two pilots
and s002–s017. All 108 expected scans have complete derivatives. Of those,
106 pass the primary scan-motion gate, and 16 participants have all four
passing FRO/CTR scans required for the main contrast. This is motion
eligibility, not a completed anatomical or downstream signal release. No
ds001927 denoising, deconvolution or PID outcome analysis has started.

## Why the queue stopped

At 05:15 Zurich time, s017 finished fMRIPrep successfully, but the queue
misread the QC command's exit status. The QC command returns 0 when every scan
passes motion, 3 for complete files with motion exclusions/review, and 2 for
incomplete/invalid files. The original queue accepted 0/1, so it incorrectly
reported incomplete derivatives when the valid exclusion returned 3.

s017's frontal baseline has maximum-axis rotation 3.4996 degrees, exceeding
the unchanged 3.15-degree limit. That scan and the participant's main FRO/CTR
contrast remain excluded. The earlier s022 control-baseline exclusion also
remains. The bug was in orchestration, not an altered motion threshold or
image-processing calculation.

## Recovery and checks

The completed s017 container identity and successful exit were verified. A
fresh QC command in a separate audit folder returned the expected 3; all six
run identities, confound hashes, motion values and inclusion decisions matched
the original receipt exactly. Missing review figures were generated, without
rerunning fMRIPrep or automatically approving anatomy.

The runner now validates the exit code against the QC receipt, accepts a
complete motion-exclusion result, and continues to stop on missing files,
crashes or inconsistent receipts. It saves successful fMRIPrep completion
before subsequent QC checks. Explicit resumption verifies the original subject
selection, completed prefix, input hashes and completion receipts; it refuses
active/unreconciled queues and uses an exclusive process lock. Eight regression
tests passed, plus the real s017 QC-command check. No production numerical
script, existing-cohort result or frozen HCP input was changed.

The pause/failure snapshots are retained in
`PROVENANCE/ds001927_preprocessing_20261126/`. The recovery receipt is
`qc_exit_recovery_20261127.json`; `cohort_queue.json` records the resumed runner
hash and remaining participants. Reproduction helpers:
`setup/run_ds001927_cohort.py`, `setup/reconcile_ds001927_qc_stop.py`, and
`setup/test_run_ds001927_cohort.py`.
