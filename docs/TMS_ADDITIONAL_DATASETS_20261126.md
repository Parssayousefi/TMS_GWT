# Additional data: access and scientific fit

Checked 26 November 2026. **OpenNeuro ds001927 is the most practical additional
resting-state fMRI dataset for this project.** It adds an informative within-person
stimulation-site comparison. A positive PID result is not guaranteed, and it
would not turn the present Aachen/SUDMEX findings into a replicated effect.

## Recommended: ds001927, Castrillon et al.

The [dataset README](https://raw.githubusercontent.com/OpenNeuroDatasets/ds001927/2.1.0/README)
describes three counterbalanced sessions in 23 healthy participants: prefrontal
(FRO), occipital (OCC) and temporoparietal control (CTR) stimulation. Each session
uses the same 1-Hz protocol and resting fMRI before and after stimulation. The
[original paper](https://doi.org/10.1126/sciadv.aay2739) motivates testing whether
information-processing responses depend on the stimulated cortical region.

The public version 2.1.0 was inspected at commit
`3239a60164f7a230e02b6a42eb36a13ff4f97f57`. The file inventory contains:

- 23 participants with all three sessions and both rest scans: **138 resting
  runs**, plus 69 T1-weighted files.
- Explicit `ses-fro`, `ses-occ`, `ses-ctr` and `task-restpre`/`task-restpost`
  labels; there is no missing active/sham allocation key to obtain.
- Resting scans plus all structural files total **22,968,325,811 bytes**
  (22.97 GB / 21.39 GiB), from anonymous S3 object listings.
- An anonymous 64-KiB imaging read returned HTTP 206 and a valid NIfTI header:
  600 volumes, TR 1.25 seconds. All resting JSON sidecars also specify TR 1.25.
- [CC0 metadata](https://raw.githubusercontent.com/OpenNeuroDatasets/ds001927/2.1.0/dataset_description.json).

This is a **site-control design, not sham control**. CTR is real stimulation
at another site. It supports statements about site specificity more directly
than an uncontrolled pre/post contrast. Its 1-Hz protocol and prefrontal target
also differ from Aachen iTBS and SUDMEX 5-Hz rTMS; identical effect directions
should not be assumed.

A narrow prospective question would be: **Does prefrontal stimulation change
information integration differently from stimulation of a control site?**
For a small analysis, specify the FRO-minus-CTR difference in pre/post change
as the main contrast for information storage and redundancy-corrected integration
(ΦR), with Holm correction across the two endpoints. OCC and the other PID
modes should have explicit secondary families. This is a proposed design for
unseen data, not a description of an observed finding. The storage definition
must be fixed explicitly: six-mode PhiID storage is not interchangeable with
the sum of regional active information storage.

Before extraction, freeze signal routing, anatomical/coverage QC, exclusions,
physical lag choices, motion models and normalization. The new TR is 1.25 s;
one sample cannot silently be treated as the same 2-s lag used in the current
cohorts. Retain separate cohort estimates and avoid pooling raw magnitudes from
different acquisitions. A first-subject preprocessing check should precede
launching the full cohort. A new cohort tests generalizability or site specificity;
it cannot repair an unsupported shared claim in the two existing datasets.

## Other candidates

| Dataset | Verified status | Recommendation |
| --- | --- | --- |
| [ds008064: depression rTMS plus psychotherapy](https://github.com/OpenNeuroDatasets/ds008064/tree/1.0.0) | Still 28 listed participants, 23 matched rest pairs. Current commit is unchanged from the prior audit. Participants table contains ID, MADRS pre/post, age and sex; no stimulation allocation. Public imaging sample access succeeds. | Potentially useful sham-controlled replication **after obtaining the allocation and psychotherapy labels**. Do not analyse pooled active/sham scans as a TMS effect. |
| [Philadelphia Theta Burst Project, ds008893](https://github.com/OpenNeuroDatasets/ds008893) | Paper describes a 100-person dataset, but the current public tree and participant table contain **only one participant**. There are 24 resting echo files: four sessions × two runs × three echoes, not 24 independent scans or participants. No T1w files found. | Wait for the actual cohort release. Do not plan a group analysis from the paper's advertised sample size. |
| [PDC: ECT/ketamine](https://www.humanconnectome.org/study/crhd-perturbation-treatment-resistant-depression-connectome-fast-acting-therapies/document/pdc-release-10) | Controlled NDA/CCF access; public completeness audit previously found 33 ECT course-end pairs, 60 ketamine pairs and 17 healthy repeat pairs, before QC. Imaging and behavioral releases exist. | Worth pursuing as a separate clinical-treatment project if access is available. It does not provide TMS replication; treatment allocation is not randomized. |
| [Concurrent TMS-fMRI, ds005498](https://pmc.ncbi.nlm.nih.gov/articles/PMC12218990/) | Paper describes 152 participants and stimulation at 11 sites, with concurrent pulse acquisitions. Full current file inventory was not audited in this review. | Scientifically interesting for perturbation propagation, but needs a pulse-aware analysis. It is not a direct substitute for the current offline pre/post resting pipeline. |
| [Moffa prefrontal TBS EEG](https://pmc.ncbi.nlm.nih.gov/articles/PMC9681877/) | Previously verified open release: 24 people, repeated iTBS/cTBS and sham visits. | Useful if changing modality deliberately. The EEG artifact, reference, volume-conduction and lag choices require a separate validated workflow. Lower priority for completing the current fMRI project. |

The [Philadelphia paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13641012/)
and the [currently released participant table](https://raw.githubusercontent.com/OpenNeuroDatasets/ds008893/main/participants.tsv)
were checked separately. The current release is pinned at
`6126be09093f4ea218bf8735bf16032302775315`; its sample header confirms TR 1 s.
Release incompleteness, rather than an inaccessible download endpoint, is the
immediate obstacle.

## What was downloaded

The initial audit downloaded only public metadata, inventories and 64-KiB file
samples for ds001927, ds008893 and ds008064. Subsequently, the user requested
the full ds001927 release. **The download is complete** in `data/ds001927/`:
433 files and 22,968,601,962 bytes including metadata and the deposited study
code. Git metadata hashes and annex content digests all matched the pinned
version 2.1.0 manifest; SHA-256 hashes were also recorded for every local file.
All 207 NIfTI headers were checked, including 600 volumes at TR 1.25 s for each
of the 138 resting runs. The study code was downloaded but not executed.
The public aggregate receipt is
`PROVENANCE/additional_datasets_20261126/download_complete.json`.

Following download, the new input review confirmed retained skull tissue and
identified a damaged sub-s022 control-session T1. A six-run fMRIPrep pilot for
sub-s001 finished successfully; all six runs passed motion QC and alignment
review, with inferior signal coverage retained as a limitation. The second
pilot, sub-s022, is running; downstream inference remains pending. See the
[preprocessing plan](DS001927_PREPROCESSING_PLAN_20261126.md). No restricted
access was obtained and no authors were contacted. Download checksums establish
file integrity, not acceptable motion, registration or signal coverage.

Reproduce the check with `setup/audit_additional_tms_candidates.py`.
The local audit is `output/additional_dataset_audit_20261126/`; it includes a
rest/T1 download manifest, HTTP results, source hashes and full metadata. A small
aggregate access receipt is published in `PROVENANCE/additional_datasets_20261126/`.
The earlier PDC and EEG inventory remains documented in
[the previous candidate review](TMS_PUBLICATION_OPTIONS_20261125.md).
