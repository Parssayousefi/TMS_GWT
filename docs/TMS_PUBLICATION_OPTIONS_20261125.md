# Publication options and public-data feasibility — 25 November 2026

PDC is the strongest new scientific fit if approved data access is already
available or can be obtained promptly. OpenNeuro ds008064 is genuinely
downloadable, but its missing treatment allocation prevents the intended
active-versus-sham analysis. For the present deadline, retain fMRI as the main
route. An EEG switch would constitute a separate analysis pipeline, rather than
a quick extension of the completed analysis.

This is a prospective recommendation for additional work, not a claim that any
unanalysed dataset contains a positive result. No full imaging/EEG collection was
downloaded, no restricted access was obtained, and no authors were contacted.

## Fit to the special issue

The [Frontiers call](https://www.frontiersin.org/research-topics/80360/psychedelic-medicine-and-neuromodulation-from-bench-to-bedside-advancing-real-world-implementation-of-next-generation-psychiatric-treatments)
accepts neuroimaging work on individual as well as combined interventions. It
does not require a psychedelic-versus-TMS experiment. The listed submission
deadline is 16 January 2027, 52 days from this review. A
[Brief Research Report](https://www.frontiersin.org/journals/psychiatry/for-authors/article-types)
allows 4,000 words and four figures/tables combined. A narrow exploratory paper
is therefore a reasonable format; interpretable controls and uncertainty still
matter more than obtaining a nominally significant result.

## What was actually checked

Public release documentation, repository metadata, participant tables, scan
sidecars, and trial registration were inspected. The cached audit can be rerun:

```powershell
& 'C:\Miniforge3\envs\gwt\python.exe' 'D:\projects\TMS_GWT\setup\audit_public_candidate_metadata.py' --cached
```

Source metadata and summaries are in
`output/candidate_dataset_audit_20261125/`. `anonymous_download_checks.json`
records successful anonymous HTTP 206 reads of the first 64 KiB of a ds008064
NIfTI and a TBS EEG ZIP. These verify actual file access, but are not full-file
integrity or imaging-quality checks. All sample counts below precede QC.

| Candidate | Verified usable design inventory | Access status | Main limitation |
|---|---|---|---|
| PDC fMRI | 33 ECT baseline/course-end pairs; 60 ketamine baseline/fourth-infusion pairs; 17 healthy repeat-scan pairs | Public manual and inventory; imaging requires approved NDA/CCF access | Treatments were not randomly assigned; post-treatment timing differs |
| ds008064 fMRI | 28 listed participants, 26 pre and 24 post resting scans, **23 matched resting pairs** | Public metadata and sample imaging bytes downloaded | No active/sham labels found in the release |
| Moffa prefrontal TBS EEG | 24 healthy participants, 119 visits; iTBS twice, cTBS twice, sham once per participant except one missing visit | Public files and condition labels verified | New EEG workflow; full block archive about 180 GB |
| Pallavicini naturalistic DMT EEG | 35 participants; DMT plus eyes-open/eyes-closed baselines; subjective scales | Public Zenodo record, approximately 1.5 GB; scale file downloaded | Naturalistic, no placebo; does not match TBS acquisition/eyes condition |

## PDC: preferred new clinical study, conditional on access

The [release page](https://www.humanconnectome.org/study/crhd-perturbation-treatment-resistant-depression-connectome-fast-acting-therapies/document/pdc-release-10)
links NDA Collection 2844 and states that acquisition follows HCP-Aging. It
contains minimally preprocessed imaging and clinical/behavioral measures.
The [public completeness inventory](https://www.humanconnectome.org/storage/app/media/documentation/PDC1.0/PDC_1.0_Release_Completeness.csv)
gives these resting-scan counts:

| Group | Baseline | Intermediate | Course end / repeat | Later follow-up |
|---|---:|---:|---:|---:|
| ECT | 44 | — | 33 | 21 |
| Ketamine | 67 | 49 after first infusion | 60 after fourth infusion | 45 |
| Healthy controls | 39 | — | 17 | — |

The [manual](https://www.humanconnectome.org/storage/app/media/documentation/PDC1.0/PDC_1.0_Release_Manual.pdf)
describes clinically assigned treatments. ECT follow-up is within a week after
the index course; ketamine scans are approximately 24 hours after the first and
fourth infusions. Healthy controls repeat after approximately two weeks. They
are not an untreated depressed group or a sham-ECT arm. The ketamine endpoint
therefore measures post-treatment change, not the acute altered state.

Access is controlled through NDA/CCF and institutional data-use approval. The
public files do not establish that the user's account is authorized. Access
instructions are linked on the release page; historical processing estimates
should not be treated as a current delivery promise. Obtain a resting-state and
clinical subset rather than the multi-terabyte complete release.

The proposed question is: **Do ECT and ketamine produce shared or distinct
changes in brain information processing after treatment?** A working title is
*Information storage and integration after electroconvulsive therapy and
ketamine treatment for depression*.

Suggested bounded analysis:

1. Use one course-end endpoint per treatment: ECT visit 2, ketamine visit 3.
   Keep early ketamine and later follow-ups as separate secondary analyses.
2. Use mean active information storage and redundancy-corrected integration
   Phi_R as two main measures. Treat raw Phi, individual atoms, workspace roles,
   and network localization as supporting analyses. These choices are proposed
   now, before inspecting PDC outcomes; they do not retroactively redefine the
   existing TMS analysis.
3. Estimate each patient's within-person change and compare with repeat-scan
   controls. If testing ECT-control, ketamine-control, and ECT-ketamine for both
   main measures, use one six-test Holm family. Report estimates and confidence
   intervals alongside corrected p values.
4. Model baseline neural values and motion; consider age, clinical severity,
   medications, and treatment selection with a parsimonious, documented model.
   Adjustment cannot turn the treatment-group comparison into a randomized
   experiment. Clinical change associations are a separate secondary family.
5. Compare raw and conditional-null-normalized information measures. Confirm
   scan length, lag in seconds, parcellation and denoising compatibility before
   interpreting differences from our existing cohorts.

Similar directions are not by themselves proof of a common mechanism, and a
nonsignificant ECT-ketamine contrast does not establish equivalence. The clinical
efficacy of ECT does not guarantee a PhiID effect. Also, ECT is not the only
neuromodulation intervention with demonstrated antidepressant efficacy:
[sham-controlled rTMS evidence](https://pubmed.ncbi.nlm.nih.gov/20439832/) exists.

## ds008064: promising replication once allocation is supplied

The [version 1.0.0 repository](https://github.com/OpenNeuroDatasets/ds008064/tree/1.0.0)
was inspected at commit `d7b64e9a048559df210f9923d57a814ed9370379`.
The participant table contains only ID, pre/post MADRS, age, and sex at birth.
No allocation fields were found in the 50 resting scan JSON files or additional
clinical/session tables in the repository tree. The missing resting scans are
pre-only D03/D05/D24, post-only D18, and no resting scan at either visit for D16.
Thus the paired imaging sample is 23, not 28.

The MRI protocol identifier links to
[NCT03289923](https://clinicaltrials.gov/study/NCT03289923), a randomized study of
six-week active/sham left-DLPFC rTMS with psychotherapy. The release and protocol
describe psychotherapy in somewhat different terms, so both stimulation
allocation and the participant-level psychotherapy protocol need clarification.
Do not assume an unverified factorial design.

A sample NIfTI header confirms 200 volumes with TR = 3 seconds. This does not
automatically match our 2-second information lag. The data are suitable for a
TMS analysis once identifiers can be joined to allocation. A pooled pre/post
analysis without labels would mix stimulation, sham, psychotherapy and time.

TRD is not a fatal confound for an active-sham comparison within TRD. An
independently defined healthy HCP workspace can be used as a reference ROI, but
its use does not validate that workspace's function in depression. Whole-brain
outcomes should carry the main inference; avoid normative claims based solely
on deviation from healthy maps. Acquisition harmonization alone also does not
make synergy estimates directly comparable.

## What Mediano's papers change

Interpreting “media group” as Pedro Mediano's group, three sources are especially
relevant:

- [Mediano et al., PNAS 2025](https://doi.org/10.1073/pnas.2423297122) provides the
  PhiID taxonomy: information storage, transfer and integration are distinct
  processes. A single global Phi value does not identify which changed.
- [Liardi et al., PLOS Computational Biology 2025](https://doi.org/10.1371/journal.pcbi.1013629)
  introduces NuMIT, showing why information components must be interpreted
  relative to total information and a specified null model.
- [Blackburne et al., October 2026 preprint, version 2](https://www.biorxiv.org/content/10.64898/2026.10.01.755915v2)
  already applies information decomposition to five psychedelic M/EEG datasets
  (80 participants), including ketamine. It reports widespread changes in
  information storage/transfer, while normalized persistent synergy is not a
  universal significant increase. It also extends NuMIT to PhiID. Its code is
  public at [infodelia](https://github.com/gblackburne/infodelia), but the paper
  directs readers to original study authors for data access. The official
  bioRxiv JATS full text and version metadata were downloaded for this review.

The latest preprint supplies a direct PhiID-normalization precedent that was
missing from our earlier note citing only the 2025 PID paper. It does not
validate every choice in our fMRI implementation. Nor are its acute M/EEG
results interchangeable with post-treatment fMRI or Yousefi's workspace roles.
Simply repeating that psychedelics affect synergy would have limited novelty.
Testing information-process changes under neuromodulation is more useful.

## EEG alternative: attractive design, not the quickest deadline route

The [Moffa et al. dataset paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9681877/)
describes 24 healthy participants with repeated iTBS/cTBS and sham visits.
There are two pre and three post eyes-open resting blocks per visit. Resting
blocks follow the TMS-EEG probe blocks, so their timing should not simply be
labelled as the probe's 2/15/30-minute time points. Analyse the rest blocks,
with an EEG-specific artifact and reference strategy. The
[public block archive](https://plus.figshare.com/articles/dataset/EEG_files_separated_per_blocks_from_the_study_effects_of_TBS_of_the_DLPFC_with_TMS-EEG_/19297701)
and condition table were verified. Most of its approximately 180 GB includes
more than the desired resting epochs.

A strong bounded question is **whether information-dynamic changes after
prefrontal TBS are reproducible and sham-specific**. Average post-block metrics
within each visit, subtract the mean pre-block metric, and retain the two active
visits for repeatability analysis. Average repeated visit changes only when
forming each participant's condition contrast. Repeated blocks and visits are
not independent participants. Published
[TMS-evoked analyses](https://doi.org/10.1111/ejn.70121) found no significant main
effects in the tested components/bands, so a positive information result should
not be presumed.

The [Pallavicini DMT record](https://zenodo.org/records/3992359) is immediately
accessible and includes subjective scales. It is naturalistic and lacks
placebo; its DMT eyes-closed comparison also differs from the eyes-open TBS
protocol. Comparing within-dataset response profiles is possible, but raw
cross-dataset synergy differences would conflate intervention and acquisition.
Sensor mixing, reference, montage, epoch length and lag require dedicated EEG
choices. Do not reuse fMRI deconvolution or a volumetric workspace mask on EEG.

## Recommendation for the current manuscript

Prioritize PDC if NDA/CCF access is already active or promptly feasible. Obtain
ds008064 allocation before treating it as a sham-controlled replication.
Otherwise finish a focused fMRI Brief Research Report, rather than tying the
deadline to a new multimodal pipeline. The TBS EEG study is a good separate
option, especially for reliability; adding DMT should answer a specific
information-process question rather than merely broaden the paper.

The completed TMS data support a cautious narrative about information amount
versus organization: Aachen's averaged raw Phi increase passes its global
family correction, with borderline normalized support, but has no sham control
and is not robust to every corrected model. SUDMEX does not retain the FPN
gateway effect after baseline adjustment, and the all-follow-up raw synergy
changes do not yield a stable normalized counterpart across sample choices.
This is not evidence that TMS reproduces the psychedelic workspace signature.

As agreed, retain the two-week randomized SUDMEX comparison for intervention
claims and the all-post average as a longitudinal sensitivity analysis.
Existing outputs are documented in `TMS_EXPLORATORY_FINDINGS_20261125.md` and
`TMS_ALL_POST_AVERAGE_20261125.md`; no new neural results were generated during
this candidate review.
