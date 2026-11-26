# TMS processing run started 25 November 2026

User authorization: continue processing the TMS data after connecting the VPN.
MATLAB R2025a (25.1.0.2943329) can reach the licence server; MATLAB, Parallel
Computing, Signal Processing and Wavelet licences passed the initial check.

## Inputs and scope

- Aachen iTBS: 60 scans from 15 participants after the existing exclusions.
- SUDMEX TMS: 133 scans from 44 participants after the existing motion QC
  decisions (124 INCLUDE, 9 TRUNCATE; 9 other scans excluded).
- Continue from the completed fMRIPrep and motion QC stages. No new exclusions
  or analytical parameter changes are introduced by this run.
- Run the copied production CONN input preparation, CONN denoising,
  parcellation, rsHRF and PhiID scripts. Group contrasts and hypothesis tests
  are outside this execution run. No preregistration has been submitted by
  this operation.
- No psychedelic or HCP analyses are rerun. The manuscript's proposed use of
  psychedelic data has changed to citation only; this run uses TMS inputs.

## Execution wrapper

`setup/run_tms_processing.py` is a new orchestration wrapper. It calls the
production scripts without editing them. MATLAB runtime copies change only
`DATASET` and/or `N_WORKERS`; source and runtime hashes and the settings are
saved next to each runtime copy.

The four-scan Aachen subject 01 pilot completed CONN and parcellation, producing
four finite 295-by-454 arrays with no constant regional signals. The bulk run
resumes these outputs, processes four participants per batch using eight CONN
workers and up to four Python workers, then uses twelve workers for rsHRF and
PhiID. These are execution settings, not analytical parameter changes.

The wrapper checks exact expected scan membership, dimensions, finite signals,
metadata presence and matrix symmetry. It stops on validation failures rather
than silently excluding failed scans. It validates the resolved paths of
regenerable CONN scratch before removing it, after regional arrays have passed
validation. It retains manifests, CONN project metadata, logs and array hashes.

## Pilot observations relevant to methods

- The production script explicitly sets `Setup.voxelresolution = 1`. CONN 25b
  therefore samples the fMRIPrep images onto its analysis-mask grid. Pilot
  cleaned images have dimensions 91-by-109-by-91, 2 mm voxels and LAS storage,
  whereas inputs have dimensions 97-by-115-by-97 and RAS storage. This is
  affine-coordinate resampling, not a new nonlinear template registration.
  Evidence: `conn_process.m` uses the mask as `Vref` and an empty `sfile`;
  `conn_create_vol.m` samples source images using their affine coordinates.
- The canonical NLin2009cAsym atlas is resampled by the unchanged parcellation
  script onto that analysis grid. All 454 labels remain, with at least 42
  voxels per label in the pilot geometry. The statement in the parcellation
  docstring implying that cleaned images retain the original voxel lattice
  is imprecise; voxel lattice and anatomical template are distinct.
- The directory is named `tools/spm12/spm`, but its `Contents.m` and execution
  log identify **SPM25, version 25.01.02 (23 January 2025)**. The installed
  copied version is used unchanged. Methods should report the actual version,
  not infer it from the folder name.
- CONN's exported NIfTI temporal spacing is zero in the pilot header. The
  manifest and regional metadata retain the correct TR of 2 s; the production
  rsHRF dataset lookup also uses 2 s.

## Run records

`logs/tms_pipeline_20261125/` contains the pilot log, saved manifest, geometry
check, array validation, scratch cleanup record and runtime provenance.
`logs/tms_pipeline_20261125/bulk/status.json` records the current stage or final
status; `bulk_runner.log` records stage transitions. Each subject batch has its
own directory under `bulk/` with logs, manifest, validation and runtime files.

This document records the launch and validated pilot, not successful completion
of the entire run. The final status and output validation must be checked.

## HRF verification and orchestration correction

On 25 November, after all 60 Aachen scans completed CONN and parcellation,
rsHRF was run on those scans using a runtime copy of the unchanged production
script, with `N_WORKERS=4`. Before this verification run, there were no Aachen
deconvolution arrays in this workspace and the bulk runner had not reached its
HRF stage. This establishes the successful new run, not the outcome of an
unidentified earlier MATLAB execution.

The scan loop completed in 1.7 minutes: 60 successes, zero errors and zero
skips. Each input has only 294 or 295 timepoints by 454 regions, so this runtime
is not comparable to voxelwise preprocessing. All 60 output arrays, MAT HRF
estimates and YAML sidecars were inspected. Shapes, finite values, nonconstant
signals, saved settings (TR 2 s, gamma order 3, threefold temporal upsampling,
threshold 1 SD, 24 s HRF, 4-8 s lag search and 50 Wiener iterations) and input
hashes passed. Outputs are neither copied inputs nor simple z-scores, and
output hashes differ across all scans.

An independent spot check reran HRF fitting and Wiener deconvolution on regions
1, 227 and 454 of three scans, including a 294-timepoint scan. Maximum absolute
deconvolution differences were zero; reconstructed HRF differences were at
most 1.12e-16. Evidence is under
`logs/tms_pipeline_20261125/hrf_verification/`: `hrf.log`, runtime provenance,
pre-run `inputs.json`, `aachen_itbs_output_audit.csv/json` and
`recompute_sample.m/json/log`.

Two of 27,240 regional HRFs have a 16 s time-to-peak, outside the descriptive
3-10 s band: sub-03 pre120 region 175 and sub-05 post120 region 440 (one-based
indices). These are flagged for scientific QC, not automatically removed or
treated as an implementation failure. Numerical replay verifies execution,
not physiological ground truth or the appropriateness of all HRF assumptions.

The new orchestration wrapper originally relocated the HRF script without its
sibling helper `parsave_hrf.m`. This could leave a saved NPY but fail before
saving HRF metadata. The wrapper now copies that helper byte-for-byte next to
the runtime script and records its hash. The identical helper was also staged
next to the future bulk HRF runtime because the already-running Python process
had loaded the earlier wrapper. New wrapper executions additionally require
both HRF sidecars during output validation. These changes affect orchestration
and validation only; no production processing algorithm was edited. The active
bulk process continues its original scan/parameter plan and will skip the
already-computed Aachen arrays when it reaches rsHRF.

## Clinical data preparation and interface

Five public source workbooks were obtained from Zenodo record 10409461,
version 4.2: CCQ-General, CCQ-Now, BIS-11, clinical craving VAS and demographics.
All published MD5 checksums match. URLs and SHA-256 hashes are recorded in
`data/sudmex_tms/clinical/source/zenodo-10409461/download_manifest.json`.

New preparation code `setup/prepare_sudmex_clinical.py` joins by source rid and
visit, verifies uniqueness and BIDS demographics/arm agreement, preserves
published scores, and maps source T1-4 to MRI t14. All 53 clinical participants
and 154 visits remain in `data/master_clinical_sessions.csv`. The existing
44-person MRI cohort alone enters `data/master_clinical_effects.csv`, which
contains T0, T1 and T1-minus-T0 values in original scale units. Forty MRI pairs
pass existing motion QC (22 active, 18 sham); processing completion is a
separate condition. Clinical identification of the three participants absent
from BIDS does not reinstate them in the MRI cohort.

Four source assessments require review: CCQ-General 001 T1 q22=0, 008 T1-4
q17=11, 020 T1 q11 missing; BIS-11 013 T1 q12=5. Published totals are never
silently repaired. The default new loader masks affected instrument changes,
with explicit opt-in for reported-score sensitivity analyses. BIS subscales
sum to its published total in all 153 assessments. One later BIS assessment
is absent. We have not reconstructed or psychometrically revalidated totals
from items. Response-label inconsistencies in the BIS dictionary are preserved.

The new `pipeline/scripts/utils/clinical_data.py` accepts the master CSV or
the `Clinical changes` worksheet and performs a one-to-one join with neural
changes using `dataset`, `subject_id`, `comparison`. This is a TMS interface;
the copied ASC/drug analysis is not repurposed by renaming clinical measures.
The seven clinical measures and input paths are documented in
`pipeline/configs/clinical_outcomes.yaml`. No primary clinical endpoint is
selected and no hypothesis test/correlation is run by this preparation.

The Excel workbook in `outputs/tms_clinical_20261125/` was built with the bundled
artifact-tool runtime. Score inputs link to the clinical sessions sheet and
changes are Excel formulas. All three tabs were rendered and visually checked.
Eight automated checks passed, including equality of all workbook values and
CSV values, valid zeros versus missingness, source flags, rejection of duplicate
identities/wrong contrasts/incorrect changes, and joins independent of row
order with randomized-arm agreement. Validation counts and source issues are
saved beside the workbook. The builder is retained at
`tmp/clinical_workbook_20261125/build.mjs`.

## TMS publication notebook (25 November 2026)

Created `output/jupyter-notebook/TMS_Global_Workspace_Publication.ipynb` from
the Jupyter skill experiment scaffold, then populated it with 29 cells (15
executable). It mirrors the relevant structure of the copied publication
notebook without its drug classification, ASC or cross-drug Mantel analyses.
The original notebook and original GWT project were not edited.

New helper `pipeline/scripts/analysis/tms_publication.py` implements the TMS
analysis interface. Participation coefficients and within-workspace rank-based
roles agree exactly with the original notebook functions on the same numerical
fixtures. The eight-module atlas mapping is independently checked against the
label table used by that notebook. The HCP workspace mask is loaded unchanged;
no HCP or psychedelic computations are rerun.

The draft specification contains 25 neural outcomes: five unique-pair global
PhiID means, two mean workspace PCs, total gateway/broadcaster counts, and both
role counts for each of eight networks. Missing workspace values remain
missing, tied ranks are neither role, and pairwise global changes require a
matching valid-ROI set. A complete TMS scan requires all five matrices plus
the completed summary sidecar. Shape, symmetry, missing-edge structure,
metadata identity/coverage and the Phi_R identity are checked before use.

The proposed primary endpoint is SUDMEX FPN gateway-count change. The draft
uses a two-sided Welch comparison of active and sham participant changes and
reports unstandardized confidence intervals plus between-arm Hedges' g. The
secondary neural FDR families are explicitly separate: other network gateways
(7), network broadcasters (8), mean PCs (2), global PhiID (5) and total roles
(2). A baseline- and delta-FD-adjusted HC3 model is a sensitivity analysis.
These are new design-specific statistical choices, not a claim that the old
paired drug tests can be used unchanged or that a preregistration exists.

Aachen pairs each intensity with its own baseline and averages both paired
changes within each person, preserving n=15 independent units. Its one-sample
tests and paired g_z concern pre/post change without a sham control. Intensity
results are descriptive. All per-dataset inference requires the complete
planned primary cohort, rather than an opportunistic subset of completed jobs.
Later SUDMEX original-arm summaries are descriptive because exposure becomes
open label and follow-up samples change.

The workbook/CSV clinical interface is integrated. Seven clinical treatment
effects form one BH family. Five neural measures (FPN/DAN gateways, the two
workspace PCs and global Phi_R) crossed with seven clinical measures form the
35-test exploratory association family. OLS with HC3 uncertainty estimates
clinical change per SD of neural change, adjusting for randomized arm, both
baseline measures and delta FD. Source flags and per-model complete-case counts
remain visible. No psychedelic phenomenology is substituted or renamed.

Thirteen validation tests passed: exact metric parity, known PC values, role/tie
accounting, unique-edge means, missing-workspace behavior, Welch/paired tests
against SciPy, BH values, within-person intensity pooling, pending-cohort
behavior, synthetic clinical joins and adjusted models, absence of drug-analysis
code, handling of zero/unestimable comparisons, regional pairing/coverage, and a synthetic figure layout check
(several properties share one test). Synthetic fixtures never enter empirical
publication outputs.

The notebook executed all 15 code cells successfully using the gwt environment.
Its initial snapshot at `output/tms_publication/20261125T113925_916849Z/`
validated 12 actual Aachen scans and reported the empirical group analyses as
pending. It saved input hashes, configuration, versions, pair flow and scan
coverage. This records a successful execution with available inputs, not a
completed analysis of either full cohort. Each later Run All writes a new
timestamped output directory. The checked notebook uses the installed gwt kernel.

## Early Aachen PhiID execution and timing evidence

To make the first TMS results available before SUDMEX denoising completes, the
60 completed Aachen deconvolved scans entered a separate execution of the
unchanged production PhiID script with `N_WORKERS=6`. Runtime copy, hash,
60-file input inventory and log are in
`logs/tms_pipeline_20261125/aachen_phiid/`. At launch, no SUDMEX deconvolved
inputs existed. The main bulk process remains in SUDMEX CONN/parcellation; when
it later reaches PhiID, the original existence guard skips completed Aachen
outputs. This changes scheduling and worker count, not estimator parameters.

A 1,000-pair timing benchmark (after 10 warmup pairs) used the unchanged
PhiIDFull function on an Aachen deconvolved scan. It took 3.70 s, extrapolating
to 381 s per 102,831-pair scan. The first full parallel scans took 6.4-6.5 min,
with later scans around 5.8-5.9 min. The benchmark saves timing only, not
inferential results. Evidence: `phiid_benchmark.m/json/log` in the run-log root.
Runtime estimates are conditional on continued license access and compute
availability; successful completion still requires the final output checks.

## Completion-triggered notebook execution

At the user's request to keep processing and analyse each completed dataset,
`setup/watch_tms_publication.py` monitors the existing production run. It requires
every QC-eligible scan in a dataset to have regional data, all HRF sidecars, all
five PhiID arrays and the completed PhiID summary. It then audits the saved HRF
outputs and executes the publication notebook with the existing specification.
No estimator, preprocessing setting, hypothesis or multiplicity family changes.

The monitor preserves an executed notebook and a SHA-256 execution receipt in
each completed dataset's timestamped output folder. Its persistent state avoids
repeating a finished dataset on restart. Pipeline errors or notebook errors are
reported rather than silently marked complete. Logs and HRF audits are kept in
`logs/tms_pipeline_20261125/publication_watch/`. The saved-input audit checks
numerical integrity; it is not a claim that deconvolution perfectly recovers
latent neuronal activity.

The public bioRxiv text was rechecked on 25 November for the comparison:
FPN gateway increases were +4.35 (LSD LAM), +4.60 (LSD LPM), +4.52 (psilocybin)
and +4.03 (mescaline) regions. All four increased DAN gateways and mean workspace
redundancy PC; synergy PC increased in three of the four conditions. These are
published aggregate findings, not reanalysis of private psychedelic data.
Source: https://www.biorxiv.org/content/10.64898/2026.09.28.755172v1.full .
The future TMS comparison is descriptive across different study designs, not a
formal drug-versus-stimulation hypothesis test.

## Early SUDMEX deconvolution

While the bulk process was denoising subjects 016/017/018/019, the first 42
validated SUDMEX regional scans entered a separate copy of the unchanged HRF
script with four workers. The other 60 inputs were completed Aachen scans and
were skipped by the production guard. The runtime copy and unchanged helper,
hashes, 102-file preflight inventory, complete log, QC summary and 42-scan audit
are in `logs/tms_pipeline_20261125/sudmex_early/`.

All 42 SUDMEX scans completed successfully, with zero errors. The scan loop took
1.7 minutes; process startup is additional. Every output was finite, differed
from both raw and merely standardized input, and had valid MAT/YAML sidecars.
Input hashes were unchanged. Estimated time-to-peak values were 4.67-8.00 s.
This advances ready scans through the same ordered chain without changing
analytical settings. No other HRF writer was active. Early SUDMEX PhiID must
wait until the separate Aachen PhiID process finishes to avoid overlapping
discovery of unfinished Aachen inputs.

## First complete Aachen notebook result

The separate Aachen PhiID run completed 60/60 scans with zero errors in 59.4
minutes. All 60 summaries reported no zero-variance regions and a successful
nonnegative-Phi_R check. At 12:23:58 UTC the monitor triggered the full-cohort
notebook; execution finished at 12:24:14 UTC. Outputs and the archived executed
notebook are in `output/tms_publication/20261125T122407_467376Z/`.

All 15 participants contributed both intensity pairs, with matching regional
coverage. The 25 estimates, confidence intervals, tests and BH-adjusted values
were independently rechecked against SciPy and statsmodels. Figure layouts
were visually inspected. FPN gateway change averaged -0.0667 regions (95% CI
-2.1957 to 2.0624; p=0.9474). No network gateway or broadcaster effect and neither
mean PC passed its specified correction. Global Phi increased by 0.001792 nats
(95% CI 0.000515 to 0.003070; global-family q=0.04689); Phi_R q=0.05731.
These are the first inspected complete Aachen results under specification v0.1,
not preregistered findings. No hypothesis or correction family was changed.

After Aachen exited successfully, early SUDMEX PhiID started with eight workers
on the 42 completed SUDMEX deconvolved scans; all 60 Aachen scans are skipped.
Its runtime settings, 102-input hash inventory and log are in `sudmex_early/`.
The bulk process was still preparing subjects 020/021/022/023, leaving the
remaining seven denoising batches before its own global HRF/PhiID stages.

## MATLAB graphics-helper resource cleanup

Inspection during processing identified 15 orphaned MATLAB R2025a CEF renderer
processes consuming substantial CPU after their parent processes had exited.
Only those individually identified renderers were stopped, after rechecking
their executable path, creation timestamp, renderer command line and absent
parent PID. No live MATLAB numerical process was stopped. The inspection and
cleanup records are `matlab_window_processes.json` and
`orphan_renderer_cleanup.json` in the run-log root.

New helper `setup/inspect_orphan_matlab_renderers.ps1` supports the same bounded
check for later job exits: read-only by default, with `-Apply` to stop verified
orphans and save a timestamped audit. It excludes active renderers, browser
parents, numerical workers, other executables and processes created before this
run. This is process cleanup only; no processing code or numerical setting changes.

## Completed early SUDMEX stages

Four sequential early executions advanced 90 ready SUDMEX scans through the
unchanged HRF and PhiID scripts while the bulk run continued CONN/parcellation.
Each used four HRF workers and eight PhiID workers. A run started only after
the preceding PhiID process exited successfully and while the bulk process was
still several subject batches away from its own HRF/PhiID stages. Completed
outputs were skipped by the original production guards.

| Run-log folder | New SUDMEX scans | HRF scan loop (min) | PhiID scan loop (min) | Errors |
|---|---:|---:|---:|---:|
| `sudmex_early` | 42 | 1.7 | 36.8 | 0 |
| `sudmex_early_02` | 21 | 0.3 | 12.8 | 0 |
| `sudmex_early_03` | 14 | 0.4 | 9.7 | 0 |
| `sudmex_early_04` | 13 | 0.2 | 8.6 | 0 |

Process/pool startup is additional to these logged scan-loop times. All runtime
copies, dependency hashes, preflight input hashes and per-execution HRF quality
CSVs are retained in those folders. Repeated technical audits confirmed all 90
SUDMEX HRF outputs and sidecars, unchanged source arrays and nontrivial
deconvolution. Their estimated HRF time-to-peak values ranged from 4.67 to 8 s.
The aggregate timing record is `early_sudmex_execution_summary.json`.

The final early process exited at approximately 13:51 UTC while the bulk run
was denoising subjects 034/036/037/039. No further early execution is scheduled;
the bulk run handles the remaining 43 SUDMEX scans. No SUDMEX group hypothesis
tests have been inspected at this point; its notebook inference still waits for
the full planned cohort.

## Recovery from a transient Windows status-file lock

At 14:12:07 UTC, the bulk runner stopped before launching CONN for subjects
046/047/048/049: atomic replacement of `bulk/status.json` raised Windows error
5 (access denied). Writing the failure status immediately afterward succeeded,
consistent with a transient file-sharing lock. This was an orchestration/status
write failure, not a failed scientific calculation. The notebook monitor
reported the failure and stopped rather than declaring the cohort complete.
At that point 115 SUDMEX regional scans, 90 SUDMEX PhiID scans, and all 60
Aachen scans remained verified on disk.

`setup/run_tms_processing.py::write_json` now retries only `PermissionError`
from atomic replacement, at most 20 attempts with capped backoff (under 9 s of
total waiting). It preserves the previous valid JSON and still raises persistent
errors. `setup/watch_tms_publication.py` uses the same writer and retries
transient reads of the bulk status file. Two regression checks passed: recovery
after two simulated sharing violations, and bounded failure with preservation
of the previous state. Test: `setup/test_tms_status_write.py`.

The failure record, relevant log tail, previous preparation log, and pre-fix
orchestration source are retained in `logs/tms_pipeline_20261125/status_lock_recovery/`.
The bulk runner was resumed with the same settings, appending to its existing
log and skipping already verified scans. No production MATLAB/Python processing
script, neural statistic, hypothesis or correction family was changed.

## Full-cohort HRF integrity audit

All 193 regional scans reached deconvolution. The final bulk HRF execution
processed the remaining 43 SUDMEX scans and skipped the completed 150 inputs.
Audits at 14:47 UTC verified all 60 Aachen and 133 SUDMEX arrays and their
MAT/YAML sidecars. Every output was finite and differed from both raw and
merely standardized input. No input ROI had zero variance. All 193 original
regional-array hashes matched their saved pre-deconvolution records, and all
193 source hashes were distinct.

Evidence is in `logs/tms_pipeline_20261125/final_validation/`: the combined
input manifest, both dataset audit CSV/JSON files, a 193-row consolidated HRF
quality CSV, and the six source-report paths/hashes. Each production HRF quality
CSV covers only that execution's newly processed scans; the consolidated report
restores complete coverage without modifying the per-execution originals.

Four scans each have one fitted regional HRF peak at 16 s, outside the audit's
descriptive 3-10 s range: Aachen 03/pre120 ROI 175 and 05/post120 ROI 440;
SUDMEX 047/t14 ROI 59 and 052/t1 ROI 58 (one-based ROI indices). These are
recorded quality flags, not newly introduced exclusions. The other regions in
those scans satisfy the range; no parameter or mathematical-integrity check
failed. Numerical integrity does not establish perfect recovery of neuronal
activity.

## Completed processing and first full-cohort results

The bulk runner completed at 14:59:03 UTC on 25 November 2026 (15:59 Zurich),
with 193 verified scans through PhiID. Its final PhiID job completed 43 scans,
skipped 150 previously verified outputs, and reported zero errors. Both the
runner and the notebook monitor exited with code 0. The legacy MATLAB log's
suggested HCP command was not executed; the frozen HCP workspace was only read.

The monitor executed the full notebook at 14:59 UTC, after all 133 SUDMEX scans
were ready. The final run is `output/tms_publication/20261125T145928_921789Z/`,
including the archived executed notebook, input provenance, execution receipt,
CSV tables and 11 figures in both PNG and PDF. All 15 code cells completed
without errors. All 193 scans have 454 valid regions and all 235 fixed workspace
regions. The primary paired cohorts are 40 SUDMEX participants (22 active,
18 sham) and 15 Aachen participants, each with both intensity visits.

`setup/verify_tms_publication_results.py` independently checked 25 SUDMEX and
25 Aachen neural tests, seven clinical treatment tests, 25 sensitivity models
and 35 clinical association models. It used SciPy t tests, statsmodels intervals
and BH correction, and a direct NumPy HC3 covariance calculation. The saved
`independent_statistics_check.json` records success, table hashes, unchanged
Aachen results, and unchanged analysis code/clinical workbook since the first
Aachen execution. Only the notebook's HRF-audit narrative was updated before
the SUDMEX run; no analysis code, endpoint or multiplicity family was changed.

The unadjusted SUDMEX FPN gateway difference in change was +3.1414 regions
(95% CI 0.0202 to 6.2626; p = 0.04861; the sole draft primary endpoint).
Mean change was -0.1364 in active and -3.2778 in sham, so this is a relative
preservation pattern rather than an observed increase within the active arm.
The already specified baseline-neural plus motion sensitivity model estimated
+1.0851 (95% CI -2.1301 to 4.3003; p = 0.49806). This attenuation cannot be
attributed to motion alone: baseline neural value and motion were adjusted
together. No secondary SUDMEX neural endpoint survived its declared FDR family.

Clinical active-minus-sham changes were -2.3206 VAS points for craving
(95% CI -4.2186 to -0.4225; q = 0.03154; n = 40) and -12.2807 BIS-total points
(95% CI -20.8564 to -3.7051; q = 0.01472; n = 39). Cognitive and motor BIS
subscales also survived the seven-outcome clinical family; CCQ-General,
CCQ-Now and nonplanning did not. None of the 35 covariate-adjusted
brain-clinical associations survived their single FDR family (minimum q =
0.96036). Clinical improvement does not establish a workspace mechanism.

Aachen's results are unchanged from its earlier completed run: FPN gateway
change -0.0667 regions (95% CI -2.1957 to 2.0624; p = 0.9474), with no supported
gateway, broadcaster or participation effect. An exploratory global Phi
increase of 0.001792 nats survived its five-measure family (q = 0.04689).
The dataset has no sham control. Phi_R did not meet q < 0.05 (q = 0.05731).

The published Yousefi paper reports classic-psychedelic FPN gateway increases
of approximately 4.0-4.6 regions, increased DAN gateways and increased mean
redundancy participation. The SUDMEX unadjusted FPN direction overlaps that
pattern, but its sensitivity and arm-specific changes make it fragile;
the broader pattern is not clearly reproduced in these TMS analyses. This
is a descriptive cross-study comparison, not a formal test that TMS and drugs
differ. No private psychedelic matrices or effect maps were reanalysed.
The notebook records the specification as not preregistered.

Published comparator:
https://www.biorxiv.org/content/10.64898/2026.09.28.755172v1.full

## Subsequent user-requested diagnostics

After these results were reported, the user requested separation of baseline
and motion adjustment and an assessment of Liardi et al.'s NuMIT method.
The post-hoc plan and full findings are in
`SUDMEX_FPN_DIAGNOSTIC_PLAN_20261125.json` and
`SUDMEX_FPN_DIAGNOSTICS_20261125.md`. Outputs are in
`output/tms_diagnostics/20261125T152110_809161Z/`.
The completed primary notebook and all production inputs remain unchanged.
The NuMIT work was a bounded method pilot, not a full normalized treatment
analysis or a rerun of any private cohort.

## Completed exploratory extension and review

The bounded exploratory screen, delayed-active SUDMEX analysis, full TMS-only
normalization branch, and family-wise comparisons are now complete. The report
is `docs/TMS_EXPLORATORY_FINDINGS_20261125.md`; its companion is
`output/jupyter-notebook/TMS_Exploratory_Findings.ipynb`. All seven notebook code
cells executed without errors on 25 November 2026 at 16:17 UTC.

Final family-wise tables and independent verification are in
`output/tms_familywise/20261125T161051_102285Z/`. The normalized branch is in
`output/tms_numit/20261125T155252_384543Z/`. The global Phi result passes
within-global-family max-T correction in raw and normalized Aachen analyses;
Phi_R does not. The effect is the mean of two matched visit changes per
participant, not 30 independent observations. Its estimate is stable to
measured motion adjustment, but that HC3 sensitivity does not retain
five-measure Holm support. The report retains all these distinctions.

Plans, numerical conditioning notes, and the review of one undefined secondary
HC3 model are retained in the corresponding TMS_*_20261125 provenance files.
Independent verification covered 977 estimable neural contrasts, 100 delayed
active contrasts, seven clinical models, and the permutation/CDF calculations.
The original completed publication outputs and all production processing code
remain unchanged; no psychedelic or HCP analyses were run.

## Baseline versus all available post visits

At the user's request, a separate all-post-average analysis was completed in
`output/tms_all_post_average/20261125T162521_698094Z/`. It uses equal visit
weights within each participant and equal participant weights: 41 SUDMEX
baselines plus 88 follow-ups, and all 60 Aachen scans. A sensitivity retains
the original 40-person acute SUDMEX cohort. Aachen is algebraically unchanged.
Normalization was extended to 49 additional scans using existing, unchanged
reference grids and baseline pools. The plan is
`TMS_ALL_POST_AVERAGE_PLAN_20261125.json`; findings and interpretation are in
`docs/TMS_ALL_POST_AVERAGE_20261125.md`.

Independent checks reproduced all 150 paired comparisons, 300 original-arm
contrasts, 14 clinical models, participant weighting, source hashes, and the
full main max-T distributions using direct SciPy statistics. The separate
temporal-information audit checked the Phi identity in all 189 used scans
without adding hypothesis tests. The updated exploratory notebook executed
all ten code cells without errors at 16:32 UTC. Earlier results remain intact.
