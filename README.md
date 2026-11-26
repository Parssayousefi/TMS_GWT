# TMS_GWT

Information decomposition and workspace analyses of two open resting-state fMRI
datasets following prefrontal TMS. The manuscript cites the earlier psychedelic
paper; its private participant data and effect maps are not included.

## Current results 26 November 2026

**Working analysis release, not a finalized submission.** The publication goal
is a PID treatment-effect story supported across both datasets. Haemodynamic
associations and fingerprint/repeatability results are not the selected narrative.
Nonsignificant changes are not interpreted as stability or equivalence.

The latest DAN audit independently reconstructs the scan-level counts. The
normalized Aachen reduction is −0.467 gateways (unadjusted t p=.0382;
Holm4 p=.1527; exact sign-flip p=.0566; maxT4 p=.1621). The normalized SUDMEX
baseline-and-motion-adjusted active–sham estimate is +0.400 gateways (p=.422).
This is an exploratory Aachen signal, without a replicated reduction.
See [the latest audit and analysis status](docs/TMS_CURRENT_ANALYSIS_STATUS_20261126.md).

The deeper PID follow-up is also complete: six processing modes in seven fixed
pair groups, raw and separately normalized, plus regional change patterns tested
against templates from the other cohort. No shared treatment effect was
established. The smallest unadjusted SUDMEX baseline-and-motion p-value among the
42 mode outcomes is .347 raw and .350 normalized; all cross-dataset Holm6 and
Holm42 values are 1. All results and diagnostic checks are in the dashboard.
[The additional-data review](docs/TMS_ADDITIONAL_DATASETS_20261126.md) identifies
ds001927 as a downloadable 23-person, three-site TMS fMRI cohort. Metadata and
file headers were verified. Its full pinned release is now downloading to
`data/ds001927/`, with per-file checksum verification; it has not been analysed.

The corrected Aachen analysis is complete: 15 participants, 60 scans, all five
information measures and 20 participation/role outcomes, four spatial scopes,
fresh conditional-null normalization, and four GWT-related expectations.
The earlier significant global Φ finding is superseded: corrected change
+0.001349 nats, 95% CI −0.001525 to +0.004224, maximum-t p5=.7441.
No GWT expectation survived its four-outcome family. Redundancy participation
increased nominally (p=.0295), but Holm4 p=.1178; the amount of whole-brain
persistent redundancy was essentially unchanged (p=.9543).

Exploratory haemodynamic associations are present. Narrower estimated HRFs
were associated with greater Φ changes across the whole brain (Holm9 p=.00849)
and outside the workspace (p=.000477). These are same-signal, model-dependent
associations and do not establish a vascular mechanism, a uniform TMS effect,
or validation of the superseded global increase.

- [Interactive analysis dashboard](dashboard/index.html) covers full results,
  correction families, clinical outcomes, targeted flow and diagnostics.
- [Corrected Aachen working draft](manuscript/aachen_brief_report/corrected_20261126/manuscript.md)
  and [supplement](manuscript/aachen_brief_report/corrected_20261126/supplement.md).
- [Manuscript DOCX](manuscript/aachen_brief_report/corrected_20261126/Aachen_Brief_Research_Report_Corrected.docx)
  and [supplement DOCX](manuscript/aachen_brief_report/corrected_20261126/Supplementary_Material_Corrected.docx).
- [Supplementary tables and scalar outputs](manuscript/aachen_brief_report/corrected_20261126/TMS_Corrected_Supplementary_Data.zip).
- [Corrected analysis report](docs/TMS_CORRECTED_ANALYSES_20261126.md).

The two-dataset haemodynamic draft in `manuscript/tms_cross_dataset_report/20261126`
is retained as an abandoned framing, not the current submission draft. Its
underlying analyses remain available in the dashboard for a complete record.

The requested control-workspace image is Supplementary Figure S1, with a
section citing Luppi and Yousefi. Dice overlaps are .776/.741 for matched SUDMEX
sham scans and .743 for the corrected Aachen paired-baseline mean. These are
descriptive partial overlaps, not evidence of identical masks or an intervention
effect. The fixed 235-region HCP reference was not changed.

## Input routing and history

Aachen anatomy was already skull-stripped; forced repeat extraction had clipped
brain tissue. The corrected branch uses `--skull-strip-t1w skip` and all other
production processing definitions are unchanged. Valid current Aachen inputs:
`reruns/aachen_skullstrip_skip_20261125/pipeline/results/`.
Original Aachen outputs under `pipeline/results/` are historical and must not
be mixed with the corrected branch. SUDMEX inputs remain in `pipeline/results/`.
The separate SUDMEX anatomical audit found no comparable extraction failure;
inferior signal coverage remains a limitation in both datasets.

The historical proposal and 25 November exploratory reports are retained for
traceability and do not describe the current corrected Aachen findings. Current
numerical results are in `output/corrected_tms_analysis/20261126T112054_284547Z`;
the current local notebook is `output/jupyter-notebook/TMS_Corrected_Analyses.ipynb`.

## Repository contents and reproduction

Git tracks analysis code, configuration, small group-level references,
manuscript documents, summary tables and the dashboard. Raw imaging, individual
time series, large intermediate matrices, licenses, credentials, local logs,
third-party software installations and private psychedelic results are excluded.
The supplementary ZIP includes deidentified public-dataset scalar outputs.
Open the dashboard locally in a browser; it needs no server or external service.

## Sharing with collaborators

The public repository is at https://github.com/Parssayousefi/TMS_GWT.
Anyone with the link can view or download it without an invitation.
For write access, the owner can invite colleagues through **Settings → Collaborators → Add people**.
See [GitHub's invitation instructions](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/repository-access-and-collaboration/inviting-collaborators-to-a-personal-repository).

For a results-only review, send `dashboard/index.html`; the interactive dashboard
and control figure are embedded and work offline. To retain the individual source
CSV downloads, send the whole `dashboard/` folder. GitHub displays HTML source;
download the repository ZIP, extract it and open `dashboard/index.html` in a browser.

The executed notebook is in [notebooks/TMS_Corrected_Analyses.ipynb](notebooks/TMS_Corrected_Analyses.ipynb).
Its stored results can be reviewed without running it. Rerunning requires the
documented local analysis inputs and updating the explicit project-root path.

## Reproducing the analyses

With the documented local inputs and the `gwt` Python environment available:

```powershell
python setup/analyse_corrected_tms.py
python setup/normalize_corrected_aachen.py --analysis <new-output-folder>
python setup/extend_corrected_tms_checks.py --analysis <new-output-folder>
python setup/run_corrected_aachen_sensitivities.py --analysis <new-output-folder>
```

Publication-build scripts point to the released analysis folder for reproducible
reporting. Update their explicit run constants only when intentionally releasing
a new analysis. They must never fall back to historical Aachen inputs.
The pipeline requires MATLAB/CONN/rsHRF and original data for production stages;
the statistical summaries alone cannot reconstruct participant imaging.
Fresh selected-pair MATLAB verification for the targeted flow extension was
unavailable because of a license-server error. Independent Python identities
and comparisons with the already computed MATLAB matrices passed.

## Rules

1. **The psychedelic cohorts are frozen.** No stage is ever re-run for `hcp`, `basel_lam`, `basel_lpm`, `missouri_psilo`, `monash_psilo`, `london_lsd`, `berlin_ketamine` or `cambridge_propofol*`. Their results were copied byte for byte from `D:\projects\GWT` and are inputs here, not outputs.
2. **TMS data go through the exact scripts that produced those results**, in the same order (table below).
3. **Only per-dataset run settings change**, as they did for every psychedelic cohort: a dataset config in `pipeline/configs/datasets/` registered in `dataset_configs.yaml`, `DATASET` at the top of `step01b_conn_denoise.m`, and a TR entry in the `dataset_tr` map of `step02_hrf_deconvolution_matlab.m`.
4. **Processing code is not edited.** Anything that must change is recorded in `PROVENANCE/`.

## Production chain

Same scripts and order as the psychedelic cohorts (Monash and Missouri were the most recent runs).

| Step | Script (under `pipeline/scripts/`) | Runs in | Output |
| --- | --- | --- | --- |
| 1. fMRIPrep 23.2.1 | `preprocessing/run_fmriprep.ps1` (reads `fmriprep_config_loader.py`) | Docker, `nipreps/fmriprep:23.2.1`; MNI152NLin2009cAsym res-2, `--fs-no-reconall`, no distortion correction | `data/<dataset>/derivatives/fmriprep` |
| 2. Motion QC (Luppi criteria) | `qc/fmriprep_motion_qc.py` | Python | `results/qc/<dataset>/exclusions.csv` |
| 3. CONN inputs | `preprocessing/01a_prepare_conn_inputs.py` | Python | `pipeline/data/conn_input/<dataset>` |
| 4. CONN denoising (5 WM + 5 CSF aCompCor, 12 motion, 0.008 to 0.09 Hz, no GSR) | `preprocessing/step01b_conn_denoise.m` | MATLAB R2025a + CONN 25b + SPM12 | cleaned 4D BOLD (temporary) |
| 5. Parcellation, Schaefer-400 + Tian-54 | `preprocessing/01c_parcellate_conn_output.py` | Python | `results/timeseries/*_timeseries.npy` |
| 6. HRF deconvolution | `preprocessing/step02_hrf_deconvolution_matlab.m` | MATLAB + rsHRF v2.2 | `*_timeseries_deconv.npy` |
| 7. ΦID, Gaussian MMI, lag 1 TR | `analysis/step04_compute_phiid_matlab.m` (Luppi's `PhiIDFull.m`) | MATLAB | `results/phiid_matrices/` |
| 8. Workspace and gateway/broadcaster roles | `analysis/05_group_average_workspace.py`, `analysis/06_identify_gateways_broadcasters.py` | Python | `results/group_analysis/<dataset>/` |
| 9. Participant-level gateway counts and statistics | reference code: `pipeline/Global_Workspace_Complete_Pipeline_publication.ipynb`; robustness variants: `tier1_robustness/` | Python | `output/` |

Steps 2 to 7 pick up whatever is new on disk or take `--dataset`, so with no psychedelic timeseries in this folder they only ever process TMS scans. The HCP workspace used in step 8 is the copied `pipeline/results/group_analysis/hcp/rest/hcp_rest_workspace_mask.npy`.

## Layout

```text
TMS_GWT/
  README.md, CLAUDE.md, license.txt (FreeSurfer)
  docs/            proposal and figures
  data/            TMS datasets only: sudmex_tms (ds003037), aachen_itbs (ds001832)
  pipeline/        copied GWT pipeline: scripts, configs, atlases, tools, reference notebook
    results/       phiid_matrices, group_analysis, qc: frozen psychedelic + HCP results; TMS results are added here
  output/          publication-notebook-2026-09-26: frozen psychedelic result tables and figures
  tier1_robustness/  tier-1 code, notebook and final tables (Part XII)
  numit-phiid/     NuMIT-PhiID code, notes and tables
  setup/           copy_from_gwt.py, verify_copy.py
  PROVENANCE/      manifest with SHA-256 of every copied file, path rewrites, environment, verification report
  tmp/             TemplateFlow cache and fMRIPrep work dirs
  logs/
```

## What was copied from GWT

`setup/copy_from_gwt.py` lists every item. In short: production code, configs, atlases, the tools shipped in the repo (SPM12, rsHRF, npy-matlab, phyid), `PhiIDFull.m`, the ΦID matrices of every GWT cohort, the group-analysis outputs, QC exclusions, the publication notebook and its outputs, the tier-1 code, notebook and tables, and the NuMIT-PhiID code and tables.

Not copied, because nothing here re-runs the psychedelic cohorts: their raw and fMRIPrep data, CONN intermediates, regional timeseries, the tier-1 copied inputs and derived arrays, the NuMIT arrays, and the manuscripts.

Every copied file is byte-identical to its source except for one rule, applied only to code and config files (`.py .m .sh .ps1 .yaml .yml .json`): absolute paths into `D:\projects\GWT` now point into `D:\projects\TMS_GWT`. Each rewritten line is in `PROVENANCE/PATCHES_paths.tsv`, and `setup/verify_copy.py` proves that the rule is the only difference. Notebooks, Markdown and all result files are untouched, because they are the record of the psychedelic runs.

A side effect is intended: the psychedelic dataset configs now point at `data/<cohort>` folders that do not exist here, so an accidental preprocessing run on a psychedelic cohort fails instead of reading or writing anything in GWT.

## External requirements

- Docker Desktop running, with `nipreps/fmriprep:23.2.1` pulled (`docker pull nipreps/fmriprep:23.2.1`).
- MATLAB R2025a with its licence server reachable (CONN, rsHRF and ΦID steps). It was offline in October 2026.
- CONN 25b at `C:\Users\YousefiP\software\conn25b\conn`, the path hardcoded in `step01b_conn_denoise.m`.
- The `gwt` conda environment (`C:\Miniforge3\envs\gwt`); package versions are in `PROVENANCE/environment_gwt_*.txt`.
- Disk: about 250 GB for TMS fMRIPrep derivatives, plus about 8 GB of work space per subject while fMRIPrep runs.

## TMS datasets

Both are registered in `pipeline/configs/dataset_configs.yaml`, with TR 2.0 s added to the rsHRF script's TR map. Those two edits to copied files are recorded in `PROVENANCE/POST_COPY_EDITS.tsv` with a diff. `setup/check_tms_registration.py` checks the registration.

| Dataset ID | Source | Processed | Conditions (session names) | Notes |
| --- | --- | --- | --- | --- |
| `sudmex_tms` | ds003037, CC0 | 44 subjects (24 active, 20 sham), 142 rest runs | t0 (control), t1, t14, t2, t3, t4 | Baseline-only anatomy via `data/sudmex_tms/bids_filter.json`; duplicate ses-t5 excluded |
| `aachen_itbs` | ds001832, CC0 | 15 of 16 subjects, 60 rest runs (sub-02 excluded: 36-slice scans, slice timing undocumented) | pre90, post90, pre120, post120 | Sessions renamed from `ses-0`…`ses-4` by `setup/rename_aachen_sessions.py`; contrasts pair each post scan with its own pre scan |

Details, exclusions and scan parameters: `data/<dataset>/PREPROCESSING_NOTES.md`. Download log: `logs/download_openneuro.log`.

## Clinical input for TMS analyses

Published SUDMEX craving and impulsivity scores were downloaded from Zenodo
version 4.2 (record 10409461) and checksum-verified on 25 November 2026.
`data/master_clinical_sessions.csv` retains 154 visits from 53 participants.
`data/master_clinical_effects.csv` contains baseline, two-week and change scores
for the existing 44-person MRI cohort. Forty have usable MRI pairs after motion
QC; this flag does not imply that subsequent processing is complete.

The reviewable workbook is
`outputs/tms_clinical_20261125/TMS_clinical_outcomes.xlsx`. Its three tabs contain
clinical changes, all clinical visits, and a dictionary with source-quality
issues. Change means T1 minus T0 in the original scale units; negative values
mean lower scores. Later open-label visits remain in the sessions table.

Use `pipeline/scripts/utils/clinical_data.py` to load either the workbook or
CSV and join by `dataset`, `subject_id`, and `comparison`. The copied ASC/drug
analysis is not a clinical loader. `pipeline/configs/clinical_outcomes.yaml`
documents the interface. The loader masks flagged outcomes by default, with
explicit `include_flagged=True` available for a declared sensitivity analysis.
The original reported totals remain unchanged on disk. No clinical hypothesis
tests have been run by this preparation step.

```python
from pathlib import Path
import sys
sys.path.insert(0, str(Path("D:/projects/TMS_GWT/pipeline/scripts/utils")))
from clinical_data import load_clinical_effects, merge_clinical_neural

clinical = load_clinical_effects(
    "D:/projects/TMS_GWT/outputs/tms_clinical_20261125/TMS_clinical_outcomes.xlsx",
    metric="craving_vas",
)
# neural_changes must have one row per dataset/subject_id/comparison.
# joined = merge_clinical_neural(neural_changes, metric="craving_vas")
```

Rebuild the source tables with `setup/prepare_sudmex_clinical.py`; the source
downloads and checksums are managed by `setup/download_sudmex_clinical.py`.
Validation is in `setup/test_clinical_data.py`. Processing and HRF verification
records are described in `PROVENANCE/TMS_PROCESSING_20261125.md`.

## TMS publication notebook

Open [TMS_Global_Workspace_Publication.ipynb](output/jupyter-notebook/TMS_Global_Workspace_Publication.ipynb)
with the **Python (gwt)** kernel. The 29-cell notebook follows the original GWT
publication structure but contains only the TMS analyses: global PhiID, mean
workspace participation, participant gateway/broadcaster counts, SUDMEX
active-versus-sham changes, Aachen paired changes, craving/impulsivity effects
and associations, motion sensitivity and descriptive follow-up.

The notebook reads the frozen 235-region HCP workspace without recalculating it.
It does not load psychedelic participant matrices or effect maps. Its proposed
primary endpoint is SUDMEX FPN gateway-count change; this is an analysis draft,
not a submitted preregistration. SUDMEX uses between-arm change differences;
Aachen averages its two intensity-specific changes within each of 15 people.
The published psychedelic findings enter only a qualitative comparison table.

Each run validates available TMS inputs and saves a new timestamped directory
under `output/tms_publication/`. Inference waits until the expected primary
pairs are processed (currently 40 SUDMEX participants and 15 Aachen participants
with both intensity visits). Incomplete input is visibly pending. The notebook
does not start preprocessing; select **Run All** again when new PhiID outputs
are available.

Implementation: `pipeline/scripts/analysis/tms_publication.py`. Rebuild the
notebook's source cells with `setup/build_tms_publication_notebook.py`; doing so
clears its cell outputs. Execute with
`C:\Miniforge3\envs\gwt\python.exe setup/execute_tms_publication_notebook.py`.
Numerical and design checks are in `setup/test_tms_publication.py`.

While an existing TMS processing run is active, use
`C:\Miniforge3\envs\gwt\python.exe -u setup/watch_tms_publication.py` to execute the
notebook automatically once each full dataset is ready. The monitor audits HRF
sidecars first and archives the executed notebook in that run's output folder.
It does not start preprocessing. Its progress, execution logs and restart state
are in `logs/tms_pipeline_20261125/publication_watch/`. Use `--status-only` for a
read-only readiness check. Local Jupyter communication requires the same host
permissions as the verified standalone notebook execution.

Start of processing, once Docker Desktop is running:

```text
.\pipeline\scripts\preprocessing\run_fmriprep.ps1 -Dataset sudmex_tms -Parallel 4
```
