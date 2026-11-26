# GWT Analysis Pipeline - Multi-Dataset Support

## Project Overview
This project implements the Integrated Information Decomposition (ΦID) framework from Luppi et al. (2024) to analyze psychedelic resting-state fMRI data. The pipeline is **dataset-agnostic** and supports multiple datasets with different BIDS structures.

## Reference Papers
- **Luppi et al. (2024)** — "A synergistic workspace for human consciousness revealed by Integrated Information Decomposition" (eLife, 12:e88173). Local copy: `docs/Luppi_2023.md`
- **Luppi et al. (2022)** — "A synergistic core for human brain evolution and cognition" (Nature Neuroscience). Local copy: `docs/luppi_2022.md`

## Key Concepts
- **Synergistic Information**: Information available only when multiple brain regions are considered together
- **Redundant Information**: Information equally available from multiple sources
- **Integrated Information (ΦR)**: Revised measure of consciousness-related information integration
- **Synergistic Workspace**: Brain regions with high synergistic interactions (DMN as gateways, FPN as broadcasters)

## Supported Datasets

### HCP (Human Connectome Project — 100 Unrelated)
- **Path**: `D:\projects\GWT\data\hcp`
- **Structure**: No sessions (single REST1_LR run per subject)
- **Subjects**: 100 healthy controls
- **Conditions**: REST only (baseline for workspace definition)
- **TR**: 0.72s
- **Volumes Removed**: 10 first
- **Run Selection**: first (REST1_LR only, per Luppi methodology)
- **Config**: `configs/datasets/hcp.yaml`
- **Role**: Defines the canonical synergistic workspace applied to all drug datasets

### London LSD Dataset
- **Path**: `D:\projects\GWT\data\london_lsd\derivatives\fmriprep`
- **Structure**: Has sessions (`sub-XX/ses-Y/func/`)
- **Subjects**: 20 participants (sub-01 through sub-20, excluding sub-08)
- **Conditions**: LSD and PCB (placebo)
- **TR**: 2.0s
- **Runs**: 3 runs per condition, **using run-3 only** (specific strategy)
- **Config**: `configs/datasets/london_lsd.yaml`

### Basel LAM Dataset
- **Path**: `D:\projects\GWT\data\basel_lam\derivatives\fmriprep`
- **Structure**: Has sessions (`sub-XX/ses-Y/func/`)
- **Subjects**: 28 subjects (001-029, excluding 007 and 023)
- **Conditions**: LSD, MDMA, amph (amphetamine), placebo
- **TR**: 1.8s
- **Runs**: First run per session
- **Config**: `configs/datasets/basel_lam.yaml`

### Basel LPM Dataset
- **Path**: `D:\projects\GWT\data\basel_lpm\derivatives\fmriprep`
- **Structure**: Has sessions (`sub-XX/ses-Y/func/`)
- **Subjects**: 33 subjects (001-034, excluding 023)
- **Conditions**: LSD, psilocybin, mescaline, placebo
- **TR**: 1.8s
- **Runs**: First run per session
- **Config**: `configs/datasets/basel_lpm.yaml`

### Berlin Ketamine Dataset
- **Path**: `D:\projects\GWT\data\berlin_ketamine`
- **Structure**: Has sessions (condition determined by session name)
- **Subjects**: 17 participants
- **Conditions**: Ketamine and placebo
- **TR**: 3.4s
- **Runs**: First run per session
- **Config**: `configs/datasets/berlin_ketamine.yaml`

### Missouri Psilocybin Dataset
- **Path**: `D:\projects\GWT\data\missouri_psilo`
- **Structure**: Has sessions (`sub-PX/ses-Y/func/`)
- **Subjects**: 7 main subjects (P1-P7; replication subjects excluded)
- **Conditions**: Baseline (control), psilocybin, methylphenidate (active control)
- **TR**: 1.761s
- **Runs**: First run per session
- **Note**: Multi-echo acquisition (5 echoes); fMRIPrep native optimal combination, combined BOLD only (as Monash). Provenance: `data/missouri_psilo/PREPROCESSING_NOTES.md`
- **Config**: `configs/datasets/missouri_psilo.yaml`

## Dataset System Architecture

### Configuration System
```
configs/
├── dataset_configs.yaml          # Registry of all datasets
├── datasets/
│   ├── hcp.yaml                 # HCP 100 unrelated (workspace baseline)
│   ├── london_lsd.yaml          # London LSD
│   ├── basel_lam.yaml           # Basel LAM
│   ├── basel_lpm.yaml           # Basel LPM
│   ├── berlin_ketamine.yaml     # Berlin Ketamine
│   └── missouri_psilo.yaml      # Missouri Psilocybin
└── paths.yaml                    # Global output paths
```

### Core Utilities
```
scripts/utils/
dataset_handler.py           # Universal data discovery & path resolution
analysis_loader.py           # Dataset-aware group-level data loading
qc_loader.py                 # QC exclusion loading
load_roi_mapping.py          # ROI-to-network mapping loader
create_roi_network_mapping.py # Generate ROI mapping CSV
```

## Project Structure
```
GWT/
CLAUDE.md                    # This documentation file
README.md                    # Project README
docs/
  Luppi_2023.md              # Luppi et al. 2024 eLife (SAPHIRE architecture)
  luppi_2022.md              # Luppi et al. 2022 Nat Neurosci (synergistic core)
  DEBUG_luppi_comparison.md  # HCP workspace comparison investigation
  luppi_supps/               # Luppi supplementary data & code
    elife-88173-fig2-data1-v1.xlsx   # Published HCP matrices (454 ROIs)
    software/MMI_PhiID/              # Original MATLAB PhiID code
scripts/
  preprocessing/
    01_extract_timeseries_luppi_method.py    # Timeseries extraction (Python)
    step02_hrf_deconvolution_matlab.m        # HRF deconvolution (MATLAB/rsHRF)
    step04_compute_phiid_matlab.m            # PhiID (MATLAB, optional)
    fmriprep_config_loader.py                # fMRIPrep config management
  analysis/
    04_compute_phiid_revised.py              # PhiID computation (Python/phyid)
    04b_compute_phiid_parallel.py            # Parallel PhiID computation
    05_group_average_workspace.py            # Group-level workspace definition
    06_identify_gateways_broadcasters.py     # Gateway/broadcaster classification
  qc/
    fmriprep_motion_qc.py                   # Motion QC analysis
    generate_supplementary_material.py       # Supplementary figures/tables
  data_acquisition/                          # BIDS conversion & data setup scripts
  test_against_luppi/                        # Validation vs Luppi reference
  utils/
    dataset_handler.py           # Universal data discovery
    analysis_loader.py           # Dataset-aware loading
    qc_loader.py                 # QC exclusion loading
    load_roi_mapping.py          # ROI-to-network mapping
    create_roi_network_mapping.py
configs/
  dataset_configs.yaml             # Dataset registry
  datasets/                        # Per-dataset configs (6 datasets)
  paths.yaml                       # Global output paths
atlases/
  schaefer_400_tian_54_combined.nii.gz  # CANONICAL combined atlas (454 ROIs)
  schaefer_400.nii.gz                   # Schaefer 400 cortical only
  tian_subcortex_s4_3T.nii.gz           # Tian S4 54 subcortical
  roi_to_network_mapping.csv            # ROI-to-network mapping (CANONICAL)
tools/
  integrated-info-decomp/          # phyid library (Python PhiID)
  rsHRF/                           # rsHRF toolbox v2.2 (MATLAB)
  spm12/                           # SPM12 (MATLAB)
  npy-matlab/                      # NumPy .npy reader for MATLAB
results/
  timeseries/                      # All datasets share output directories
  phiid_matrices/
  group_analysis/
  networks/
  statistics/
logs/
```

## Quick Start

### Running the Pipeline

**For any dataset:**
```bash
# Discover available scans
python scripts/utils/dataset_handler.py --dataset london_lsd --discover
```

**Per-dataset CONN extraction (production):**
```bash
# 1a — prepare CONN inputs (Python, ~5 min/scan)
python scripts/preprocessing/01a_prepare_conn_inputs.py --dataset hcp

# 1b — run CONN denoise (MATLAB, ~8 min/scan; edit DATASET at top of script)
# In MATLAB: cd('D:\projects\GWT\pipeline'); run('scripts/preprocessing/01b_conn_denoise.m')

# 1c — parcellate to production-format timeseries (Python, ~1 min/scan)
python scripts/preprocessing/01c_parcellate_conn_output.py --dataset hcp
```

Optional filters (subject/session/task) on stages 1a and 1c:
```bash
python scripts/preprocessing/01a_prepare_conn_inputs.py --dataset basel_lam --subjects 001 --sessions 1 --tasks LSD
```

**Legacy fMRIPrep-aCompCor pipeline (deprecated, kept for reference):**
```bash
# Single command per scan (Python, ~1-2 min/scan)
python scripts/preprocessing/01_extract_timeseries_luppi_method.py --dataset hcp
```

## Analysis Pipeline

### Phase 1: Data Preprocessing (CONN-based, three sub-stages)

**Production extraction uses CONN V25.b's native aCompCor** (validated at n=20 vs Luppi 2024 — see `scripts/test_against_luppi/WORKSPACE_DISCREPANCY.md`). The legacy script `01_extract_timeseries_luppi_method.py` (fMRIPrep-precomputed aCompCor + nilearn) is preserved for reference but **deprecated**.

**1a. Prepare CONN Inputs** (`01a_prepare_conn_inputs.py`, Python)
- Dataset-agnostic via `dataset_handler.discover_all_scans()`. Handles all registered datasets uniformly (HCP, London LSD, Basel LAM, Basel LPM, Berlin Ketamine, Cambridge Propofol, Monash PsiConnect, Missouri Psilocybin).
- Per-scan input dirs at `data/conn_input/{dataset}/{scan_id}/` where scan_id = `sub-XXX[_ses-Y]_task-Z_run-N`.
- For each scan: trims BOLD per QC + uniform trim (first N only, no end trim), decompresses BOLD + T1 to plain `.nii`, extracts 6-column motion.txt from fMRIPrep confounds.
- **aCompCor is NOT bridged from fMRIPrep** — only motion + T1 cross over. The 5 WM + 5 CSF aCompCor components are re-computed natively by CONN in stage 1b from the T1 `structural_segment` tissue masks, matching Luppi 2022/2024 exactly.
- Honours `results/qc/{dataset}/exclusions.csv` via `qc_loader.QCExclusions` — EXCLUDE → skip, TRUNCATE → `select_volume_range`.
- Writes `data/conn_input/{dataset}/manifest.json` covering all valid scans.

**1b. CONN Denoise** (`01b_conn_denoise.m`, MATLAB)
- Edit `DATASET = 'hcp'` (etc.) at top of script, then `run('scripts/preprocessing/01b_conn_denoise.m')`.
- Reads manifest, runs CONN V25.b setup (structural_segment + ROI extraction) then denoising.
- Confound matrix: 5 WM aCompCor + 5 CSF aCompCor + 6 motion + 6 derivatives. Linear detrend, RegBP, bandpass 0.008–0.09 Hz, no GSR, no ART, no condition effect, no despiking (paper-3 HCP methodology applied uniformly).
- CONN-cleaned 4D BOLDs copied to `results/timeseries_intermediate/{dataset}/{scan_id}_conn_cleaned.nii`.
- Two-pass batch (Setup → Denoising) per CONN HCP example to avoid step-5 ordering bug.
- Wall time ~8 min/scan; serial fallback for PhiID via `parfor (..., max_workers)` if Parallel Toolbox license unavailable.

**1c. Parcellate** (`01c_parcellate_conn_output.py`, Python)
- Loads each cleaned 4D, applies `atlases/schaefer_400_tian_54_combined.nii.gz` via `NiftiLabelsMasker(standardize=False, detrend=False, low_pass=None, high_pass=None)` — pure spatial averaging since CONN already denoised.
- Output filename via `dataset_handler.get_output_filename()` = production format `dataset-{name}_sub-XXX[_ses-Y]_task-Z_run-N_timeseries.npy` written to `results/timeseries/`.
- `_metadata.yaml` sidecar matches production schema (TR, n_timepoints, denoising config).
- Auto-cleans intermediate cleaned 4D (~4.3 GB) and per-scan input scratch (~5 GB) after successful parcellation. Use `--keep-cleaned` or `--keep-inputs` to override.

Per-dataset workflow:
```bash
python scripts/preprocessing/01a_prepare_conn_inputs.py --dataset hcp
# then in MATLAB: edit DATASET in 01b_conn_denoise.m and run it
python scripts/preprocessing/01c_parcellate_conn_output.py --dataset hcp
```

**2. HRF Deconvolution** (`step02_hrf_deconvolution_matlab.m`)
- **MATLAB script** using rsHRF toolbox v2.2 (Wu et al. 2013)
- 3-stage pipeline: event detection (1 SD threshold), HRF estimation (Gamma basis order 3, 24s length), iterative Wiener deconvolution (50 iterations)
- **CRITICAL**: Required for accurate information-theoretic analysis
- Reads `*_timeseries.npy`, outputs `*_timeseries_deconv.npy` (via npy-matlab)
- Prerequisites: SPM12, rsHRF toolbox, npy-matlab (all in `tools/`)

### Phase 2: Integrated Information Decomposition

**3. Compute PhiID Measures** (`04_compute_phiid_revised.py`)
- **Dataset-aware**: Handles flexible naming and sessions
- Calculate time-delayed mutual information (TDMI) between all ROI pairs
- Decompose using Gaussian MMI-PID solver from phyid library
- **Primary synergy measure**: `sts` atom (persistent synergy / double-synergy), matching Luppi supplementary code
- Non-finite values are filtered before averaging (`isfinite` guard, matching MATLAB `PhiIDFull.m`)
- Measures saved per ROI pair:
  - `synergy`: persistent synergy (`sts` atom) -- used for workspace definition
  - `redundancy`: persistent redundancy (`rtr` atom)
  - `synergy_total`: Eq. 5 total synergy (TDMI - max individual MIs) -- saved for comparison
  - `phi`: Eq. 6 integrated information
  - `phi_r`: Eq. 7 revised integrated information (Phi + Red)
- Output: `sub-{subject}_[ses-{session}_]task-{task}_run-{run}_{measure}.npy`

**Why `sts` / `rtr` alone, not `sts+stx+sty+str` / `rtr+rtx+rty+rts`?**
The eLife paper's `Red()` function definition `min{I(X_past;[X_fut,Y_fut]), I(Y_past;[X_fut,Y_fut])}` is the MMI *source-side* redundancy and in the PhiID lattice equals the sum `rtr+rtx+rty+rts` (coarse source, all targets). Luppi's published "redundancy matrix" is NOT this sum — it is the **persistent (double-) redundancy atom `rtr`** alone (redundant in both source and target), and synergy is the **persistent synergy atom `sts`** alone. This is explicit in Luppi's own supplementary README (`scripts/test_against_luppi/luppi_supps/software/MMI_PhiID/README.md`, lines 41-64: `synergy_mat(row,col) = atoms.sts; redundancy_mat(row,col) = atoms.rtr`). Our pipeline matches the supplementary code verbatim — `step04_compute_phiid_matlab.m` calls `PhiIDFull.m` directly and assigns `syn_mat = atoms.sts`, `red_mat = atoms.rtr`, `phi_r = phi + atoms.rtr`. The sum form is computed as a diagnostic (`synergy_total`) only. The atom choice has been empirically tested against Luppi's Fig-2 matrices (see `scripts/test_against_luppi/diagnostic_compare_luppi.py`): `sts` gives a stronger fit than `str+stx+sty+sts`, so the single-atom choice is both definition-matching and empirically better. Do NOT "fix" this by summing atoms.

**4. Group-Level Workspace Analysis** (`05_group_average_workspace.py`)
- Compute group-average synergy and redundancy matrices
- Identify synergistic workspace regions (high synergy rank)
- Calculate node strengths for synergy/redundancy networks
- Supports `--task` argument for analyzing specific conditions
```bash
# Control condition (default)
python scripts/analysis/05_group_average_workspace.py --dataset london_lsd

# Drug condition
python scripts/analysis/05_group_average_workspace.py --dataset london_lsd --task LSD
```

### Workspace Definition (CRITICAL)

**The synergistic workspace is defined ONCE from the HCP dataset (100 healthy controls)
and applied to ALL drug datasets.** This follows Luppi et al. (2024) methodology.

- Workspace source: `results/group_analysis/hcp/rest/hcp_rest_workspace_mask.npy`
- All downstream scripts (06, 06b, 07, 08) default to `--workspace_dataset hcp`
- Do NOT define workspace from drug datasets — this would be circular reasoning
- To override (rare): pass `--workspace_dataset <other_dataset>` explicitly

### Phase 3: Network Analysis

**5. Gateway/Broadcaster Classification** (`06_identify_gateways_broadcasters.py`)
- Calculate participation coefficients for synergy/redundancy networks
- Classify workspace regions:
  - **Gateways** (DMN): High synergy participation coefficient
  - **Broadcasters** (FPN): High redundancy participation coefficient

**Workspace Dataset (`--workspace_dataset`):**

All downstream scripts default to `--workspace_dataset hcp`, loading the canonical HCP workspace. This ensures the workspace is defined from healthy resting-state data and applied consistently across all drug analyses.

```bash
# Uses HCP workspace by default — no flag needed
python scripts/analysis/06_identify_gateways_broadcasters.py --dataset london_lsd

# Drug condition analysis (HCP workspace is default)
python scripts/analysis/06_identify_gateways_broadcasters.py \
    --dataset london_lsd \
    --task LSD \
    --control_task PCB

# Override: use dataset-specific workspace (rare)
python scripts/analysis/06_identify_gateways_broadcasters.py \
    --dataset london_lsd \
    --workspace_dataset london_lsd
```

## Key Implementation Details

### Dataset Handler Utility

**Core Functions:**
```python
from dataset_handler import load_dataset_handler

# Initialize handler
handler = load_dataset_handler('london_lsd')

# Auto-discover subjects (respects include/exclude in config)
subjects = handler.get_subjects()

# Get sessions for a subject (returns [None] if no sessions)
sessions = handler.get_sessions('01')

# Get tasks for subject/session
tasks = handler.get_tasks('01', session='2')

# Select run based on dataset strategy (last, specific, first)
run = handler.select_run('01', 'MDMA', session='2')

# Get file paths
func_file = handler.get_functional_file('01', 'MDMA', 2, '2')
confounds = handler.get_confounds_file('01', 'MDMA', 2, '2')

# Generate standardized output filename
filename = handler.get_output_filename(
    subject='01', task='MDMA', run=2, session='2',
    suffix='timeseries', extension='npy'
)

# Discover ALL available scans
scans = handler.discover_all_scans()
# Returns list of dicts: [{'subject': '01', 'session': '2', 'task': 'MDMA',
#                          'run': 2, 'func_file': '...', 'confounds_file': '...'}]
```

### Analysis Loader Utility

**For group-level analyses:**
```python
from analysis_loader import (
    load_group_matrices,
    load_group_all_measures,
    load_baseline_group_data,
    load_drug_comparison_data,
    get_control_task_name,
    discover_processed_scans
)

# Discover what's been processed
processed = discover_processed_scans('basel_mmm', task_filter='MDMA')

# Load single measure for all subjects
synergy_mats, subjects = load_group_matrices(
    dataset='basel_mmm',
    task='MDMA',
    measure='synergy'
)
# Returns: (n_subjects, 454, 454), list of subject IDs

# Load all measures at once
all_mats, subjects = load_group_all_measures(
    dataset='basel_mmm',
    task='MDMA'
)
# Returns: {'synergy': array, 'redundancy': array, 'phi': array, 'phi_r': array}

# Get control task automatically
control_task = get_control_task_name('basel_mmm')  # Returns 'placebo'

# Load baseline data
control_mats, subjects, control_task = load_baseline_group_data('basel_mmm')

# Load drug vs control comparison
control_mats, drug_mats, control_subj, drug_subj = load_drug_comparison_data(
    dataset='basel_mmm',
    drug_task='MDMA'
)
```

### ROI-to-Network Mapping (CRITICAL)

**IMPORTANT**: Always use the centralized ROI-to-network mapping system.

**Canonical Source**: `atlases/roi_to_network_mapping.csv`
- **456 lines total**: 1 header + 455 data rows
- **ROI 0**: Background (Network_ID = -1, not used in analysis)
- **ROI 1-400**: Schaefer 400 cortical parcels (7 Yeo networks)
- **ROI 401-454**: Tian 54 subcortical structures

**CSV Columns:**
| Column | Purpose | Example Values |
|--------|---------|----------------|
| `ROI` | ROI index (1-454) | 1, 2, ... 454 |
| `Network_ID` | Network ID for PC computation (8 networks) | 0-7 |
| `Network_Name` | Network name for PC computation | Visual, Default, Subcortical |
| `Structure_Name` | Detailed structure name for labeling | Hippocampus, Thalamus, 7Networks_LH_Vis_1 |
| `Atlas` | Source atlas | Schaefer400, Tian54 |
| `Hemisphere` | Hemisphere | LH, RH |
| `ROI_Label` | Original atlas label | HIP-head-m1-rh, 7Networks_LH_Vis_1 |

**CRITICAL - Column Usage:**
- **For participation coefficient (PC) computation**: Use `Network_ID` and `Network_Name` (8 networks)
- **For detailed result labeling**: Use `Structure_Name` (e.g., "Hippocampus" instead of "Subcortical")

This design keeps PC computation consistent with Luppi et al. (8 networks) while allowing detailed labeling of subcortical structures.

### ROI Indexing Convention (CRITICAL)

**Atlas ROIs are 1-indexed (1-454), Python arrays are 0-indexed (0-453).**

| Context | Index Range | Example |
|---------|-------------|---------|
| Atlas CSV (`roi_to_network_mapping.csv`) | 1-454 | ROI 1 = first Visual parcel |
| Python arrays (matrices, rsn_labels) | 0-453 | `array[0]` = first Visual parcel |
| CSV output files (ROI column) | **Use 1-454** | Match atlas convention |

**Rule**: When writing CSV files with an ROI column, use `np.arange(1, N_ROIS + 1)` not `np.arange(N_ROIS)`:
```python
# CORRECT - matches atlas convention
summary_data = {'ROI': np.arange(1, N_ROIS + 1)}  # 1, 2, 3, ... 454

# WRONG - off-by-one from atlas
summary_data = {'ROI': np.arange(N_ROIS)}  # 0, 1, 2, ... 453
```

**Loading Function:**
```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / 'utils'))
from load_roi_mapping import load_roi_network_mapping

BASE_DIR = Path(__file__).parent.parent.parent
rsn_labels, rsn_names, roi_df = load_roi_network_mapping(BASE_DIR)

# Returns:
# - rsn_labels: np.ndarray (454,) with Network_ID for each brain ROI
#   (indices 0-453 map to ROI 1-454)
# - rsn_names: list of 8 network names
# - roi_df: pd.DataFrame with full mapping (455 rows including background)

# For detailed labeling of results:
structure_names = roi_df[roi_df['ROI'] > 0]['Structure_Name'].values  # (454,)
```

**Network ID Mapping (for PC computation):**
- 0: Visual
- 1: Somatomotor
- 2: DorsAttn (Dorsal Attention)
- 3: SalVentAttn (Salience/Ventral Attention)
- 4: Limbic
- 5: Cont (Control/Frontoparietal)
- 6: Default (Default Mode Network)
- 7: Subcortical

**Structure Names (for detailed labeling):**
- Cortical ROIs: Same as ROI_Label (e.g., 7Networks_LH_Vis_1)
- Subcortical ROIs: Hippocampus, Thalamus, Striatum, Amygdala, Globus_Pallidus

**Plot-only label columns** (added by `scripts/utils/add_schaefer17_subgroups.py`):
- `Network17_ID` — Schaefer 17-net parcel ID (1..400). Empty for subcortical.
- `Network17_Subgroup` — 17-net subgroup token (VisCent, VisPeri, SomMotA/B, DorsAttnA/B, SalVentAttnA/B, LimbicA/B, ContA/B/C, DefaultA/B/C, TempPar). Empty for subcortical.
- `Substructure` — fine subcortical label (HIP / AMY / NAc / PUT / CAU / GP / vTHA / dTHA), parsed from Tian S4 `ROI_Label`. Empty for cortical.

For axis labels and visualization ONLY. Do NOT use `Network17_*` for participation-coefficient module partition — PC computation is locked to the 8-module partition (`Network_ID 0..7`) for Luppi compliance.

**7-net → 17-net is not 1:1.** The augmentor maps each 7-net parcel to its modal-overlap 17-net parcel. ~14% of cortical parcels at network boundaries fall outside their parent 7-net group in the 17-net partition — e.g. occipito-temporal Visual parcels (`7Networks_LH_Vis_1`) get reassigned to `DorsAttnA_TempOcc`, medial-temporal Visual parcels go to `DefaultC_PHC`. This is the genuine Schaefer 7 vs 17 disagreement (the two parcellations are generated independently in Schaefer 2018), not a bug. The plot will show those parcels in the 17-net block they spatially overlap with.

### File Naming Convention

**Without sessions (London LSD):**
```
sub-{subject}_task-{task}_run-{run}_{suffix}.{ext}

Examples:
- sub-01_task-LSD_run-3_timeseries.npy
- sub-01_task-PCB_run-3_timeseries_deconv.npy
- sub-01_task-LSD_run-3_synergy.npy
```

**With sessions:**
```
sub-{subject}_ses-{session}_task-{task}_run-{run}_{suffix}.{ext}

Examples:
- sub-01_ses-2_task-MDMA_run-2_timeseries.npy
- sub-01_ses-2_task-MDMA_run-2_timeseries_deconv.npy
- sub-16_ses-3_task-MDMA_run-2_phi_r.npy
```

## Adding a New Dataset

### Step 1: Create Dataset Config
Create `configs/datasets/my_dataset.yaml`:
```yaml
dataset_id: "my_dataset"
name: "My Dataset Name"

paths:
  bids_root: "/path/to/bids/root"
  fmriprep_dir: "/path/to/fmriprep/derivatives"

bids_structure:
  has_sessions: true  # or false

subjects:
  auto_discover: true  # or list explicitly

tasks:
  drug:
    bids_task_name: "drug"
    description: "Drug condition"
    is_control: false
  placebo:
    bids_task_name: "placebo"
    description: "Placebo"
    is_control: true

run_selection:
  strategy: "last"  # Options: "last", "first", "specific", "all"

file_patterns:
  bold: "sub-{subject}_[ses-{session}_]task-{task}_run-{run}_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz"
  confounds: "sub-{subject}_[ses-{session}_]task-{task}_run-{run}_desc-confounds_timeseries.tsv"
```

### Step 2: Register Dataset
Add to `configs/dataset_configs.yaml`:
```yaml
datasets:
  my_dataset:
    name: "My Dataset"
    config_file: "configs/datasets/my_dataset.yaml"
    description: "Description"
```

### Step 3: Test & Run
```bash
# Test discovery
python scripts/utils/dataset_handler.py --dataset my_dataset --discover


```

## Key Python Dependencies
- `pybids`: BIDS dataset handling
- `nilearn`: Neuroimaging data manipulation
- `jpype1` + JIDT: Information theory calculations (via phyid)
- `bct-py`: Brain Connectivity Toolbox for network metrics
- `statsmodels`: Statistical analysis
- `scipy`: Scientific computing
- `numpy`, `pandas`: Data handling
- `matplotlib`, `seaborn`: Visualization
- `pyyaml`: Configuration file handling


## Critical Parameters

**Preprocessing (configurable per dataset):**
- **TR**: Dataset-specific (HCP: 0.72s, Basel: 1.8s, London: 2.0s, Berlin: 3.4s, Missouri: 1.761s)
- **Volumes Removed**: First N only, no end-of-scan trim (Luppi 2024). HCP: 10 first; all others: 5 first; 0 last in all cases.
- **Bandpass**: 0.008-0.09 Hz
- **Spatial Smoothing**: None (ROI averaging only; Luppi et al. 2024 do not report voxel-level smoothing)
- **CompCor**: 5 WM + 5 CSF regenerated inside CONN (stage 1b) from the T1 `structural_segment` — aCompCor is NOT bridged from fMRIPrep confounds. Matches Luppi 2022 Nat Neurosci Methods exactly.
- **Global Signal Regression**: False

**HRF Deconvolution (MATLAB/rsHRF):**
- **Method**: Wiener deconvolution (rsHRF toolbox v2.2, Wu et al. 2013)
- **Iterations**: 50 (rsHRF default)
- **HRF Basis**: Gamma functions, order 3
- **HRF Length**: 24 seconds
- **Event Threshold**: 1 SD
- **Smoothing/Low-pass**: rsHRF defaults (auto from TR)

**Analysis:**
- **Time Delay (τ)**: 1 TR for TDMI calculations
- **Atlas**: Schaefer 400 + Tian 54 (454 ROIs), loaded from `atlases/schaefer_400_tian_54_combined.nii.gz`
- **PhiID Method**: Gaussian MMI-PID (phyid library)
- **Primary Synergy**: `sts` atom (persistent synergy), matching Luppi supplementary code
- **Group Averaging**: **Arithmetic mean** across subjects (`np.nanmean(..., axis=0)`) — `scripts/analysis/05_group_average_workspace.py:131-134`. NOT Fisher z, NOT median. PhiID atoms are MI-based information-theoretic quantities (bounded below by 0, not correlations); Fisher z / atanh are theoretically inappropriate here. Arithmetic mean is the Luppi-compliant choice. Do NOT "improve" this by adding a z-transform.
- **Zero-variance ROI handling**: Pairs involving a zero-variance ROI are stored as **NaN** (not 0) throughout the pipeline, because synergy/redundancy are typically positive and a 0 would pull node-strength rankings toward the low end and bias gateway/broadcaster classification. The MATLAB driver initialises all 5 output matrices with `nan(n_rois)` at `scripts/analysis/step04_compute_phiid_matlab.m:132-137` and leaves skipped pairs as NaN. All downstream Python aggregation is NaN-aware: subject-level `np.nanmean`, node strength via `np.nansum` with all-NaN rows → NaN strength, and workspace/PC ranks use `scipy.stats.rankdata(nan_policy='omit')` so excluded ROIs propagate as NaN rather than getting a spurious rank. Zero-variance ROIs should be rare (<5% of 454) in properly preprocessed data — the step04 `_summary.yaml` records `n_zero_variance_rois` per scan so higher rates flag signal dropout or registration failures.
- **Workspace Source**: **HCP-only**. The synergistic workspace is defined exclusively from HCP (100 healthy controls, produced by `05_group_average_workspace.py --dataset hcp` → `results/group_analysis/hcp/rest/hcp_rest_workspace_mask.npy`). All drug-dataset analyses project this HCP mask onto their data via `--workspace_dataset hcp` (default in `06_identify_gateways_broadcasters.py:811`). Drug datasets NEVER redefine the workspace — doing so would make drug-vs-placebo comparisons circular. The only exception is `--pooled` mode, which uses a pooled-baseline workspace and must be opted into explicitly.
- **Workspace → gateway/broadcaster operationalization (two-step, binary sign split)**: Exact Luppi et al. (2024) Methods implementation. Step 1 and Step 2 use different node metrics and the split is a pure sign-of-rank-difference test — no percentile thresholds anywhere.
  - **Step 1 — Workspace identification** (`05_group_average_workspace.py`): For each of the 454 ROIs, compute **node strength** (sum of edges) separately on the group-average synergy and redundancy matrices (`compute_node_strength`, lines 45-61). Rank the 454 strength values for each measure with `scipy.stats.rankdata(..., nan_policy='omit')`, gradient = `synergy_rank − redundancy_rank`, workspace = `gradient > 0` (`identify_workspace`, lines 64-103; sign test at lines 99-101). Binary, no threshold.
  - **Step 2 — Gateway vs broadcaster** (`06_identify_gateways_broadcasters.py`): Within workspace ROIs only, compute **participation coefficient** (Rubinov & Sporns 2010: P_i = 1 − Σ_s (κ_is / k_i)²) on the group-average synergy and redundancy **connectivity matrices** — NOT node strengths (`compute_participation_coefficient`, lines 95-150). The within-workspace gate is at line 184: `valid_mask = workspace_mask & np.isfinite(synergy_pc) & np.isfinite(redundancy_pc)`; PC ranks are taken over `synergy_pc[valid_mask]` / `redundancy_pc[valid_mask]` (NOT the full 454-ROI vector — brain-wide ranking would mix non-workspace ROIs into the distribution and bias the split). Gateway = `pc_gradient > 0`, broadcaster = `pc_gradient < 0` (`classify_gateways_broadcasters`, lines 152-211). Binary, no threshold.
  - **Module partition**: Exactly **8 modules** — 7 Yeo networks (Visual, Somatomotor, DorsAttn, SalVentAttn, Limbic, Cont, Default) + 1 Subcortical module (Tian S4). Source: `atlases/roi_to_network_mapping.csv` loaded via `scripts/utils/load_roi_mapping.py:54-70` (drops the background row, returns `Network_ID ∈ {0..7}`). The PC denominator in `compute_participation_coefficient` iterates over `np.unique(module_labels)` so module count is derived from the CSV, not hardcoded.
  - **Do NOT change to**: brain-wide PC ranking (breaks within-workspace gate), node-strength-only classification (collapses step 2 into step 1), any percentile threshold (breaks the binary sign split), or a different partition granularity such as 17 Yeo networks (Luppi uses the 7+1=8 partition).
- **Significance Threshold**: F>9 for NBS
- **Permutations**: 5,000 for statistical testing

## Expected Outcomes

### Primary Analyses
1. **Synergistic Workspace Changes**: How drugs affect the brain's synergistic core
2. **DMN Gateway Disruption**: Changes in information gathering
3. **FPN Broadcasting Alterations**: Effects on information dissemination
4. **Integrated Information**: Overall ΦR changes (consciousness-related measure)

### Key Comparisons
- Drug vs Control (within-subject, across datasets)
- Network-specific effects (DMN, FPN, sensory networks)
- Cross-dataset comparisons (future enhancement)
- All comparisons use HCP-defined workspace (`--workspace_dataset hcp` default)

## Computational Considerations
- **Storage**: ~10GB per subject for intermediate files
- **Memory Requirements**: Peak ~32GB during ΦID computation, available memory is 150 gb
- **Critical Bottleneck**: Pairwise TDMI calculations (454×454 ROI pairs)


## Version Control
- Use Git for all scripts
- Exclude large data files (`.gitignore`)
- Track analysis parameters in config files
- Document any deviations from standard pipeline
## Troubleshooting

**Dataset not found:**
```bash
python -c "from scripts.utils.dataset_handler import load_dataset_handler;
           handler = load_dataset_handler('london_lsd'); print(handler)"
```

**No scans discovered:**
- Check fmriprep_dir path in dataset config
- Verify BIDS naming conventions
- Check has_sessions setting

**Wrong run selected:**
- Review run_selection strategy in dataset config
- Use `--discover` to see what runs exist
