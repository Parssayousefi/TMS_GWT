%% Stage 1b CAMBRIDGE-PROPOFOL VARIANT — CONN denoise WITH ART.
%
% Verbatim fork of step01b_conn_denoise.m, hardcoded to cambridge_propofol
% only, with two additions:
%   1. batch.Setup.preprocessing.steps adds 'functional_art' (CONN-bundled
%      Artifact Detection Tool) with Luppi 2024 / DOC-cohort conservative
%      thresholds: FD > 0.5 mm, global-signal Z > 3.
%   2. batch.Denoising.confounds adds 'scrubbing' (per-volume outlier
%      indicator regressors auto-emitted by CONN's ART step). Outliers are
%      REGRESSED OUT, not censored — preserves volume count for PhiID.
%
% Every other parameter (CONN V25.b, aCompCor 5+5, motion + 6 derivatives,
% bandpass 0.008-0.09 Hz, RegBP, no GSR, no smoothing, no despiking, single
% 'rest' condition, two-pass Setup -> Denoising, output index 2 -> d*.nii)
% is bit-identical to step01b_conn_denoise.m.
%
% Output OVERWRITES the existing no-ART cambridge results at
%   results/timeseries_intermediate/cambridge_propofol/
% The CONN project file keeps an _art suffix
% (conn_cambridge_propofol_art.mat) so anyone inspecting the project after
% the fact knows ART was applied. Downstream stages (01c, step02, step04,
% 05, 06) need NO changes — they keep using --dataset cambridge_propofol.
%
% Prerequisites:
%   - CONN V25.b at C:\Users\YousefiP\software\conn25b\conn\
%   - SPM12 at D:\projects\TMS_GWT\pipeline\tools\spm12\spm\
%   - data/conn_input/cambridge_propofol/manifest.json from
%     01a_prepare_conn_inputs.py (unchanged from non-ART run)
%
% Usage:
%   >> cd('D:\projects\TMS_GWT\pipeline')
%   >> run('scripts/preprocessing/step01b_conn_denoise_cambridge_art.m')

%% ============================================================
%% USER CONFIG
%% ============================================================
DATASET = 'cambridge_propofol';   % HARDCODED — this script is propofol-only
assert(strcmp(DATASET, 'cambridge_propofol'), ...
       'This script is propofol-only. Use step01b_conn_denoise.m for other datasets.');

% Number of parallel CONN workers (per-dataset, scan-level parallelism via
% CONN's batch.parallel feature). Each worker is an independent MATLAB+CONN
% process spawned by 'Background process' profile — no PCT license needed.
%   0   = serial (original behaviour)
%   >0  = parallel; CONN chunks scans across N workers via conn_jobmanager.
% Per-worker RAM peak ~8-10 GB (SPM segmentation + aCompCor). On a 150 GB
% Windows workstation, 12 is aggressive-but-safe; drop to 8 under memory
% pressure. CONN auto-caps N to n_scans, so small datasets are fine.
N_WORKERS = 15;

% Skip scans whose final output ({scan_id}_conn_cleaned.nii) already exists
% in results/timeseries_intermediate/cambridge_propofol/. Set to false to
% force full re-run (e.g. after parameter changes).
SKIP_DONE = false;

%% ============================================================
%% PATH SETUP
%% ============================================================
project_dir = 'D:\projects\TMS_GWT\pipeline';
conn_dir    = 'C:\Users\YousefiP\software\conn25b\conn';
spm_dir     = fullfile(project_dir, 'tools', 'spm12', 'spm');

if exist(spm_dir, 'dir')
    addpath(spm_dir);
else
    error('SPM12 not found at %s', spm_dir);
end
if exist(conn_dir, 'dir')
    addpath(conn_dir);
else
    error('CONN not found at %s', conn_dir);
end

%% ============================================================
%% PARALLEL PROFILE (OS-dependent)
%% ============================================================
if ispc
    parallel_profile = 'Background process (Windows)';
else
    parallel_profile = 'Background process (Unix,Mac)';
end
if N_WORKERS > 0
    fprintf('Parallel mode: N_WORKERS=%d, profile=''%s'' (projected peak RAM ~%d GB)\n', ...
            N_WORKERS, parallel_profile, N_WORKERS * 9);
else
    fprintf('Parallel mode: SERIAL (N_WORKERS=0)\n');
end

%% ============================================================
%% RESOLVE DATASET (single, hardcoded — no auto-discover)
%% ============================================================
conn_input_root = fullfile(project_dir, 'data', 'conn_input');
datasets = {DATASET};
fprintf('============================================================\n');
fprintf('CONN denoising — single dataset (ART variant): %s\n', DATASET);
fprintf('CONN: %s\n', conn_dir);
fprintf('SPM:  %s\n', spm_dir);
fprintf('ART thresholds: FD > 0.5 mm, GS Z > 3 (Luppi 2024 DOC convention)\n');
fprintf('============================================================\n');

per_dataset_results = cell(numel(datasets), 4);  % {dataset, n_found, n_total, n_skipped}

for ds_idx = 1:numel(datasets)
    current_dataset = datasets{ds_idx};
    fprintf('\n');
    fprintf('############################################################\n');
    fprintf('# DATASET %d/%d: %s\n', ds_idx, numel(datasets), current_dataset);
    fprintf('############################################################\n');

%% ============================================================
%% LOAD MANIFEST
%% ============================================================
manifest_path = fullfile(project_dir, 'data', 'conn_input', current_dataset, 'manifest.json');
if ~exist(manifest_path, 'file')
    fprintf('  SKIP: manifest.json not found at %s\n', manifest_path);
    per_dataset_results{ds_idx, 1} = current_dataset;
    per_dataset_results{ds_idx, 2} = 0;
    per_dataset_results{ds_idx, 3} = 0;
    per_dataset_results{ds_idx, 4} = 0;
    continue;
end
fprintf('Loading manifest: %s\n', manifest_path);
manifest = jsondecode(fileread(manifest_path));

if ~isfield(manifest, 'scans') || isempty(manifest.scans)
    fprintf('  SKIP: manifest.json contains no scans\n');
    per_dataset_results{ds_idx, 1} = current_dataset;
    per_dataset_results{ds_idx, 2} = 0;
    per_dataset_results{ds_idx, 3} = 0;
    per_dataset_results{ds_idx, 4} = 0;
    continue;
end

scans_struct = manifest.scans;
if ~isstruct(scans_struct)
    fprintf('  SKIP: manifest.scans is not a struct array (check JSON)\n');
    per_dataset_results{ds_idx, 1} = current_dataset;
    per_dataset_results{ds_idx, 2} = 0;
    per_dataset_results{ds_idx, 3} = 0;
    per_dataset_results{ds_idx, 4} = 0;
    continue;
end
n_scans = numel(scans_struct);
fprintf('Scans in manifest: %d\n', n_scans);
for s = 1:n_scans
    fprintf('  [%d] %s  (%d volumes, TR=%.3fs)\n', s, scans_struct(s).scan_id, ...
            scans_struct(s).n_volumes, scans_struct(s).tr);
end

% TRs must be uniform within a dataset (asserted by the dataset config). We
% use the first scan's TR as the project-level RT.
unique_trs = unique([scans_struct.tr]);
if numel(unique_trs) > 1
    warning('Manifest contains mixed TRs %s — using first scan TR (%.3fs) for the CONN project.', ...
            mat2str(unique_trs), scans_struct(1).tr);
end

%% ============================================================
%% PREPARE OUTPUT DIRECTORIES
%% ============================================================
% ART variant OVERWRITES the existing non-ART cambridge results. Project
% filename keeps an _art tag for provenance, but the cleaned 4D NIfTIs land
% under the canonical cambridge_propofol path so 01c (and all downstream)
% works unchanged.
intermediate_root = fullfile(project_dir, 'results', 'timeseries_intermediate', ...
                             'cambridge_propofol');
project_subdir = fullfile(intermediate_root, 'conn_project');
cleaned_dir = intermediate_root;
if ~exist(project_subdir, 'dir'), mkdir(project_subdir); end
if ~exist(cleaned_dir, 'dir'), mkdir(cleaned_dir); end

project_filename = fullfile(project_subdir, sprintf('conn_%s_art.mat', current_dataset));
fprintf('CONN project file: %s\n', project_filename);

%% ============================================================
%% SKIP_DONE: FILTER OUT ALREADY-COMPLETED SCANS
%% ============================================================
n_manifest = n_scans;
n_skipped  = 0;
if SKIP_DONE
    keep_mask = true(1, n_manifest);
    for s = 1:n_manifest
        done_path = fullfile(cleaned_dir, sprintf('%s_conn_cleaned.nii', scans_struct(s).scan_id));
        if exist(done_path, 'file')
            fprintf('  SKIP_DONE: %s (output exists)\n', scans_struct(s).scan_id);
            keep_mask(s) = false;
        end
    end
    n_skipped = sum(~keep_mask);
    if n_skipped > 0
        scans_struct = scans_struct(keep_mask);
        n_scans = numel(scans_struct);
        fprintf('Filtered: %d scans to process (%d skipped as already done)\n', n_scans, n_skipped);
    end
    if n_scans == 0
        fprintf('\n');
        fprintf('------------------------------------------------------------\n');
        fprintf('Dataset %s: ALL %d scans already done — skipping CONN batch\n', current_dataset, n_manifest);
        fprintf('------------------------------------------------------------\n');
        per_dataset_results{ds_idx, 1} = current_dataset;
        per_dataset_results{ds_idx, 2} = n_manifest;
        per_dataset_results{ds_idx, 3} = n_manifest;
        per_dataset_results{ds_idx, 4} = n_skipped;
        continue;
    end
end

% Fresh project per run (matches batch.Setup.isnew=1).
if exist(project_filename, 'file')
    fprintf('Removing existing CONN project: %s\n', project_filename);
    delete(project_filename);
    [proj_dir, proj_base, ~] = fileparts(project_filename);
    sibling_dir = fullfile(proj_dir, proj_base);
    if exist(sibling_dir, 'dir')
        fprintf('Removing existing CONN project dir: %s\n', sibling_dir);
        rmdir(sibling_dir, 's');
    end
end

%% ============================================================
%% BUILD CONN BATCH STRUCTURE
%% ============================================================
% Design choice: each scan = one CONN "subject" with one session. This
% duplicates SPM segment work for multi-scan-per-subject datasets (e.g. one T1
% segmented N times for N scans) but keeps the input/output file naming clean
% and the failure mode per-scan rather than per-subject. The wasted compute
% (~20 s per duplicate scan) is a small fraction of CONN's per-scan cost (~8
% min) and not worth optimising for the first production version.
clear batch;

batch.filename = project_filename;

batch.Setup.isnew     = 1;
batch.Setup.nsubjects = n_scans;
batch.Setup.RT        = scans_struct(1).tr;
batch.Setup.acquisitiontype = 1;                 % continuous
batch.Setup.voxelresolution = 1;                 % 2 mm isotropic analysis space (matches CONN HCP example)

batch.Setup.functionals = repmat({{}}, [n_scans, 1]);
batch.Setup.structurals = cell(n_scans, 1);
for s = 1:n_scans
    batch.Setup.functionals{s}{1} = {scans_struct(s).bold_trimmed};
    batch.Setup.structurals{s}    = scans_struct(s).T1w;
end

% Single 'rest' condition spanning each scan in full (paper-3 HCP convention;
% applied uniformly across all datasets — we treat each drug condition as a
% rest-style continuous block for denoising purposes).
batch.Setup.conditions.names = {'rest'};
for s = 1:n_scans
    batch.Setup.conditions.onsets{1}{s}{1}    = 0;
    batch.Setup.conditions.durations{1}{s}{1} = inf;
end

% Realignment covariate from fMRIPrep motion params (so CONN doesn't recompute)
batch.Setup.covariates.names = {'realignment'};
batch.Setup.covariates.files{1} = repmat({{}}, [n_scans, 1]);
for s = 1:n_scans
    batch.Setup.covariates.files{1}{s}{1} = {scans_struct(s).motion_txt};
end

% Output index 2 -> "d*.nii" prefix files (confound-corrected 4D)
batch.Setup.outputfiles = [0, 1, 0, 0, 0, 0];

% analyses=[1,2] is required to trigger ROI extraction (Step 5 of CONN Setup).
% We do not consume the seed-to-voxel / ROI-to-ROI connectivity outputs;
% setting analyses=[] would skip Step 5 entirely and break Denoising.
batch.Setup.analyses  = [1, 2];
batch.Setup.overwrite = 'Yes';
batch.Setup.done      = 1;

% PROPOFOL/DOC METHODOLOGY (Luppi 2024, page 17/18):
%   - structural_segment: CONN segments T1 -> GM/WM/CSF for aCompCor
%   - functional_art: identify outlier volumes (FD > 0.5 mm OR global-signal
%     |Z| > 3 — CONN's "conservative" preset, matching Luppi's DOC/anaesthesia
%     description). ART writes per-scan outlier indicator regressors that
%     CONN auto-imports as a 'scrubbing' covariate; the Denoising pass below
%     regresses them out (it does NOT remove volumes — preserves volume count
%     for PhiID downstream).
% IMPORTANT: art_thresholds(1) = global-signal Z, art_thresholds(2) = FD mm
% (see CONN conn_batch.m:343-363). [3 0.5] is the "conservative" preset
% matching Luppi 2024 anaesthesia/DOC convention. DO NOT reverse.
batch.Setup.preprocessing.steps          = {'structural_segment', 'functional_art'};
batch.Setup.preprocessing.art_thresholds = [3, 0.5];   % [GS_Z, FD_mm] — conservative preset

%% ============================================================
%% RUN CONN — PASS 1: SETUP (segmentation + ROI extraction + ART)
%% ============================================================
% CONN's official HCP example splits Setup and Denoising into two separate
% conn_batch() calls. Combining them into one batch causes the
% Setup.Denoising.update step to fire before ROI extraction has completed,
% producing a "Not ready to process step conn_process_5" error.
fprintf('\n============================================================\n');
fprintf('Calling conn_batch — PASS 1: Setup (segmentation + ROI + ART)\n');
fprintf('Per-scan wall time: ~10-30 minutes (segmentation + import + ART)\n');
fprintf('============================================================\n');

if N_WORKERS > 0
    batch.parallel.N               = N_WORKERS;
    batch.parallel.profile         = parallel_profile;
    batch.parallel.immediatereturn = 0;   % block until all workers done before PASS 2
    fprintf('  PASS 1 parallel: N=%d (%d scans chunked across workers)\n', N_WORKERS, n_scans);
end

tic_total = tic;
conn_batch(batch);

%% ============================================================
%% RUN CONN — PASS 2: DENOISING
%% ============================================================
clear batch;
batch.filename = project_filename;

batch.Denoising.filter      = [0.008, 0.09];  % paper-3 / paper-5 bandpass
batch.Denoising.detrending  = 1;              % linear
batch.Denoising.despiking   = 0;              % no despiking
batch.Denoising.regbp       = 1;              % RegBP: regression then bandpass

% Confound matrix: 5 WM aCompCor + 5 CSF aCompCor + (6 motion + 6 derivatives)
%                  + scrubbing (one binary indicator per ART-flagged volume).
% scrubbing: dimensions=inf -> regress every outlier indicator column;
% deriv=0 because outlier indicators are binary spikes (derivatives would
% smear the spike across two TRs). Per Luppi 2024 page 18 the artefacts
% identified by ART are listed as confound regressors, not as removed volumes
% — preserves volume count for PhiID array shapes.
batch.Denoising.confounds.names      = {'White Matter', 'CSF', 'realignment', 'scrubbing'};
batch.Denoising.confounds.dimensions = {5, 5, 6, inf};
batch.Denoising.confounds.deriv      = {0, 0, 1, 0};   % motion: raw + 1st deriv -> 6+6=12 cols
batch.Denoising.confounds.power      = {1, 1, 1, 1};

batch.Denoising.overwrite = 'Yes';
batch.Denoising.done      = 1;

fprintf('\n============================================================\n');
fprintf('Calling conn_batch — PASS 2: Denoising (with scrubbing regressors)\n');
fprintf('============================================================\n');

if N_WORKERS > 0
    batch.parallel.N               = N_WORKERS;
    batch.parallel.profile         = parallel_profile;
    batch.parallel.immediatereturn = 0;   % block until denoising done before copyfile stage
    fprintf('  PASS 2 parallel: N=%d\n', N_WORKERS);
end

conn_batch(batch);
elapsed_total = toc(tic_total);

fprintf('\n============================================================\n');
fprintf('CONN batch finished in %.1f minutes\n', elapsed_total / 60);
fprintf('============================================================\n');

%% ============================================================
%% DISCOVER AND COPY DENOISED 4D OUTPUTS
%% ============================================================
% CONN writes confound-corrected timeseries with prefix 'd' next to the source
% functional file. So `dbold_trimmed.nii` should appear in each scan's
% data/conn_input/{dataset}/{scan_id}/ directory after CONN finishes. We copy
% these to results/timeseries_intermediate/cambridge_propofol/{scan_id}_conn_cleaned.nii
% for the parcellate stage.

fprintf('\nSearching for CONN denoised outputs...\n');
n_found = 0;
for s = 1:n_scans
    scan_id = scans_struct(s).scan_id;
    bold_path = scans_struct(s).bold_trimmed;
    [bold_dir, bold_base, bold_ext] = fileparts(bold_path);

    candidates = {
        fullfile(bold_dir, ['d' bold_base bold_ext]), ...
        fullfile(bold_dir, 'dbold_trimmed.nii'), ...
    };

    found = '';
    for c = 1:numel(candidates)
        if exist(candidates{c}, 'file')
            found = candidates{c};
            break;
        end
    end

    % Fallback: scan the conn project directory for niftiDATA files
    if isempty(found)
        glob_pattern = fullfile(project_subdir, '**', sprintf('*Subject%03d*.nii', s));
        d = dir(glob_pattern);
        if ~isempty(d)
            [~, idx] = max([d.datenum]);
            found = fullfile(d(idx).folder, d(idx).name);
        end
    end

    if ~isempty(found)
        dst = fullfile(cleaned_dir, sprintf('%s_conn_cleaned.nii', scan_id));
        copyfile(found, dst);
        info = dir(dst);
        fprintf('  %s: copied %s -> %s (%.1f MB)\n', ...
                scan_id, found, dst, info.bytes / 1e6);
        n_found = n_found + 1;
    else
        fprintf('  %s: NO denoised output found! (looked in %s)\n', scan_id, bold_dir);
    end
end

%% ============================================================
%% PER-SCAN ART OUTLIER SUMMARY (CSV for Methods/Supplementary)
%% ============================================================
% CONN's ART step writes art_regression_outliers_*.mat per scan into the
% scan's INPUT directory (data/conn_input/{dataset}/{scan_id}/), NOT into
% the conn project's preprocessing dir. Each .mat contains an R matrix with
% one column per flagged outlier volume; we tabulate the count and outlier
% fraction so the supplementary table can quote per-scan censoring
% intensity. IMPORTANT: 01c auto-cleans the conn_input scan dirs after
% parcellation, so this CSV must be written by 01b before 01c runs.
fprintf('\nReading ART outlier summaries...\n');
art_summary = cell(n_scans, 3);   % {scan_id, n_outliers, fraction_outliers}
for s = 1:n_scans
    scan_id = scans_struct(s).scan_id;
    art_summary{s, 1} = scan_id;
    art_summary{s, 2} = NaN;
    art_summary{s, 3} = NaN;
    bold_path = scans_struct(s).bold_trimmed;
    [bold_dir, ~, ~] = fileparts(bold_path);
    art_files = dir(fullfile(bold_dir, 'art_regression_outliers_*.mat'));
    if isempty(art_files)
        continue;
    end
    art_data = load(fullfile(art_files(1).folder, art_files(1).name));
    if isfield(art_data, 'R')
        n_out = size(art_data.R, 2);
    else
        n_out = 0;
    end
    art_summary{s, 2} = n_out;
    art_summary{s, 3} = n_out / scans_struct(s).n_volumes;
    fprintf('  %s: %d outlier volumes (%.1f%% of %d)\n', ...
            scan_id, n_out, 100 * art_summary{s, 3}, scans_struct(s).n_volumes);
end

csv_path = fullfile(intermediate_root, 'art_outlier_summary.csv');
fid = fopen(csv_path, 'w');
fprintf(fid, 'scan_id,n_outliers,fraction_outliers\n');
for s = 1:n_scans
    n_out = art_summary{s, 2};
    frac  = art_summary{s, 3};
    if isnan(n_out)
        fprintf(fid, '%s,,\n', art_summary{s, 1});
    else
        fprintf(fid, '%s,%d,%.4f\n', art_summary{s, 1}, n_out, frac);
    end
end
fclose(fid);
fprintf('ART outlier summary written to: %s\n', csv_path);

fprintf('\n');
fprintf('------------------------------------------------------------\n');
fprintf('Dataset %s done: %d/%d scans denoised by CONN (%d previously done, skipped)\n', ...
        current_dataset, n_found, n_scans, n_skipped);
fprintf('Cleaned 4D BOLDs at: %s\n', cleaned_dir);
fprintf('------------------------------------------------------------\n');

per_dataset_results{ds_idx, 1} = current_dataset;
per_dataset_results{ds_idx, 2} = n_found + n_skipped;
per_dataset_results{ds_idx, 3} = n_manifest;
per_dataset_results{ds_idx, 4} = n_skipped;

end  % end of dataset loop

%% ============================================================
%% FINAL SUMMARY
%% ============================================================
fprintf('\n');
fprintf('============================================================\n');
fprintf('Stage 1b (cambridge ART variant) complete\n');
total_found = 0;
total_scans = 0;
total_skipped = 0;
for k = 1:size(per_dataset_results, 1)
    ds = per_dataset_results{k, 1};
    f  = per_dataset_results{k, 2};
    n  = per_dataset_results{k, 3};
    sk = per_dataset_results{k, 4};
    if isempty(ds), ds = '<unknown>'; end
    if isempty(f), f = 0; end
    if isempty(n), n = 0; end
    if isempty(sk), sk = 0; end
    if sk > 0
        fprintf('  %-20s  %d/%d scans denoised (%d skipped)\n', ds, f, n, sk);
    else
        fprintf('  %-20s  %d/%d scans denoised\n', ds, f, n);
    end
    total_found   = total_found + f;
    total_scans   = total_scans + n;
    total_skipped = total_skipped + sk;
end
if total_skipped > 0
    fprintf('  %-20s  %d/%d scans denoised (%d skipped)\n', 'TOTAL', total_found, total_scans, total_skipped);
else
    fprintf('  %-20s  %d/%d scans denoised\n', 'TOTAL', total_found, total_scans);
end
fprintf('============================================================\n');
fprintf('\nNext step:\n');
fprintf('  python scripts/preprocessing/01c_parcellate_conn_output.py --dataset cambridge_propofol\n');
