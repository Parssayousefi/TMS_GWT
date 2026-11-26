%% HRF Deconvolution using rsHRF toolbox (Wu et al. 2013)
%
% Processes all *_timeseries.npy files -> *_timeseries_deconv.npy
% Follows Luppi et al. 2024 methodology (rsHRF toolbox v2.5)
%
% 3-stage pipeline:
%   Stage 1: Event detection (threshold crossing at 1 SD)
%   Stage 2: HRF estimation (Gamma basis functions, order 3)
%   Stage 3: Iterative Wiener deconvolution (50 iterations, default smoothing & low-pass)
%
% All math performed by rsHRF toolbox functions.
%
% Prerequisites:
%   - SPM12 in tools/spm12/
%   - npy-matlab in tools/npy-matlab/
%   - rsHRF toolbox in tools/rsHRF/ (already present)
%   - Signal Processing Toolbox (resample)
%   - Wavelet Toolbox (wavedec, wnoisest)
%
% Usage:
%   >> run('scripts/preprocessing/02_hrf_deconvolution_matlab.m')

%% User config
% Parallel workers for outer scan-level parfor. Each worker is an independent
% MATLAB process loading SPM + rsHRF + Wavelet + Signal Processing toolboxes
% (~3-5 GB baseline) plus tiny per-scan data (~10 MB). On a 150 GB box, 12
% workers is a safe default. Requires Parallel Computing Toolbox license.
N_WORKERS = 24;

%% Setup paths
project_dir = 'D:\projects\TMS_GWT\pipeline';
addpath(fullfile(project_dir, 'tools', 'rsHRF'));
addpath(fullfile(project_dir, 'tools', 'spm12', 'spm'));
addpath(genpath(fullfile(project_dir, 'tools', 'npy-matlab')));

ts_dir = fullfile(project_dir, 'results', 'timeseries');

% TR lookup. missouri_psilo added 2026-08-26: the old exclusion note
% ("insufficient CompCor components") referred to the deprecated legacy
% fMRIPrep-aCompCor extraction path and never applied to this CONN chain,
% which regenerates aCompCor inside CONN.
% Full list restored 2026-08-10: a cambridge-only map silently skipped every
% other dataset. Already-deconvolved scans are skipped by the exist() guard
% below, so listing them here cannot recompute anything.
% sudmex_tms and aachen_itbs added 2026-11-05 in TMS_GWT (TR 2.0 s from their
% BIDS sidecars); see PROVENANCE/POST_COPY_EDITS.tsv.
dataset_tr = containers.Map( ...
    {'hcp', 'basel_lam', 'basel_lpm', 'london_lsd', 'berlin_ketamine', ...
     'cambridge_propofol', 'cambridge_propofol_uniform', 'monash_psilo', ...
     'missouri_psilo', 'sudmex_tms', 'aachen_itbs'}, ...
    {0.72, 1.8, 1.8, 2.0, 3.4, 2.0, 2.0, 0.91, 1.761, 2.0, 2.0});

%% Discover files
all_files = dir(fullfile(ts_dir, 'dataset-*_timeseries.npy'));

% Exclude _deconv files (missouri_psilo filter removed 2026-08-26 - dataset
% activated; see TR-map comment above)
mask = ~contains({all_files.name}, '_deconv');
all_files = all_files(mask);
n_files = length(all_files);

fprintf('============================================================\n');
fprintf('HRF Deconvolution — rsHRF toolbox (Wu et al. 2013)\n');
fprintf('============================================================\n');
fprintf('Found %d timeseries files\n', n_files);
fprintf('Output directory: %s\n', ts_dir);
fprintf('Workers: %d\n', N_WORKERS);
fprintf('============================================================\n\n');

%% Pre-resolve per-file TRs (avoids containers.Map broadcast in parfor) and
%% filter out scans whose deconvolved output already exists. Any filename we
%% can't parse or map to a known dataset is skipped here rather than inside
%% the parallel loop — keeps the pool fed with only valid work.
todo_names   = cell(n_files, 1);
todo_folders = cell(n_files, 1);
todo_TRs     = zeros(n_files, 1);
n_todo = 0;
for f = 1:n_files
    fname = all_files(f).name;
    out_name = strrep(fname, '_timeseries.npy', '_timeseries_deconv.npy');
    out_path = fullfile(ts_dir, out_name);
    if exist(out_path, 'file')
        fprintf('SKIP: %s (already computed)\n', fname);
        continue;
    end
    tok = regexp(fname, 'dataset-(\w+?)_sub-', 'tokens');
    if isempty(tok)
        fprintf('SKIP: %s (cannot parse dataset name)\n', fname);
        continue;
    end
    ds = tok{1}{1};
    if ~isKey(dataset_tr, ds)
        fprintf('SKIP: %s (unknown dataset: %s)\n', fname, ds);
        continue;
    end
    n_todo = n_todo + 1;
    todo_names{n_todo}   = fname;
    todo_folders{n_todo} = all_files(f).folder;
    todo_TRs(n_todo)     = dataset_tr(ds);
end
% Trim pre-allocated cells to actual size
todo_names   = todo_names(1:n_todo);
todo_folders = todo_folders(1:n_todo);
todo_TRs     = todo_TRs(1:n_todo);
n_pre_skip   = n_files - n_todo;

fprintf('\n%d files to process, %d already done / skipped\n\n', n_todo, n_pre_skip);

if n_todo == 0
    fprintf('Nothing to do.\n');
    return;
end

%% Start parallel pool (addpath above already ran so workers inherit paths)
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local', N_WORKERS);
elseif pool.NumWorkers ~= N_WORKERS
    delete(pool);
    pool = parpool('local', N_WORKERS);
end
fprintf('Parallel pool ready: %d workers\n\n', pool.NumWorkers);

%% Process files in parallel
% Per-iteration results: 0=success, 2=error (pre-filter already handled skips)
statuses = zeros(n_todo, 1);
errors   = cell(n_todo, 1);
qc_rows  = cell(n_todo, 1);
tic_total = tic;

parfor k = 1:n_todo
    fname      = todo_names{k};
    input_path = fullfile(todo_folders{k}, fname);
    TR         = todo_TRs(k);

    % Output path: replace _timeseries.npy with _timeseries_deconv.npy
    out_name = strrep(fname, '_timeseries.npy', '_timeseries_deconv.npy');
    out_path = fullfile(ts_dir, out_name);

    try
        % Load timeseries (.npy format)
        bold_sig = readNPY(input_path);  % (timepoints x n_rois)
        [n_tp, n_rois] = size(bold_sig);
        fprintf('[%d/%d] Processing: %s\n', k, n_todo, fname);
        fprintf('  TR=%.3fs, %d timepoints, %d ROIs\n', TR, n_tp, n_rois);

        tic_scan = tic;

        % ---- Preprocessing (as per rsHRF demo) ----
        % Capture stats before mean-centering (for metadata)
        mean_original = mean(bold_sig(:));
        std_original = std(bold_sig(:));
        % Mean-center (bandpass already applied in step 01)
        bold_sig = bold_sig - mean(bold_sig);
        sigma = std(bold_sig);

        % Detect zero-variance ROIs (e.g. signal dropout)
        zero_var_mask = sigma < 1e-10;
        n_zero_var = sum(zero_var_mask);
        if n_zero_var > 0
            fprintf('  WARNING: %d zero-variance ROIs detected, will be set to zero in output\n', n_zero_var);
        end

        % ---- rsHRF parameters ----
        % Built in one struct() call so parfor classifies it as temporary.
        dt_val = TR / 3;
        para = struct( ...
            'TR', TR, ...
            'name', 'Gamma functions', ...
            'order', 3, ...
            'T', 3, ...               % 3x temporal upsampling
            'T0', 1, ...              % reference time bin
            'dt', dt_val, ...         % fine-scale time resolution
            'AR_lag', 1, ...          % AR(1) noise autocorrelation
            'thr', 1, ...             % 1 SD event detection threshold
            'len', 24, ...            % HRF length in seconds
            'localK', 2, ...          % local peak detection window
            'lag', fix(4/dt_val):fix(8/dt_val) ...  % 4-8s onset search
        );

        % ---- Stage 1+2: HRF estimation (valid ROIs only) ----
        % rsHRF_estimation_temporal_basis handles:
        %   - Event detection (rsHRF_find_event_vector)
        %   - Basis function generation (spm_Gpdf for Gamma)
        %   - HRF fitting with lag search (GLS with AR(1))
        valid_idx = find(~zero_var_mask);
        [beta_hrf, bf, event_bold] = rsHRF_estimation_temporal_basis(...
            bold_sig(:, valid_idx), para, [], 0);  % flag_parfor=0

        % Reconstruct HRF waveforms from basis coefficients
        hrfa_valid = bf * beta_hrf(1:size(bf, 2), :);  % (microtime_samples x n_valid)

        % Map back to full ROI space
        hrfa = zeros(size(hrfa_valid, 1), n_rois);
        hrfa(:, valid_idx) = hrfa_valid;

        % Extract HRF parameters: Height, Time-to-Peak, FWHM
        PARA = zeros(3, n_rois);
        for hh = 1:length(valid_idx)
            roi_idx = valid_idx(hh);
            PARA(:, roi_idx) = rsHRF_get_HRF_parameters(hrfa(:, roi_idx), para.dt);
        end

        % ---- Resample HRF from microtime to TR grid ----
        if para.T > 1
            hrfa_TR = resample(hrfa, 1, para.T);
        else
            hrfa_TR = hrfa;
        end

        % ---- Stage 3: Iterative Wiener deconvolution ----
        % Following rsHRF_deconv_job.m lines 103-113:
        %   - Z-score the input signal
        %   - Normalize HRF by signal SD
        %   - Use default smoothing and low-pass filtering (auto from TR)
        %   - 50 iterations (toolbox default)
        zdata = zscore(bold_sig(:, valid_idx));
        data_deconv = zeros(size(bold_sig));

        for jj = 1:length(valid_idx)
            roi_idx = valid_idx(jj);
            hrf = hrfa_TR(:, roi_idx);
            data_deconv(:, roi_idx) = rsHRF_iterative_wiener_deconv(...
                zdata(:, jj), ...          % z-scored denoised signal
                hrf ./ sigma(roi_idx), ... % HRF normalized by signal SD
                'TR', para.TR, ...         % TR for smoothing/low-pass defaults
                'MaxIter', 50);            % toolbox default
        end

        elapsed_scan = toc(tic_scan);

        % ---- Save deconvolved timeseries as .npy ----
        writeNPY(data_deconv, out_path);

        % ---- Save HRF estimation info as .mat ----
        mat_name = strrep(fname, '_timeseries.npy', '_hrf.mat');
        mat_path = fullfile(ts_dir, mat_name);
        parsave_hrf(mat_path, para, beta_hrf, hrfa, PARA, event_bold);

        % ---- Save metadata YAML (compatible with downstream pipeline) ----
        yaml_name = strrep(fname, '_timeseries.npy', '_timeseries_deconv_metadata.yaml');
        yaml_path = fullfile(ts_dir, yaml_name);
        fid_yaml = fopen(yaml_path, 'w');
        fprintf(fid_yaml, 'method: rsHRF_Wu2013_gamma\n');
        fprintf(fid_yaml, 'toolbox: rsHRF_matlab_v2.5\n');
        fprintf(fid_yaml, 'basis_function: Gamma functions\n');
        fprintf(fid_yaml, 'basis_order: 3\n');
        fprintf(fid_yaml, 'tr: %.3f\n', TR);
        fprintf(fid_yaml, 'thr: 1\n');
        fprintf(fid_yaml, 'hrf_len_s: 24\n');
        fprintf(fid_yaml, 'onset_search_s: 4-8\n');
        fprintf(fid_yaml, 'wiener_max_iter: 50\n');
        fprintf(fid_yaml, 'deconv_internal_smoothing: default\n');
        fprintf(fid_yaml, 'deconv_internal_low_pass: default\n');
        fprintf(fid_yaml, 'original_shape:\n- %d\n- %d\n', n_tp, n_rois);
        fprintf(fid_yaml, 'deconv_shape:\n- %d\n- %d\n', size(data_deconv));
        fprintf(fid_yaml, 'mean_original: %.15e\n', mean_original);
        fprintf(fid_yaml, 'std_original: %.15e\n', std_original);
        fprintf(fid_yaml, 'mean_deconv: %.15e\n', mean(data_deconv(:)));
        fprintf(fid_yaml, 'std_deconv: %.15e\n', std(data_deconv(:)));
        fprintf(fid_yaml, 'mean_ttp_s: %.4f\n', mean(PARA(2,:)));
        fprintf(fid_yaml, 'mean_height: %.4f\n', mean(PARA(1,:)));
        fprintf(fid_yaml, 'mean_fwhm_s: %.4f\n', mean(PARA(3,:)));
        fprintf(fid_yaml, 'input_file: %s\n', input_path);
        fprintf(fid_yaml, 'output_file: %s\n', out_path);
        fprintf(fid_yaml, 'note: >-\n');
        fprintf(fid_yaml, '  HRF deconvolution via rsHRF MATLAB toolbox (Wu et al. 2013),\n');
        fprintf(fid_yaml, '  same toolbox cited by Luppi et al. 2024 (doi:10.7554/eLife.88173),\n');
        fprintf(fid_yaml, '  with default smoothing and low-pass filtering (auto from TR and rest mode).\n');
        fclose(fid_yaml);

        fprintf('  Saved: %s\n', out_name);
        fprintf('  Mean TTP=%.2fs, Mean Height=%.4f, Elapsed=%.1fs\n', ...
            mean(PARA(2, :)), mean(PARA(1, :)), elapsed_scan);

        % ---- Collect QC metrics for summary CSV ----
        % Parse subject/session/task from filename
        sub_tok = regexp(fname, 'sub-([^_]+)', 'tokens');
        ses_tok = regexp(fname, 'ses-([^_]+)', 'tokens');
        task_tok = regexp(fname, 'task-([^_]+)', 'tokens');
        sub_id = sub_tok{1}{1};
        if ~isempty(ses_tok); ses_id = ses_tok{1}{1}; else; ses_id = ''; end
        if ~isempty(task_tok); task_id = task_tok{1}{1}; else; task_id = ''; end

        % Event counts per ROI
        n_events = cellfun(@numel, event_bold);
        % % of ROIs with TTP in physiological range [3, 10]s
        pct_valid = 100 * sum(PARA(2,:) >= 3 & PARA(2,:) <= 10) / n_rois;

        % Re-parse dataset name for QC row (pre-resolved TR broadcast doesn't
        % carry the dataset label; cheap regex, local to the worker)
        ds_tok = regexp(fname, 'dataset-(\w+?)_sub-', 'tokens');
        ds = ds_tok{1}{1};

        qc_rows{k} = {fname, ds, sub_id, ses_id, task_id, TR, ...
            mean(PARA(2,:)), median(PARA(2,:)), std(PARA(2,:)), ...
            mean(PARA(1,:)), std(PARA(1,:)), ...
            mean(PARA(3,:)), std(PARA(3,:)), ...
            mean(n_events), pct_valid, elapsed_scan};

        statuses(k) = 0;  % success (redundant — default — but explicit)

    catch ME
        fprintf('  ERROR: %s\n    %s\n', fname, ME.message);
        errors{k}   = sprintf('%s: %s', fname, ME.message);
        statuses(k) = 2;
    end
end

%% Summary
elapsed_total = toc(tic_total);

% Aggregate from per-iter arrays and drop empty slots from failed iterations
n_success = sum(statuses == 0);
n_error   = sum(statuses == 2);
n_skip    = n_pre_skip;
errors    = errors(~cellfun('isempty', errors));
qc_rows   = qc_rows(~cellfun('isempty', qc_rows));

fprintf('\n============================================================\n');
fprintf('HRF DECONVOLUTION COMPLETE\n');
fprintf('============================================================\n');
fprintf('  Success: %d\n', n_success);
fprintf('  Skipped: %d (already computed / unparseable / unknown dataset)\n', n_skip);
fprintf('  Errors:  %d\n', n_error);
fprintf('  Total:   %d files\n', n_files);
fprintf('  Wall clock: %.1f minutes\n', elapsed_total / 60);

if ~isempty(errors)
    fprintf('\nError details:\n');
    for e = 1:length(errors)
        fprintf('  - %s\n', errors{e});
    end
end

%% Write QC summary CSV
if ~isempty(qc_rows)
    qc_path = fullfile(ts_dir, 'hrf_deconvolution_quality.csv');
    fid = fopen(qc_path, 'w');
    fprintf(fid, 'filename,dataset,subject,session,task,TR,');
    fprintf(fid, 'mean_ttp_s,median_ttp_s,sd_ttp_s,');
    fprintf(fid, 'mean_height,sd_height,');
    fprintf(fid, 'mean_fwhm_s,sd_fwhm_s,');
    fprintf(fid, 'mean_n_events,pct_valid_hrf,elapsed_s\n');
    for r = 1:length(qc_rows)
        row = qc_rows{r};
        fprintf(fid, '%s,%s,%s,%s,%s,%.3f,', row{1:6});
        fprintf(fid, '%.3f,%.3f,%.3f,', row{7:9});
        fprintf(fid, '%.4f,%.4f,', row{10:11});
        fprintf(fid, '%.3f,%.3f,', row{12:13});
        fprintf(fid, '%.1f,%.1f,%.1f\n', row{14:16});
    end
    fclose(fid);
    fprintf('\nQC summary saved: %s\n', qc_path);
end

fprintf('\nNext step: PhiID computation (04_compute_phiid_revised.py)\n');
fprintf('============================================================\n');
