%% Compute PhiID using Luppi's original MATLAB PhiIDFull.m (parallel version)
%
%
% For each deconvolved timeseries file in results/timeseries/:
%   1. Loads the .npy timeseries (timepoints x 454 ROIs)
%   2. Computes all 16 PhiID atoms for every ROI pair using PhiIDFull.m
%      (PhiIDFull normalizes to unit variance internally — no pre-z-scoring)
%   3. Saves 5 matrices + 1 summary YAML per subject
%
% Output files (compatible with downstream Python scripts):
%   {prefix}_synergy.npy        — persistent synergy (atoms.sts)
%   {prefix}_redundancy.npy     — persistent redundancy (atoms.rtr)
%   {prefix}_synergy_total.npy  — Eq. 5 total synergy
%   {prefix}_phi.npy            — Eq. 6: TDMI - I(X_past;X_fut) - I(Y_past;Y_fut)
%   {prefix}_phi_r.npy          — Eq. 7: Phi + Red
%   {prefix}_summary.yaml       — metadata
%
% Prerequisites:
%   - Luppi PhiIDFull.m in docs/luppi_supps/software/MMI_PhiID/
%   - JIDT (infodynamics.jar) in MMI_PhiID/private/
%   - npy-matlab in tools/npy-matlab/
%
% Usage:
%   >> cd('D:\projects\TMS_GWT\pipeline')
%   >> run('scripts/analysis/step04_compute_phiid_matlab.m')
%
% Parallelization:
%   Uses parfor across subjects with 30 workers.
%   Each worker independently calls PhiIDFull.m (JIDT auto-loads per worker).

%% Configuration
N_WORKERS = 24;

%% Setup paths (must be done before parpool so workers inherit them)
project_dir = 'D:\projects\TMS_GWT\pipeline';
phiid_path  = fullfile(project_dir, 'docs', 'luppi_supps', 'software', 'MMI_PhiID');
npy_path    = fullfile(project_dir, 'tools', 'npy-matlab');

addpath(phiid_path);
addpath(genpath(npy_path));

ts_dir  = fullfile(project_dir, 'results', 'timeseries');
out_dir = fullfile(project_dir, 'results', 'phiid_matrices');

if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% Discover deconvolved timeseries files
all_files = dir(fullfile(ts_dir, 'dataset-*_timeseries_deconv.npy'));
mask = contains({all_files.name}, '_timeseries_deconv.npy');
all_files = all_files(mask);
n_files = length(all_files);

fprintf('============================================================\n');
fprintf('PhiID Computation — Luppi MATLAB PhiIDFull.m (PARALLEL)\n');
fprintf('============================================================\n');
fprintf('Found %d deconvolved timeseries files\n', n_files);
fprintf('Workers: %d\n', N_WORKERS);
fprintf('Input:  %s\n', ts_dir);
fprintf('Output: %s\n', out_dir);
fprintf('============================================================\n\n');

%% Filter out already-computed files
todo_idx = [];
for f = 1:n_files
    prefix = strrep(all_files(f).name, '_timeseries_deconv.npy', '');
    phi_r_path = fullfile(out_dir, [prefix '_phi_r.npy']);
    if exist(phi_r_path, 'file')
        fprintf('SKIP: %s (already computed)\n', all_files(f).name);
    else
        todo_idx(end+1) = f; %#ok<AGROW>
    end
end

n_todo = length(todo_idx);
fprintf('\n%d files to process, %d already done\n\n', n_todo, n_files - n_todo);

if n_todo == 0
    fprintf('Nothing to do.\n');
    return;
end

% Extract file info into simple arrays for parfor (no struct indexing)
todo_names   = cell(n_todo, 1);
todo_folders = cell(n_todo, 1);
for k = 1:n_todo
    todo_names{k}   = all_files(todo_idx(k)).name;
    todo_folders{k}  = all_files(todo_idx(k)).folder;
end

%% Start parallel pool
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local', N_WORKERS);
elseif pool.NumWorkers ~= N_WORKERS
    delete(pool);
    pool = parpool('local', N_WORKERS);
end
fprintf('Parallel pool ready: %d workers\n\n', pool.NumWorkers);

%% Process files in parallel
% Collect results per file: 0=success, 1=error, status message
statuses  = zeros(n_todo, 1);
messages  = cell(n_todo, 1);

tic_total = tic;

parfor k = 1:n_todo
    fname      = todo_names{k};
    input_path = fullfile(todo_folders{k}, fname);
    prefix     = strrep(fname, '_timeseries_deconv.npy', '');

    try
        fprintf('[%d/%d] START: %s\n', k, n_todo, fname);
        tic_scan = tic;

        % Load timeseries (timepoints x n_rois)
        bold = readNPY(input_path);
        [n_tp, n_rois] = size(bold);

        % No z-scoring: PhiIDFull normalizes to unit variance internally
        % (matching Luppi's usage in README: PhiIDFull([my_BOLD(row,:); my_BOLD(col,:)]))

        % Detect zero-variance ROIs (singular covariance breaks PhiIDFull)
        valid_rois = find(std(bold) > 0);
        n_invalid = n_rois - length(valid_rois);
        if n_invalid > 0
            fprintf('  [%s] WARNING: %d zero-variance ROIs, skipping pairs\n', prefix, n_invalid);
        end
        valid_set = false(1, n_rois);
        valid_set(valid_rois) = true;

        % Initialize output matrices with NaN (not zero). Zero-variance ROI pairs
        % stay NaN so downstream Python uses nansum/nanmean to exclude them — zero
        % would bias rankings because synergy/redundancy are typically positive.
        syn_mat       = nan(n_rois);
        red_mat       = nan(n_rois);
        syn_total_mat = nan(n_rois);
        phi_mat       = nan(n_rois);
        phi_r_mat     = nan(n_rois);

        % Compute PhiID for each ROI pair
        for i = 1:n_rois
            if mod(i, 100) == 0
                fprintf('  [%s] ROI %d/%d\n', prefix, i, n_rois);
            end

            for j = (i+1):n_rois
                % Skip pairs involving zero-variance ROIs
                if ~valid_set(i) || ~valid_set(j)
                    continue;
                end

                % PhiIDFull expects D x T matrix (2 x n_timepoints)
                X = [bold(:,i)'; bold(:,j)'];
                atoms = PhiIDFull(X, 1);  % tau = 1

                % Persistent synergy & redundancy (direct atoms)
                syn_mat(i,j) = atoms.sts;
                syn_mat(j,i) = atoms.sts;
                red_mat(i,j) = atoms.rtr;
                red_mat(j,i) = atoms.rtr;

                % Total synergy (all source-synergy atoms)
                syn_total = atoms.str + atoms.stx + atoms.sty + atoms.sts;
                syn_total_mat(i,j) = syn_total;
                syn_total_mat(j,i) = syn_total;

                % Phi (Eq. 6) = TDMI - I(X_past;X_fut) - I(Y_past;Y_fut)
                % TDMI = sum of all 16 atoms (M matrix row 16)
                tdmi = atoms.rtr + atoms.rtx + atoms.rty + atoms.rts ...
                     + atoms.xtr + atoms.xtx + atoms.xty + atoms.xts ...
                     + atoms.ytr + atoms.ytx + atoms.yty + atoms.yts ...
                     + atoms.str + atoms.stx + atoms.sty + atoms.sts;
                % I(X_past;X_fut) = rtr + rtx + xtr + xtx (M matrix row 8)
                I_xta = atoms.rtr + atoms.rtx + atoms.xtr + atoms.xtx;
                % I(Y_past;Y_fut) = rtr + rty + ytr + yty (M matrix row 11)
                I_ytb = atoms.rtr + atoms.rty + atoms.ytr + atoms.yty;
                phi = tdmi - I_xta - I_ytb;
                phi_mat(i,j) = phi;
                phi_mat(j,i) = phi;

                % Eq. 7: Phi_R = Phi + Red
                phi_r = phi + atoms.rtr;
                phi_r_mat(i,j) = phi_r;
                phi_r_mat(j,i) = phi_r;
            end
        end

        elapsed_scan = toc(tic_scan);

        % ---- Save all 5 matrices as .npy ----
        writeNPY(syn_mat,       fullfile(out_dir, [prefix '_synergy.npy']));
        writeNPY(red_mat,       fullfile(out_dir, [prefix '_redundancy.npy']));
        writeNPY(syn_total_mat, fullfile(out_dir, [prefix '_synergy_total.npy']));
        writeNPY(phi_mat,       fullfile(out_dir, [prefix '_phi.npy']));
        writeNPY(phi_r_mat,     fullfile(out_dir, [prefix '_phi_r.npy']));

        % ---- Save summary YAML ----
        triu_idx = triu(true(n_rois), 1);
        yaml_path = fullfile(out_dir, [prefix '_summary.yaml']);
        fid = fopen(yaml_path, 'w');

        ds_tok = regexp(prefix, 'dataset-(\w+?)_sub-', 'tokens', 'once');
        if ~isempty(ds_tok); fprintf(fid, 'dataset: %s\n', ds_tok{1}); end

        sub_tok = regexp(prefix, 'sub-([^_]+)', 'tokens', 'once');
        if ~isempty(sub_tok); fprintf(fid, 'subject: ''%s''\n', sub_tok{1}); end

        ses_tok = regexp(prefix, 'ses-([^_]+)', 'tokens', 'once');
        if ~isempty(ses_tok)
            fprintf(fid, 'session: ''%s''\n', ses_tok{1});
        else
            fprintf(fid, 'session: null\n');
        end

        task_tok = regexp(prefix, 'task-([^_]+)', 'tokens', 'once');
        if ~isempty(task_tok); fprintf(fid, 'task: %s\n', task_tok{1}); end

        run_tok = regexp(prefix, 'run-(\d+)', 'tokens', 'once');
        if ~isempty(run_tok); fprintf(fid, 'run: %s\n', run_tok{1}); end

        fprintf(fid, 'tau: 1\n');
        fprintf(fid, 'method: gaussian_matlab\n');
        fprintf(fid, 'implementation: PhiIDFull.m (Luppi/Mediano MATLAB)\n');
        fprintf(fid, 'n_timepoints: %d\n', n_tp);
        fprintf(fid, 'n_rois: %d\n', n_rois);
        fprintf(fid, 'elapsed_minutes: %.2f\n', elapsed_scan / 60);
        % NaN-aware summary stats — skip zero-variance pairs, not include them as 0
        fprintf(fid, 'n_zero_variance_rois: %d\n', n_invalid);
        fprintf(fid, 'mean_redundancy: %.15e\n', mean(red_mat(triu_idx), 'omitnan'));
        fprintf(fid, 'mean_synergy: %.15e\n', mean(syn_mat(triu_idx), 'omitnan'));
        fprintf(fid, 'mean_synergy_total: %.15e\n', mean(syn_total_mat(triu_idx), 'omitnan'));
        fprintf(fid, 'mean_phi: %.15e\n', mean(phi_mat(triu_idx), 'omitnan'));
        fprintf(fid, 'mean_phi_r: %.15e\n', mean(phi_r_mat(triu_idx), 'omitnan'));
        phi_finite = phi_mat(triu_idx);
        phi_finite = phi_finite(isfinite(phi_finite));
        fprintf(fid, 'negative_phi_percentage: %.4f\n', ...
            100 * sum(phi_finite < 0) / max(length(phi_finite), 1));
        phi_r_finite = phi_r_mat(triu_idx);
        phi_r_finite = phi_r_finite(isfinite(phi_r_finite));
        fprintf(fid, 'positive_phi_r_check: %s\n', ...
            mat2str(all(phi_r_finite >= 0)));
        fclose(fid);

        fprintf('[%d/%d] DONE: %s (%.1f min)\n', k, n_todo, prefix, elapsed_scan / 60);
        messages{k} = sprintf('OK in %.1f min', elapsed_scan / 60);

    catch ME
        fprintf('[%d/%d] ERROR: %s — %s\n', k, n_todo, fname, ME.message);
        statuses(k) = 1;
        messages{k} = ME.message;
    end
end

%% Summary
elapsed_total = toc(tic_total);
n_success = sum(statuses == 0);
n_error   = sum(statuses == 1);

fprintf('\n============================================================\n');
fprintf('PHIID COMPUTATION COMPLETE (MATLAB, PARALLEL)\n');
fprintf('============================================================\n');
fprintf('  Success: %d\n', n_success);
fprintf('  Skipped: %d (already computed)\n', n_files - n_todo);
fprintf('  Errors:  %d\n', n_error);
fprintf('  Total:   %d files\n', n_files);
fprintf('  Wall clock: %.1f minutes\n', elapsed_total / 60);

if n_error > 0
    fprintf('\nError details:\n');
    for k = 1:n_todo
        if statuses(k) == 1
            fprintf('  - %s: %s\n', todo_names{k}, messages{k});
        end
    end
end

fprintf('\nNext step: python scripts/analysis/05_group_average_workspace.py --dataset hcp\n');
fprintf('============================================================\n');
