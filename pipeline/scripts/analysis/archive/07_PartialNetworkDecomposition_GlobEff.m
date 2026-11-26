%% Partial Network Decomposition (PND) — drug vs placebo, per subject
%
% Adapts Luppi et al. (2024) PND to the GWT pipeline. For every paired
% (drug, placebo) PhiID matrix in results/phiid_matrices/, binarises both
% networks at matched density, runs the PND decomposition of shortest-path
% efficiency, and saves per-subject atom matrices + a degree-preserving
% null distribution.
%
% The PND maths itself (PartialNetworkDecomposition_GlobEff, distance_bin,
% isequalntol) is preserved verbatim from the Luppi/Mediano 2023 release at
% the bottom of this file as local functions.
%
% Output per subject-pair-measure:
%   {prefix}_pnd_red.npy            — redundant efficiency (454x454)
%   {prefix}_pnd_unx_drug.npy       — unique to drug network
%   {prefix}_pnd_uny_placebo.npy    — unique to placebo network
%   {prefix}_pnd_syn.npy            — synergistic (requires both)
%   {prefix}_pnd_dominant_mode.npy  — M matrix (-1/0/1/2/3)
%   {prefix}_pnd_null_distribution.npy  — N_NULLS x 4 atom proportions
%   {prefix}_pnd_summary.yaml
%
% Usage:
%   >> cd('D:\projects\TMS_GWT\pipeline')
%   >> run('scripts/analysis/07_PartialNetworkDecomposition_GlobEff.m')

%% ============================================================
%% SECTION A — SCRIPT DRIVER (loads inputs, calls PND, saves)
%% ============================================================

%% Configuration
project_dir = 'D:\projects\TMS_GWT\pipeline';
N_WORKERS   = 24;
DENSITY     = 0.10;
N_NULLS     = 1000;
MEASURES    = {'synergy', 'redundancy'};
SEED        = 42;
DRY_RUN     = false;     % if true, process only 1 subject per (dataset,drug,measure) and N_NULLS=10

% Condition lives in the SESSION name (task is always "rest" in this pipeline).
% Each entry: {drug_session, placebo_session}
DATASET_PAIRS = struct();
DATASET_PAIRS.london_lsd      = { {'LSD',          'PCB'     } };
DATASET_PAIRS.berlin_ketamine = { {'ketamine',     'placebo' } };
DATASET_PAIRS.basel_lam       = { {'lsd',          'placebo' }, ...
                                  {'mdma',         'placebo' }, ...
                                  {'amph',         'placebo' } };
DATASET_PAIRS.basel_lpm       = { {'lsd',          'placebo' }, ...
                                  {'psilocybin',   'placebo' }, ...
                                  {'mescaline',    'placebo' } };
DATASET_PAIRS.missouri_psilo  = { {'psilocybin',   'baseline'}, ...
                                  {'methylphenidate','baseline'} };

DATASETS = fieldnames(DATASET_PAIRS);

if DRY_RUN
    N_NULLS = 10;
    fprintf('*** DRY_RUN active: 1 subject per group, N_NULLS=%d ***\n', N_NULLS);
end

%% Paths (must be set before parpool so workers inherit them)
npy_path = fullfile(project_dir, 'tools', 'npy-matlab');
addpath(genpath(npy_path));

ts_dir   = fullfile(project_dir, 'results', 'phiid_matrices');
out_root = fullfile(project_dir, 'results', 'pnd');
if ~exist(out_root, 'dir'); mkdir(out_root); end

%% Build the work list: one entry per (dataset, drug_ses, placebo_ses, subject, measure, session_run_key)
work = struct('dataset', {}, 'drug', {}, 'placebo', {}, 'subject', {}, ...
              'measure', {}, 'drug_file', {}, 'placebo_file', {}, ...
              'prefix', {}, 'out_dir', {});

fprintf('============================================================\n');
fprintf('PND — Drug vs Placebo (Luppi 2024 methodology)\n');
fprintf('============================================================\n');

for d = 1:numel(DATASETS)
    ds   = DATASETS{d};
    pairs = DATASET_PAIRS.(ds);
    for p = 1:numel(pairs)
        drug_ses    = pairs{p}{1};
        placebo_ses = pairs{p}{2};
        for mi = 1:numel(MEASURES)
            measure = MEASURES{mi};
            % Find all drug-session matrices for this measure
            drug_glob = sprintf('dataset-%s_sub-*_ses-%s_task-*_run-*_%s.npy', ds, drug_ses, measure);
            drug_files = dir(fullfile(ts_dir, drug_glob));

            n_paired = 0;
            for f = 1:numel(drug_files)
                drug_name = drug_files(f).name;
                % Parse subject + run from drug filename
                sub_tok = regexp(drug_name, 'sub-([^_]+)', 'tokens', 'once');
                run_tok = regexp(drug_name, 'run-(\d+)',   'tokens', 'once');
                if isempty(sub_tok) || isempty(run_tok); continue; end
                subject = sub_tok{1};
                run_id  = run_tok{1};

                % Construct placebo counterpart — same subject, same run, placebo session
                placebo_name = sprintf('dataset-%s_sub-%s_ses-%s_task-rest_run-%s_%s.npy', ...
                                       ds, subject, placebo_ses, run_id, measure);
                placebo_path = fullfile(ts_dir, placebo_name);
                if ~exist(placebo_path, 'file')
                    % Try discovering any placebo run for this subject as fallback
                    pc_glob  = sprintf('dataset-%s_sub-%s_ses-%s_task-*_run-*_%s.npy', ...
                                       ds, subject, placebo_ses, measure);
                    pc_files = dir(fullfile(ts_dir, pc_glob));
                    if isempty(pc_files)
                        fprintf('  [%s/%s vs %s/%s] UNPAIRED: sub-%s (no placebo)\n', ...
                                ds, measure, drug_ses, placebo_ses, subject);
                        continue;
                    end
                    placebo_path = fullfile(pc_files(1).folder, pc_files(1).name);
                end

                drug_path = fullfile(drug_files(f).folder, drug_name);
                out_dir   = fullfile(out_root, ds, measure);
                if ~exist(out_dir, 'dir'); mkdir(out_dir); end
                prefix = sprintf('dataset-%s_sub-%s_task-%s-vs-%s_run-%s', ...
                                 ds, subject, drug_ses, placebo_ses, run_id);

                % Skip if already done
                if exist(fullfile(out_dir, [prefix '_pnd_syn.npy']), 'file')
                    fprintf('  SKIP: %s (%s) — already computed\n', prefix, measure);
                    continue;
                end

                entry.dataset      = ds;
                entry.drug         = drug_ses;
                entry.placebo      = placebo_ses;
                entry.subject      = subject;
                entry.measure      = measure;
                entry.drug_file    = drug_path;
                entry.placebo_file = placebo_path;
                entry.prefix       = prefix;
                entry.out_dir      = out_dir;
                work(end+1) = entry; %#ok<SAGROW>
                n_paired = n_paired + 1;

                if DRY_RUN && n_paired >= 1; break; end
            end
            fprintf('  [%s] %s vs %s (%s): %d paired subjects queued\n', ...
                    ds, drug_ses, placebo_ses, measure, n_paired);
        end
    end
end

n_todo = numel(work);
fprintf('\nTotal work items: %d\n', n_todo);
if n_todo == 0
    fprintf('Nothing to do.\n');
    return;
end

%% Parallel pool
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool('local', N_WORKERS);
elseif pool.NumWorkers ~= N_WORKERS
    delete(pool);
    pool = parpool('local', N_WORKERS);
end
fprintf('Parallel pool ready: %d workers\n\n', pool.NumWorkers);

%% Process in parallel
statuses = zeros(n_todo, 1);
messages = cell(n_todo, 1);

tic_total = tic;

parfor k = 1:n_todo
    entry  = work(k);
    prefix = entry.prefix;
    out_dir = entry.out_dir;

    try
        tic_k = tic;

        M_drug    = readNPY(entry.drug_file);
        M_placebo = readNPY(entry.placebo_file);
        M_drug(~isfinite(M_drug)) = 0;
        M_placebo(~isfinite(M_placebo)) = 0;
        M_drug(M_drug < 0)       = 0;
        M_placebo(M_placebo < 0) = 0;

        net_drug    = binarize_top_density(M_drug,    DENSITY);
        net_placebo = binarize_top_density(M_placebo, DENSITY);

        % ---- CORE PND CALL (unchanged maths) ----
        [PND, Mmode] = PartialNetworkDecomposition_GlobEff(net_drug, net_placebo);

        % Scalar summary: proportion of joint-network efficiency per atom
        [p_red, p_unx, p_uny, p_syn] = atom_proportions(PND);

        % Degree-preserving nulls
        null_props = zeros(N_NULLS, 4);
        for n = 1:N_NULLS
            r_drug    = rewire_degree_preserving(net_drug,    SEED + 1000*k + n);
            r_placebo = rewire_degree_preserving(net_placebo, SEED + 1000*k + n + N_NULLS);
            PND_n = PartialNetworkDecomposition_GlobEff(r_drug, r_placebo);
            [rr, ru, uu, ss] = atom_proportions(PND_n);
            null_props(n,:) = [rr, ru, uu, ss];
        end
        pval_red = mean(null_props(:,1) >= p_red);
        pval_unx = mean(null_props(:,2) >= p_unx);
        pval_uny = mean(null_props(:,3) >= p_uny);
        pval_syn = mean(null_props(:,4) >= p_syn);

        % Save matrices
        writeNPY(PND.red,    fullfile(out_dir, [prefix '_pnd_red.npy']));
        writeNPY(PND.unx,    fullfile(out_dir, [prefix '_pnd_unx_drug.npy']));
        writeNPY(PND.uny,    fullfile(out_dir, [prefix '_pnd_uny_placebo.npy']));
        writeNPY(PND.syn,    fullfile(out_dir, [prefix '_pnd_syn.npy']));
        writeNPY(Mmode,      fullfile(out_dir, [prefix '_pnd_dominant_mode.npy']));
        writeNPY(null_props, fullfile(out_dir, [prefix '_pnd_null_distribution.npy']));

        % YAML summary
        yaml_path = fullfile(out_dir, [prefix '_pnd_summary.yaml']);
        fid = fopen(yaml_path, 'w');
        fprintf(fid, 'dataset: %s\n', entry.dataset);
        fprintf(fid, 'subject: ''%s''\n', entry.subject);
        fprintf(fid, 'drug_session: %s\n', entry.drug);
        fprintf(fid, 'placebo_session: %s\n', entry.placebo);
        fprintf(fid, 'measure: %s\n', entry.measure);
        fprintf(fid, 'density: %.4f\n', DENSITY);
        fprintf(fid, 'n_nulls: %d\n', N_NULLS);
        fprintf(fid, 'n_rois: %d\n', size(net_drug, 1));
        fprintf(fid, 'n_edges_drug: %d\n',    sum(net_drug(:))/2);
        fprintf(fid, 'n_edges_placebo: %d\n', sum(net_placebo(:))/2);
        fprintf(fid, 'prop_red: %.15e\n', p_red);
        fprintf(fid, 'prop_unx_drug: %.15e\n', p_unx);
        fprintf(fid, 'prop_uny_placebo: %.15e\n', p_uny);
        fprintf(fid, 'prop_syn: %.15e\n', p_syn);
        fprintf(fid, 'pval_red: %.6f\n', pval_red);
        fprintf(fid, 'pval_unx_drug: %.6f\n', pval_unx);
        fprintf(fid, 'pval_uny_placebo: %.6f\n', pval_uny);
        fprintf(fid, 'pval_syn: %.6f\n', pval_syn);
        fprintf(fid, 'elapsed_minutes: %.2f\n', toc(tic_k)/60);
        fclose(fid);

        fprintf('[%d/%d] DONE: %s (%s) %.1f min\n', k, n_todo, prefix, entry.measure, toc(tic_k)/60);
        messages{k} = sprintf('OK in %.1f min', toc(tic_k)/60);

    catch ME
        fprintf('[%d/%d] ERROR: %s — %s\n', k, n_todo, prefix, ME.message);
        statuses(k) = 1;
        messages{k} = ME.message;
    end
end

%% Summary
elapsed_total = toc(tic_total);
n_success = sum(statuses == 0);
n_error   = sum(statuses == 1);

fprintf('\n============================================================\n');
fprintf('PND COMPUTATION COMPLETE\n');
fprintf('============================================================\n');
fprintf('  Success: %d\n', n_success);
fprintf('  Errors:  %d\n', n_error);
fprintf('  Total:   %d work items\n', n_todo);
fprintf('  Wall clock: %.1f minutes\n', elapsed_total/60);

if n_error > 0
    fprintf('\nError details:\n');
    for k = 1:n_todo
        if statuses(k) == 1
            fprintf('  - %s: %s\n', work(k).prefix, messages{k});
        end
    end
end

fprintf('\nOutputs: %s\n', out_root);
fprintf('============================================================\n');


%% ============================================================
%% SECTION B — CORE PND MATHS (unchanged, now a local function)
%% ============================================================
function [PND, M] = PartialNetworkDecomposition_GlobEff(netx, nety)

%%PARTIALNETWORKDECOMPOSITION Decomposes global efficiency between two networks
%
%    PND = PARTIALNETWORKDECOMPOSITION(NETX, NETY) computes the partial network
%    decomposition of global efficiency (i.e. inverse shortest path length) for
%    two binary undirected networks NETX and NETY. Returns a structure with the
%    redundant, unique, and synergistic contributions corresponding to both
%    networks.
%
%    [PND, M] = PARTIALNETWORKDECOMPOSITION(NETX, NETY) returns a matrix M
%    representing the dominant 'mode' of each pair of nodes, where M(i,j) = -1
%    for disconnected pairs, 0 for redundant pairs, 1 (resp.  2) for pairs with
%    unique information from network X (resp. Y), and 3 for synergistic pairs.
%
% This function includes code from the Brain Connectivity Toolbox by Rubinov et
% al: https://sites.google.com/site/bctnet/
%
% Pedro Mediano and Andrea Luppi, January 2023

%% Parameter checks
if ~(ismatrix(netx) && ismatrix(nety) && ...
    size(netx,1) == size(netx,2) && all(size(netx) == size(nety)))
  error("Networks must be square matrices of the same size");
end
if ~(issymmetric(netx) && issymmetric(nety))
  error("Both matrices must be symmetric (i.e. undirected networks).");
end
isxbinary = islogical(netx) || length(unique(netx(:))) <= 2;
isybinary = islogical(nety) || length(unique(nety(:))) <= 2;
if ~(isxbinary && isybinary)
  error("Both inputs must be binary networks (i.e. false/true or 0/1 edges).");
end

% Binarise networks and convert to double
netx = double(netx ~= 0);
nety = double(nety ~= 0);


%% Compute efficiencies in both single networks and their union
Ex  = 1./distance_bin(netx);
Ey  = 1./distance_bin(nety);
Exy = 1./distance_bin(netx | nety);

% Turn diagonal entries to NaN's
D = length(netx);
Ex(1:D+1:end)  = nan;
Ey(1:D+1:end)  = nan;
Exy(1:D+1:end) = nan;


%% Compare these efficiencies to compute PND atoms
PND = [];
PND.red = min(Ex, Ey);
PND.unx = Ex - PND.red;
PND.uny = Ey - PND.red;
PND.syn = Exy - max(Ex, Ey);

% Check that the atoms sum up to the right things
assert(isequalntol( Ex, PND.red + PND.unx));
assert(isequalntol( Ey, PND.red + PND.uny));
assert(isequalntol(Exy, PND.red + PND.unx + PND.uny + PND.syn));


%% If requested, compute the dominant mode of each path
if nargout > 1
  M = -1*ones(size(netx));
  M(PND.syn > 0) = 3;
  M(PND.uny > 0 & M < 0) = 2;
  M(PND.unx > 0 & M < 0) = 1;
  M(PND.red > 0 & M < 0) = 0;
  M(1:length(M)+1:end) = nan;

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Small utility function to test whether two arrays containing NaNs are equal
% within a certain tolerance level
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = isequalntol(A, B, tol)
  if nargin < 3 || isempty(tol)
    tol = 1e-12;
  end
  nan_equals = isequal(isnan(A), isnan(B));
  if nan_equals
    num_equals = all(abs(A(~isnan(A)) - B(~isnan(B))) < tol, 'all');
  else
    num_equals = false;
  end
  tf = nan_equals && num_equals;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function taken verbatim from the Brain Connectivity Toolbox on Jan 20th 2023
%
%                   https://sites.google.com/site/bctnet/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D=distance_bin(A)
%DISTANCE_BIN       Distance matrix
%
%   D = distance_bin(A);
%
%   The distance matrix contains lengths of shortest paths between all
%   pairs of nodes. An entry (u,v) represents the length of shortest path
%   from node u to node v. The average shortest path length is the
%   characteristic path length of the network.
%
%   Input:      A,      binary directed/undirected connection matrix
%
%   Output:     D,      distance matrix
%
%   Notes:
%       Lengths between disconnected nodes are set to Inf.
%       Lengths on the main diagonal are set to 0.
%
%   Algorithm: Algebraic shortest paths.
%
%
%   Mika Rubinov, U Cambridge
%   Jonathan Clayden, UCL
%   2007-2013

% Modification history:
% 2007: Original (MR)
% 2013: Bug fix, enforce zero distance for self-connections (JC)

A=double(A~=0);                 %binarize and convert to double format

l=1;                            %path length
Lpath=A;                        %matrix of paths l
D=A;                            %distance matrix

Idx=true;
while any(Idx(:))
    l=l+1;
    Lpath=Lpath*A;
    Idx=(Lpath~=0)&(D==0);
    D(Idx)=l;
end

D(~D)=inf;                      %assign inf to disconnected nodes
D(1:length(A)+1:end)=0;         %clear diagonal
end


%% ============================================================
%% SECTION E — helpers added for pipeline integration
%% ============================================================

function B = binarize_top_density(W, density)
% Keep the top `density` fraction of off-diagonal edges (upper-triangle
% entries) and return a symmetric logical matrix. Ties broken arbitrarily
% by sort; exact density may deviate by at most 1 edge due to integer rounding.
    N = size(W, 1);
    W(1:N+1:end) = 0;                          % zero diagonal
    W(~isfinite(W)) = 0;
    triu_mask = triu(true(N), 1);
    vals = W(triu_mask);
    n_possible = numel(vals);
    n_keep = max(1, round(density * n_possible));
    if n_keep >= n_possible
        thr = -inf;
    else
        sorted = sort(vals, 'descend');
        thr = sorted(n_keep);
    end
    B_upper = false(N);
    B_upper(triu_mask) = vals >= thr;
    B = B_upper | B_upper';                    % symmetrize
end


function R = rewire_degree_preserving(net, seed)
% Maslov-Sneppen double-edge swap on a binary undirected graph.
% Preserves the degree sequence exactly. Attempts 10 * n_edges swaps.
    rng(seed);
    R = logical(net);
    N = size(R, 1);
    R(1:N+1:end) = false;

    [i_idx, j_idx] = find(triu(R, 1));
    m = numel(i_idx);
    if m < 2
        return;
    end

    n_attempts = 10 * m;
    for t = 1:n_attempts
        a = randi(m);
        b = randi(m);
        if a == b; continue; end
        u1 = i_idx(a); v1 = j_idx(a);
        u2 = i_idx(b); v2 = j_idx(b);
        % avoid self-loops and existing edges in the two proposed new pairs
        if u1 == u2 || u1 == v2 || v1 == u2 || v1 == v2
            continue;
        end
        if rand < 0.5
            nu1 = u1; nv1 = v2; nu2 = u2; nv2 = v1;
        else
            nu1 = u1; nv1 = u2; nu2 = v1; nv2 = v2;
        end
        if R(nu1, nv1) || R(nu2, nv2)
            continue;
        end
        % Apply swap
        R(u1, v1) = false; R(v1, u1) = false;
        R(u2, v2) = false; R(v2, u2) = false;
        R(nu1, nv1) = true; R(nv1, nu1) = true;
        R(nu2, nv2) = true; R(nv2, nu2) = true;
        i_idx(a) = min(nu1, nv1); j_idx(a) = max(nu1, nv1);
        i_idx(b) = min(nu2, nv2); j_idx(b) = max(nu2, nv2);
    end
end


function [pr, pux, puy, ps] = atom_proportions(PND)
% Scalar proportion of joint-network efficiency per atom.
    total = nansum(PND.red(:)) + nansum(PND.unx(:)) + nansum(PND.uny(:)) + nansum(PND.syn(:));
    if total == 0
        pr = 0; pux = 0; puy = 0; ps = 0;
        return;
    end
    pr  = nansum(PND.red(:)) / total;
    pux = nansum(PND.unx(:)) / total;
    puy = nansum(PND.uny(:)) / total;
    ps  = nansum(PND.syn(:)) / total;
end
