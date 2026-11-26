"""TMS-only publication analyses using the frozen GWT metric definitions.

No preprocessing is launched here. Only TMS matrices, the fixed HCP workspace,
atlas labels, existing motion QC and published SUDMEX clinical scores are read.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import yaml

NETWORKS = ('VIS', 'SMN', 'DAN', 'VAN', 'LIM', 'FPN', 'DMN', 'SUB')
DATASETS = ('sudmex_tms', 'aachen_itbs')
MATRIX_NAMES = ('synergy', 'redundancy', 'synergy_total', 'phi', 'phi_r')
PRIMARY = 'gateways_FPN'
KEYS = ['dataset', 'subject_id', 'comparison']
METRIC_INFO = {
    **{'global_' + name: ('global_phiid', name.replace('_', ' ').title()) for name in MATRIX_NAMES},
    'mean_synergy_pc': ('workspace_pc', 'Mean workspace synergy PC'),
    'mean_redundancy_pc': ('workspace_pc', 'Mean workspace redundancy PC'),
    'gateways_total': ('total_roles', 'Total workspace gateways'),
    'broadcasters_total': ('total_roles', 'Total workspace broadcasters'),
    **{'gateways_' + n: ('primary' if n == 'FPN' else 'other_network_gateways', n + ' gateways') for n in NETWORKS},
    **{'broadcasters_' + n: ('network_broadcasters', n + ' broadcasters') for n in NETWORKS},
}
METRICS = list(METRIC_INFO)
PAPER_URL = 'https://www.biorxiv.org/content/10.64898/2026.09.28.755172v1.full'


def find_root(start=None):
    start = Path(start or Path.cwd()).resolve()
    for candidate in (start, *start.parents):
        if (candidate / 'data/sudmex_tms').is_dir() and (candidate / 'pipeline/results/phiid_matrices').is_dir():
            return candidate
    raise FileNotFoundError('Run inside TMS_GWT, or pass its root explicitly.')


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            digest.update(block)
    return digest.hexdigest()


def load_reference(root):
    """Read the fixed 235-region mask; never rederive it from HCP or TMS data."""
    root = Path(root)
    atlas = root / 'pipeline/atlases/roi_to_network_mapping.csv'
    label_file = root / 'pipeline/atlases/Schaefer2018_400Parcels_7Networks_order_Tian_Subcortex_S4_label.txt'
    mask_file = root / 'pipeline/results/group_analysis/hcp/rest/hcp_rest_workspace_mask.npy'
    frame = pd.read_csv(atlas)
    frame = frame[frame.ROI > 0].sort_values('ROI').reset_index(drop=True)
    if frame.ROI.tolist() != list(range(1, 455)):
        raise ValueError('Atlas ROI order must be 1..454')
    modules = frame.Network_ID.to_numpy(dtype=int)
    # Independent cross-check against the exact atlas mapping used by Yousefi.
    original = np.full(454, -1)
    yeo = {'Vis': 0, 'SomMot': 1, 'DorsAttn': 2, 'SalVentAttn': 3, 'Limbic': 4, 'Cont': 5, 'Default': 6}
    lines = label_file.read_text().splitlines()
    for name, numbers in zip(lines[::2], lines[1::2]):
        index = int(numbers.split()[0])
        roi, module = (index - 54, yeo[name.split('_')[2]]) if name.startswith('7Networks_') else (index + 400, 7)
        original[roi - 1] = module
    if not np.array_equal(modules, original) or set(modules) != set(range(8)):
        raise ValueError('Atlas module mapping differs from the publication notebook')
    workspace = np.load(mask_file).astype(bool)
    if workspace.shape != (454,) or int(workspace.sum()) != 235:
        raise ValueError('Unexpected fixed HCP workspace')
    return {'workspace': workspace, 'modules': modules, 'labels': frame,
            'input_files': [atlas, label_file, mask_file]}


def participation_coefficient(matrix, modules):
    """Same NaN-aware, diagonal-excluded eight-module PC as publication cell 2."""
    matrix = np.array(matrix, dtype=float, copy=True)
    np.fill_diagonal(matrix, 0)
    total = np.nansum(matrix, axis=1)
    safe = np.where(total == 0, np.inf, total)
    ss = np.zeros(len(matrix))
    for module in np.unique(modules):
        ss += (np.nansum(matrix[:, modules == module], axis=1) / safe) ** 2
    result = 1 - ss
    result[(total == 0) | (~np.isfinite(matrix).any(axis=1))] = np.nan
    return result


def classify_roles(syn_pc, red_pc, workspace):
    valid = workspace & np.isfinite(syn_pc) & np.isfinite(red_pc)
    idx = np.flatnonzero(valid)
    gradient = stats.rankdata(syn_pc[idx], method='average') - stats.rankdata(red_pc[idx], method='average')
    gateway = np.zeros(len(workspace), bool)
    broadcaster = np.zeros(len(workspace), bool)
    gateway[idx[gradient > 0]] = True
    broadcaster[idx[gradient < 0]] = True
    score = np.full(len(workspace), np.nan)
    score[idx] = gradient
    return gateway, broadcaster, score


def metric_values(matrices, reference):
    workspace, modules = reference['workspace'], reference['modules']
    upper = np.triu_indices(len(workspace), 1)
    row = {}
    for name in MATRIX_NAMES:
        values = matrices[name][upper]
        row['global_' + name] = float(np.nanmean(values)) if np.isfinite(values).any() else np.nan
    syn = participation_coefficient(matrices['synergy'], modules)
    red = participation_coefficient(matrices['redundancy'], modules)
    gateway, broadcaster, score = classify_roles(syn, red, workspace)
    valid_ws = workspace & np.isfinite(syn) & np.isfinite(red)
    complete = bool(valid_ws.sum() == workspace.sum())
    row.update(n_valid_workspace=int(valid_ws.sum()), full_workspace=complete,
               n_valid_pairs=int(np.isfinite(matrices['synergy'][upper]).sum()))
    # Match the publication's complete-workspace PC/count analysis. Missing
    # regions must not appear as fewer gateways or an average over fewer nodes.
    row['mean_synergy_pc'] = float(syn[workspace].mean()) if complete else np.nan
    row['mean_redundancy_pc'] = float(red[workspace].mean()) if complete else np.nan
    row['gateways_total'] = int(gateway.sum()) if complete else np.nan
    row['broadcasters_total'] = int(broadcaster.sum()) if complete else np.nan
    row['ties_total'] = int(np.sum(score[workspace] == 0)) if complete else np.nan
    for i, network in enumerate(NETWORKS):
        row['gateways_' + network] = int((gateway & (modules == i)).sum()) if complete else np.nan
        row['broadcasters_' + network] = int((broadcaster & (modules == i)).sum()) if complete else np.nan
    profiles = pd.DataFrame({'roi': np.flatnonzero(workspace) + 1,
                             'network': [NETWORKS[m] for m in modules[workspace]],
                             'synergy_pc': syn[workspace], 'redundancy_pc': red[workspace],
                             'rank_score': score[workspace], 'gateway': gateway[workspace], 'broadcaster': broadcaster[workspace]})
    return row, profiles


def inventory(root):
    """Expected scans come from existing QC, including excluded/missing scans."""
    root = Path(root)
    phi = root / 'pipeline/results/phiid_matrices'
    ts = root / 'pipeline/results/timeseries'
    participants = pd.read_csv(root / 'data/sudmex_tms/participants.tsv', sep='\t')
    arm_lookup = participants.set_index('participant_id').group.map({1: 'sham', 2: 'active'}).to_dict()
    records, files = [], [root / 'data/sudmex_tms/participants.tsv']
    for dataset in DATASETS:
        qc_file = root / f'pipeline/results/qc/{dataset}/exclusions.csv'
        files.append(qc_file)
        qc = pd.read_csv(qc_file, dtype={'subject': str, 'session': str})
        if qc.duplicated(['subject', 'session', 'task']).any():
            raise ValueError('Multiple runs require an explicit run selection, not silent averaging')
        for row in qc.itertuples():
            sid = row.subject.zfill(3 if dataset == 'sudmex_tms' else 2)
            prefix = f'dataset-{dataset}_sub-{sid}_ses-{row.session}_task-rest_run-1'
            summary_file = phi / (prefix + '_summary.yaml')
            summary = {}
            if summary_file.is_file():
                summary = yaml.safe_load(summary_file.read_text(encoding='utf-8')) or {}
            ready = all((phi / (prefix + '_' + m + '.npy')).is_file() for m in MATRIX_NAMES)
            ready = ready and {'n_rois', 'n_timepoints', 'method', 'positive_phi_r_check'}.issubset(summary)
            arm = arm_lookup.get('sub-' + sid) if dataset == 'sudmex_tms' else 'iTBS'
            if dataset == 'sudmex_tms' and arm is None:
                raise ValueError('Missing SUDMEX randomized arm: ' + sid)
            records.append({'dataset': dataset, 'subject_id': dataset + '_' + sid, 'subject': sid,
                            'session': row.session, 'arm': arm, 'qc_status': row.status,
                            'eligible': row.status in ('INCLUDE', 'TRUNCATE'), 'prefix': prefix,
                            'mean_fd': row.usable_mean_fd_mm, 'regional_ready': (ts / (prefix + '_timeseries.npy')).is_file(),
                            'hrf_ready': all((ts / (prefix + s)).is_file() for s in ('_timeseries_deconv.npy', '_hrf.mat', '_timeseries_deconv_metadata.yaml')),
                            'phiid_ready': bool(ready)})
    return pd.DataFrame(records), files


def extract_scan_metrics(root, scan_inventory, reference):
    root = Path(root)
    rows, regions, files = [], [], []
    phi = root / 'pipeline/results/phiid_matrices'
    for scan in scan_inventory[scan_inventory.eligible & scan_inventory.phiid_ready].itertuples():
        if scan.dataset not in DATASETS or not scan.prefix.startswith('dataset-' + scan.dataset + '_'):
            raise ValueError('Non-TMS input refused')
        paths = {m: phi / (scan.prefix + '_' + m + '.npy') for m in MATRIX_NAMES}
        summary_path = phi / (scan.prefix + '_summary.yaml')
        summary = yaml.safe_load(summary_path.read_text(encoding='utf-8'))
        if summary['n_rois'] != 454 or summary['n_timepoints'] < 240 or summary.get('tau') != 1 or summary['method'] != 'gaussian_matlab':
            raise ValueError('Unexpected PhiID production parameters: ' + scan.prefix)
        matrices = {m: np.load(p) for m, p in paths.items()}
        valid_mask = None
        for name, matrix in matrices.items():
            if matrix.shape != (454, 454) or np.isinf(matrix).any() or not np.allclose(matrix, matrix.T, equal_nan=True):
                raise ValueError('Malformed/asymmetric matrix: ' + str(paths[name]))
            if not np.isnan(np.diag(matrix)).all():
                raise ValueError('Production PhiID matrices must have NaN diagonal')
            finite = np.isfinite(matrix)
            if valid_mask is not None and not np.array_equal(finite, valid_mask):
                raise ValueError('Inconsistent missing pairs across PhiID matrices')
            valid_mask = finite
        valid_rois = valid_mask.any(axis=1)
        if valid_rois.sum() < 2 or summary.get('n_zero_variance_rois') != int((~valid_rois).sum()):
            raise ValueError('Invalid ROI coverage or inconsistent summary: ' + scan.prefix)
        if str(summary.get('subject')) != scan.subject or summary.get('dataset') != scan.dataset or summary.get('session') != scan.session:
            raise ValueError('PhiID metadata identity mismatch: ' + scan.prefix)
        expected = np.outer(valid_rois, valid_rois) & ~np.eye(454, dtype=bool)
        if not np.array_equal(valid_mask, expected):
            raise ValueError('Unexpected isolated missing edges: ' + scan.prefix)
        if not np.allclose(matrices['phi_r'], matrices['phi'] + matrices['redundancy'], equal_nan=True, atol=1e-10):
            raise ValueError('Phi_R identity failed: ' + scan.prefix)
        values, profile = metric_values(matrices, reference)
        identity = {k: getattr(scan, k) for k in ('dataset', 'subject_id', 'session', 'arm', 'mean_fd', 'prefix')}
        rows.append({**identity, 'n_timepoints': summary['n_timepoints'], 'n_valid_rois': int(valid_rois.sum()),
                     'valid_roi_signature': ''.join('1' if v else '0' for v in valid_rois), **values})
        for key in ('dataset', 'subject_id', 'session', 'arm'):
            profile[key] = identity[key]
        regions.append(profile)
        files.extend([*paths.values(), summary_path])
    columns = ['dataset', 'subject_id', 'session', 'arm', 'mean_fd', 'prefix', 'n_timepoints', 'n_valid_rois',
               'valid_roi_signature', 'n_valid_workspace', 'full_workspace', 'n_valid_pairs', 'ties_total', *METRICS]
    region_columns = ['roi', 'network', 'synergy_pc', 'redundancy_pc', 'rank_score', 'gateway', 'broadcaster', 'dataset', 'subject_id', 'session', 'arm']
    return pd.DataFrame(rows, columns=columns), (pd.concat(regions, ignore_index=True) if regions else pd.DataFrame(columns=region_columns)), files


def paired_changes(scan_inventory, scans):
    """One record per expected QC-eligible participant/contrast; missing stays NaN."""
    plans = {'sudmex_tms': [('t1_minus_t0', 't0', 't1')],
             'aachen_itbs': [('post90_minus_pre90', 'pre90', 'post90'), ('post120_minus_pre120', 'pre120', 'post120')]}
    available = scans.set_index(['dataset', 'subject_id', 'session'])
    rows = []
    for dataset, contrasts in plans.items():
        good = scan_inventory[(scan_inventory.dataset == dataset) & scan_inventory.eligible]
        for comparison, pre, post in contrasts:
            subjects = sorted(set(good.loc[good.session == pre, 'subject_id']) & set(good.loc[good.session == post, 'subject_id']))
            for sid in subjects:
                row = {'dataset': dataset, 'subject_id': sid, 'comparison': comparison, 'baseline_session': pre, 'followup_session': post,
                       'arm': good.loc[good.subject_id == sid, 'arm'].iloc[0]}
                keys = [(dataset, sid, pre), (dataset, sid, post)]
                ready = all(key in available.index for key in keys)
                row['pair_ready'] = ready
                a, b = (available.loc[k] for k in keys) if ready else (None, None)
                same_rois = ready and a.valid_roi_signature == b.valid_roi_signature
                row['same_valid_rois'] = bool(same_rois)
                row['delta_fd'] = float(b.mean_fd - a.mean_fd) if ready else np.nan
                for metric in METRICS:
                    row[metric + '_baseline'] = a[metric] if ready else np.nan
                    row[metric + '_followup'] = b[metric] if ready else np.nan
                    # Coverage changes can imitate a global-mean effect. Keep raw
                    # values, but do not test a change over different ROI pairs.
                    row[metric + '_delta'] = (b[metric] - a[metric]) if ready and same_rois else np.nan
                rows.append(row)
    return pd.DataFrame(rows)


def pool_aachen(pairs):
    source = pairs[pairs.dataset == 'aachen_itbs']
    rows = []
    numeric = ['delta_fd'] + [m + suffix for m in METRICS for suffix in ('_baseline', '_followup', '_delta')]
    for sid, block in source.groupby('subject_id'):
        if set(block.comparison) != {'post90_minus_pre90', 'post120_minus_pre120'}:
            continue
        item = {'dataset': 'aachen_itbs', 'subject_id': sid, 'arm': 'iTBS', 'comparison': 'mean_of_two_paired_changes',
                'pair_ready': bool(block.pair_ready.all()), 'same_valid_rois': bool(block.same_valid_rois.all())}
        # Do not let pandas silently average one visit when the other is missing.
        for metric in numeric:
            values = block[metric].to_numpy(dtype=float)
            item[metric] = float(values.mean()) if np.isfinite(values).all() else np.nan
        rows.append(item)
    return pd.DataFrame(rows)


def bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    result = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        order = np.argsort(p[ok])
        corrected = p[ok][order] * ok.sum() / np.arange(1, ok.sum() + 1)
        corrected = np.minimum.accumulate(corrected[::-1])[::-1].clip(0, 1)
        values = np.empty(ok.sum())
        values[order] = corrected
        result[ok] = values
    return result


def welch_summary(active, sham):
    active, sham = np.asarray(active, float), np.asarray(sham, float)
    active, sham = active[np.isfinite(active)], sham[np.isfinite(sham)]
    result = {'n_active': len(active), 'n_sham': len(sham), 'estimate': np.nan, 'ci_low': np.nan, 'ci_high': np.nan,
              'p': np.nan, 't': np.nan, 'df': np.nan, 'effect_size': np.nan, 'effect_size_type': 'Hedges_g_between_change_scores'}
    if min(len(active), len(sham)) < 3:
        return {**result, 'status': 'insufficient_observations'}
    n1, n0 = len(active), len(sham)
    v1, v0 = active.var(ddof=1), sham.var(ddof=1)
    estimate = float(active.mean() - sham.mean())
    se2 = v1 / n1 + v0 / n0
    result['estimate'] = estimate
    if se2 == 0:
        return {**result, 'status': 'zero_variance'}
    df = se2**2 / ((v1 / n1)**2 / (n1 - 1) + (v0 / n0)**2 / (n0 - 1))
    se = np.sqrt(se2)
    t = estimate / se
    ci = stats.t.ppf(.975, df) * se
    pooled = np.sqrt(((n1 - 1)*v1 + (n0 - 1)*v0) / (n1 + n0 - 2))
    correction = 1 - 3 / (4 * (n1 + n0 - 2) - 1)
    return {**result, 'estimate': estimate, 'ci_low': estimate-ci, 'ci_high': estimate+ci, 't': t, 'df': df,
            'p': 2*stats.t.sf(abs(t), df), 'effect_size': correction*estimate/pooled,
            'mean_change_active': float(active.mean()), 'mean_change_sham': float(sham.mean()), 'status': 'ok'}


def paired_summary(changes):
    changes = np.asarray(changes, float)
    changes = changes[np.isfinite(changes)]
    n = len(changes)
    result = {'n': n, 'estimate': float(changes.mean()) if n else np.nan, 'ci_low': np.nan, 'ci_high': np.nan,
              't': np.nan, 'df': n - 1, 'p': np.nan, 'effect_size': np.nan, 'effect_size_type': 'Hedges_g_z'}
    if n < 3 or changes.std(ddof=1) == 0:
        return {**result, 'status': 'insufficient_observations_or_zero_variance'}
    se = changes.std(ddof=1) / np.sqrt(n)
    t = result['estimate'] / se
    ci = stats.t.ppf(.975, n - 1) * se
    return {**result, 't': t, 'p': 2*stats.t.sf(abs(t), n - 1), 'ci_low': result['estimate']-ci,
            'ci_high': result['estimate']+ci, 'effect_size': (1-3/(4*(n-1)-1))*result['estimate']/changes.std(ddof=1), 'status': 'ok'}


def correct_families(table):
    table = table.copy()
    table['q_family'] = np.nan
    if len(table):
        for _, block in table.groupby('family'):
            table.loc[block.index, 'q_family'] = bh(block.p)
    return table


def neural_tests(pairs, between_arms):
    rows = []
    if pairs.empty or not pairs.pair_ready.all():
        return pd.DataFrame()
    if pairs.subject_id.duplicated().any():
        raise ValueError('Participants must occur once; average Aachen visits before inference')
    for metric in METRICS:
        delta = metric + '_delta'
        result = welch_summary(pairs.loc[pairs.arm == 'active', delta], pairs.loc[pairs.arm == 'sham', delta]) if between_arms else paired_summary(pairs[delta])
        rows.append({'metric': metric, 'family': METRIC_INFO[metric][0], 'label': METRIC_INFO[metric][1], **result})
    return correct_families(pd.DataFrame(rows))


def fit_hc3(frame, outcome, predictor, covariates, standardize_predictor=False):
    cols = list(dict.fromkeys([outcome, predictor, *covariates]))
    data = frame[cols].replace([np.inf, -np.inf], np.nan).dropna()
    result = {'n': len(data), 'estimate': np.nan, 'ci_low': np.nan, 'ci_high': np.nan, 'p': np.nan,
              'model': outcome + ' ~ ' + ' + '.join([predictor, *covariates])}
    if data.empty or data[predictor].nunique() < 2:
        return {**result, 'status': 'constant_or_missing_predictor'}
    x = data[[predictor, *covariates]].astype(float).copy()
    dropped = [c for c in covariates if x[c].nunique() < 2]
    x = x.drop(columns=dropped)
    if standardize_predictor:
        x[predictor] = (x[predictor] - x[predictor].mean()) / x[predictor].std(ddof=1)
    x = sm.add_constant(x, has_constant='add')
    if len(data) <= x.shape[1] + 2 or np.linalg.matrix_rank(x) < x.shape[1]:
        return {**result, 'status': 'insufficient_observations_or_rank_deficient'}
    fit = sm.OLS(data[outcome].astype(float), x).fit(cov_type='HC3', use_t=True)
    lo, hi = fit.conf_int().loc[predictor]
    return {**result, 'estimate': fit.params[predictor], 'ci_low': lo, 'ci_high': hi, 'p': fit.pvalues[predictor],
            'dropped_constant_covariates': ','.join(dropped), 'status': 'ok'}


def motion_sensitivity(pairs):
    if pairs.empty or not pairs.pair_ready.all():
        return pd.DataFrame()
    frame = pairs.copy()
    frame['active'] = frame.arm.eq('active').astype(int)
    rows = []
    for metric in METRICS:
        result = fit_hc3(frame, metric + '_delta', 'active', [metric + '_baseline', 'delta_fd'])
        rows.append({'metric': metric, 'family': METRIC_INFO[metric][0], **result})
    return correct_families(pd.DataFrame(rows))


def clinical_analyses(root, sudmex_pairs, source_path, neural_metrics, include_flagged=False):
    if sudmex_pairs.empty or not sudmex_pairs.pair_ready.all():
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    utils = Path(root) / 'pipeline/scripts/utils'
    if str(utils) not in sys.path:
        sys.path.insert(0, str(utils))
    from clinical_data import METRICS as CLINICAL_METRICS, merge_clinical_neural
    # The validated clinical table supplies the same explicit t0/t1 labels.
    frame = merge_clinical_neural(sudmex_pairs.drop(columns=['baseline_session', 'followup_session']),
                                  source_path, include_flagged=include_flagged)
    frame['active'] = frame.arm.eq('active').astype(int)
    clinical_rows, association_rows = [], []
    for clinical in CLINICAL_METRICS:
        delta = clinical + '_delta'
        result = welch_summary(frame.loc[frame.arm == 'active', delta], frame.loc[frame.arm == 'sham', delta])
        clinical_rows.append({'metric': clinical, 'family': 'clinical_outcomes', **result})
        for neural in neural_metrics:
            result = fit_hc3(frame, delta, neural + '_delta',
                             ['active', clinical + '_t0', neural + '_baseline', 'delta_fd'], standardize_predictor=True)
            association_rows.append({'clinical_metric': clinical, 'neural_metric': neural,
                                     'family': 'clinical_neural_associations', 'units': 'clinical change per 1 SD neural change', **result})
    return frame, correct_families(pd.DataFrame(clinical_rows)), correct_families(pd.DataFrame(association_rows))


def regional_changes(regions, pairs):
    if regions.empty or pairs.empty or not pairs.pair_ready.all():
        return pd.DataFrame()
    rows = []
    for pair in pairs.itertuples():
        if not pair.same_valid_rois:
            continue
        one = regions[(regions.dataset == pair.dataset) & (regions.subject_id == pair.subject_id)]
        a = one[one.session == pair.baseline_session].set_index('roi')
        b = one[one.session == pair.followup_session].set_index('roi')
        if not a.index.equals(b.index):
            raise ValueError('Regional identities differ across paired sessions')
        delta = b[['synergy_pc', 'redundancy_pc', 'rank_score']] - a[['synergy_pc', 'redundancy_pc', 'rank_score']]
        delta['network'] = a.network
        delta['subject_id'], delta['arm'], delta['dataset'], delta['comparison'] = pair.subject_id, pair.arm, pair.dataset, pair.comparison
        rows.append(delta.reset_index())
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def comparison_to_paper(sudmex, aachen, alpha=.05):
    # Published qualitative statements only. No psychedelic arrays or effect maps.
    statements = {
        'gateways_FPN': 'Classic psychedelics increased FPN gateway counts.',
        'gateways_DAN': 'Classic psychedelics increased DAN gateway counts (exploratory network family).',
        'mean_redundancy_pc': 'Mean workspace redundancy PC increased under classic psychedelics.',
    }
    rows = []
    for metric, statement in statements.items():
        row = {'metric': metric, 'Yousefi_published_finding': statement, 'source': PAPER_URL}
        for name, table in [('SUDMEX', sudmex), ('Aachen', aachen)]:
            if table.empty:
                row[name] = 'Pending complete cohort'
            else:
                result = table.set_index('metric').loc[metric]
                if not np.isfinite(result.estimate):
                    row[name] = 'Not estimable'
                    continue
                direction = 'higher' if result.estimate > 0 else ('lower' if result.estimate < 0 else 'zero estimated change')
                row[name] = f'{direction}; effect not established' if not np.isfinite(result.q_family) or result.q_family >= alpha else direction
        rows.append(row)
    return pd.DataFrame(rows)


def write_provenance(root, out, inputs, config):
    root, out = Path(root).resolve(), Path(out).resolve()
    allowed = root / 'output/tms_publication'
    if allowed not in out.parents:
        raise ValueError('Publication outputs must stay in output/tms_publication/<run>')
    out.mkdir(parents=True, exist_ok=True)
    paths = sorted({Path(p).resolve() for p in inputs}, key=str)
    record = {'executed_utc': datetime.now(timezone.utc).isoformat(), 'python': platform.python_version(),
              'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': __import__('scipy').__version__,
              'statsmodels': __import__('statsmodels').__version__, 'config': config,
              'inputs': [{'path': str(p.relative_to(root)), 'sha256': sha256(p)} for p in paths]}
    (out / 'provenance.json').write_text(json.dumps(record, indent=2), encoding='utf-8')
    return record
