"""Post-hoc baseline/motion diagnostics; completed primary analysis stays intact."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'output/tms_publication/20261125T145928_921789Z'
MODELS = {
    'unadjusted_HC3': [],
    'baseline_only': ['fpn_baseline'],
    'motion_change_only': ['fd_change'],
    'baseline_and_motion_change': ['fpn_baseline', 'fd_change'],
    'baseline_and_both_motion_levels': ['fpn_baseline', 'fd_baseline', 'fd_change'],
}


def fit(frame, covariates, outcome='fpn_change'):
    columns = ['active', *covariates]
    data = frame[[outcome, *columns]].replace([np.inf, -np.inf], np.nan).dropna()
    x = sm.add_constant(data[columns].astype(float), has_constant='add')
    assert np.linalg.matrix_rank(x) == x.shape[1]
    model = sm.OLS(data[outcome].astype(float), x).fit(cov_type='HC3', use_t=True)
    low, high = model.conf_int().loc['active']
    result = {'n': len(data), 'estimate': model.params['active'], 'se_HC3': model.bse['active'],
              'ci_low': low, 'ci_high': high, 'p': model.pvalues['active'],
              'residual_df': model.df_resid, 'r_squared': model.rsquared,
              'formula': outcome + ' ~ active' + ''.join(' + ' + c for c in covariates)}
    return result, model


def bootstrap(frame, count, seed):
    rng = np.random.default_rng(seed)
    indices = [np.flatnonzero(frame.active.to_numpy() == arm) for arm in [0, 1]]
    models = list(MODELS.items())[:4]
    designs = {name: np.column_stack([np.ones(len(frame)), frame[['active', *cov]].to_numpy(float)])
               for name, cov in models}
    y = frame.fpn_change.to_numpy(float)
    estimates = np.empty((count, len(models)))
    for i in range(count):
        sampled = np.concatenate([rng.choice(idx, size=len(idx), replace=True) for idx in indices])
        for j, (name, _) in enumerate(models):
            estimates[i, j] = np.linalg.lstsq(designs[name][sampled], y[sampled], rcond=None)[0][1]
    samples = pd.DataFrame(estimates, columns=[name for name, _ in models])
    rows = []
    for name in samples:
        low, high = samples[name].quantile([.025, .975])
        rows.append({'model': name, 'replicates': count, 'seed': seed,
                     'bootstrap_ci_low': low, 'bootstrap_ci_high': high})
    return samples, pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bootstrap', type=int, default=5000)
    args = parser.parse_args()
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out = ROOT / 'output/tms_diagnostics' / stamp
    out.mkdir(parents=True, exist_ok=False)
    paths = [SOURCE / 'participant_paired_changes.csv', SOURCE / 'participant_scan_metrics.csv',
             SOURCE / 'clinical_neural_participant_changes.csv', SOURCE / 'primary_FPN_result.csv',
             SOURCE / 'motion_baseline_adjusted_neural_effects.csv',
             ROOT / 'PROVENANCE/SUDMEX_FPN_DIAGNOSTIC_PLAN_20261125.json']
    hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    pairs = pd.read_csv(paths[0])
    frame = pairs.loc[pairs.dataset.eq('sudmex_tms')].copy().reset_index(drop=True)
    assert len(frame) == 40 and frame.pair_ready.all() and frame.same_valid_rois.all()
    assert frame.arm.value_counts().to_dict() == {'active': 22, 'sham': 18}
    assert not frame.subject_id.duplicated().any()
    frame['active'] = frame.arm.eq('active').astype(int)
    frame = frame.rename(columns={'gateways_FPN_baseline': 'fpn_baseline',
                                  'gateways_FPN_followup': 'fpn_followup',
                                  'gateways_FPN_delta': 'fpn_change', 'delta_fd': 'fd_change'})
    scans = pd.read_csv(paths[1]).query("dataset == 'sudmex_tms'")
    for session, tag in [('t0', 'baseline'), ('t1', 'followup')]:
        scan = scans.loc[scans.session.eq(session)].set_index('subject_id')
        frame['fd_' + tag] = frame.subject_id.map(scan.mean_fd)
        frame['timepoints_' + tag] = frame.subject_id.map(scan.n_timepoints)
    assert frame[['fpn_baseline', 'fpn_followup', 'fd_baseline', 'fd_followup']].notna().all().all()
    np.testing.assert_allclose(frame.fpn_change, frame.fpn_followup-frame.fpn_baseline)
    np.testing.assert_allclose(frame.fd_change, frame.fd_followup-frame.fd_baseline)
    frame.to_csv(out/'participants.csv', index=False)

    descriptors = []
    for variable in ['fpn_baseline', 'fpn_followup', 'fpn_change', 'fd_baseline', 'fd_followup', 'fd_change']:
        a, b = [frame.loc[frame.arm.eq(arm), variable] for arm in ['active', 'sham']]
        pooled = np.sqrt(((len(a)-1)*a.var() + (len(b)-1)*b.var())/(len(a)+len(b)-2))
        descriptors.append({'variable': variable, 'active_mean': a.mean(), 'active_sd': a.std(),
                            'sham_mean': b.mean(), 'sham_sd': b.std(), 'active_minus_sham': a.mean()-b.mean(),
                            'standardized_difference': (a.mean()-b.mean())/pooled})
    pd.DataFrame(descriptors).to_csv(out/'arm_descriptives.csv', index=False)

    rows, coefficients, models = [], [], {}
    for name, covariates in MODELS.items():
        result, model = fit(frame, covariates)
        rows.append({'model': name, **result})
        models[name] = model
        for term in model.params.index:
            lo, hi = model.conf_int().loc[term]
            coefficients.append({'model': name, 'term': term, 'estimate': model.params[term],
                                 'ci_low': lo, 'ci_high': hi, 'p': model.pvalues[term]})
    table = pd.DataFrame(rows)
    table.to_csv(out/'fpn_model_comparison.csv', index=False)
    pd.DataFrame(coefficients).to_csv(out/'model_coefficients.csv', index=False)
    original = pd.read_csv(paths[3]).iloc[0]
    original_adjusted = pd.read_csv(paths[4]).query("metric == 'gateways_FPN'").iloc[0]
    joint = table.set_index('model').loc['baseline_and_motion_change']
    for c in ['estimate', 'ci_low', 'ci_high', 'p']:
        np.testing.assert_allclose(joint[c], original_adjusted[c], rtol=1e-9, atol=1e-10)
    ancova, _ = fit(frame, ['fpn_baseline', 'fd_change'], outcome='fpn_followup')
    for c in ['estimate', 'ci_low', 'ci_high', 'p']:
        np.testing.assert_allclose(joint[c], ancova[c], rtol=1e-9, atol=1e-10)

    decomposition = []
    for variable in ['fpn_baseline', 'fd_change']:
        imbalance = frame.loc[frame.active.eq(1), variable].mean()-frame.loc[frame.active.eq(0), variable].mean()
        slope = models['baseline_and_motion_change'].params[variable]
        decomposition.append({'covariate': variable, 'active_minus_sham': imbalance,
                              'joint_model_slope': slope, 'contribution_to_raw_minus_adjusted': imbalance*slope})
    difference = table.iloc[0].estimate-joint.estimate
    np.testing.assert_allclose(sum(r['contribution_to_raw_minus_adjusted'] for r in decomposition), difference)
    pd.DataFrame(decomposition).to_csv(out/'coefficient_attenuation_decomposition.csv', index=False)

    correlations = []
    for arm, block in [('all', frame), *list(frame.groupby('arm'))]:
        for x, y in [('fpn_baseline', 'fpn_change'), ('fpn_baseline', 'fpn_followup'),
                     ('fd_change', 'fpn_change'), ('fd_baseline', 'fpn_baseline')]:
            r, p = stats.pearsonr(block[x], block[y])
            correlations.append({'arm': arm, 'n': len(block), 'x': x, 'y': y, 'r': r, 'p_descriptive': p})
    pd.DataFrame(correlations).to_csv(out/'descriptive_correlations.csv', index=False)

    loo = []
    for i, participant in frame.iterrows():
        sample = frame.drop(index=i)
        for name, cov in MODELS.items():
            result, _ = fit(sample, cov)
            loo.append({'omitted_subject': participant.subject_id, 'arm': participant.arm, 'model': name, **result})
        a = sample.loc[sample.active.eq(1), 'fpn_change']
        b = sample.loc[sample.active.eq(0), 'fpn_change']
        loo.append({'omitted_subject': participant.subject_id, 'arm': participant.arm,
                    'model': 'original_Welch', 'n': len(sample), 'estimate': a.mean()-b.mean(),
                    'p': stats.ttest_ind(a, b, equal_var=False).pvalue})
    loo = pd.DataFrame(loo)
    loo.to_csv(out/'leave_one_out.csv', index=False)
    loo_summary = loo.groupby('model').agg(estimate_min=('estimate', 'min'), estimate_max=('estimate', 'max'),
                                         p_min=('p', 'min'), p_max=('p', 'max'),
                                         count_p_below_05=('p', lambda s: int((s<.05).sum())))
    loo_summary.to_csv(out/'leave_one_out_summary.csv')
    influence = models['baseline_and_motion_change'].get_influence()
    pd.DataFrame({'subject_id': frame.subject_id, 'arm': frame.arm,
                  'leverage': influence.hat_matrix_diag, 'cooks_distance': influence.cooks_distance[0],
                  'studentized_residual': influence.resid_studentized_external}).to_csv(out/'influence.csv', index=False)
    samples, bootstrap_table = bootstrap(frame, args.bootstrap, 2026112502)
    samples.to_csv(out/'bootstrap_coefficients.csv', index=False)
    bootstrap_table.to_csv(out/'bootstrap_intervals.csv', index=False)

    clinical = pd.read_csv(paths[2])
    clinical['active'] = clinical.arm.eq('active').astype(int)
    clinical_rows = []
    names = ['craving_vas', 'craving_ccq_general', 'craving_ccq_now', 'impulsivity_bis_total',
             'impulsivity_bis_cognitive', 'impulsivity_bis_motor', 'impulsivity_bis_nonplanning']
    for metric in names:
        result, _ = fit(clinical, [metric+'_t0'], outcome=metric+'_delta')
        clinical_rows.append({'metric': metric, **result})
    clinical_table = pd.DataFrame(clinical_rows)
    clinical_table['q_family'] = multipletests(clinical_table.p, method='fdr_bh')[1]
    clinical_table.to_csv(out/'clinical_baseline_adjusted.csv', index=False)

    colors = {'active': '#6A3D9A', 'sham': '#7F8791'}
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 3, figsize=(13.3, 4.2))
    for arm, block in frame.groupby('arm'):
        for row in block.itertuples():
            axes[0].plot([0, 1], [row.fpn_baseline, row.fpn_followup], color=colors[arm], alpha=.15, lw=.8)
        means = [block.fpn_baseline.mean(), block.fpn_followup.mean()]
        axes[0].plot([0, 1], means, 'o-', color=colors[arm], lw=3, label=arm)
        axes[1].scatter(block.fpn_baseline, block.fpn_followup, color=colors[arm], alpha=.75, label=arm)
        axes[2].scatter(block.fd_change, block.fpn_change, color=colors[arm], alpha=.75)
    axes[0].set_xticks([0, 1]); axes[0].set_xticklabels(['Baseline', 'Two weeks'])
    axes[0].set(ylabel='FPN gateway count', title='Individual visits and arm means')
    axes[0].legend(frameon=False)
    axes[1].plot([0, 30], [0, 30], '--', color='#B5BCC6', lw=1)
    axes[1].set(xlabel='Baseline FPN gateways', ylabel='Follow-up FPN gateways', title='Baseline and follow-up')
    axes[2].axhline(0, color='#B5BCC6', lw=1); axes[2].axvline(0, color='#B5BCC6', lw=1)
    axes[2].set(xlabel='Change in mean FD (mm)', ylabel='Change in FPN gateways', title='Motion change and gateway change')
    fig.tight_layout()
    for ext in ['png', 'pdf']: fig.savefig(out/('baseline_motion_diagnostics.'+ext), dpi=180, bbox_inches='tight')
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4.2))
    ax.errorbar(table.estimate, np.arange(len(table)),
                xerr=np.array([table.estimate-table.ci_low, table.ci_high-table.estimate]),
                fmt='o', capsize=3, color='#6A3D9A')
    ax.axvline(0, color='#B5BCC6', lw=1)
    ax.set_yticks(np.arange(len(table)))
    ax.set_yticklabels(['Unadjusted', 'Baseline FPN only', 'Motion change only',
                        'Baseline FPN + motion change', 'Baseline FPN + both motion levels'])
    ax.invert_yaxis()
    ax.set(xlabel='Active minus sham effect (gateway regions; HC3 95% CI)',
           title='Same participants and uncertainty method in every model')
    fig.tight_layout()
    for ext in ['png', 'pdf']: fig.savefig(out/('adjustment_comparison.'+ext), dpi=180, bbox_inches='tight')
    plt.close(fig)

    for relative, digest in hashes.items():
        assert hashlib.sha256((ROOT/relative).read_bytes()).hexdigest() == digest
    record = {'created_utc': datetime.now(timezone.utc).isoformat(), 'type': 'post_hoc_diagnostic',
              'source_run': str(SOURCE), 'n': len(frame), 'original_welch_p': original.p,
              'original_welch_estimate': original.estimate, 'input_sha256': hashes,
              'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'checks': ['Original adjusted result reproduced', 'ANCOVA/change coefficient equivalence',
                         'Exact coefficient-attenuation identity', 'No source files changed'],
              'no_new_exclusions': True, 'bootstrap_replicates': args.bootstrap}
    (out/'provenance.json').write_text(json.dumps(record, indent=2)+'\n')
    print('OUTPUT:', out)
    print(table[['model','estimate','ci_low','ci_high','p']].to_string(index=False))
    print('\nAlgebraic attenuation (not causal):\n', pd.DataFrame(decomposition).to_string(index=False))
    print('\nLeave-one-out:\n', loo_summary.to_string())
    print('\nClinical baseline adjustment:\n', clinical_table[['metric','n','estimate','ci_low','ci_high','p','q_family']].to_string(index=False))


if __name__ == '__main__':
    main()
