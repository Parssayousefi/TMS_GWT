"""Post-hoc family-wise inference, retaining all original endpoints and tests."""
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'output/tms_publication/20261125T145928_921789Z'
ROUND1 = ROOT/'output/tms_exploration/20261125T154020_735179Z'
FOLLOWUP = ROOT/'output/tms_exploration/20261125T155441_285503Z'
NUMIT = ROOT/'output/tms_numit/20261125T155252_384543Z'
DIAG = ROOT/'output/tms_diagnostics/20261125T152110_809161Z'
sys.path.insert(0,str(ROOT/'pipeline/scripts/analysis'))
import tms_publication as tms
sys.path.insert(0,str(ROOT/'setup'))
from explore_tms_round1 import regression


def holm(values):
    p = np.asarray(values,dtype=float)
    result = np.full(p.shape,np.nan)
    good = np.isfinite(p)
    # An unestimable member remains in the declared family as p=1 for the
    # correction burden, but receives no reported adjusted p of its own.
    if good.any(): result[good] = multipletests(np.where(good,p,1.),method='holm')[1][good]
    return result


def exact_max_t(values):
    """All joint sign flips; single-step max-|t| simultaneous inference."""
    values = np.asarray(values,float)
    n,k = values.shape
    assert np.isfinite(values).all() and n<=20
    usable = values.std(axis=0,ddof=1)>1e-13
    data = values[:,usable]
    n_perm = 2**n
    observed = abs(data.mean(0)/(data.std(0,ddof=1)/np.sqrt(n)))
    ss = (data**2).sum(0)
    distribution = np.empty((n_perm,data.shape[1]))
    for first in range(0,n_perm,2048):
        labels = np.arange(first,min(first+2048,n_perm))
        signs = 2*((labels[:,None]>>np.arange(n))&1)-1
        means = signs@data/n
        variance = np.maximum((ss[None]-n*means**2)/(n-1),0)
        with np.errstate(divide='ignore',invalid='ignore'):
            distribution[first:first+len(labels)] = abs(means/np.sqrt(variance/n))
    raw = (distribution>=observed[None]-1e-10).mean(0)
    maximum = distribution.max(1)
    adjusted = (maximum[:,None]>=observed[None]-1e-10).mean(0)
    raw_full,adj_full = np.full(k,np.nan),np.full(k,np.nan)
    raw_full[usable],adj_full[usable] = raw,adjusted
    return raw_full,adj_full


def family(metric):
    if metric in tms.METRIC_INFO: return tms.METRIC_INFO[metric][0]
    if metric.startswith('network_'): return 'network_continuous_24'
    if metric.startswith('anatomy_'): return 'subcortical_14'
    raise ValueError(metric)


def check_max_t():
    rng = np.random.default_rng(412)
    x = rng.normal(size=(8,3)) + [.2,.1,-.1]
    p,pa = exact_max_t(x)
    all_t = []
    for label in range(2**len(x)):
        signs = np.array([1 if label&(1<<i) else -1 for i in range(len(x))])
        all_t.append(abs(stats.ttest_1samp(x*signs[:,None],0,axis=0).statistic))
    observed = abs(stats.ttest_1samp(x,0,axis=0).statistic)
    all_t = np.asarray(all_t)
    np.testing.assert_allclose(p,(all_t>=observed[None]-1e-10).mean(0))
    np.testing.assert_allclose(pa,(all_t.max(1)[:,None]>=observed[None]-1e-10).mean(0))
    _,reordered = exact_max_t(x[:,[2,0,1]])
    np.testing.assert_array_equal(reordered,pa[[2,0,1]])
    assert np.all(pa>=p)


def main():
    check_max_t()
    out = ROOT/'output/tms_familywise'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out.mkdir(parents=True,exist_ok=False)
    input_files = []
    def read(path):
        input_files.append(path)
        return pd.read_csv(path)
    result_tables, motion_rows = [], []
    for branch,location in [('raw',SOURCE),('numit',NUMIT)]:
        pairs = read(location/'participant_paired_changes.csv')
        aachen = pairs.loc[pairs.dataset.eq('aachen_itbs')].groupby('subject_id').mean(numeric_only=True)
        if branch=='raw':
            estimates = read(SOURCE/'aachen_neural_results.csv').copy()
            sud = read(ROUND1/'neural_exploration.csv')
            sud = sud.loc[sud.dataset.eq('sudmex_tms')&sud.metric.isin(tms.METRICS)].copy()
        else:
            all_results = read(NUMIT/'neural_contrasts.csv')
            estimates = all_results.loc[all_results.dataset.eq('aachen_itbs')&all_results.contrast.eq('pooled_change')].copy()
            sud = all_results.loc[all_results.dataset.eq('sudmex_tms')&all_results.contrast.eq('baseline_adjusted_active_sham')].copy()
        estimates['family'] = estimates.metric.map(family)
        estimates['p_holm_family'] = estimates.groupby('family').p.transform(holm)
        estimates['p_holm_all25'] = holm(estimates.p)
        estimates['p_sign_flip'] = np.nan
        estimates['p_maxT_family'] = np.nan
        for name,block in estimates.groupby('family'):
            values = aachen[[m+'_delta' for m in block.metric]].to_numpy()
            raw,adjusted = exact_max_t(values)
            estimates.loc[block.index,'p_sign_flip'] = raw
            estimates.loc[block.index,'p_maxT_family'] = adjusted
        _,all25 = exact_max_t(aachen[[m+'_delta' for m in estimates.metric]].to_numpy())
        estimates['p_maxT_all25'] = all25
        estimates['dataset'],estimates['contrast'],estimates['branch'] = 'aachen_itbs','pooled_change',branch
        for m in [m for m in tms.METRICS if m.startswith('global_')]:
            data = pd.DataFrame({'outcome':aachen[m+'_delta'],'intercept':1.,'delta_fd':aachen.delta_fd})
            motion_rows.append({'branch':branch,'metric':m,**regression(data,'outcome',['intercept','delta_fd'],'intercept')})
        sud['family'] = sud.metric.map(family)
        sud['p_holm_family'] = sud.groupby('family').p.transform(holm)
        sud['p_holm_all25'] = holm(sud.p)
        sud['branch'] = branch
        result_tables.extend([estimates,sud])
    main_results = pd.concat(result_tables,ignore_index=True)
    main_results.to_csv(out/'main_comparisons_familywise.csv',index=False)
    motion = pd.DataFrame(motion_rows)
    motion['p_holm_global5'] = motion.groupby('branch').p.transform(holm)
    motion.to_csv(out/'aachen_motion_global_familywise.csv',index=False)
    # Secondary exploration: all three nonpooled intensity contrasts belong
    # to the same question; do not split them to lower the correction burden.
    extended = read(ROUND1/'neural_exploration.csv')
    validity_path = ROOT/'output/tms_familywise/20261125T160334_991348Z/numerically_unestimable_contrasts.csv'
    if validity_path.is_file():
        invalid = read(validity_path)
        extended['reported_original_p'] = extended.p
        for row in invalid.loc[invalid.reference_seed.eq('raw')].itertuples():
            match = extended.dataset.eq(row.dataset)&extended.contrast.eq(row.contrast)&extended.metric.eq(row.metric)
            extended.loc[match,['p','ci_low','ci_high']] = np.nan
            extended.loc[match,'status'] = 'undefined_HC3_unit_leverage'
    extended['family'] = extended.metric.map(family)
    extended['contrast_group'] = extended.contrast.map({
        'baseline_adjusted_active_sham':'main', 'pooled_change':'main',
        '90_change':'intensity_followups','120_change':'intensity_followups',
        'intensity_120_minus_90':'intensity_followups',
        'intensity_adjusted_order_baseline_motion':'adjusted_intensity_sensitivity'})
    extended['p_holm_family'] = extended.groupby(['dataset','contrast_group','family']).p.transform(holm)
    extended.to_csv(out/'extended_exploration_familywise.csv',index=False)
    clinical = read(DIAG/'clinical_baseline_adjusted.csv')
    clinical['p_holm_clinical7'] = holm(clinical.p)
    clinical.to_csv(out/'clinical_familywise.csv',index=False)
    association = read(SOURCE/'clinical_neural_associations.csv')
    association['p_holm_association35'] = holm(association.p)
    association.to_csv(out/'clinical_neural_associations_familywise.csv',index=False)
    for filename in ['clinical_baseline_moderation.csv','multivariate_profiles.csv']:
        frame = read(ROUND1/filename)
        frame['p_holm_family'] = holm(frame.p)
        frame.to_csv(out/filename,index=False)
    delayed = read(FOLLOWUP/'delayed_active_results.csv')
    delayed['family'] = delayed.metric.map(family)
    delayed['p_holm_family_across_both_contrasts'] = delayed.groupby(['model','family']).p.transform(holm)
    delayed.to_csv(out/'delayed_active_familywise.csv',index=False)
    input_files += [ROOT/'PROVENANCE/TMS_FAMILYWISE_PLAN_20261125.json',Path(__file__)]
    record = {'status':'post_hoc_familywise_complete','created_utc':datetime.now(timezone.utc).isoformat(),
              'maxT_n_sign_flips':32768,'method':'Exact two-sided single-step max-absolute-t within pre-existing original metric families',
              'independent_maxT_scipy_exhaustive_check':'passed','post_hoc':True,
              'scope':'Family-wise, not study-wide; original and normalized branches are reported together, not selected by significance',
              'input_sha256':{str(p.relative_to(ROOT)):tms.sha256(p) for p in input_files}}
    (out/'provenance.json').write_text(json.dumps(record,indent=2)+'\n')
    print('OUTPUT:',out)
    print(main_results.loc[main_results.family.eq('global_phiid')&main_results.dataset.eq('aachen_itbs'),
                           ['branch','metric','estimate','p','p_holm_family','p_sign_flip','p_maxT_family','p_maxT_all25']].to_string(index=False))


if __name__=='__main__': main()
