"""Fixed-HCP-workspace audit and TMS-only scope sensitivity.

Production preprocessing, PhiID estimators, frozen reference masks, and historical
results are read-only. New summaries are derived from existing pairwise matrices.
"""
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'pipeline/scripts/analysis'), str(ROOT/'pipeline/scripts/utils'), str(ROOT/'setup')]
import tms_publication as tms
import tms_numit as tn
from analyse_tms_familywise import exact_max_t, holm, check_max_t
from analyse_tms_all_post_average import hc3, paired_family_tests
from clinical_data import merge_clinical_neural, METRICS as CLINICAL

RAW = ROOT/'output/tms_publication/20261125T145928_921789Z'
NUMIT = ROOT/'output/tms_numit/20261125T155252_384543Z'
ALLPOST = ROOT/'output/tms_all_post_average/20261125T162521_698094Z'
KEY = ['dataset', 'subject_id', 'session']
SEEDS = [2026112511, 2026112512]
SCOPES = ['whole_brain', 'within_workspace', 'outside_workspace', 'between']
MEASURES = list(tms.MATRIX_NAMES)
ROLE_METRICS = [m for m in tms.METRICS if not m.startswith('global_')]
INFO_METRICS = [s+'__'+m for s in SCOPES for m in MEASURES]
METRICS = INFO_METRICS + ROLE_METRICS
MAIN_METRICS = [s+'__'+m for s in SCOPES[:2] for m in MEASURES] + ROLE_METRICS
WS25 = ['within_workspace__'+m for m in MEASURES] + ROLE_METRICS
HASHES = {}


def track(path):
    path = Path(path)
    HASHES[str(path.relative_to(ROOT))] = tms.sha256(path)
    return path


def read(path):
    return pd.read_csv(track(path))


def save_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, default=lambda x: x.item() if isinstance(x, np.generic) else str(x))+'\n', encoding='utf-8')


def edge_scopes(workspace):
    workspace = np.asarray(workspace, bool)
    i, j = np.triu_indices(len(workspace), 1)
    masks = {'whole_brain': np.ones(len(i), bool),
             'within_workspace': workspace[i] & workspace[j],
             'outside_workspace': ~workspace[i] & ~workspace[j],
             'between': workspace[i] ^ workspace[j]}
    assert np.all(sum(masks[s].astype(int) for s in SCOPES[1:]) == 1)
    return np.column_stack((i, j)), masks


def scoped_means(values, masks):
    values = np.asarray(values, float)
    if not np.isfinite(values).all():
        raise ValueError('Incomplete edge coverage: do not silently change the sampled graph')
    return {s: float(values[mask].mean()) for s, mask in masks.items()}


def audit_reference(out):
    ref = tms.load_reference(ROOT)
    for path in ref['input_files']: track(path)
    track(ROOT/'pipeline/Global_Workspace_Complete_Pipeline_publication.ipynb')
    edges, masks = edge_scopes(ref['workspace'])
    assert {s: int(m.sum()) for s,m in masks.items()} == {
        'whole_brain': 102831, 'within_workspace': 27495,
        'outside_workspace': 23871, 'between': 51465}
    pd.DataFrame({'roi': np.arange(1,455), 'workspace': ref['workspace'],
                  'network': [tms.NETWORKS[m] for m in ref['modules']]}).to_csv(out/'fixed_hcp_workspace.csv', index=False)
    save_json(out/'analysis_specification.json', {
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'scope': 'TMS only; no frozen HCP or psychedelic reanalysis',
        'workspace': 'The unchanged 235-region HCP boolean mask; both endpoints must be in the mask',
        'edge_counts': {s: int(m.sum()) for s,m in masks.items()},
        'unchanged_GWT_outcomes': 'PCs use all 454 destinations, then rank and summarize only the fixed 235 regions',
        'new_outcome': 'Mean Phi/PhiR within workspace is an added summary, not an outcome reported in Yousefi',
        'Luppi_distinction': 'Luppi used edgewise NBS and workspace overlap, not this scalar mean',
        'aachen': 'Average the two paired post-pre differences per person; n=15',
        'sudmex': 'Randomized t0/t1 contrast with baseline and motion sensitivity; later follow-up averages separate',
        'families': 'Five information metrics within each scope; additionally 10 whole/workspace metrics, workspace25 and expanded30',
        'spatial_diagnostics': 'All four scopes jointly corrected over 20 endpoints; outside-minus-workspace directly tested over five',
        'normalization': 'Use the existing shared baseline-derived conditional references; do not relearn a prior after changing mask',
        'exploratory': True})
    return ref, edges, masks


def extract_raw(out, ref, edges, masks):
    inventory, paths = tms.inventory(ROOT)
    scans, regions, files = tms.extract_scan_metrics(ROOT, inventory, ref)
    for p in paths+files: track(p)
    old = read(RAW/'participant_scan_metrics.csv').set_index(KEY).sort_index()
    check = scans.set_index(KEY).sort_index()
    assert old.index.equals(check.index) and len(check) == 193
    np.testing.assert_allclose(check[tms.METRICS], old[tms.METRICS], rtol=1e-12, atol=1e-12)
    assert scans.n_valid_rois.eq(454).all() and scans.n_valid_workspace.eq(235).all()
    for row in scans.itertuples():
        for m in MEASURES:
            values = np.load(ROOT/'pipeline/results/phiid_matrices'/(row.prefix+'_'+m+'.npy'))[edges[:,0],edges[:,1]]
            for s,value in scoped_means(values,masks).items(): scans.loc[row.Index, s+'__'+m] = value
    scans.to_csv(out/'raw_scan_metrics.csv', index=False)
    regions.to_csv(out/'raw_regional_metrics.csv', index=False)
    print('Raw scan audit: 193 scans; all original 25 metrics reproduced.', flush=True)
    return scans


def extract_normalized(out, raw, edges, masks):
    existing = pd.concat([read(NUMIT/'scan_metrics_all_seeds.csv'),
                          read(ALLPOST/'new_normalized_scan_metrics_all_seeds.csv')], ignore_index=True)
    existing.reference_seed = existing.reference_seed.astype(str)
    assert not existing.duplicated(KEY+['reference_seed']).any()
    needed = set(read(ALLPOST/'normalized_scan_metrics_used.csv').prefix)
    existing = existing[existing.prefix.isin(needed)].copy()
    assert len(needed)==189 and len(existing)==567
    lookup = existing.set_index(['prefix','reference_seed'])
    references = {}
    for dataset in tms.DATASETS:
        for seed in SEEDS:
            folder = NUMIT/'references'/dataset/str(seed)
            grid = np.load(track(folder/'tdmi_grid.npy'))
            tables = {m: np.load(track(folder/(m+'.npy')), mmap_mode='r') for m in MEASURES}
            references[dataset,seed] = grid,tables
    rows=[]
    for n,scan in enumerate(raw.loc[raw.prefix.isin(needed)].itertuples(),1):
        cache = NUMIT/'cache'/(scan.prefix+'_tdmi.npy')
        if cache.is_file(): tdmi = np.load(track(cache))
        else:
            ts = np.load(track(ROOT/'pipeline/results/timeseries'/(scan.prefix+'_timeseries_deconv.npy')))
            tdmi = tn.pg.tdmi_from_sigma4(tn.pg.lagged_sigma_from_timeseries(ts,1,edges))
        added={str(seed):{} for seed in SEEDS}; added['combined']={}
        for m in MEASURES:
            values=np.load(ROOT/'pipeline/results/phiid_matrices'/(scan.prefix+'_'+m+'.npy'))[edges[:,0],edges[:,1]]
            qs=[]
            for seed in SEEDS:
                grid,tables=references[scan.dataset,seed]
                q=tn.conditional_cdf(grid,tables[m],tdmi,values); qs.append(q)
                added[str(seed)].update({s+'__'+m:v for s,v in scoped_means(q,masks).items()})
            added['combined'].update({s+'__'+m:v for s,v in scoped_means((qs[0]+qs[1])/2,masks).items()})
        for seed,new in added.items():
            item=lookup.loc[(scan.prefix,seed)].to_dict()
            for m in MEASURES:
                np.testing.assert_allclose(new['whole_brain__'+m],item['global_'+m],rtol=1e-12,atol=1e-12)
            rows.append({**item,'prefix':scan.prefix,'reference_seed':seed,**new})
        if n%20==0 or n==189: print(f'Normalized workspace extraction: {n}/189 scans',flush=True)
    frame=pd.DataFrame(rows)
    frame.to_csv(out/'normalized_scan_metrics_all_seeds.csv',index=False)
    frame[frame.reference_seed.eq('combined')].to_csv(out/'normalized_scan_metrics.csv',index=False)
    return frame


def make_pairs(scans, plan):
    rows=[]
    assert not scans.duplicated(KEY).any()
    for p in plan.itertuples():
        block=scans.loc[scans.dataset.eq(p.dataset)&scans.subject_id.eq(p.subject_id)].set_index('session')
        a,b=block.loc[p.baseline_session],block.loc[p.followup_session]
        assert a.valid_roi_signature==b.valid_roi_signature
        row={k:getattr(p,k) for k in ['dataset','subject_id','arm','comparison','baseline_session','followup_session']}
        row['delta_fd']=b.mean_fd-a.mean_fd
        for m in METRICS:
            row.update({m+'_baseline':a[m],m+'_followup':b[m],m+'_delta':b[m]-a[m]})
        rows.append(row)
    return pd.DataFrame(rows)


def aggregate_aachen(pairs):
    a=pairs.loc[pairs.dataset.eq('aachen_itbs')]
    assert a.groupby('subject_id').size().eq(2).all()
    return a.groupby('subject_id').mean(numeric_only=True).sort_index()


def test_aachen(frame, metrics=MAIN_METRICS):
    rows=[{'metric':m, 'family':m.split('__')[0] if '__' in m else tms.METRIC_INFO[m][0],
           **tms.paired_summary(frame[m+'_delta'])} for m in metrics]
    table=pd.DataFrame(rows)
    for _,block in table.groupby('family'):
        p,pc=exact_max_t(frame[[m+'_delta' for m in block.metric]])
        table.loc[block.index,'p_sign_flip']=p
        table.loc[block.index,'p_maxT_family']=pc
        table.loc[block.index,'p_holm_family']=holm(block.p)
    table['p_holm_all']=holm(table.p)
    table['p_maxT_all']=exact_max_t(frame[[m+'_delta' for m in metrics]])[1]
    info=[m for m in metrics if m.startswith(('whole_brain__','within_workspace__'))]
    table.loc[table.metric.isin(info),'p_maxT_whole_and_workspace10']=exact_max_t(frame[[m+'_delta' for m in info]])[1]
    table.loc[table.metric.isin(info),'p_holm_whole_and_workspace10']=holm(table.loc[table.metric.isin(info),'p'])
    if set(WS25).issubset(metrics):
        idx=[table.index[table.metric.eq(m)][0] for m in WS25]
        table.loc[idx,'p_maxT_workspace25']=exact_max_t(frame[[m+'_delta' for m in WS25]])[1]
    return table


def test_sudmex(pairs, metrics=MAIN_METRICS):
    frame=pairs.loc[pairs.dataset.eq('sudmex_tms')].copy()
    frame['active']=frame.arm.eq('active').astype(float);frame['intercept']=1.
    rows=[]
    for m in metrics:
        family=m.split('__')[0] if '__' in m else tms.METRIC_INFO[m][0]
        rows.append({'model':'unadjusted_change_difference','metric':m,'family':family,
                     **tms.welch_summary(frame.loc[frame.active.eq(1),m+'_delta'],frame.loc[frame.active.eq(0),m+'_delta'])})
        for model,covs in [('baseline_adjusted',['intercept','active',m+'_baseline']),
                           ('baseline_motion_adjusted',['intercept','active',m+'_baseline','delta_fd'])]:
            rows.append({'model':model,'metric':m,'family':family,**hc3(frame,m+'_followup',covs,'active')})
    result=pd.DataFrame(rows)
    result['p_holm_family']=result.groupby(['model','family']).p.transform(holm)
    result['p_holm_all30']=result.groupby('model').p.transform(holm)
    for model,group in result.groupby('model'):
        selected=group.metric.str.contains('__')
        result.loc[group.index[selected],'p_holm_whole_and_workspace10']=holm(group.loc[selected,'p'])
    return result


def robustness(out, pairs_by_branch):
    robust=[]; loo=[]; intensity=[]; spatial=[]
    order={}
    for line in track(ROOT/'data/aachen_itbs/derivatives/README.txt').read_text().splitlines():
        found=re.match(r'^\s*(\d+)\s+([1342/]+)\s+',line)
        if found: order['aachen_itbs_'+found[1].zfill(2)]=int(found[2].startswith('3'))
    for branch,pairs in pairs_by_branch.items():
        pooled=aggregate_aachen(pairs)
        a90=pairs.loc[pairs.comparison.eq('post90_minus_pre90')].set_index('subject_id').sort_index()
        a120=pairs.loc[pairs.comparison.eq('post120_minus_pre120')].set_index('subject_id').sort_index()
        for m in [s+'__'+v for s in SCOPES[:2] for v in MEASURES]:
            family=m.split('__')[0];x=pooled[m+'_delta'].to_numpy()
            frame=pd.DataFrame({'y':x,'intercept':1.,'delta_fd':pooled.delta_fd})
            robust.append({'branch':branch,'metric':m,'family':family,'check':'motion_HC3',**hc3(frame,'y',['intercept','delta_fd'],'intercept')})
            d90=a90[m+'_delta'];d120=a120[m+'_delta']
            for contrast,values in [('90_change',d90),('120_change',d120),('120_minus_90',d120-d90)]:
                intensity.append({'branch':branch,'metric':m,'family':family,'contrast':contrast,**tms.paired_summary(values)})
            frame=pd.DataFrame({'y':d120-d90,'intercept':1.,'order':np.array([order[s] for s in a90.index])-.5,
                                'baseline':a120[m+'_baseline']-a90[m+'_baseline'],'fd':a120.delta_fd-a90.delta_fd})
            intensity.append({'branch':branch,'metric':m,'family':family,'contrast':'120_minus_90_order_baseline_motion',
                              **hc3(frame,'y',['intercept','order','baseline','fd'],'intercept')})
            if m.endswith(('__phi','__phi_r')):
                rng=np.random.default_rng(2026112526)
                boots=x[rng.integers(len(x),size=(10000,len(x)))].mean(1)
                lo,hi=np.quantile(boots,[.025,.975])
                robust.extend([{'branch':branch,'metric':m,'family':family,'check':'bootstrap10000','estimate':x.mean(),'ci_low':lo,'ci_high':hi},
                               {'branch':branch,'metric':m,'family':family,'check':'wilcoxon','estimate':np.median(x),'p':stats.wilcoxon(x,method='exact').pvalue}])
                for k,subject in enumerate(pooled.index):
                    loo.append({'branch':branch,'metric':m,'omitted_subject':subject,**tms.paired_summary(np.delete(x,k))})
        d=pd.DataFrame({m+'_delta':pooled['outside_workspace__'+m+'_delta']-pooled['within_workspace__'+m+'_delta'] for m in MEASURES})
        p,pc=exact_max_t(d)
        for k,m in enumerate(MEASURES):
            spatial.append({'branch':branch,'metric':m,'contrast':'outside_minus_within_workspace',
                            **tms.paired_summary(d[m+'_delta']),'p_sign_flip':p[k],'p_maxT_5':pc[k]})
    robust=pd.DataFrame(robust)
    chosen=robust.check.eq('motion_HC3')
    robust.loc[chosen,'p_holm_scope5']=robust[chosen].groupby(['branch','family']).p.transform(holm)
    robust.loc[chosen,'p_holm_whole_and_workspace10']=robust[chosen].groupby('branch').p.transform(holm)
    robust.to_csv(out/'aachen_robustness.csv',index=False)
    pd.DataFrame(loo).to_csv(out/'aachen_leave_one_out.csv',index=False)
    intensity=pd.DataFrame(intensity)
    intensity['contrast_family']=np.where(intensity.contrast.str.endswith('motion'),'adjusted_intensity','three_intensity_contrasts')
    intensity['p_holm_family']=intensity.groupby(['branch','family','contrast_family']).p.transform(holm)
    intensity.to_csv(out/'aachen_intensity_results.csv',index=False)
    spatial=pd.DataFrame(spatial);spatial['p_holm5']=spatial.groupby('branch').p.transform(holm)
    spatial.to_csv(out/'aachen_spatial_interactions.csv',index=False)


def all_post(out, scans_by_branch):
    plan=read(ALLPOST/'raw_participant_changes.csv')
    rows=[]; arms=[];changes=[]
    families={m:'within_workspace5' for m in ['within_workspace__'+m for m in MEASURES]}
    for branch,scans in scans_by_branch.items():
        look=scans.set_index(KEY);records=[]
        for p in plan.loc[plan.dataset.eq('sudmex_tms')].itertuples():
            pre=look.loc[(p.dataset,p.subject_id,'t0')]
            post=look.loc[[(p.dataset,p.subject_id,s) for s in p.post_sessions.split('|')]]
            item={k:getattr(p,k) for k in ['dataset','subject_id','arm','n_pre','n_post','post_sessions','delta_fd']}
            assert post.valid_roi_signature.eq(pre.valid_roi_signature).all()
            for m in MAIN_METRICS: item.update({m+'_baseline':pre[m],m+'_followup':post[m].mean(),m+'_delta':post[m].mean()-pre[m]})
            records.append(item)
        frame=pd.DataFrame(records);frame['branch']=branch;changes.append(frame)
        res=paired_family_tests(frame,list(families),families,2026112527)
        res['branch']=branch;rows.append(res)
        arm=test_sudmex(frame);arm['branch']=branch;arms.append(arm)
    pd.concat(changes,ignore_index=True).to_csv(out/'sudmex_all_post_participant_changes.csv',index=False)
    pd.concat(rows,ignore_index=True).to_csv(out/'sudmex_all_post_workspace_results.csv',index=False)
    pd.concat(arms,ignore_index=True).to_csv(out/'sudmex_all_post_arm_results.csv',index=False)


def clinical(out,pairs_by_branch):
    track(ROOT/'data/master_clinical_effects.csv')
    results=[]
    neural=['gateways_FPN','gateways_DAN','mean_synergy_pc','mean_redundancy_pc',
            'within_workspace__phi_r','within_workspace__phi','whole_brain__phi_r','whole_brain__phi']
    for branch,pairs in pairs_by_branch.items():
        sud=pairs.loc[pairs.dataset.eq('sudmex_tms')].drop(columns=['baseline_session','followup_session'])
        frame=merge_clinical_neural(sud).assign(active=lambda d:d.arm.eq('active').astype(float),intercept=1.)
        for cm in CLINICAL:
            for nm in neural:
                d=frame.copy();d['neural_change_z']=(d[nm+'_delta']-d[nm+'_delta'].mean())/d[nm+'_delta'].std(ddof=1)
                result=hc3(d,cm+'_delta',['intercept','neural_change_z','active',cm+'_t0',nm+'_baseline','delta_fd'],'neural_change_z')
                results.append({'branch':branch,'clinical_metric':cm,'neural_metric':nm,**result})
    result=pd.DataFrame(results);result['p_holm56']=result.groupby('branch').p.transform(holm)
    # Matching the original five-neural-outcome question, replacing only its PhiR scope.
    selected=result.neural_metric.isin(neural[:5])
    result.loc[selected,'p_holm_workspace35']=result[selected].groupby('branch').p.transform(holm)
    result.to_csv(out/'clinical_neural_associations.csv',index=False)


def main():
    started=time.perf_counter();check_max_t()
    out=ROOT/'output/tms_workspace_scope'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out.mkdir(parents=True);print('OUTPUT:',out,flush=True)
    ref,edges,masks=audit_reference(out)
    raw=extract_raw(out,ref,edges,masks)
    normalized=extract_normalized(out,raw,edges,masks)
    plan=read(RAW/'participant_paired_changes.csv')
    main_tables=[];sud_tables=[];diagnostic=[];seed_tables=[]
    pairs_by_branch={};scans_by_branch={'raw':raw,'normalized':normalized[normalized.reference_seed.eq('combined')]}
    for branch,scans in scans_by_branch.items():
        pairs=make_pairs(scans,plan);pairs.to_csv(out/(branch+'_paired_changes.csv'),index=False)
        pairs_by_branch[branch]=pairs
        pooled=aggregate_aachen(pairs);pooled.to_csv(out/(branch+'_aachen_pooled_changes.csv'))
        result=test_aachen(pooled);result['branch']=branch;main_tables.append(result)
        result=test_sudmex(pairs);result['branch']=branch;sud_tables.append(result)
        result=test_aachen(pooled,INFO_METRICS);result['branch']=branch;diagnostic.append(result)
    pd.concat(main_tables,ignore_index=True).to_csv(out/'aachen_main_results.csv',index=False)
    pd.concat(sud_tables,ignore_index=True).to_csv(out/'sudmex_randomized_results.csv',index=False)
    pd.concat(diagnostic,ignore_index=True).to_csv(out/'aachen_all_spatial20_results.csv',index=False)
    for seed,scans in normalized.groupby('reference_seed'):
        if seed=='combined':continue
        result=test_aachen(aggregate_aachen(make_pairs(scans,plan)),['within_workspace__'+m for m in MEASURES])
        result['reference_seed']=seed;seed_tables.append(result)
    pd.concat(seed_tables,ignore_index=True).to_csv(out/'normalization_seed_stability.csv',index=False)
    robustness(out,pairs_by_branch);all_post(out,scans_by_branch);clinical(out,pairs_by_branch)
    for path in [Path(__file__),ROOT/'setup/analyse_tms_familywise.py',ROOT/'setup/analyse_tms_all_post_average.py',
                 ROOT/'pipeline/scripts/analysis/tms_publication.py',ROOT/'pipeline/scripts/analysis/tms_numit.py']:track(path)
    for relative,digest in HASHES.items():assert tms.sha256(ROOT/relative)==digest
    save_json(out/'provenance.json',{'status':'complete','created_utc':datetime.now(timezone.utc).isoformat(),
                                   'elapsed_seconds':time.perf_counter()-started,'input_sha256':HASHES,
                                   'raw_scans':193,'normalized_scans':189,'aachen_participants':15,'sudmex_randomized_participants':40,
                                   'sudmex_all_followup_participants':41,'mask_unchanged':True,'old_metrics_reproduced':True,
                                   'raw_units':'nats','normalized_units':'conditional CDF units; multiply by 100 for percentile points',
                                   'uncertainty':'All confidence intervals marginal; normalized reference treated as fixed; comparisons exploratory'})
    (ROOT/'output/tms_workspace_scope/LATEST.txt').write_text(out.name+'\n')
    print(pd.concat(main_tables).query("metric.str.startswith('within_workspace')",engine='python')[
        ['branch','metric','estimate','ci_low','ci_high','p','p_maxT_family','p_holm_family','p_maxT_all']].to_string(index=False),flush=True)
    print('COMPLETE:',out,flush=True)
    return out


if __name__=='__main__':
    with threadpool_limits(limits=1):main()
