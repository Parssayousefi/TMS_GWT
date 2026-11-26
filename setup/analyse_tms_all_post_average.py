"""Pre versus all available post visits, with explicit participant weighting."""
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'output/tms_publication/20261125T145928_921789Z'
NUMIT=ROOT/'output/tms_numit/20261125T155252_384543Z'
sys.path.insert(0,str(ROOT/'pipeline/scripts/analysis'))
sys.path.insert(0,str(ROOT/'setup'))
import tms_publication as tms
import tms_numit as tn
from analyse_tms_familywise import holm

KEY=['dataset','subject_id','session']
SEEDS=[2026112511,2026112512]
POST=['t1','t14','t2','t3','t4']
CLINICAL=['craving_vas','craving_ccq_general','craving_ccq_now','impulsivity_bis_total',
          'impulsivity_bis_cognitive','impulsivity_bis_motor','impulsivity_bis_nonplanning']


def save_json(path,value):
    path.write_text(json.dumps(value,indent=2,default=str)+'\n',encoding='utf-8')


def aggregate(scans):
    rows,contributions=[],[]
    assert not scans.duplicated(KEY).any()
    for (dataset,subject),frame in scans.groupby(['dataset','subject_id']):
        if dataset=='sudmex_tms':
            if 't0' not in set(frame.session): continue
            pre=frame.loc[frame.session.eq('t0')]
            post=frame.loc[frame.session.isin(POST)]
            if not len(post): continue
        elif dataset=='aachen_itbs':
            assert set(frame.session)=={'pre90','pre120','post90','post120'}
            pre=frame.loc[frame.session.str.startswith('pre')]
            post=frame.loc[frame.session.str.startswith('post')]
        else: raise ValueError('Only TMS inputs permitted')
        selected=pd.concat([pre,post])
        assert selected.arm.nunique()==1 and selected.valid_roi_signature.nunique()==1
        assert selected[tms.METRICS].notna().all().all()
        arm=selected.arm.iloc[0]
        item={'dataset':dataset,'subject_id':subject,'arm':arm,'n_pre':len(pre),'n_post':len(post),
              'pre_sessions':'|'.join(sorted(pre.session)),
              'post_sessions':'|'.join(sorted(post.session,key=lambda s: POST.index(s) if s in POST else s)),
              'mean_pre_fd':pre.mean_fd.mean(),'mean_post_fd':post.mean_fd.mean(),
              'delta_fd':post.mean_fd.mean()-pre.mean_fd.mean()}
        for metric in tms.METRICS:
            item[metric+'_baseline']=float(pre[metric].mean())
            item[metric+'_followup']=float(post[metric].mean())
            item[metric+'_delta']=item[metric+'_followup']-item[metric+'_baseline']
        rows.append(item)
        for phase,block in [('pre',pre),('post',post)]:
            for scan in block.itertuples():
                contributions.append({**{k:getattr(scan,k) for k in [*KEY,'arm','prefix','mean_fd','n_timepoints']},
                                      'phase':phase,'within_person_weight':1/len(block)})
    result=pd.DataFrame(rows)
    assert not result.duplicated(['dataset','subject_id']).any()
    return result,pd.DataFrame(contributions)


def extend_normalization(out,needed,hashes):
    old_path=NUMIT/'participant_scan_metrics.csv'
    hashes[str(old_path.relative_to(ROOT))]=tms.sha256(old_path)
    existing=pd.read_csv(old_path)
    existing_keys=set(existing.prefix)
    missing=needed.loc[~needed.prefix.isin(existing_keys)].copy()
    edges=np.column_stack(np.triu_indices(454,1))
    reference=tms.load_reference(ROOT)
    refs={}
    for dataset in missing.dataset.unique():
        for seed in SEEDS:
            folder=NUMIT/'references'/dataset/str(seed)
            paths=[folder/'tdmi_grid.npy',*[folder/(m+'.npy') for m in tms.MATRIX_NAMES]]
            for path in paths: hashes[str(path.relative_to(ROOT))]=tms.sha256(path)
            refs[(dataset,seed)]=(np.load(paths[0]),{m:np.load(folder/(m+'.npy'),mmap_mode='r') for m in tms.MATRIX_NAMES})
    new,quality=[],[]
    sample=np.random.default_rng(2026112520).choice(len(edges),64,replace=False)
    for i,scan in enumerate(missing.itertuples()):
        path=ROOT/'pipeline/results/timeseries'/(scan.prefix+'_timeseries_deconv.npy')
        hashes[str(path.relative_to(ROOT))]=tms.sha256(path)
        ts=np.load(path)
        assert ts.shape==(scan.n_timepoints,454) and np.isfinite(ts).all()
        sigma=tn.pg.lagged_sigma_from_timeseries(ts,1,edges)
        tdmi=tn.pg.tdmi_from_sigma4(sigma)
        assert np.isfinite(tdmi).all() and (tdmi>0).all()
        check=tn.stable_phiid_atoms(sigma[sample])
        qs={str(seed):{} for seed in SEEDS}; qs['combined']={}
        for m in tms.MATRIX_NAMES:
            path=ROOT/'pipeline/results/phiid_matrices'/(scan.prefix+'_'+m+'.npy')
            hashes[str(path.relative_to(ROOT))]=tms.sha256(path)
            raw=np.load(path)
            observed=raw[edges[:,0],edges[:,1]]
            np.testing.assert_allclose(check[m],observed[sample],rtol=1e-7,atol=2e-9)
            for seed in SEEDS:
                grid,tables=refs[(scan.dataset,seed)]
                qs[str(seed)][m]=tn.conditional_cdf(grid,tables[m],tdmi,observed)
            qs['combined'][m]=(qs[str(SEEDS[0])][m]+qs[str(SEEDS[1])][m])/2
        for seed,measures in qs.items():
            matrices={}
            for m,values in measures.items():
                matrix=np.full((454,454),np.nan)
                matrix[edges[:,0],edges[:,1]]=matrix[edges[:,1],edges[:,0]]=values
                matrices[m]=matrix
            metrics,_=tms.metric_values(matrices,reference)
            new.append({**{k:getattr(scan,k) for k in [*KEY,'arm','prefix','mean_fd','n_timepoints','n_valid_rois','valid_roi_signature']},
                        'reference_seed':seed,**metrics})
        quality.append({**{k:getattr(scan,k) for k in KEY},'tdmi_min':tdmi.min(),'tdmi_max':tdmi.max(),
                        **{'mean_seed_difference_'+m:float(np.mean(abs(qs[str(SEEDS[0])][m]-qs[str(SEEDS[1])][m]))) for m in tms.MATRIX_NAMES}})
        if (i+1)%10==0 or i+1==len(missing): print(f'Normalized additional follow-up scans: {i+1}/{len(missing)}',flush=True)
    new=pd.DataFrame(new)
    new.to_csv(out/'new_normalized_scan_metrics_all_seeds.csv',index=False)
    pd.DataFrame(quality).to_csv(out/'new_normalization_quality.csv',index=False)
    combined=pd.concat([existing,new.loc[new.reference_seed.eq('combined')]],ignore_index=True)
    combined=combined.loc[combined.prefix.isin(needed.prefix)].copy()
    assert len(combined)==len(needed) and not combined.duplicated(KEY).any()
    combined.to_csv(out/'normalized_scan_metrics_used.csv',index=False)
    return combined


def paired_family_tests(frame,metrics,families,seed):
    data=frame[[m+'_delta' for m in metrics]].to_numpy(float)
    n=len(data)
    usable=np.isfinite(data).all(0)&(data.std(axis=0,ddof=1)>1e-13)
    x=data[:,usable]
    used_metrics=np.array(metrics)[usable]
    observed=abs(x.mean(0)/(x.std(0,ddof=1)/np.sqrt(n)))
    ss=(x*x).sum(0)
    exact=n<=20
    total=2**n if exact else 200000
    rng=np.random.default_rng(seed)
    counts=np.zeros(len(x.T),dtype=int)
    family_counts=np.zeros(len(x.T),dtype=int)
    all_counts=np.zeros(len(x.T),dtype=int)
    groups={f:np.flatnonzero([families[m]==f for m in used_metrics]) for f in set(families.values())}
    for first in range(0,total,4096):
        batch=min(4096,total-first)
        if exact:
            signs=2*((np.arange(first,first+batch)[:,None]>>np.arange(n))&1)-1
        else:
            signs=rng.integers(0,2,size=(batch,n),dtype=np.int8)*2-1
            if first==0: signs[0]=1
        means=signs@x/n
        variance=np.maximum((ss[None]-n*means**2)/(n-1),0)
        with np.errstate(divide='ignore',invalid='ignore'):
            simulated=abs(means/np.sqrt(variance/n))
        simulated[(variance==0)&(means==0)]=0
        assert not np.isnan(simulated).any()
        counts+=(simulated>=observed[None]-1e-10).sum(0)
        all_counts+=(simulated.max(1)[:,None]>=observed[None]-1e-10).sum(0)
        for indices in groups.values():
            if len(indices):
                family_counts[indices]+=(simulated[:,indices].max(1)[:,None]>=observed[indices][None]-1e-10).sum(0)
    rows=[]; lookup={m:i for i,m in enumerate(used_metrics)}
    for m in metrics:
        values=frame[m+'_delta'].to_numpy(float)
        result=tms.paired_summary(values) if m in lookup else {'n':n,'estimate':values.mean(),'ci_low':np.nan,'ci_high':np.nan,'p':np.nan,'status':'constant'}
        k=lookup.get(m)
        rows.append({'metric':m,'family':families[m],**result,
                     'p_sign_flip':counts[k]/total if k is not None else np.nan,
                     'p_maxT_family':family_counts[k]/total if k is not None else np.nan,
                     'p_maxT_all_metrics':all_counts[k]/total if k is not None else np.nan,
                     'sign_flip_method':'exact' if exact else 'Monte_Carlo_including_identity',
                     'sign_vectors':total,'seed':seed})
    table=pd.DataFrame(rows)
    table['p_holm_family']=table.groupby('family').p.transform(holm)
    return table


def hc3(frame,outcome,predictors,term):
    data=frame[[outcome,*predictors]].dropna()
    keep=[c for c in predictors if c=='intercept' or data[c].nunique()>1]
    if term not in keep: return {'n':len(data),'status':'constant_predictor','p':np.nan}
    x=data[keep].astype(float)
    if np.linalg.matrix_rank(x)<x.shape[1] or len(data)<=x.shape[1]+2:
        return {'n':len(data),'status':'rank_or_sample_limit','p':np.nan}
    model=sm.OLS(data[outcome],x).fit()
    if np.max(model.get_influence().hat_matrix_diag)>=1-1e-10:
        return {'n':len(data),'status':'undefined_HC3_unit_leverage','p':np.nan}
    model=sm.OLS(data[outcome],x).fit(cov_type='HC3',use_t=True)
    low,high=model.conf_int().loc[term]
    return {'n':len(data),'estimate':model.params[term],'ci_low':low,'ci_high':high,
            'p':model.pvalues[term],'se':model.bse[term],'status':'ok'}


def arm_tests(pairs,metrics,families):
    pairs=pairs.assign(intercept=1.,active=pairs.arm.eq('active').astype(float))
    rows=[]
    for m in metrics:
        a=pairs.loc[pairs.active.eq(1),m+'_delta']; b=pairs.loc[pairs.active.eq(0),m+'_delta']
        rows.append({'model':'unadjusted_change_difference','metric':m,'family':families[m],**tms.welch_summary(a,b)})
        for model,covariates in [('baseline_adjusted',['intercept','active',m+'_baseline']),
                                  ('baseline_motion_adjusted',['intercept','active',m+'_baseline','delta_fd'])]:
            rows.append({'model':model,'metric':m,'family':families[m],
                         **hc3(pairs,m+'_followup',covariates,'active')})
    table=pd.DataFrame(rows)
    table['p_holm_family']=table.groupby(['model','family']).p.transform(holm)
    return table


def clinical_aggregate(neural,contributions,out):
    clinical=pd.read_csv(ROOT/'data/master_clinical_sessions.csv').set_index(KEY)
    assert not clinical.index.duplicated().any()
    rows,coverage=[],[]
    for person in neural.loc[neural.dataset.eq('sudmex_tms')].itertuples():
        item={k:getattr(person,k) for k in ['dataset','subject_id','arm','n_post','post_sessions','delta_fd']}
        sessions=person.post_sessions.split('|')
        for m in CLINICAL:
            values=[];bad=[]
            for s in ['t0',*sessions]:
                key=(person.dataset,person.subject_id,s)
                if key not in clinical.index:
                    bad.append(s+':missing_visit'); values.append(np.nan); continue
                row=clinical.loc[key]
                valid=row[m+'_qc']=='ok' and pd.notna(row[m])
                values.append(float(row[m]) if valid else np.nan)
                if not valid: bad.append(s+':'+str(row[m+'_qc']))
            complete=np.isfinite(values).all()
            item[m+'_baseline']=values[0]
            item[m+'_followup']=float(np.mean(values[1:])) if complete else np.nan
            item[m+'_delta']=item[m+'_followup']-item[m+'_baseline']
            coverage.append({'subject_id':person.subject_id,'metric':m,'n_post_required':len(sessions),
                             'complete_matched_window':bool(complete),'issues':'|'.join(bad)})
        rows.append(item)
    frame=pd.DataFrame(rows)
    frame.to_csv(out/'clinical_all_post_participant_changes.csv',index=False)
    pd.DataFrame(coverage).to_csv(out/'clinical_matched_visit_coverage.csv',index=False)
    rows=[]
    for m in CLINICAL:
        d=frame.dropna(subset=[m+'_baseline',m+'_followup']).copy()
        rows.append({'contrast':'overall_pre_post','metric':m,**tms.paired_summary(d[m+'_delta'])})
        d=d.assign(intercept=1.,active=d.arm.eq('active').astype(float))
        rows.append({'contrast':'original_arm_baseline_adjusted','metric':m,
                     **hc3(d,m+'_followup',['intercept','active',m+'_baseline'],'active')})
    results=pd.DataFrame(rows)
    results['p_holm_clinical7']=results.groupby('contrast').p.transform(holm)
    results.to_csv(out/'clinical_all_post_results.csv',index=False)


def main():
    start=time.perf_counter()
    out=ROOT/'output/tms_all_post_average'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out.mkdir(parents=True,exist_ok=False)
    print('OUTPUT:',out,flush=True)
    paths=[SOURCE/'participant_scan_metrics.csv',SOURCE/'participant_paired_changes.csv',
           ROOT/'data/master_clinical_sessions.csv',ROOT/'PROVENANCE/TMS_ALL_POST_AVERAGE_PLAN_20261125.json']
    hashes={str(p.relative_to(ROOT)):tms.sha256(p) for p in paths}
    scans=pd.read_csv(paths[0])
    raw,contributions=aggregate(scans)
    assert raw.groupby('dataset').size().to_dict()=={'aachen_itbs':15,'sudmex_tms':41}
    contributions.to_csv(out/'scan_contributions.csv',index=False)
    raw.to_csv(out/'raw_participant_changes.csv',index=False)
    needed=scans.loc[scans.prefix.isin(contributions.prefix)].copy()
    needed.to_csv(out/'raw_scan_metrics_used.csv',index=False)
    assert len(needed)==189
    normalized_scans=extend_normalization(out,needed,hashes)
    normalized,_=aggregate(normalized_scans)
    normalized.to_csv(out/'normalized_participant_changes.csv',index=False)
    old_pairs=pd.read_csv(paths[1])
    matched=set(old_pairs.loc[old_pairs.dataset.eq('sudmex_tms'),'subject_id'])
    families={m:tms.METRIC_INFO[m][0] for m in tms.METRICS}
    main_tables,arm_tables,descriptive=[],[],[]
    for branch,pairs in [('raw',raw),('normalized',normalized)]:
        for dataset,block in pairs.groupby('dataset'):
            result=paired_family_tests(block,tms.METRICS,families,2026112521)
            result['branch'],result['dataset'],result['cohort']=branch,dataset,'all_eligible'
            main_tables.append(result)
            if dataset=='sudmex_tms':
                result=paired_family_tests(block.loc[block.subject_id.isin(matched)],tms.METRICS,families,2026112522)
                result['branch'],result['dataset'],result['cohort']=branch,dataset,'original_acute_40'
                main_tables.append(result)
                for cohort,d in [('all_eligible',block),('original_acute_40',block.loc[block.subject_id.isin(matched)])]:
                    arm_result=arm_tests(d,tms.METRICS,families)
                    arm_result['branch'],arm_result['cohort']=branch,cohort
                    arm_tables.append(arm_result)
                for arm,d in block.groupby('arm'):
                    for m in tms.METRICS:
                        values=d[m+'_delta']
                        ci=stats.t.interval(.95,len(values)-1,loc=values.mean(),scale=stats.sem(values))
                        descriptive.append({'branch':branch,'arm':arm,'metric':m,'n':len(d),
                                            'mean_pre':d[m+'_baseline'].mean(),'mean_post':d[m+'_followup'].mean(),
                                            'mean_change':values.mean(),'ci_low':ci[0],'ci_high':ci[1],
                                            'interpretation':'Descriptive within-original-arm interval; no independent efficacy claim'})
            print('Completed family tests:',branch,dataset,flush=True)
    main_results=pd.concat(main_tables,ignore_index=True)
    arms=pd.concat(arm_tables,ignore_index=True)
    main_results.to_csv(out/'neural_pre_post_results.csv',index=False)
    arms.to_csv(out/'original_arm_contrasts.csv',index=False)
    pd.DataFrame(descriptive).to_csv(out/'within_original_arm_descriptive.csv',index=False)
    clinical_aggregate(raw,contributions,out)
    # Exact matched-window identity, independently of test implementation.
    aachen=raw.loc[raw.dataset.eq('aachen_itbs')].set_index('subject_id')
    original=old_pairs.loc[old_pairs.dataset.eq('aachen_itbs')].groupby('subject_id').mean(numeric_only=True)
    max_error=0.
    for m in tms.METRICS:
        max_error=max(max_error,float(abs(aachen[m+'_delta']-original[m+'_delta']).max()))
        np.testing.assert_allclose(aachen[m+'_delta'],original[m+'_delta'],atol=1e-12,rtol=1e-10)
    for relative,digest in hashes.items(): assert tms.sha256(ROOT/relative)==digest
    record={'status':'analysis_complete','created_utc':datetime.now(timezone.utc).isoformat(),
            'elapsed_seconds':time.perf_counter()-start,'subjects':raw.groupby('dataset').size().to_dict(),
            'contributing_scans':contributions.groupby(['dataset','phase']).size().to_dict(),
            'aachen_pooling_max_difference':max_error,'main_metrics':len(tms.METRICS),
            'input_sha256':hashes,'source_sha256':{str(p.relative_to(ROOT)):tms.sha256(p) for p in [Path(__file__),ROOT/'pipeline/scripts/analysis/tms_numit.py',ROOT/'pipeline/scripts/analysis/tms_publication.py']},
            'warnings':['SUDMEX original-sham follow-ups mix sham and subsequent active treatment',
                        'Visit horizons and number of follow-ups differ across individuals; no missing-data imputation',
                        'Whole-cohort pre/post is not a randomized treatment contrast',
                        'Normalized metrics inherit the previously documented null-model/finite-record limitations']}
    record['contributing_scans']={'/'.join(k):int(v) for k,v in record['contributing_scans'].items()}
    save_json(out/'provenance.json',record)
    print(main_results.query("cohort=='all_eligible' and metric in ['global_phi','global_phi_r','gateways_FPN']")[[
        'branch','dataset','metric','n','estimate','ci_low','ci_high','p','p_maxT_family','p_holm_family']].to_string(index=False))
    print('COMPLETE:',out,flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=1): main()
