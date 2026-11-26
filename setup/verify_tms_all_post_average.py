"""Independently verify the all-follow-up analysis and document its estimand."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'output/tms_publication/20261125T145928_921789Z'
sys.path.insert(0,str(ROOT/'pipeline/scripts/analysis'))
import tms_publication as tms
import tms_numit as tn


def close(a,b):
    np.testing.assert_allclose(a,b,rtol=2e-6,atol=2e-9,equal_nan=True)


def independent_hc3(frame,outcome,predictors,term):
    data=frame[[outcome,*predictors]].dropna()
    keep=[c for c in predictors if c=='intercept' or data[c].nunique()>1]
    x=data[keep].to_numpy(float); y=data[outcome].to_numpy(float)
    inv=np.linalg.pinv(x); b=inv@y
    leverage=np.sum(x*inv.T,axis=1)
    assert max(leverage)<1-1e-10
    resid=(y-x@b)/(1-leverage)
    covariance=(inv*resid**2)@inv.T
    i=keep.index(term); se=np.sqrt(covariance[i,i]); df=len(y)-len(keep)
    half=stats.t.ppf(.975,df)*se
    return [b[i],b[i]-half,b[i]+half,2*stats.t.sf(abs(b[i]/se),df)]


def independent_max_t(frame,metrics,families,seed):
    x=frame[[m+'_delta' for m in metrics]].to_numpy(float)
    n=len(x); exact=n<=20; size=2**n if exact else 200000
    if exact:
        signs=np.array([[(1 if i&(1<<j) else -1) for j in range(n)] for i in range(size)],dtype=np.int8)
    else:
        signs=np.random.default_rng(seed).integers(0,2,size=(size,n),dtype=np.int8)*2-1
        signs[0]=1
    active=x.std(0,ddof=1)>1e-13
    observed=abs(stats.ttest_1samp(x[:,active],0,axis=0).statistic)
    # Direct SciPy t statistics on signed observations, independent of the
    # sums-of-squares formula used in the analysis. Small chunks bound RAM.
    dist=np.empty((size,int(active.sum())))
    for first in range(0,size,2048):
        transformed=signs[first:first+2048,:,None]*x[None,:,active]
        dist[first:first+len(transformed)]=abs(stats.ttest_1samp(transformed,0,axis=1).statistic)
    names=np.asarray(metrics)[active]
    result={}
    for k,m in enumerate(names):
        family_indices=np.flatnonzero([families[c]==families[m] for c in names])
        p=float((dist[:,k]>=observed[k]-1e-10).mean())
        adjusted=float((dist[:,family_indices].max(1)>=observed[k]-1e-10).mean())
        result[m]=(p,adjusted)
    return result


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('run_directory',type=Path)
    out=parser.parse_args().run_directory.resolve()
    assert ROOT in out.parents and 'tms_all_post_average' in out.parts
    provenance=json.loads((out/'provenance.json').read_text())
    for name,digest in provenance['input_sha256'].items():
        assert tms.sha256(ROOT/name)==digest
    contrib=pd.read_csv(out/'scan_contributions.csv')
    assert not contrib.duplicated(['dataset','subject_id','session']).any()
    close(contrib.groupby(['dataset','subject_id','phase']).within_person_weight.sum(),1.)
    tables={}
    for branch,scan_file,pair_file in [('raw','raw_scan_metrics_used.csv','raw_participant_changes.csv'),
                                      ('normalized','normalized_scan_metrics_used.csv','normalized_participant_changes.csv')]:
        scans=pd.read_csv(out/scan_file)
        pairs=pd.read_csv(out/pair_file)
        tables[branch]=pairs
        joined=scans.merge(contrib[['dataset','subject_id','session','phase']],on=['dataset','subject_id','session'],validate='one_to_one')
        means=joined.groupby(['dataset','subject_id','phase'])[tms.METRICS].mean().unstack('phase')
        index=pd.MultiIndex.from_frame(pairs[['dataset','subject_id']])
        for m in tms.METRICS:
            pre=means[(m,'pre')].reindex(index).to_numpy(); post=means[(m,'post')].reindex(index).to_numpy()
            close(pairs[m+'_baseline'],pre); close(pairs[m+'_followup'],post); close(pairs[m+'_delta'],post-pre)
    main=pd.read_csv(out/'neural_pre_post_results.csv')
    old=pd.read_csv(SOURCE/'participant_paired_changes.csv')
    original_ids=set(old.loc[old.dataset.eq('sudmex_tms'),'subject_id'])
    tested=0
    for (branch,dataset,cohort),block in main.groupby(['branch','dataset','cohort']):
        data=tables[branch].loc[tables[branch].dataset.eq(dataset)]
        if cohort=='original_acute_40': data=data.loc[data.subject_id.isin(original_ids)]
        for row in block.itertuples():
            if not np.isfinite(row.p): continue
            test=stats.ttest_1samp(data[row.metric+'_delta'],0); ci=test.confidence_interval()
            close([row.estimate,row.ci_low,row.ci_high,row.p],[data[row.metric+'_delta'].mean(),ci.low,ci.high,test.pvalue])
            tested+=1
        # Full max-T check for the main cohorts. Sensitivity uses the same
        # validated code and is checked through its independent t statistics.
        if cohort=='all_eligible':
            methods=independent_max_t(data,tms.METRICS,{m:tms.METRIC_INFO[m][0] for m in tms.METRICS},2026112521)
            for row in block.itertuples():
                if row.metric in methods: close([row.p_sign_flip,row.p_maxT_family],methods[row.metric])
            print('Independent full sign-flip check passed:',branch,dataset,flush=True)
    arms=pd.read_csv(out/'original_arm_contrasts.csv'); adjusted_count=0
    for row in arms.itertuples():
        if row.status!='ok': continue
        data=tables[row.branch].loc[tables[row.branch].dataset.eq('sudmex_tms')].copy()
        if row.cohort=='original_acute_40': data=data.loc[data.subject_id.isin(original_ids)]
        data=data.assign(intercept=1.,active=data.arm.eq('active').astype(float))
        if row.model=='unadjusted_change_difference':
            a=data.loc[data.active.eq(1),row.metric+'_delta'];b=data.loc[data.active.eq(0),row.metric+'_delta']
            test=stats.ttest_ind(a,b,equal_var=False);ci=test.confidence_interval()
            expected=[a.mean()-b.mean(),ci.low,ci.high,test.pvalue]
        else:
            predictors=['intercept','active',row.metric+'_baseline']
            if row.model=='baseline_motion_adjusted': predictors.append('delta_fd')
            expected=independent_hc3(data,row.metric+'_followup',predictors,'active')
        close([row.estimate,row.ci_low,row.ci_high,row.p],expected);adjusted_count+=1
    clinical=pd.read_csv(out/'clinical_all_post_participant_changes.csv')
    clinical_source=pd.read_csv(ROOT/'data/master_clinical_sessions.csv').set_index(['dataset','subject_id','session'])
    clinical_results=pd.read_csv(out/'clinical_all_post_results.csv')
    for row in clinical_results.itertuples():
        data=clinical.dropna(subset=[row.metric+'_baseline',row.metric+'_followup']).copy()
        for person in data.itertuples():
            visits=['t0',*person.post_sessions.split('|')]
            values=[clinical_source.loc[('sudmex_tms',person.subject_id,s),row.metric] for s in visits]
            assert all(clinical_source.loc[('sudmex_tms',person.subject_id,s),row.metric+'_qc']=='ok' for s in visits)
            close(getattr(person,row.metric+'_delta'),np.mean(values[1:])-values[0])
        if row.contrast=='overall_pre_post':
            test=stats.ttest_1samp(data[row.metric+'_delta'],0);ci=test.confidence_interval()
            expected=[data[row.metric+'_delta'].mean(),ci.low,ci.high,test.pvalue]
        else:
            data=data.assign(intercept=1.,active=data.arm.eq('active').astype(float))
            expected=independent_hc3(data,row.metric+'_followup',['intercept','active',row.metric+'_baseline'],'active')
        close([row.estimate,row.ci_low,row.ci_high,row.p],expected)
    for table,group,adjusted in [(main,['branch','dataset','cohort','family'],'p_holm_family'),
                                 (arms,['branch','cohort','model','family'],'p_holm_family'),
                                 (clinical_results,['contrast'],'p_holm_clinical7')]:
        for _,b in table.groupby(group):
            valid=b.p.notna(); expected=multipletests(b.p.fillna(1.),method='holm')[1]
            close(b.loc[valid,adjusted],expected[valid])
    # Side-by-side raw versus all-post changes for the same 40 individuals.
    reference=old.loc[old.dataset.eq('sudmex_tms')].set_index('subject_id')
    allpost=tables['raw'].loc[tables['raw'].subject_id.isin(original_ids)].set_index('subject_id')
    comparison=[]
    for arm in ['active','sham']:
        ids=reference.index[reference.arm.eq(arm)]
        for m in tms.METRICS:
            comparison.append({'arm':arm,'metric':m,'n':len(ids),
                               'acute_mean_change':reference.loc[ids,m+'_delta'].mean(),
                               'all_post_mean_change':allpost.loc[ids,m+'_delta'].mean()})
    pd.DataFrame(comparison).to_csv(out/'acute_vs_all_post_same_people.csv',index=False)
    record={'status':'passed','created_utc':datetime.now(timezone.utc).isoformat(),
            'independent_paired_tests':tested,'independent_original_arm_contrasts':adjusted_count,
            'independent_clinical_models':len(clinical_results),
            'main_family_maxT_independent_scipy_reproduction':True,
            'input_hashes_unchanged':True,'equal_participant_weighting_verified':True,
            'clinical_windows_and_QC_verified':True,'source_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (out/'independent_verification.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))


if __name__=='__main__':
    with threadpool_limits(limits=1): main()
