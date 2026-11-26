"""Test shared spatial patterns of CHANGE, using external cohort templates.

Six PhiID modes; one cohort supplies a mean-centered regional effect template
for testing in the other. No stability or haemodynamic analysis. Participant
tests, rather than parcel-wise p values, retain the proper sampling unit.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import sys
import hashlib

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'setup'),str(ROOT/'pipeline/scripts/analysis')]
import explore_pid_processing_modes as mode_run
from analyse_tms_familywise import exact_max_t

SOURCE=ROOT/'output/pid_processing_modes/20261126T130022_906563Z'
MODES=list(mode_run.MODES)
HASHES={}


def track(p):
    HASHES[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    return p


def save_json(p,obj):
    p.write_text(json.dumps(obj,indent=2,default=str)+'\n',encoding='utf-8')


def paired_fit(x):
    n=len(x); se=stats.sem(x); mean=x.mean(); ci=mean+np.array([-1,1])*stats.t.ppf(.975,n-1)*se
    return dict(n=n,estimate=mean,ci_low=ci[0],ci_high=ci[1],p_t=stats.ttest_1samp(x,0).pvalue)


def hc3(y,baseline,active,fd,adjust=True):
    x=pd.DataFrame({'constant':1.,'active':np.asarray(active,dtype=float),'baseline':baseline})
    if adjust:x['delta_fd']=fd
    fit=sm.OLS(y,x).fit(cov_type='HC3',use_t=True)
    # Direct sandwich verification, distinct from statsmodels.
    a=x.to_numpy();pinv=np.linalg.pinv(a);beta=pinv@y;res=y-a@beta
    lev=np.sum(a*pinv.T,axis=1);cov=(pinv*(res/(1-lev))**2)@pinv.T
    np.testing.assert_allclose(beta,fit.params,atol=1e-12)
    np.testing.assert_allclose(cov,fit.cov_params(),atol=1e-12)
    ci=fit.conf_int().loc['active']
    return dict(n=len(y),estimate=fit.params.active,ci_low=ci.iloc[0],ci_high=ci.iloc[1],p=fit.pvalues.active)


def main():
    out=ROOT/'output/pid_change_patterns'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out.mkdir(parents=True);print('OUTPUT:',out,flush=True)
    spec={'created_utc':datetime.now(timezone.utc).isoformat(),'status':'Exploratory analysis, not preregistered; specified before regional outcome extraction.',
          'question':'Do two independent cohorts share regional redistribution of PhiID processing modes after stimulation?',
          'modes':mode_run.MODES,'primary_features':'All 454 regions; each feature is the mean of the mode over its 453 incident pairs.',
          'primary_templates':'Aachen mean of participant-averaged two-visit changes; SUDMEX baseline-plus-motion-adjusted active coefficient at each region. Center template spatially and scale absolute weights to sum to one.',
          'validation':'Use only the OTHER dataset to construct each template. Project individual pre/post regional vectors onto that external template. Aachen paired sign flips; SUDMEX followup projected score ~ active + baseline projected score + deltaFD, HC3/t.',
          'inference':'Two-sided cohort p values; replication requires BOTH projected effects positive. Conjunction p=max(two cohort p), or 1 if either is nonpositive. Holm6 for all six modes. Aachen also exact joint maxT6. No parcel-independence assumption.',
          'sensitivities':'Cortical 400-node feature subset, still all-destination node means; Holm6 within branch and Holm12 across both branches. Baseline-only SUDMEX scalar model. Every omission in the full-region primary, external templates kept fixed.',
          'limits':'Spatial redistribution, not global increase; no inference of stability or causal direction. No pooling of the uncontrolled Aachen contrast with the randomized contrast. External templates avoid within-cohort outcome selection but this follow-up is exploratory.'}
    save_json(out/'analysis_specification.json',spec)
    track(Path(__file__));track(Path(mode_run.__file__));track(Path(mode_run.tn.pg.__file__))
    prior=pd.read_csv(track(SOURCE/'scan_modes_and_atoms.csv')).set_index(['subject_id','session'])
    scans=pd.read_csv(track(mode_run.CORE/'raw_scan_metrics.csv'))
    scans=scans.loc[[(r.subject_id,r.session) in prior.index for r in scans.itertuples()]].reset_index(drop=True)
    assert len(scans)==140
    plan=pd.read_csv(track(mode_run.CORE/'raw_paired_changes.csv'))
    edges=np.array(np.triu_indices(454,k=1)).T;i,j=edges.T
    nodes={};checks=[]
    for k,r in enumerate(scans.itertuples()):
        base=mode_run.CORRECTED if r.dataset=='aachen_itbs' else ROOT
        ts=np.load(track(base/'pipeline/results/timeseries'/f'{r.prefix}_timeseries_deconv.npy'))
        values=mode_run.modes(mode_run.tn.pg.phiid_atoms(mode_run.tn.pg.lagged_sigma_from_timeseries(ts,1,edges)))
        arr=np.stack([(np.bincount(i,weights=values[m],minlength=454)+np.bincount(j,weights=values[m],minlength=454))/453 for m in MODES])
        expected=prior.loc[(r.subject_id,r.session),['whole_brain__'+m for m in MODES]].to_numpy(float)
        np.testing.assert_allclose(arr.mean(1),expected,atol=1e-12)
        for roi in [0,225,453]:
            direct=np.array([values[m][(i==roi)|(j==roi)].mean() for m in MODES])
            np.testing.assert_allclose(arr[:,roi],direct,atol=1e-12)
        nodes[(r.subject_id,r.session)]=arr
        checks.append({'prefix':r.prefix,'max_scalar_difference':float(abs(arr.mean(1)-expected).max())})
        if (k+1)%35==0:print(f'Regional extraction {k+1}/140',flush=True)
    pp={}
    for dataset in ['aachen_itbs','sudmex_tms']:
        group=plan.loc[plan.dataset.eq(dataset)]
        ids=sorted(group.subject_id.unique());pre=[];post=[];fd=[];arms=[]
        for subject in ids:
            rows=group.loc[group.subject_id.eq(subject)]
            pre.append(np.mean([nodes[(subject,s)] for s in rows.baseline_session],axis=0))
            post.append(np.mean([nodes[(subject,s)] for s in rows.followup_session],axis=0))
            fd.append(rows.delta_fd.mean());arms.append(rows.arm.iloc[0])
        pp[dataset]={'pre':np.stack(pre),'post':np.stack(post),'fd':np.array(fd),'active':np.array(arms)=='active','ids':ids}
    aa,sud=pp['aachen_itbs'],pp['sudmex_tms']
    np.savez_compressed(out/'regional_modes.npz',**{f'{d}_{field}':v[field] for d,v in pp.items() for field in ['pre','post','fd','active']})
    np.testing.assert_equal(aa['pre'].shape,(15,6,454));np.testing.assert_equal(sud['pre'].shape,(40,6,454))
    # Templates are estimated separately in each cohort and never tuned to cross-cohort test p values.
    aa_template=(aa['post']-aa['pre']).mean(0)
    sud_template=np.empty((6,454))
    for m in range(6):
        for roi in range(454):
            x=np.column_stack([np.ones(40),sud['active'],sud['pre'][:,m,roi],sud['fd']])
            sud_template[m,roi]=np.linalg.lstsq(x,sud['post'][:,m,roi],rcond=None)[0][1]
    results=[];templates=[];scores=[];loo=[];pattern_correlations=[]
    for branch,features in [('all_regions',np.arange(454)),('cortical_nodes',np.arange(400))]:
        at=aa_template[:,features].copy();st=sud_template[:,features].copy()
        at-=at.mean(1,keepdims=True);st-=st.mean(1,keepdims=True)
        at/=np.abs(at).sum(1,keepdims=True);st/=np.abs(st).sum(1,keepdims=True)
        assert np.isfinite(at).all() and np.isfinite(st).all()
        for donor,t in [('aachen',at),('sudmex',st)]:
            for k,m in enumerate(MODES):
                templates.extend({'branch':branch,'donor':donor,'mode':m,'roi':int(roi)+1,'weight':float(w)} for roi,w in zip(features,t[k]))
        # Both receiver scores use an externally fitted donor template.
        apre=np.einsum('nmr,mr->nm',aa['pre'][:,:,features],st)
        apost=np.einsum('nmr,mr->nm',aa['post'][:,:,features],st)
        spre=np.einsum('nmr,mr->nm',sud['pre'][:,:,features],at)
        spost=np.einsum('nmr,mr->nm',sud['post'][:,:,features],at)
        ad=apost-apre
        ap,amax=exact_max_t(ad)
        for k,m in enumerate(MODES):
            pattern_correlations.append({'branch':branch,'mode':m,'pearson_r':stats.pearsonr(at[k],st[k]).statistic,
                                         'spearman_rho':stats.spearmanr(at[k],st[k]).statistic,'inference':'Descriptive only; no parcel-wise p value.'})
            results.append(dict(branch=branch,dataset='aachen_itbs',mode=m,model='paired_external_template',
                                p=ap[k],p_maxT6=amax[k],**paired_fit(ad[:,k])))
            for model,adj in [('baseline_motion_adjusted',True),('baseline_adjusted',False)]:
                results.append(dict(branch=branch,dataset='sudmex_tms',mode=m,model=model,
                                    **hc3(spost[:,k],spre[:,k],sud['active'],sud['fd'],adj)))
            for d,info,b,f in [('aachen_itbs',aa,apre,apost),('sudmex_tms',sud,spre,spost)]:
                for n,subject in enumerate(info['ids']):
                    scores.append(dict(branch=branch,dataset=d,subject_id=subject,mode=m,pre_score=b[n,k],post_score=f[n,k],delta=f[n,k]-b[n,k]))
        if branch=='all_regions':
            for omit in range(15):
                x=np.delete(ad,omit,axis=0);ps,pm=exact_max_t(x)
                for k,m in enumerate(MODES):loo.append(dict(dataset='aachen_itbs',omission=f'omit_{omit+1:02d}',mode=m,estimate=x[:,k].mean(),p=ps[k],p_maxT6=pm[k]))
            for omit in range(40):
                keep=np.arange(40)!=omit
                for k,m in enumerate(MODES):
                    loo.append(dict(dataset='sudmex_tms',omission=f'omit_{omit+1:02d}',mode=m,
                                    **hc3(spost[keep,k],spre[keep,k],sud['active'][keep],sud['fd'][keep])))
    result=pd.DataFrame(results)
    result['p_holm6']=result.groupby(['branch','dataset','model']).p.transform(lambda p:multipletests(p,method='holm')[1])
    result['p_holm12']=result.groupby(['dataset','model']).p.transform(lambda p:multipletests(p,method='holm')[1])
    result.to_csv(out/'external_template_tests.csv',index=False)
    conjunction=[]
    for branch in ['all_regions','cortical_nodes']:
        for m in MODES:
            a=result.loc[result.branch.eq(branch)&result['mode'].eq(m)&result.dataset.eq('aachen_itbs')].iloc[0]
            b=result.loc[result.branch.eq(branch)&result['mode'].eq(m)&result.dataset.eq('sudmex_tms')&result.model.eq('baseline_motion_adjusted')].iloc[0]
            positive=a.estimate>0 and b.estimate>0
            conjunction.append(dict(branch=branch,mode=m,both_positive=bool(positive),aachen_estimate=a.estimate,sudmex_estimate=b.estimate,
                                    p_aachen=a.p,p_sudmex=b.p,p=max(a.p,b.p) if positive else 1.))
    cc=pd.DataFrame(conjunction);cc['p_holm6']=cc.groupby('branch').p.transform(lambda p:multipletests(p,method='holm')[1]);cc['p_holm12']=multipletests(cc.p,method='holm')[1]
    cc.to_csv(out/'cross_dataset_conjunction.csv',index=False)
    pd.DataFrame(pattern_correlations).to_csv(out/'pattern_correlations.csv',index=False)
    pd.DataFrame(templates).to_csv(out/'external_template_weights.csv',index=False)
    pd.DataFrame(scores).to_csv(out/'participant_projected_scores.csv',index=False)
    ll=pd.DataFrame(loo);ll['p_holm6']=ll.groupby(['dataset','omission']).p.transform(lambda p:multipletests(p,method='holm')[1]);ll.to_csv(out/'leave_one_out.csv',index=False)
    pd.DataFrame(checks).to_csv(out/'extraction_checks.csv',index=False)
    save_json(out/'source_sha256.json',HASHES)
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in HASHES.items())
    save_json(out/'verification.json',{'passed':True,'scans':140,'mode_node_means_match_validated_previous_extraction':True,
                                      'hc3_matches_independent_sandwich':True,'sources_unchanged':True,'leave_one_out_rows':len(ll),
                                      'primary_conjunction_holm6_passes':int(((cc.branch=='all_regions')&(cc.p_holm6<.05)).sum())})
    (out.parent/'LATEST.txt').write_text(out.name+'\n',encoding='utf-8')
    print(cc.to_string(index=False),flush=True)
    print(pd.DataFrame(pattern_correlations).to_string(index=False),flush=True)
    print('COMPLETE',out,flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
