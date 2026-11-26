"""Bounded target-based PhiID follow-up; reads TMS-only frozen time series.

No processing code is modified. Definitions and caveats were recorded before
extracting targeted outcomes in PROVENANCE/TARGET_FLOW_WORKLOG_20261125.md.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import sys
import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import savemat
import statsmodels.api as sm
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'pipeline/scripts/analysis'))
import tms_numit as tn
sys.path.insert(0,str(ROOT/'setup'))
from analyse_tms_familywise import exact_max_t,holm

METRICS=['net_te','net_unique','net_synergistic','persistent_synergy','persistent_redundancy']
OLD=['synergy','redundancy','synergy_total','phi','phi_r']
VARIANTS={'lag1_full':(1,False),'lag2_full':(2,False),'lag1_equal_length':(1,True)}
SCOPES=['left','left_minus_right']
SEED=2026112519
HASHES={}

def track(path):
    key=str(path.relative_to(ROOT))
    if key not in HASHES:HASHES[key]=hashlib.sha256(path.read_bytes()).hexdigest()
    return path

def directed(atoms):
    out={name:atoms[name] for name in tn.pg.ATOM_NAMES}
    out['te_out']=atoms['xtr']+atoms['xty']+atoms['str']+atoms['sty']
    out['te_in']=atoms['ytr']+atoms['ytx']+atoms['str']+atoms['stx']
    out['net_te']=out['te_out']-out['te_in']
    out['net_unique']=(atoms['xtr']+atoms['xty'])-(atoms['ytr']+atoms['ytx'])
    out['net_synergistic']=atoms['sty']-atoms['stx']
    out['persistent_synergy']=atoms['sts']
    out['persistent_redundancy']=atoms['rtr']
    np.testing.assert_allclose(out['net_te'],out['net_unique']+out['net_synergistic'],atol=1e-12)
    return out

def cmi(sigma,x,y,z):
    def ld(index):
        sign,val=np.linalg.slogdet(sigma[:,index,:][:,:,index])
        assert (sign>0).all()
        return val
    return .5*(ld([x,z])+ld([y,z])-ld([z])-ld([x,y,z]))

def numerical_unit_checks():
    rng=np.random.default_rng(SEED)
    ts=np.zeros((30000,2));noise=rng.normal(size=ts.shape)
    for t in range(1,len(ts)):
        ts[t,0]=.55*ts[t-1,0]+noise[t,0]
        ts[t,1]=.65*ts[t-1,0]+.3*ts[t-1,1]+noise[t,1]
    sig=tn.pg.lagged_sigma_from_timeseries(ts[1000:],1,np.array([[0,1]]))
    a=directed(tn.pg.phiid_atoms(sig))
    swapped=directed(tn.pg.phiid_atoms(sig[:,[1,0,3,2],:][:,:,[1,0,3,2]]))
    assert a['te_out'][0]>.1 and a['te_in'][0]<.001
    for name in METRICS[:3]:np.testing.assert_allclose(a[name],-swapped[name],atol=1e-12)
    for name in METRICS[3:]:np.testing.assert_allclose(a[name],swapped[name],atol=1e-12)
    np.testing.assert_allclose(a['te_out'],cmi(sig,0,3,1),atol=1e-12)
    np.testing.assert_allclose(a['te_in'],cmi(sig,1,2,0),atol=1e-12)
    return {'synthetic_x_to_y_te':float(a['te_out'][0]),'synthetic_y_to_x_te':float(a['te_in'][0]),
            'source_swap_and_component_identity':'passed'}

def get_inputs():
    source=ROOT/'output/tms_workspace_scope'/(ROOT/'output/tms_workspace_scope/LATEST.txt').read_text().strip()
    pairs=pd.read_csv(track(source/'raw_paired_changes.csv'))
    scans=pd.read_csv(track(source/'raw_scan_metrics.csv'))
    w1=pd.read_csv(track(ROOT/'output/aachen_target_anatomy/subject_region_weights.csv'))
    w1['definition']='deposited_individual_provisional'
    w2=pd.read_csv(track(ROOT/'output/sudmex_target_anatomy/subject_region_weights.csv'))
    weights=pd.concat([w1,w2],ignore_index=True)
    assert not weights.duplicated(['definition','subject_id','side','roi']).any()
    assert pairs.loc[pairs.dataset.eq('sudmex_tms')].comparison.eq('t1_minus_t0').all()
    return pairs,scans,weights

def extract(out):
    pairs,scans,weights=get_inputs()
    rows=[];checks=[];matlab_inputs=[];matlab_expected=[];matlab_meta=[]
    for sindex,(subject,wb) in enumerate(weights.groupby('subject_id')):
        p=pairs.loc[pairs.subject_id.eq(subject)]
        sessions=set(p.baseline_session)|set(p.followup_session)
        sb=scans.loc[scans.subject_id.eq(subject)&scans.session.isin(sessions)]
        series={r.session:np.load(track(ROOT/'pipeline/results/timeseries'/(r.prefix+'_timeseries_deconv.npy'))) for r in sb.itertuples()}
        equal=min(len(x) for x in series.values())
        assert all(x.shape[1]==454 and np.isfinite(x).all() for x in series.values())
        for r in sb.itertuples():
            for variant,(lag,trim) in VARIANTS.items():
                ts=series[r.session][:equal] if trim else series[r.session]
                for definition,block in wb.groupby('definition'):
                    sideweights={side:b.sort_values('roi').weight.to_numpy() for side,b in block.groupby('side')}
                    seeds=np.unique(np.concatenate([np.flatnonzero(w>0) for w in sideweights.values()]))
                    remote=np.setdiff1d(np.arange(454),seeds)
                    pr=np.array([(i,j) for i in seeds for j in remote],int)
                    sigma=tn.pg.lagged_sigma_from_timeseries(ts,lag,pr)
                    atoms=tn.pg.phiid_atoms(sigma);values=directed(atoms)
                    assert all(np.isfinite(v).all() for v in values.values())
                    error=max(float(np.max(abs(values['te_out']-cmi(sigma,0,3,1)))),
                              float(np.max(abs(values['te_in']-cmi(sigma,1,2,0)))))
                    assert error<2e-10
                    check={'prefix':r.prefix,'definition':definition,'variant':variant,'n_pairs':len(pr),'max_cmi_error':error}
                    if variant=='lag1_full':
                        for name in OLD:
                            old=np.load(track(ROOT/'pipeline/results/phiid_matrices'/(r.prefix+'_'+name+'.npy')))
                            err=float(np.max(abs(atoms[name]-old[pr[:,0],pr[:,1]])))
                            check['max_saved_'+name+'_error']=err
                            assert err<2e-8,(r.prefix,name,err)
                    checks.append(check)
                    for side,w in sideweights.items():
                        assert abs(w.sum()-1)<1e-10
                        # Common destinations for both sides. Seeds have anatomical weights.
                        pw=w[pr[:,0]]/len(remote)
                        assert abs(pw.sum()-1)<1e-10
                        record={'dataset':r.dataset,'subject_id':subject,'session':r.session,'arm':r.arm,
                                'definition':definition,'variant':variant,'side':side,'mean_fd':r.mean_fd,
                                'n_timepoints':len(ts),'n_seed_parcels':int((w>0).sum()),'n_remote_parcels':len(remote)}
                        record.update({name:float(pw@v) for name,v in values.items()})
                        rows.append(record)
                    # Bounded check in the original MATLAB code, fixed selection independent of effects.
                    if subject in ['aachen_itbs_01','sudmex_tms_001'] and variant!='lag1_equal_length':
                        for idx in np.linspace(0,len(pr)-1,3,dtype=int):
                            i,j=pr[idx]
                            matlab_inputs.append(np.array([ts[:-lag,i],ts[:-lag,j],ts[lag:,i],ts[lag:,j]]))
                            matlab_expected.append([atoms[name][idx] for name in tn.pg.ATOM_NAMES])
                            matlab_meta.append({'prefix':r.prefix,'lag':lag,'x_roi':int(i+1),'y_roi':int(j+1)})
        if (sindex+1)%5==0:print('Extracted subjects:',sindex+1,flush=True)
    pd.DataFrame(rows).to_csv(out/'scan_metrics.csv',index=False)
    pd.DataFrame(checks).to_csv(out/'numerical_scan_checks.csv',index=False)
    cells=np.empty((len(matlab_inputs),1),dtype=object)
    for i,v in enumerate(matlab_inputs):cells[i,0]=v
    savemat(out/'matlab_validation_input.mat',{'vectors':cells,'expected':np.asarray(matlab_expected),'atom_names':np.array(tn.pg.ATOM_NAMES,dtype=object)})
    pd.DataFrame(matlab_meta).to_csv(out/'matlab_validation_pairs.csv',index=False)
    return make_pairs(pd.DataFrame(rows),pairs)

def make_pairs(scans,sourcepairs):
    names=METRICS+['te_out','te_in']+tn.pg.ATOM_NAMES
    rows=[]
    for (definition,variant,subject),block in scans.groupby(['definition','variant','subject_id']):
        for p in sourcepairs.loc[sourcepairs.subject_id.eq(subject)].itertuples():
            bs=block.loc[block.session.eq(p.baseline_session)].set_index('side')
            ps=block.loc[block.session.eq(p.followup_session)].set_index('side')
            for scope in SCOPES:
                if 'l' not in bs.index or (scope=='left_minus_right' and 'r' not in bs.index):continue
                row={'dataset':p.dataset,'subject_id':subject,'arm':p.arm,'comparison':p.comparison,
                     'definition':definition,'variant':variant,'scope':scope,'delta_fd':p.delta_fd}
                for name in names:
                    pre=bs.loc['l',name];post=ps.loc['l',name]
                    if scope=='left_minus_right':pre-=bs.loc['r',name];post-=ps.loc['r',name]
                    row[name+'_baseline']=pre;row[name+'_followup']=post;row[name+'_delta']=post-pre
                rows.append(row)
    return pd.DataFrame(rows)

def ancova(y,baseline,arms,fd=None):
    """HC3 ANCOVA for many common arm assignments, residualized covariates."""
    y=np.asarray(y,float);baseline=np.asarray(baseline,float)
    z=np.column_stack([np.ones(len(y)),baseline-baseline.mean()])
    if fd is not None:z=np.column_stack([z,np.asarray(fd)-np.mean(fd)])
    q=np.linalg.qr(z,mode='reduced')[0]
    ar=arms-(arms@q)@q.T;yr=y-q@(q.T@y)
    ss=(ar**2).sum(1)
    assert ss.min()>1e-8
    beta=ar@yr/ss
    residual=yr[None,:]-beta[:,None]*ar
    leverage=(q**2).sum(1)[None,:]+ar**2/ss[:,None]
    assert leverage.max()<1-1e-9
    se=np.sqrt(((ar*residual/(1-leverage))**2).sum(1))/ss
    return beta,se,beta/se,len(y)-z.shape[1]-1

def trial_stats(block,b=32768,motion=False):
    subjects=sorted(set.intersection(*[set(block.loc[block.scope.eq(scope),'subject_id']) for scope in SCOPES]))
    blocks={scope:block.loc[block.scope.eq(scope)].set_index('subject_id').loc[subjects] for scope in SCOPES}
    labels=blocks['left'].arm.eq('active').to_numpy(float)
    rng=np.random.default_rng(SEED)
    assignments=np.vstack([labels]+[rng.permutation(labels) for _ in range(b)])
    results=[];distributions=[]
    for scope in SCOPES:
        d=blocks[scope]
        for metric in METRICS:
            y=d[metric+'_followup'].to_numpy();baseline=d[metric+'_baseline'].to_numpy()
            fd=d.delta_fd.to_numpy() if motion else None
            beta,se,t,df=ancova(y,baseline,assignments,fd)
            # Independently verify the observed robust coefficient and two permuted fits.
            for idx in [0,1,min(7,b)]:
                x=np.column_stack([np.ones(len(y)),baseline-baseline.mean(),assignments[idx]])
                if motion:x=np.column_stack([x,fd-fd.mean()])
                fit=sm.OLS(y,x).fit(cov_type='HC3',use_t=True)
                np.testing.assert_allclose([beta[idx],se[idx],t[idx]],[fit.params[2],fit.bse[2],fit.tvalues[2]],rtol=1e-8,atol=1e-10)
            low,high=beta[0]+np.array([-1,1])*stats.t.ppf(.975,df)*se[0]
            rec={'scope':scope,'metric':metric,'n':len(y),'n_active':int(labels.sum()),'n_sham':int((1-labels).sum()),
                 'estimate':beta[0],'ci_low':low,'ci_high':high,'t':t[0],'p_hc3_t':2*stats.t.sf(abs(t[0]),df),
                 'active_pre':float(baseline[labels==1].mean()),'active_post':float(y[labels==1].mean()),
                 'sham_pre':float(baseline[labels==0].mean()),'sham_post':float(y[labels==0].mean()),
                 'active_delta':float((y-baseline)[labels==1].mean()),'sham_delta':float((y-baseline)[labels==0].mean())}
            results.append(rec);distributions.append(abs(t))
    tdist=np.array(distributions).T;observed=tdist[0];dist=tdist[1:]
    p=(1+(dist>=observed[None]-1e-10).sum(0))/(b+1)
    adj=(1+(dist.max(1)[:,None]>=observed[None]-1e-10).sum(0))/(b+1)
    adj5=np.empty(10)
    for start in [0,5]:adj5[start:start+5]=(1+(dist[:,start:start+5].max(1)[:,None]>=observed[None,start:start+5]-1e-10).sum(0))/(b+1)
    answer=pd.DataFrame(results);answer['p_permutation']=p;answer['p_maxT5']=adj5;answer['p_maxT10']=adj
    answer['p_holm10']=holm(p)
    return answer

def aachen_stats(block):
    results=[]
    for scope,b in block.groupby('scope',sort=False):
        d=b.groupby('subject_id').mean(numeric_only=True)
        values=d[[m+'_delta' for m in METRICS]].to_numpy()
        p,adjusted=exact_max_t(values)
        for j,metric in enumerate(METRICS):
            v=values[:,j];n=len(v);mean=v.mean();se=stats.sem(v);t=mean/se
            ci=mean+np.array([-1,1])*stats.t.ppf(.975,n-1)*se
            results.append({'scope':scope,'metric':metric,'n':n,'estimate':mean,'ci_low':ci[0],'ci_high':ci[1],
              't':t,'p_t':2*stats.t.sf(abs(t),n-1),'p_sign_flip':p[j],'p_maxT5':adjusted[j],
              'pre':d[metric+'_baseline'].mean(),'post':d[metric+'_followup'].mean(),
              'positive_participants':int((v>0).sum())})
    df=pd.DataFrame(results);df['p_holm10']=holm(df.p_sign_flip)
    return df

def sensitivities(pairs,out):
    loo=[];aac_adjust=[];intensity=[]
    primary=pairs.loc[pairs.variant.eq('lag1_full')]
    for (dataset,definition,scope),b in primary.groupby(['dataset','definition','scope']):
        d=b.groupby('subject_id').mean(numeric_only=True) if dataset=='aachen_itbs' else b.set_index('subject_id')
        for metric in METRICS:
            for removed in d.index:
                left=d.drop(index=removed)
                if dataset=='sudmex_tms':
                    arms=left.arm.eq('active').to_numpy(float)[None,:]
                    beta,se,t,df=ancova(left[metric+'_followup'],left[metric+'_baseline'],arms)
                    est=beta[0];p=2*stats.t.sf(abs(t[0]),df)
                else:
                    est=left[metric+'_delta'].mean();p=stats.ttest_1samp(left[metric+'_delta'],0).pvalue
                loo.append({'dataset':dataset,'definition':definition,'scope':scope,'metric':metric,'omitted':removed,'estimate':est,'p_nominal':p})
            if dataset=='aachen_itbs':
                x=np.column_stack([np.ones(len(d)),d[metric+'_baseline']-d[metric+'_baseline'].mean(),d.delta_fd])
                fit=sm.OLS(d[metric+'_delta'],x).fit(cov_type='HC3',use_t=True)
                ci=fit.conf_int().to_numpy()[0]
                aac_adjust.append({'scope':scope,'metric':metric,'n':len(d),'estimate':float(fit.params.iloc[0]),'ci_low':ci[0],'ci_high':ci[1],'p':float(fit.pvalues.iloc[0])})
        if dataset=='aachen_itbs':
            for comparison,v in b.groupby('comparison'):
                for metric in METRICS:
                    vals=v[metric+'_delta'];ci=stats.t.interval(.95,len(vals)-1,loc=vals.mean(),scale=stats.sem(vals))
                    intensity.append({'scope':scope,'comparison':comparison,'metric':metric,'n':len(vals),'estimate':vals.mean(),'ci_low':ci[0],'ci_high':ci[1],'p_nominal':stats.ttest_1samp(vals,0).pvalue})
    pd.DataFrame(loo).to_csv(out/'leave_one_out.csv',index=False)
    a=pd.DataFrame(aac_adjust);a['p_holm10']=holm(a.p);a.to_csv(out/'aachen_baseline_motion_sensitivity.csv',index=False)
    i=pd.DataFrame(intensity);i['p_holm20']=holm(i.p_nominal);i.to_csv(out/'aachen_intensity_descriptive.csv',index=False)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--permutations',type=int,default=32768);parser.add_argument('--reuse',type=Path)
    args=parser.parse_args()
    out=ROOT/'output/targeted_information_flow'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out.mkdir(parents=True)
    with threadpool_limits(limits=2):
        validation=numerical_unit_checks()
        if args.reuse:
            pairs=pd.read_csv(track(args.reuse/'participant_visit_changes.csv'))
        else:pairs=extract(out)
        pairs.to_csv(out/'participant_visit_changes.csv',index=False)
        results=[]
        for (dataset,definition,variant),block in pairs.groupby(['dataset','definition','variant']):
            if dataset=='sudmex_tms':
                models=['baseline','baseline_motion'] if variant=='lag1_full' else ['baseline']
                for model in models:
                    result=trial_stats(block,args.permutations,motion=model=='baseline_motion')
                    result['model']=model;result['dataset']=dataset;result['definition']=definition;result['variant']=variant;results.append(result)
            else:
                result=aachen_stats(block);result['model']='paired_mean';result['dataset']=dataset;result['definition']=definition;result['variant']=variant;results.append(result)
            print('Tested',dataset,definition,variant,flush=True)
        all_results=pd.concat(results,ignore_index=True);all_results.to_csv(out/'inference.csv',index=False)
        sensitivities(pairs,out)
    validation.update({'created_utc':datetime.now(timezone.utc).isoformat(),'seed':SEED,'permutations':args.permutations,
      'units':'nats; marginal 95% CIs are not multiplicity-adjusted',
      'sudmex_family':'5 measures x 2 scopes, shared-label max-|HC3 t|; Holm10 sensitivity',
      'aachen_family':'5 measures each scope; Holm10 wider sensitivity; provisional localization',
      'permutation_assumption':'Exchangeable treatment labels within imaging completers, fixed group counts; original randomization blocks not reconstructed.',
      'aachen_sign_flip_assumption':'Symmetric participant change distributions; not a sham randomization test.',
      'input_sha256':HASHES,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    (out/'provenance.json').write_text(json.dumps(validation,indent=2),encoding='utf-8')
    (out.parent/'LATEST.txt').write_text(out.name,encoding='utf-8')
    print('COMPLETE',out,flush=True)
    print(all_results.loc[all_results.variant.eq('lag1_full')&all_results.model.ne('baseline_motion'),['dataset','definition','scope','metric','n','estimate','p_maxT5','p_maxT10','p_holm10']].to_string(index=False))

if __name__=='__main__':main()
