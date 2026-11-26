"""Bounded follow-up of outside-workspace Phi; all variants retained.

TMS only. Does not alter preprocessing, the frozen HCP mask, or saved PhiID.
"""
from pathlib import Path
from datetime import datetime, timezone
import json
import sys

import numpy as np
import pandas as pd
from scipy import stats
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'setup'))
import analyse_tms_workspace_scope as base
SOURCE=ROOT/'output/tms_workspace_scope'/(ROOT/'output/tms_workspace_scope/LATEST.txt').read_text().strip()


def summarize(data, family):
    rows=[]
    p,pc=base.exact_max_t(data.to_numpy())
    for i,m in enumerate(data):
        rows.append({'metric':m,'family':family,**base.tms.paired_summary(data[m]),'p_sign_flip':p[i],'p_maxT_family':pc[i]})
    result=pd.DataFrame(rows);result['p_holm_family']=base.holm(result.p)
    return result


def pair(scans, metrics):
    rows=[]
    for (branch,subject),b in scans.groupby(['branch','subject_id']):
        b=b.set_index('session')
        for intensity in [90,120]:
            a,c=b.loc['pre'+str(intensity)],b.loc['post'+str(intensity)]
            r={'branch':branch,'subject_id':subject,'intensity':intensity,'delta_fd':c.mean_fd-a.mean_fd}
            for m in metrics:r.update({m:c[m]-a[m],m+'_baseline':a[m]})
            rows.append(r)
    return pd.DataFrame(rows)


def main():
    out=ROOT/'output/aachen_spatial_followup'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out.mkdir(parents=True);print('OUTPUT:',out,flush=True)
    hashes={}
    def track(p):hashes[str(p.relative_to(ROOT))]=base.tms.sha256(p);return p
    base.save_json(out/'analysis_scope.json',{
      'exploratory':True,'question':'What supports or limits the apparent outside-workspace Phi increase?',
      'checks':['visit and intensity consistency','participant influence','motion and baseline differences',
                'fixed-network block contributions and leave-one-network-out','TDMI minus self-information identity',
                'equal 294-volume estimation window','Phi estimated on denoised BOLD without HRF deconvolution'],
      'limits':'No sham control; same participants throughout; alternative preprocessing sensitivities do not create independent replications',
      'multiplicity':'Joint maxT within each complete follow-up question; all variants retained; previous five/20/40 tests remain reference results'})
    ref=base.tms.load_reference(ROOT)
    for p in ref['input_files']:track(p)
    edges,masks=base.edge_scopes(ref['workspace']);i,j=edges.T
    modules=ref['modules'];w=ref['workspace']
    networks={};counts=[]
    for a in range(8):
        for b in range(a,8):
            sel=masks['outside_workspace']&(((modules[i]==a)&(modules[j]==b))|((modules[i]==b)&(modules[j]==a)))
            name=base.tms.NETWORKS[a]+'_'+base.tms.NETWORKS[b]
            counts.append({'network_pair':name,'edges':int(sel.sum()),'weight_of_outside_mean':sel.sum()/23871})
            if sel.any():networks['block__'+name]=sel
    for n in range(8):
        sel=masks['outside_workspace']&(modules[i]!=n)&(modules[j]!=n)
        networks['without__'+base.tms.NETWORKS[n]]=sel
    pd.DataFrame(counts).to_csv(out/'network_pair_counts.csv',index=False)
    pd.DataFrame({'network':base.tms.NETWORKS,'workspace_nodes':[int(((modules==n)&w).sum()) for n in range(8)],
                  'outside_nodes':[int(((modules==n)&~w).sum()) for n in range(8)]}).to_csv(out/'network_node_counts.csv',index=False)
    scans=pd.read_csv(track(SOURCE/'raw_scan_metrics.csv')).query("dataset=='aachen_itbs'")
    norm_refs={}
    for seed in base.SEEDS:
        folder=base.NUMIT/'references/aachen_itbs'/str(seed)
        norm_refs[seed]=(np.load(track(folder/'tdmi_grid.npy')),np.load(track(folder/'phi.npy'),mmap_mode='r'))
    values=[];profiles=[];max_identity=0.;max_matlab=0.
    for n,scan in enumerate(scans.itertuples(),1):
        prefix=scan.prefix
        deconv=np.load(track(ROOT/'pipeline/results/timeseries'/(prefix+'_timeseries_deconv.npy')))
        denoised=np.load(track(ROOT/'pipeline/results/timeseries'/(prefix+'_timeseries.npy')))
        rawphi=np.load(track(ROOT/'pipeline/results/phiid_matrices'/(prefix+'_phi.npy')))[i,j]
        rawphir=np.load(track(ROOT/'pipeline/results/phiid_matrices'/(prefix+'_phi_r.npy')))[i,j]
        for branch,ts in [('deconvolved',deconv),('equal294',deconv[:294]),('no_deconvolution',denoised)]:
            sigma=base.tn.pg.lagged_sigma_from_timeseries(ts,1,edges)
            tdmi=base.tn.pg.tdmi_from_sigma4(sigma)
            past=ts[:-1];future=ts[1:]
            past=(past-past.mean(0))/past.std(0);future=(future-future.mean(0))/future.std(0)
            ac=(past*future).mean(0);ais=-.5*np.log1p(-ac**2)
            own=ais[i]+ais[j];phi=tdmi-own
            redundancy=np.min(np.stack([-.5*np.log1p(-sigma[:,a,b]**2) for a,b in [(0,2),(0,3),(1,2),(1,3)]]),axis=0)
            if branch=='deconvolved':
                err=float(abs(phi-rawphi).max());max_matlab=max(max_matlab,err)
                np.testing.assert_allclose(phi,rawphi,rtol=1e-7,atol=2e-9)
                np.testing.assert_allclose(phi+redundancy,rawphir,rtol=1e-7,atol=2e-9)
            r={'branch':branch,'subject_id':scan.subject_id,'session':scan.session,'mean_fd':scan.mean_fd,'n_timepoints':len(ts)}
            for scope,sel in masks.items():
                r[scope+'__phi']=phi[sel].mean();r[scope+'__phi_r']=(phi+redundancy)[sel].mean()
                r[scope+'__tdmi']=tdmi[sel].mean();r[scope+'__self_sum']=own[sel].mean();r[scope+'__redundancy']=redundancy[sel].mean()
                max_identity=max(max_identity,abs(r[scope+'__phi']-(r[scope+'__tdmi']-r[scope+'__self_sum'])))
            for metric in ['phi','phi_r','tdmi','self_sum','redundancy']:
                r['difference__'+metric]=r['outside_workspace__'+metric]-r['within_workspace__'+metric]
            values.append(r)
            if branch=='deconvolved':
                q=[]
                for grid,reference in norm_refs.values():q.append(base.tn.conditional_cdf(grid,reference,tdmi,rawphi))
                q=(q[0]+q[1])/2
                for pb,edgevals in [('raw',rawphi),('normalized',q)]:
                    profiles.append({'branch':pb,'subject_id':scan.subject_id,'session':scan.session,'mean_fd':scan.mean_fd,
                        **{k:edgevals[v].mean() for k,v in networks.items()},
                        **{s+'__phi':edgevals[v].mean() for s,v in masks.items()},
                        'difference__phi':edgevals[masks['outside_workspace']].mean()-edgevals[masks['within_workspace']].mean()})
        if n%10==0:print(f'Signal, graph and estimator checks: {n}/60 scans',flush=True)
    signals=pd.DataFrame(values);profile=pd.DataFrame(profiles)
    signals.to_csv(out/'signal_decomposition_scans.csv',index=False);profile.to_csv(out/'network_profile_scans.csv',index=False)
    signal_metrics=[c for c in signals if '__' in c];profile_metrics=[c for c in profile if '__' in c]
    sp=pair(signals,signal_metrics);pp=pair(profile,profile_metrics)
    sp.to_csv(out/'signal_participant_changes.csv',index=False);pp.to_csv(out/'network_participant_changes.csv',index=False)
    signal_tests=[];network_tests=[];visit_tests=[];adjusted=[];influence=[];nonparam=[]
    order={}
    import re
    for line in track(ROOT/'data/aachen_itbs/derivatives/README.txt').read_text().splitlines():
        m=re.match(r'^\s*(\d+)\s+([1342/]+)\s+',line)
        if m:order['aachen_itbs_'+m[1].zfill(2)]=120 if m[2].startswith('3') else 90
    for branch,block in sp.groupby('branch'):
        pooled=block.groupby('subject_id').mean(numeric_only=True)
        cols=[s+'__'+m for s in ['outside_workspace','within_workspace','difference'] for m in ['tdmi','self_sum','redundancy']]
        r=summarize(pooled[cols],'temporal_components9');r['branch']=branch;signal_tests.append(r)
    variants={}
    for branch,block in sp.groupby('branch'):
        pooled=block.groupby('subject_id').mean(numeric_only=True)
        for m in ['whole_brain__phi','within_workspace__phi','outside_workspace__phi','difference__phi']:
            variants[branch+'::'+m]=pooled[m]
    summarize(pd.DataFrame(variants),'estimator_variants12').to_csv(out/'estimator_sensitivity_results.csv',index=False)
    pd.concat(signal_tests,ignore_index=True).to_csv(out/'temporal_component_results.csv',index=False)
    for branch,block in pp.groupby('branch'):
        pooled=block.groupby('subject_id').mean(numeric_only=True)
        for family,cols in [('network_blocks',[m for m in profile_metrics if m.startswith('block__')]),
                            ('network_omissions8',[m for m in profile_metrics if m.startswith('without__')])]:
            r=summarize(pooled[cols],family);r['branch']=branch;network_tests.append(r)
        visit={}
        for m in ['outside_workspace__phi','difference__phi']:
            wide=block.pivot(index='subject_id',columns='intensity',values=m).sort_index()
            first=np.array([wide.loc[s,order[s]] for s in wide.index]);second=np.array([wide.loc[s,210-order[s]] for s in wide.index])
            for name,x in [('90',wide[90]),('120',wide[120]),('120_minus_90',wide[120]-wide[90]),
                           ('first_visit',first),('second_visit',second),('second_minus_first',second-first)]:visit[m+'::'+name]=np.asarray(x)
            x=pooled[m].to_numpy()
            for k,s in enumerate(pooled.index):influence.append({'branch':branch,'metric':m,'omitted_subject':s,**base.tms.paired_summary(np.delete(x,k))})
            rng=np.random.default_rng(2026112528);boot=x[rng.integers(15,size=(10000,15))].mean(1);lo,hi=np.quantile(boot,[.025,.975])
            nonparam.append({'branch':branch,'metric':m,'positive_changes':int((x>0).sum()),'n':15,
                             'wilcoxon_p':stats.wilcoxon(x,method='exact').pvalue,
                             'sign_test_p':stats.binomtest(int((x>0).sum()),15,.5).pvalue,
                             'bootstrap_low':lo,'bootstrap_high':hi,'visit_correlation':wide[90].corr(wide[120])})
            baseline=pooled[m+'_baseline']
            frame=pd.DataFrame({'y':x,'intercept':1.,'fd':pooled.delta_fd,'baseline_centered':baseline-baseline.mean()})
            for name,predictors in [('motion',['intercept','fd']),('baseline',['intercept','baseline_centered']),
                                    ('baseline_motion',['intercept','baseline_centered','fd'])]:
                adjusted.append({'branch':branch,'metric':m,'model':name,
                                 **base.hc3(frame,'y',predictors,'intercept')})
        r=summarize(pd.DataFrame(visit),'visit_and_intensity12');r['branch']=branch;visit_tests.append(r)
    pd.concat(network_tests,ignore_index=True).to_csv(out/'network_results.csv',index=False)
    pd.concat(visit_tests,ignore_index=True).to_csv(out/'visit_consistency_results.csv',index=False)
    adjusted=pd.DataFrame(adjusted);adjusted['p_holm6']=adjusted.groupby('branch').p.transform(base.holm)
    adjusted.to_csv(out/'baseline_motion_sensitivity.csv',index=False)
    pd.DataFrame(influence).to_csv(out/'participant_influence.csv',index=False)
    pd.DataFrame(nonparam).to_csv(out/'nonparametric_checks.csv',index=False)
    track(Path(__file__))
    for p,h in hashes.items():assert base.tms.sha256(ROOT/p)==h
    base.save_json(out/'provenance.json',{'status':'complete','input_sha256':hashes,'n':15,'scans':60,
                                        'max_matlab_phi_error':max_matlab,'max_temporal_identity_error':max_identity,
                                        'normalized_reference_unchanged':True,'no_deconvolution_sensitivity':'Gaussian Phi from denoised BOLD; same regional data and lag; raw nats only',
                                        'equal_length_sensitivity':'First 294 samples of previously deconvolved signals; HRF fitting is not repeated'})
    (ROOT/'output/aachen_spatial_followup/LATEST.txt').write_text(out.name+'\n')
    print('COMPLETE:',out,flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=1):main()
