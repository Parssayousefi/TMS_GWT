"""Independent full-design OLS/permutation and paired-inference verification."""
from pathlib import Path
import hashlib,itertools,json
import numpy as np
import pandas as pd
from scipy import stats

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output/targeted_information_flow'/(ROOT/'output/targeted_information_flow/LATEST.txt').read_text().strip()
METRICS=['net_te','net_unique','net_synergistic','persistent_synergy','persistent_redundancy']

def main():
    provenance=json.loads((OUT/'provenance.json').read_text())
    for name,digest in provenance['input_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    scans=pd.read_csv(OUT/'scan_metrics.csv');pairs=pd.read_csv(OUT/'participant_visit_changes.csv');results=pd.read_csv(OUT/'inference.csv')
    np.testing.assert_allclose(scans.net_te,scans.net_unique+scans.net_synergistic,atol=2e-15)
    for p in pairs.itertuples():
        d=scans.loc[scans.subject_id.eq(p.subject_id)&scans.definition.eq(p.definition)&scans.variant.eq(p.variant)]
        pre,post=('t0','t1') if p.dataset=='sudmex_tms' else (('pre90','post90') if '90' in p.comparison else ('pre120','post120'))
        for metric in METRICS:
            a=d.loc[d.session.eq(pre)].set_index('side')[metric];b=d.loc[d.session.eq(post)].set_index('side')[metric]
            change=b['l']-a['l']
            if p.scope=='left_minus_right':change-=b['r']-a['r']
            np.testing.assert_allclose(change,getattr(p,metric+'_delta'),atol=2e-15)
    max_t_error=0.;max_p_error=0.
    # Independent batched full-model HC3 (does not use the residualization helper).
    for definition in ['common_baseline_target','individual_baseline_target']:
        block=pairs.loc[pairs.definition.eq(definition)&pairs.variant.eq('lag1_full')]
        subjects=sorted(block.subject_id.unique())
        left=block.loc[block.scope.eq('left')].set_index('subject_id').loc[subjects]
        arm=left.arm.eq('active').to_numpy(float)
        rng=np.random.default_rng(provenance['seed']);B=provenance['permutations']
        perm=np.vstack([arm]+[rng.permutation(arm) for _ in range(B)])
        abs_t=[];observed=[]
        for scope in ['left','left_minus_right']:
            d=block.loc[block.scope.eq(scope)].set_index('subject_id').loc[subjects]
            for metric in METRICS:
                y=d[metric+'_followup'].to_numpy();base=d[metric+'_baseline'].to_numpy()
                x=np.ones((B+1,len(y),3));x[:,:,1]=base-base.mean();x[:,:,2]=perm
                inv=np.linalg.inv(np.einsum('bni,bnj->bij',x,x))
                pinv=inv@x.transpose(0,2,1)
                coef=pinv@y
                residual=y[None]-np.einsum('bni,bi->bn',x,coef)
                h=np.einsum('bni,bin->bn',x,pinv)
                se=np.sqrt(np.sum((pinv[:,2,:]*residual/(1-h))**2,axis=1))
                t=coef[:,2]/se;abs_t.append(abs(t))
                r=results.loc[results.definition.eq(definition)&results.variant.eq('lag1_full')&results.model.eq('baseline')&results.scope.eq(scope)&results.metric.eq(metric)].iloc[0]
                np.testing.assert_allclose([coef[0,2],t[0]],[r.estimate,r.t],rtol=1e-8,atol=1e-10)
                max_t_error=max(max_t_error,abs(t[0]-r.t));observed.append(r)
        dist=np.array(abs_t).T
        corrected=(1+(dist[1:].max(1)[:,None]>=dist[0][None]-1e-10).sum(0))/(B+1)
        raw=(1+(dist[1:]>=dist[0][None]-1e-10).sum(0))/(B+1)
        for j,r in enumerate(observed):
            np.testing.assert_allclose([corrected[j],raw[j]],[r.p_maxT10,r.p_permutation],atol=1e-14)
            max_p_error=max(max_p_error,abs(corrected[j]-r.p_maxT10))
    for scope in ['left','left_minus_right']:
        d=pairs.loc[pairs.dataset.eq('aachen_itbs')&pairs.variant.eq('lag1_full')&pairs.scope.eq(scope)].groupby('subject_id').mean(numeric_only=True)
        v=d[[m+'_delta' for m in METRICS]].to_numpy()
        signs=np.array(list(itertools.product([-1,1],repeat=len(v))))
        t=np.abs(stats.ttest_1samp(signs[:,:,None]*v[None,:,:],0,axis=1).statistic)
        obs=np.abs(stats.ttest_1samp(v,0,axis=0).statistic)
        p=(t.max(1)[:,None]>=obs[None]-1e-10).mean(0)
        r=results.loc[results.dataset.eq('aachen_itbs')&results.variant.eq('lag1_full')&results.scope.eq(scope)].set_index('metric').loc[METRICS]
        np.testing.assert_allclose(p,r.p_maxT5,atol=1e-14)
    numerical=pd.read_csv(OUT/'numerical_scan_checks.csv')
    report={'status':'passed','input_hashes_checked':len(provenance['input_sha256']),'scan_rows':len(scans),
      'participant_visit_rows':len(pairs),'max_OLS_t_error':max_t_error,'max_maxT10_p_error':max_p_error,
      'max_saved_matrix_error':float(numerical.filter(regex='max_saved_').max().max()),
      'max_independent_CMI_error':float(numerical.max_cmi_error.max()),
      'paired_changes':'independently reconstructed from scan rows',
      'aachen_inference':'independent SciPy enumeration of every sign pattern',
      'sudmex_inference':'full-design batched OLS and HC3 sandwich, both localization definitions; 32768 permutations each',
      'anatomical_QA':'All 16 Aachen sagittal panels and SUDMEX common target reviewed. Aachen localization remains provisional due extracortical mask portions.'}
    (OUT/'independent_verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
