"""Independent checks plus bounded method/HRF diagnostics for the spatial follow-up."""
from pathlib import Path
import sys
import json

import h5py
import numpy as np
import pandas as pd
from scipy import stats
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'setup'))
import deepen_aachen_spatial_effect as follow
from analyse_tms_workspace_scope import exact_max_t,holm,hc3,tms,tn


def main():
    out=ROOT/'output/aachen_spatial_followup'/(ROOT/'output/aachen_spatial_followup/LATEST.txt').read_text().strip()
    profile=pd.read_csv(out/'network_profile_scans.csv')
    weights=pd.read_csv(out/'network_pair_counts.csv').query('edges>0')
    computed=sum(profile['block__'+r.network_pair]*r.weight_of_outside_mean for r in weights.itertuples())
    np.testing.assert_allclose(computed,profile.outside_workspace__phi,rtol=1e-11,atol=1e-12)
    pairs=pd.read_csv(out/'signal_participant_changes.csv')
    d=pairs.groupby(['branch','subject_id']).mean(numeric_only=True)
    method={}
    for scope in ['whole_brain','within_workspace','outside_workspace','difference']:
        for measure in ['phi','phi_r']:
            metric=scope+'__'+measure
            method[metric]=d.xs('deconvolved')[metric]-d.xs('no_deconvolution')[metric]
    contrasts=follow.summarize(pd.DataFrame(method),'deconvolution_interactions8')
    contrasts.to_csv(out/'method_interactions.csv',index=False)
    # Independent slogdet calculation on denoised BOLD, without the cofactor TDMI helper.
    manifest=pd.read_csv(follow.SOURCE/'raw_scan_metrics.csv').query("dataset=='aachen_itbs'")
    ws=tms.load_reference(ROOT)['workspace'];edges=np.column_stack(np.triu_indices(454,1))
    rng=np.random.default_rng(2026112529);selected=rng.choice(len(edges),80,replace=False)
    errors=[];hrf=[];hashes={}
    signal=pd.read_csv(out/'signal_decomposition_scans.csv').set_index(['branch','subject_id','session'])
    for scan in manifest.itertuples():
        path=ROOT/'pipeline/results/timeseries'/(scan.prefix+'_timeseries.npy')
        ts=np.load(path)
        sigma=tn.pg.lagged_sigma_from_timeseries(ts,1,edges[selected])
        a=tn.stable_phiid_atoms(sigma)
        tdmi=tn.pg.tdmi_from_sigma4(sigma)
        errors.append(float(abs(a['tdmi']-tdmi).max()))
        np.testing.assert_allclose(a['tdmi'],tdmi,rtol=1e-9,atol=1e-10)
        path=ROOT/'pipeline/results/timeseries'/(scan.prefix+'_hrf.mat')
        hashes[str(path.relative_to(ROOT))]=tms.sha256(path)
        with h5py.File(path,'r') as f:params=f['PARA'][:]
        assert params.shape==(454,3) and np.isfinite(params).all()
        row={'branch':'hrf','subject_id':scan.subject_id,'session':scan.session,'mean_fd':scan.mean_fd}
        for scope,mask in [('within_workspace',ws),('outside_workspace',~ws)]:
            for k,name in enumerate(['height','peak_seconds','width_seconds']):row[scope+'__'+name]=params[mask,k].mean()
        for name in ['height','peak_seconds','width_seconds']:row['difference__'+name]=row['outside_workspace__'+name]-row['within_workspace__'+name]
        hrf.append(row)
    hrf=pd.DataFrame(hrf);metrics=[m for m in hrf if '__' in m]
    hp=follow.pair(hrf,metrics);pooled=hp.groupby('subject_id').mean(numeric_only=True)
    hrf.to_csv(out/'hrf_parameter_scans.csv',index=False)
    follow.summarize(pooled[metrics],'hrf_parameters9').to_csv(out/'hrf_parameter_tests.csv',index=False)
    association=[]
    for scope in ['outside_workspace','difference']:
        for name in ['height','peak_seconds','width_seconds']:
            x=pooled[scope+'__'+name]
            frame=pd.DataFrame({'y':d.xs('deconvolved')[scope+'__phi'],'intercept':1.,'fd':pooled.delta_fd,'hrf_change':(x-x.mean())/x.std(ddof=1)})
            association.append({'outcome':scope+'__phi','hrf_parameter':name,**hc3(frame,'y',['intercept','fd','hrf_change'],'hrf_change')})
    association=pd.DataFrame(association);association['p_holm6']=holm(association.p)
    association.to_csv(out/'hrf_phi_associations.csv',index=False)
    # Independently calculate the direct deconvolution interaction tests.
    arr=pd.DataFrame(method).to_numpy();signs=2*((np.arange(32768)[:,None]>>np.arange(15))&1)-1
    randomized=stats.ttest_1samp(signs[:,:,None]*arr[None,:,:],0,axis=1).statistic
    observed=abs(stats.ttest_1samp(arr,0,axis=0).statistic)
    adjusted=(abs(randomized).max(1)[:,None]>=observed[None]-1e-10).mean(0)
    np.testing.assert_allclose(adjusted,contrasts.p_maxT_family,atol=1e-12)
    source=json.loads((out/'provenance.json').read_text())
    for p,h in source['input_sha256'].items():assert tms.sha256(ROOT/p)==h
    for p,h in hashes.items():assert tms.sha256(ROOT/p)==h
    record={'status':'passed','max_independent_logdet_error':max(errors),'network_partition_identity':'passed',
            'independent_method_interaction_maxT':'passed','all_inputs_unchanged':True,'hrf_input_sha256':hashes,
            'script_sha256':tms.sha256(Path(__file__))}
    (out/'followup_verification.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
    print(contrasts[['metric','estimate','p','p_maxT_family']].to_string(index=False))
    print('HRF corrected minima:',pd.read_csv(out/'hrf_parameter_tests.csv').p_maxT_family.min(),association.p_holm6.min())
    print('VERIFICATION PASSED, independent logdet error',max(errors))


if __name__=='__main__':
    with threadpool_limits(limits=1):main()
