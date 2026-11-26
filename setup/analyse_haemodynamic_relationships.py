"""Bounded HRF association follow-up requested after corrected TMS review.

One new whole-brain scope joins the six existing associations: nine tests.
Participant influence and rank correlations are sensitivities, not replications.
"""
from pathlib import Path
import json
import sys
import numpy as np
import pandas as pd
from scipy import stats
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'setup'),str(ROOT/'pipeline/scripts/analysis')]
import analyse_tms_workspace_scope as ws
import analyse_corrected_tms as core
RUN=ROOT/'output/corrected_tms_analysis/20261126T112054_284547Z'
OUT=RUN/'haemodynamic_followup'


def fit(data):
    frame=pd.DataFrame({'y':data.y,'intercept':1.,'fd':data.fd,'x_z':(data.x-data.x.mean())/data.x.std(ddof=1)})
    return ws.hc3(frame,'y',['intercept','fd','x_z'],'x_z')


def main():
    OUT.mkdir(exist_ok=True)
    core.save_json(OUT/'analysis_specification.json',{'question':'Do corrected haemodynamic associations generalize to whole-brain Phi and withstand participant omission?',
        'prompted_by':'User haemodynamic-narrative question after initial six association results were known.',
        'new_tests':['whole_brain Phi with whole-brain mean HRF height','whole_brain Phi with whole-brain mean HRF peak latency','whole_brain Phi with whole-brain mean HRF width'],
        'family':'All nine associations: three HRF parameters at whole-brain, outside-workspace and outside-minus-inside scopes. Two-sided Holm9.',
        'sensitivity':'Same nine models repeated omitting each participant; unadjusted Spearman with Holm9, descriptive dependence check.',
        'interpretation':'Same-signal association, not mediation, vascular mechanism, causal effect, or validation of the superseded global Phi result.'})
    files=[RUN/'sensitivities/hrf_parameter_scans.csv',RUN/'sensitivities/signal_participant_changes.csv',RUN/'sensitivities/hrf_phi_associations.csv']
    hrf=pd.read_csv(files[0]);info=pd.read_csv(files[1]).query("branch == 'deconvolved'")
    parameters=['height','peak_seconds','width_seconds']
    for name in parameters:hrf['whole_brain__'+name]=(235*hrf['within_workspace__'+name]+219*hrf['outside_workspace__'+name])/454
    rows=[]
    for sid,block in hrf.groupby('subject_id'):
        block=block.set_index('session');r={'subject_id':sid}
        for scope in ['whole_brain','outside_workspace','difference']:
            for parameter in parameters:
                key=scope+'__'+parameter
                r[key]=np.mean([block.loc['post'+str(i),key]-block.loc['pre'+str(i),key] for i in [90,120]])
        rows.append(r)
    hp=pd.DataFrame(rows).set_index('subject_id').sort_index();ip=info.groupby('subject_id').mean(numeric_only=True).sort_index()
    assert hp.index.equals(ip.index) and len(hp)==15
    results=[];loo=[];ranks=[];observations=[]
    for scope in ['whole_brain','outside_workspace','difference']:
        for parameter in parameters:
            data=pd.DataFrame({'y':ip[scope+'__phi'],'fd':ip.delta_fd,'x':hp[scope+'__'+parameter]})
            label={'scope':scope,'parameter':parameter}
            results.append({**label,**fit(data)})
            # HRF summaries contain mathematical ties obscured by summation
            # roundoff. Fix rank ties at far below their measurement resolution.
            rho,p=stats.spearmanr(np.round(data.x,12),np.round(data.y,12));ranks.append({**label,'rho':rho,'p':p,'n':len(data),'model':'unadjusted Spearman sensitivity; rank ties at 12 decimals'})
            for sid in data.index:loo.append({**label,'omitted_subject':sid,**fit(data.drop(sid))})
            for sid,r in data.iterrows():observations.append({**label,'subject_id':sid,'hrf_change':r.x,'phi_change':r.y,'delta_fd':r.fd})
    result=pd.DataFrame(results);result['p_holm9']=ws.holm(result.p);result.to_csv(OUT/'association_results.csv',index=False)
    ranked=pd.DataFrame(ranks);ranked['p_holm9']=ws.holm(ranked.p);ranked.to_csv(OUT/'rank_sensitivity.csv',index=False)
    omitted=pd.DataFrame(loo);omitted['p_holm9']=omitted.groupby('omitted_subject').p.transform(ws.holm);omitted.to_csv(OUT/'leave_one_out.csv',index=False)
    summary=omitted.groupby(['scope','parameter']).agg(estimate_min=('estimate','min'),estimate_max=('estimate','max'),p_min=('p','min'),p_max=('p','max'),p_holm9_max=('p_holm9','max'),n_corrected_omissions=('p_holm9',lambda x:int((x<.05).sum())))
    summary.to_csv(OUT/'influence_summary.csv')
    pd.DataFrame(observations).to_csv(OUT/'participant_associations.csv',index=False)
    original=pd.read_csv(files[2]).copy();original['scope']=original.outcome.str.replace('__phi','',regex=False);original=original.rename(columns={'hrf_parameter':'parameter'}).set_index(['scope','parameter'])
    got=result.set_index(['scope','parameter']).loc[original.index]
    np.testing.assert_allclose(got[['estimate','ci_low','ci_high','p']],original[['estimate','ci_low','ci_high','p']],rtol=1e-9,atol=1e-12)
    core.save_json(OUT/'verification.json',{'status':'complete','existing_six_models_reproduced':True,'new_models':3,'influence_fits':135,
        'source_sha256':{str(p.relative_to(ROOT)):ws.tms.sha256(p) for p in [Path(__file__),*files]}})
    print(result[['scope','parameter','estimate','ci_low','ci_high','p','p_holm9']].to_string(index=False))
    print(summary.to_string())


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
