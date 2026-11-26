"""Descriptive temporal-information context for the all-follow-up contrast.

No additional hypothesis search: verify the Phi identity and summarize its
total-information and self-predictability terms using the same visit weights.
"""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'pipeline/scripts/analysis'))
import tms_numit as tn


def main():
    parser=argparse.ArgumentParser();parser.add_argument('run_directory',type=Path)
    out=parser.parse_args().run_directory.resolve()
    assert ROOT in out.parents and 'tms_all_post_average' in out.parts
    scans=pd.read_csv(out/'raw_scan_metrics_used.csv')
    contributions=pd.read_csv(out/'scan_contributions.csv')
    edges=np.column_stack(np.triu_indices(454,1))
    rows=[];maximum=0.
    for i,scan in enumerate(scans.itertuples()):
        ts=np.load(ROOT/'pipeline/results/timeseries'/(scan.prefix+'_timeseries_deconv.npy'))
        a,b=ts[:-1],ts[1:]
        a=(a-a.mean(0))/a.std(0);b=(b-b.mean(0))/b.std(0)
        ac=(a*b).mean(0); self_mi=-.5*np.log1p(-ac**2)
        tdmi=tn.pg.tdmi_from_sigma4(tn.pg.lagged_sigma_from_timeseries(ts,1,edges)).mean()
        error=abs(tdmi-2*self_mi.mean()-scan.global_phi)
        maximum=max(maximum,float(error));assert error<1e-10
        rows.append({**{k:getattr(scan,k) for k in ['dataset','subject_id','session','arm']},
                     'mean_tdmi':tdmi,'mean_self_mi':self_mi.mean(),'mean_autocorrelation':ac.mean(),
                     'global_phi':scan.global_phi})
        if (i+1)%50==0: print(f'Temporal-information identity: {i+1}/{len(scans)}',flush=True)
    table=pd.DataFrame(rows)
    table.to_csv(out/'temporal_information_scans.csv',index=False)
    table=table.merge(contributions[['dataset','subject_id','session','phase']],on=['dataset','subject_id','session'],validate='one_to_one')
    metrics=['mean_tdmi','mean_self_mi','mean_autocorrelation','global_phi']
    means=table.groupby(['dataset','subject_id','arm','phase'])[metrics].mean().unstack('phase')
    delta=pd.DataFrame({m:means[(m,'post')]-means[(m,'pre')] for m in metrics}).reset_index()
    delta.to_csv(out/'temporal_information_participant_changes.csv',index=False)
    delta.groupby(['dataset','arm'])[metrics].agg(['mean','std','count']).to_csv(out/'temporal_information_descriptive_arm_summary.csv')
    delta.groupby('dataset')[metrics].mean().to_csv(out/'temporal_information_descriptive_summary.csv')
    (out/'temporal_identity_verification.json').write_text(json.dumps({'scans':len(scans),'maximum_identity_error':maximum,
                                                                    'new_hypothesis_tests':0},indent=2)+'\n')
    print(delta.groupby('dataset')[metrics].mean().to_string())


if __name__=='__main__':
    with threadpool_limits(limits=1):main()
