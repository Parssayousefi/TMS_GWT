"""Bounded delayed-active phase analysis and selected Aachen robustness checks."""
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'output/tms_publication/20261125T145928_921789Z'
sys.path.insert(0,str(ROOT/'pipeline/scripts/analysis'))
sys.path.insert(0,str(ROOT/'setup'))
import tms_publication as tms
from explore_tms_round1 import regression, one_sample


def main():
    out = ROOT/'output/tms_exploration'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out.mkdir(parents=True,exist_ok=False)
    paths = [SOURCE/'participant_scan_metrics.csv', SOURCE/'participant_paired_changes.csv',
             ROOT/'PROVENANCE/TMS_FOLLOWUP_EXPLORATION_PLAN_20261125.json']
    hashes = {str(p.relative_to(ROOT)):tms.sha256(p) for p in paths}
    scans = pd.read_csv(paths[0])
    sham = scans.loc[scans.dataset.eq('sudmex_tms') & scans.arm.eq('sham')]
    changes = []
    for subject,frame in sham.groupby('subject_id'):
        frame = frame.set_index('session')
        if not {'t1','t14'}.issubset(frame.index): continue
        for contrast in ['delayed_active_change','delayed_active_minus_prior_sham_change']:
            if contrast.endswith('sham_change') and 't0' not in frame.index: continue
            rows = [frame.loc[s] for s in (['t1','t14'] if contrast=='delayed_active_change' else ['t0','t1','t14'])]
            assert all(r.valid_roi_signature==rows[0].valid_roi_signature for r in rows)
            item = {'subject_id':subject,'contrast':contrast}
            for metric in [*tms.METRICS,'mean_fd']:
                value = frame.loc['t14',metric]-frame.loc['t1',metric]
                if contrast.endswith('sham_change'): value -= frame.loc['t1',metric]-frame.loc['t0',metric]
                item[metric] = value
            changes.append(item)
    changes = pd.DataFrame(changes)
    changes.to_csv(out/'delayed_active_participant_changes.csv',index=False)
    records = []
    for contrast,frame in changes.groupby('contrast'):
        frame = frame.assign(intercept=1.)
        for metric in tms.METRICS:
            for model in ['unadjusted','motion_adjusted']:
                result = one_sample(frame[metric]) if model=='unadjusted' else regression(frame,metric,['intercept','mean_fd'],'intercept')
                records.append({'dataset':'sudmex_tms','contrast':contrast,'model':model,'metric':metric,**result})
    results = pd.DataFrame(records)
    results['q_block'] = results.groupby(['contrast','model']).p.transform(tms.bh)
    results['q_followup_branch'] = tms.bh(results.p)
    previous = pd.read_csv(ROOT/'output/tms_exploration/20261125T154020_735179Z/complete_new_test_inventory.csv')
    results['q_raw_explorations'] = tms.bh(np.r_[previous.p.to_numpy(),results.p.to_numpy()])[-len(results):]
    results.to_csv(out/'delayed_active_results.csv',index=False)

    pairs = pd.read_csv(paths[1])
    pooled = pairs.loc[pairs.dataset.eq('aachen_itbs')].groupby('subject_id').mean(numeric_only=True)
    rng = np.random.default_rng(2026112515)
    robust, loo = [], []
    for metric in ['global_phi','global_phi_r']:
        values = pooled[metric+'_delta'].to_numpy()
        signs = 2*((np.arange(2**len(values))[:,None]>>np.arange(len(values)))&1)-1
        p = float((abs(signs@values/len(values))>=abs(values.mean())-1e-14).mean())
        boots = values[rng.integers(len(values),size=(10000,len(values)))].mean(axis=1)
        robust.append({'metric':metric,'check':'exact_sign_flip_mean','estimate':values.mean(),'p':p})
        robust.append({'metric':metric,'check':'wilcoxon','estimate':np.median(values),'p':stats.wilcoxon(values,method='exact').pvalue})
        lo,hi = np.quantile(boots,[.025,.975])
        robust.append({'metric':metric,'check':'bootstrap_mean','estimate':values.mean(),'ci_low':lo,'ci_high':hi})
        frame = pd.DataFrame({'outcome':values,'intercept':1.,'delta_fd':pooled.delta_fd.to_numpy()})
        robust.append({'metric':metric,'check':'motion_adjusted_HC3',**regression(frame,'outcome',['intercept','delta_fd'],'intercept')})
        for i,subject in enumerate(pooled.index):
            loo.append({'metric':metric,'omitted_subject':subject,**one_sample(np.delete(values,i))})
    pd.DataFrame(robust).to_csv(out/'aachen_selected_global_robustness.csv',index=False)
    pd.DataFrame(loo).to_csv(out/'aachen_selected_global_leave_one_out.csv',index=False)
    for p,h in hashes.items(): assert tms.sha256(ROOT/p)==h
    record = {'status':'bounded_followup_exploration_complete','created_utc':datetime.now(timezone.utc).isoformat(),
              'participants_by_contrast':changes.groupby('contrast').size().to_dict(),'tests':len(results),
              'input_sha256':hashes,'script_sha256':tms.sha256(Path(__file__)),
              'aachen_sensitivity_selection':'global_phi and global_phi_r selected after initial results; not independent confirmation',
              'limitations':'Delayed active phase has no randomized concurrent control and selective retention. Period/order/expectancy confounding prevents causal efficacy interpretation.'}
    (out/'provenance.json').write_text(json.dumps(record,indent=2)+'\n')
    print('OUTPUT:',out)
    print(results.sort_values('p').head(12)[['contrast','model','metric','n','estimate','ci_low','ci_high','p','q_block','q_followup_branch']].to_string(index=False))
    print(pd.DataFrame(robust).to_string(index=False))


if __name__=='__main__': main()
