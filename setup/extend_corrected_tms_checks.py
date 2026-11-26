"""Complete four-scope motion, all-post and independent-reference sensitivities."""
from pathlib import Path
import sys
import argparse
import pandas as pd
from threadpoolctl import threadpool_limits
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'setup'),str(ROOT/'pipeline/scripts/analysis')]
import analyse_tms_workspace_scope as ws
import analyse_corrected_tms as core


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--analysis',type=Path,required=True);args=parser.parse_args();out=args.analysis.resolve()
    plan=pd.read_csv(ws.ALLPOST/'raw_participant_changes.csv').query("dataset == 'sudmex_tms'")
    allpost=[];motion=[];seeds=[]
    families={m:m.split('__')[0] if '__' in m else ws.tms.METRIC_INFO[m][0] for m in ws.METRICS}
    for branch in ['raw','normalized']:
        scans=pd.read_csv(out/(branch+'_scan_metrics.csv')).set_index(ws.KEY)
        records=[]
        for p in plan.itertuples():
            pre=scans.loc[(p.dataset,p.subject_id,'t0')]
            post=scans.loc[[(p.dataset,p.subject_id,s) for s in p.post_sessions.split('|')]]
            assert post.valid_roi_signature.eq(pre.valid_roi_signature).all()
            record={k:getattr(p,k) for k in ['dataset','subject_id','arm','n_pre','n_post','post_sessions']}
            record['delta_fd']=post.mean_fd.mean()-pre.mean_fd
            for m in ws.METRICS:record.update({m+'_baseline':pre[m],m+'_followup':post[m].mean(),m+'_delta':post[m].mean()-pre[m]})
            records.append(record)
        changes=pd.DataFrame(records);changes.to_csv(out/(branch+'_all_post_all40_changes.csv'),index=False)
        result=ws.paired_family_tests(changes,ws.METRICS,families,2026112527);result['branch']=branch;allpost.append(result)
        pairs=pd.read_csv(out/(branch+'_paired_changes.csv'));pooled=ws.aggregate_aachen(pairs)
        for m in ws.INFO_METRICS:
            frame=pd.DataFrame({'y':pooled[m+'_delta'],'intercept':1.,'delta_fd':pooled.delta_fd})
            motion.append({'branch':branch,'metric':m,'family':m.split('__')[0],**ws.hc3(frame,'y',['intercept','delta_fd'],'intercept')})
    pd.concat(allpost,ignore_index=True).to_csv(out/'sudmex_all_post_all40_results.csv',index=False)
    motion=pd.DataFrame(motion);motion['p_holm_scope5']=motion.groupby(['branch','family']).p.transform(ws.holm)
    motion['p_holm20']=motion.groupby('branch').p.transform(ws.holm);motion.to_csv(out/'aachen_all_spatial_motion.csv',index=False)
    seed_scans=pd.read_csv(out/'normalization/scan_metrics_all_seeds_scoped.csv')
    acute_plan=pd.read_csv(ws.RAW/'participant_paired_changes.csv').query("dataset == 'aachen_itbs'")
    for seed,block in seed_scans.groupby('reference_seed'):
        paired=ws.make_pairs(block,acute_plan);r=ws.test_aachen(ws.aggregate_aachen(paired),ws.METRICS);r['reference_seed']=seed;seeds.append(r)
    pd.concat(seeds,ignore_index=True).to_csv(out/'normalization_seed_tests.csv',index=False)
    core.save_json(out/'extended_checks_receipt.json',{'status':'complete','source_sha256':ws.tms.sha256(Path(__file__)),'families':'Full40 all-post endpoints; motion5 per scope and20 combined; same40 endpoints at each independent reference seed; seeds are numerical sensitivities.'})
    print('Extended checks complete',flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
