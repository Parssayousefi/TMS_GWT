"""Fresh corrected-Aachen conditional references, same frozen NuMIT method.

Use reviewed helpers via an isolated namespace with explicit input routing. No
historical Aachen null/cache or metrics are reused. SUDMEX normalization is reused
only after confirming the original input hashes, because its inputs are unchanged.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import json
import sys
import types
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'setup'),str(ROOT/'pipeline/scripts/analysis')]
import analyse_corrected_tms as core
import analyse_tms_workspace_scope as ws
import tms_publication as tms


def verify_reused_sudmex(previous):
    """Guard acute and later longitudinal inputs, cached TDMI and null tables."""
    manifest=previous/'provenance.json'
    provenance=json.loads(manifest.read_text(encoding='utf-8'))
    checked={rel:digest for rel,digest in provenance['input_sha256'].items() if 'sudmex' in rel.lower()}
    assert len(checked)>=941
    for rel,digest in checked.items():
        assert tms.sha256(ROOT/rel)==digest, 'Changed reused SUDMEX input: '+rel
    return checked


def runtime(out):
    original=ROOT/'setup/run_tms_numit_sensitivity.py'
    source=original.read_text(encoding='utf-8')
    source=source.replace("ROOT = Path(__file__).resolve().parents[1]",'ROOT = Path('+repr(str(ROOT))+')')
    source=source.replace("ROOT/'pipeline/results/timeseries'", "CORRECTED/'pipeline/results/timeseries'")
    source=source.replace("ROOT/'pipeline/results/phiid_matrices'", "CORRECTED/'pipeline/results/phiid_matrices'")
    generated=out/'runtime_numit_helpers.py';generated.write_text(source,encoding='utf-8')
    module=types.ModuleType('corrected_numit_runtime');module.__file__=str(generated)
    exec(compile(source,str(generated),'exec'),module.__dict__)
    module.CORRECTED=core.CORRECTED
    module.tms=types.SimpleNamespace(**{name:getattr(tms,name) for name in dir(tms) if not name.startswith('__')})
    module.tms.DATASETS=('aachen_itbs',)
    return module,original,generated


def main():
    p=argparse.ArgumentParser();p.add_argument('--analysis',required=True,type=Path);args=p.parse_args()
    analysis=args.analysis.resolve();out=analysis/'normalization';out.mkdir(exist_ok=True)
    assert tms.sha256(core.RELEASE)==core.RELEASE_SHA
    scans=pd.read_csv(analysis/'raw_scan_metrics.csv');aachen=scans.loc[scans.dataset.eq('aachen_itbs')].sort_values(ws.KEY).reset_index(drop=True)
    assert len(aachen)==60
    helper,original,generated=runtime(out)
    prepared=helper.prepare(out,aachen)
    for rel,digest in prepared['input_sha256'].items():assert tms.sha256(ROOT/rel)==digest
    refs=helper.build_references(out,prepared)
    interpolation=helper.validate_interpolation(out,aachen,refs)
    calibration=helper.calibrate_finite_records(out,prepared,refs)
    tables,pairs=helper.extract_metrics(out,aachen,refs)
    helper.monte_carlo_stability(out,tables,pairs)
    edges,masks=ws.edge_scopes(tms.load_reference(ROOT)['workspace'])
    # Extend the same conditional percentiles to all four fixed spatial scopes.
    for idx,row in tables.iterrows():
        seeds=helper.SEEDS if str(row.reference_seed)=='combined' else [int(row.reference_seed)]
        tdmi=np.load(out/'cache'/(row.prefix+'_tdmi.npy'))
        for measure in helper.MEASURES:
            observed=np.load(core.CORRECTED/'pipeline/results/phiid_matrices'/(row.prefix+'_'+measure+'.npy'))[edges[:,0],edges[:,1]]
            q=np.mean([helper.tn.conditional_cdf(refs[('aachen_itbs',seed)][0],refs[('aachen_itbs',seed)][1][measure],tdmi,observed) for seed in seeds],axis=0)
            for scope,value in ws.scoped_means(q,masks).items():tables.loc[idx,scope+'__'+measure]=value
    tables.to_csv(out/'scan_metrics_all_seeds_scoped.csv',index=False)
    normalized_aachen=tables.loc[tables.reference_seed.eq('combined')].drop(columns='reference_seed')
    previous=ROOT/'output/tms_workspace_scope'/(ROOT/'output/tms_workspace_scope/LATEST.txt').read_text().strip()
    old_prepared=json.loads((ws.NUMIT/'cache/prepared.json').read_text())
    checked={rel:digest for rel,digest in old_prepared['input_sha256'].items() if 'dataset-sudmex_tms_' in rel}
    assert checked and all(tms.sha256(ROOT/rel)==digest for rel,digest in checked.items())
    checked.update(verify_reused_sudmex(previous))
    sud=pd.read_csv(previous/'normalized_scan_metrics.csv').query("dataset == 'sudmex_tms'")
    assert not sud.duplicated(ws.KEY).any()
    normalized=pd.concat([sud,normalized_aachen],ignore_index=True)
    normalized.to_csv(analysis/'normalized_scan_metrics.csv',index=False)
    core.run_statistics(analysis,{'raw':scans,'normalized':normalized})
    for rel,digest in prepared['input_sha256'].items():assert tms.sha256(ROOT/rel)==digest
    core.save_json(out/'verification.json',{'completed_utc':datetime.now(timezone.utc).isoformat(),'fresh_aachen_baseline_pool':True,
        'interpolation_max_absolute_error':interpolation,'calibration':calibration,'corrected_input_sha256':prepared['input_sha256'],
        'sudmex_unchanged_input_sha256':checked,'sudmex_normalized_source':str(previous/'normalized_scan_metrics.csv'),
        'source_sha256':{str(f.relative_to(ROOT)):tms.sha256(f) for f in [Path(__file__),original,generated,previous/'normalized_scan_metrics.csv']}})
    core.save_json(analysis/'status.json',{'raw':'complete','normalized':'complete','target_flow':'separate parallel analysis','manuscript':'pending'})
    print('Corrected normalization and all statistics complete',flush=True)


if __name__=='__main__':
    with threadpool_limits(limits=2):main()
