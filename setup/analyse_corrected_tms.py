"""Recompute TMS summaries with explicit, audited corrected Aachen routing.

Historical outputs and production processing code are read-only.  Families and
contrasts follow the existing scope analysis; the full 40-endpoint view is kept.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import json
import sys
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'setup'), str(ROOT/'pipeline/scripts/analysis')]
import tms_publication as tms
import analyse_tms_workspace_scope as ws

CORRECTED = ROOT/'reruns/aachen_skullstrip_skip_20261125'
RELEASE = ROOT/'output/aachen_alignment_review/corrected_cohort/cohort_release_receipt.json'
RELEASE_SHA = 'd02401c42be2ea20c9e6de008a30c9665d92e7151467bb0aa4c873a8dfd81895'


def matrix_root(dataset):
    assert dataset in tms.DATASETS
    return CORRECTED if dataset == 'aachen_itbs' else ROOT


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, default=str)+'\n', encoding='utf-8')


def extract(out):
    assert tms.sha256(RELEASE) == RELEASE_SHA
    reference = tms.load_reference(ROOT)
    inventory, inputs = tms.inventory(ROOT)
    qcpath = CORRECTED/'pipeline/results/qc/aachen_itbs/exclusions.csv'
    qc = pd.read_csv(qcpath, dtype={'subject':str, 'session':str}).set_index(['subject','session'])
    for idx, scan in inventory.loc[inventory.dataset.eq('aachen_itbs')].iterrows():
        if (scan.subject, scan.session) not in qc.index:
            # Subject 02 remains excluded; do not promote a missing scan.
            assert not scan.eligible
            continue
        row = qc.loc[(scan.subject, scan.session)]
        inventory.loc[idx, 'mean_fd'] = row.usable_mean_fd_mm
        inventory.loc[idx, 'qc_status'] = row.status
        inventory.loc[idx, 'eligible'] = row.status in ('INCLUDE', 'TRUNCATE')
        phi = CORRECTED/'pipeline/results/phiid_matrices'
        ts = CORRECTED/'pipeline/results/timeseries'
        inventory.loc[idx,'phiid_ready'] = all((phi/(scan.prefix+'_'+m+'.npy')).exists() for m in tms.MATRIX_NAMES)
        inventory.loc[idx,'regional_ready'] = (ts/(scan.prefix+'_timeseries.npy')).exists()
        inventory.loc[idx,'hrf_ready'] = (ts/(scan.prefix+'_timeseries_deconv.npy')).exists()
    inputs += [qcpath, RELEASE, *reference['input_files']]
    scans, regions = [], []
    edges, masks = ws.edge_scopes(reference['workspace'])
    for dataset in tms.DATASETS:
        base = matrix_root(dataset)
        block, regional, files = tms.extract_scan_metrics(base, inventory.loc[inventory.dataset.eq(dataset)], reference)
        assert len(block) == (60 if dataset == 'aachen_itbs' else 133)
        for idx, scan in block.iterrows():
            for measure in tms.MATRIX_NAMES:
                path = base/'pipeline/results/phiid_matrices'/(scan.prefix+'_'+measure+'.npy')
                values = np.load(path)[edges[:,0],edges[:,1]]
                for scope, value in ws.scoped_means(values, masks).items():
                    block.loc[idx,scope+'__'+measure] = value
        scans.append(block); regions.append(regional); inputs += files
    scans = pd.concat(scans,ignore_index=True)
    assert not scans.duplicated(ws.KEY).any()
    historical = pd.read_csv(ws.RAW/'participant_scan_metrics.csv').query("dataset == 'sudmex_tms'").set_index(ws.KEY).sort_index()
    unchanged = scans.query("dataset == 'sudmex_tms'").set_index(ws.KEY).sort_index()
    np.testing.assert_allclose(unchanged[tms.METRICS],historical[tms.METRICS],rtol=1e-12,atol=1e-12)
    inventory.to_csv(out/'scan_inventory.csv',index=False)
    scans.to_csv(out/'raw_scan_metrics.csv',index=False)
    pd.concat(regions,ignore_index=True).to_csv(out/'raw_region_metrics.csv',index=False)
    save_json(out/'raw_input_sha256.json',{str(p.relative_to(ROOT)):tms.sha256(p) for p in set(inputs)})
    return scans, inventory


def analyse_branch(out, scans, plan, branch):
    pairs = ws.make_pairs(scans, plan)
    assert pairs.groupby('dataset').size().to_dict() == {'aachen_itbs':30,'sudmex_tms':40}
    pairs.to_csv(out/(branch+'_paired_changes.csv'), index=False)
    pooled = ws.aggregate_aachen(pairs)
    assert len(pooled) == 15
    pooled.to_csv(out/(branch+'_aachen_participant_changes.csv'))
    for label, metrics in [('main30',ws.MAIN_METRICS),('information20',ws.INFO_METRICS),('all40',ws.METRICS)]:
        table=ws.test_aachen(pooled,metrics);table['branch']=branch
        table.to_csv(out/(branch+'_aachen_'+label+'.csv'), index=False)
    sud=ws.test_sudmex(pairs,ws.MAIN_METRICS);sud['branch']=branch
    sud.to_csv(out/(branch+'_sudmex_acute.csv'),index=False)
    return pairs


def run_statistics(out, branches):
    plan = pd.read_csv(ws.RAW/'participant_paired_changes.csv')
    pairs={branch:analyse_branch(out,scans,plan,branch) for branch,scans in branches.items()}
    ws.robustness(out,pairs)
    ws.clinical(out,pairs)
    ws.all_post(out,branches)
    save_json(out/'statistics_source_sha256.json',{str(p.relative_to(ROOT)):tms.sha256(p) for p in [Path(__file__),Path(ws.__file__),Path(tms.__file__),ws.RAW/'participant_paired_changes.csv',*map(lambda p:ROOT/p,ws.HASHES)]})
    return pairs


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--out',type=Path)
    args=parser.parse_args()
    out=args.out or ROOT/'output/corrected_tms_analysis'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
    out.mkdir(parents=True,exist_ok=True)
    print('OUTPUT:',out,flush=True)
    specification={'created_utc':datetime.now(timezone.utc).isoformat(),'aachen_inputs':str(CORRECTED),'sudmex_inputs':str(ROOT/'pipeline/results'),
        'release_sha256':RELEASE_SHA,'scope':ws.SCOPES,'information_measures':ws.MEASURES,'role_metrics':ws.ROLE_METRICS,
        'primary_comparisons':{'aachen':'Average the two within-visit changes per person; exact joint sign flips, n15.',
                               'sudmex':'Two-week randomized active versus sham; HC3 baseline and baseline+motion adjustment.'},
        'families':'Five information measures per scope; ten whole+workspace; twenty four-scope information; thirty main endpoints; forty scope+role endpoints. Original role subfamilies retained. No selection of endpoints from corrected p values.',
        'longitudinal_sensitivity':'SUDMEX all available post scans averaged within person; later active treatment precludes interpreting original arms as continuing randomized exposure.'}
    save_json(out/'analysis_specification.json',specification)
    scans, _ = extract(out)
    run_statistics(out,{'raw':scans})
    (out.parent/'LATEST.txt').write_text(out.name+'\n')
    save_json(out/'status.json',{'raw':'complete','normalized':'pending','target_flow':'separate parallel analysis','manuscript':'pending'})
    print('Raw analysis complete',flush=True)


if __name__ == '__main__':
    with threadpool_limits(limits=2): main()
