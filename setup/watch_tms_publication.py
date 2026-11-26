"""Run the publication notebook once per completed TMS dataset.

The production pipeline must already be running. This monitor only reads its
outputs, audits saved HRF files, and executes the separate publication notebook.
It never starts or modifies preprocessing and never reads psychedelic results.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'pipeline/scripts/analysis'))
import tms_publication as tms
from run_tms_processing import write_json

LOGS = ROOT / 'logs/tms_pipeline_20261125/publication_watch'
STATE = LOGS / 'state.json'
BULK_STATE = ROOT / 'logs/tms_pipeline_20261125/bulk/status.json'
OUTPUTS = ROOT / 'output/tms_publication'
NOTEBOOK = ROOT / 'output/jupyter-notebook/TMS_Global_Workspace_Publication.ipynb'


def now():
    return datetime.now(timezone.utc).isoformat()


def emit(event, **fields):
    print(json.dumps({'utc': now(), 'event': event, **fields}), flush=True)


def save(state):
    state['updated_utc'] = now()
    write_json(STATE, state)


def progress():
    scans, _ = tms.inventory(ROOT)
    expected = scans[scans.eligible]
    report = {}
    for dataset in tms.DATASETS:
        block = expected[expected.dataset == dataset]
        if block.empty:
            raise ValueError('No expected eligible scans for ' + dataset)
        report[dataset] = {'expected': len(block), **{
            key: int(block[key].sum()) for key in ('regional_ready', 'hrf_ready', 'phiid_ready')}}
    return report


def run_logged(arguments, destination):
    with destination.open('w', encoding='utf-8') as stream:
        subprocess.run([sys.executable, '-X', 'utf8', *arguments], cwd=ROOT,
                       stdout=stream, stderr=subprocess.STDOUT, check=True)


def analyse(dataset):
    emit('dataset_ready', dataset=dataset)
    audit_dir = LOGS / 'hrf_audits' / dataset
    run_logged([str(ROOT / 'setup/audit_tms_hrf.py'), '--dataset', dataset,
                '--out', str(audit_dir)], LOGS / (dataset + '_hrf_audit.log'))
    audit = json.loads((audit_dir / (dataset + '_output_audit.json')).read_text())
    before = set(OUTPUTS.glob('*/run_status.json'))
    run_logged([str(ROOT / 'setup/execute_tms_publication_notebook.py')],
               LOGS / (dataset + '_notebook.log'))
    created = set(OUTPUTS.glob('*/run_status.json')) - before
    if len(created) != 1:
        raise RuntimeError('Expected exactly one successful notebook output directory')
    status_file = created.pop()
    result = json.loads(status_file.read_text())
    if not result['cohort_readiness'].get(dataset):
        raise RuntimeError('Notebook did not confirm the completed cohort: ' + dataset)
    archive = status_file.parent / NOTEBOOK.name
    shutil.copy2(NOTEBOOK, archive)
    record = {'completed_utc': now(), 'output_directory': str(status_file.parent),
              'executed_notebook': str(archive), 'notebook_sha256': tms.sha256(archive),
              'hrf_scans_audited': audit['scans_passed'], 'run_status': result}
    (status_file.parent / 'execution_receipt.json').write_text(
        json.dumps(record, indent=2), encoding='utf-8')
    emit('notebook_complete', dataset=dataset, **record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--status-only', action='store_true', help='Print readiness without writing or executing anything')
    parser.add_argument('--interval', type=float, default=30)
    args = parser.parse_args()
    if args.status_only:
        emit('progress', datasets=progress())
        return
    if not 5 <= args.interval <= 60:
        parser.error('--interval must be between 5 and 60 seconds')
    LOGS.mkdir(parents=True, exist_ok=True)
    state = json.loads(STATE.read_text()) if STATE.exists() else {'started_utc': now(), 'completed': {}}
    # Refuse to silently skip an output that was removed after an earlier run.
    for dataset, record in state['completed'].items():
        if not Path(record['executed_notebook']).is_file():
            raise FileNotFoundError('Missing archived notebook for ' + dataset)
    previous = None
    state['status'] = 'watching'
    state.pop('error', None)
    save(state)
    try:
        while True:
            try:
                report = progress()
            except (yaml.YAMLError, OSError) as exc:
                emit('transient_read_retry', error=str(exc))
                time.sleep(args.interval)
                continue
            if report != previous:
                emit('progress', datasets=report)
                previous = report
            state['progress'] = report
            save(state)
            for dataset, counts in report.items():
                complete = all(counts[key] == counts['expected']
                               for key in ('regional_ready', 'hrf_ready', 'phiid_ready'))
                if complete and dataset not in state['completed']:
                    state['status'] = 'executing_' + dataset
                    save(state)
                    state['completed'][dataset] = analyse(dataset)
                    state['status'] = 'watching'
                    save(state)
            if set(state['completed']) == set(tms.DATASETS):
                state['status'] = 'complete'
                save(state)
                emit('all_datasets_analysed')
                return
            if BULK_STATE.exists():
                try:
                    bulk = json.loads(BULK_STATE.read_text())
                except (OSError, json.JSONDecodeError) as exc:
                    emit('transient_bulk_state_read_retry', error=str(exc))
                    time.sleep(args.interval)
                    continue
                if bulk.get('status') in ('failed', 'error'):
                    raise RuntimeError('Bulk pipeline reported failure: ' + str(bulk))
                if bulk.get('status') == 'complete':
                    if progress() != report:
                        continue  # The final scan may have appeared during this iteration.
                    raise RuntimeError('Pipeline completed without all expected notebook inputs')
            time.sleep(args.interval)
    except Exception as exc:
        state['status'] = 'failed'
        state['error'] = repr(exc)
        save(state)
        emit('monitor_failed', error=repr(exc))
        raise


if __name__ == '__main__':
    main()
