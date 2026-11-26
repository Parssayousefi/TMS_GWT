"""Read public candidate metadata only; never download full imaging recordings."""
from __future__ import annotations

import argparse
import collections
import csv
import io
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import requests

OUT = Path(__file__).resolve().parents[1] / 'output' / 'candidate_dataset_audit_20261125'
OUT.mkdir(parents=True, exist_ok=True)
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--cached', action='store_true', help='Use saved metadata only; make no network requests.')
ARGS = parser.parse_args()


def fetch(name, url):
    try:
        if ARGS.cached:
            return name, (OUT / name).read_text(encoding='utf-8-sig')
        response = requests.get(url, timeout=40)
        response.raise_for_status()
        if len(response.content) > 20_000_000:
            raise ValueError('Metadata response exceeds 20 MB')
        (OUT / name).write_bytes(response.content)
        return name, response.text
    except Exception as exc:
        return name, {'error': str(exc)}


urls = {
    'ds008064_tree.json': 'https://api.github.com/repos/OpenNeuroDatasets/ds008064/git/trees/1.0.0?recursive=1',
    'ds008064_participants.tsv': 'https://raw.githubusercontent.com/OpenNeuroDatasets/ds008064/1.0.0/participants.tsv',
    'ds008064_description.json': 'https://raw.githubusercontent.com/OpenNeuroDatasets/ds008064/1.0.0/dataset_description.json',
    'figshare_tbs_collection.json': 'https://api.figshare.com/v2/collections/5910329',
    'figshare_tbs_articles.json': 'https://api.figshare.com/v2/collections/5910329/articles?page_size=100',
    'dmt_zenodo_record.json': 'https://zenodo.org/api/records/3992359',
    'pdc_completeness.csv': 'https://www.humanconnectome.org/storage/app/media/documentation/PDC1.0/PDC_1.0_Release_Completeness.csv',
}
with ThreadPoolExecutor(max_workers=5) as pool:
    results = dict(pool.map(lambda item: fetch(*item), urls.items()))

summary = {'sources': urls, 'errors': {k: v for k, v in results.items() if isinstance(v, dict)}}
tree_text = results.get('ds008064_tree.json')
if isinstance(tree_text, str):
    tree = json.loads(tree_text)
    paths = [x['path'] for x in tree['tree'] if x['type'] == 'blob']
    rests = [p for p in paths if 'task-rest' in p and p.endswith('.nii.gz')]
    session_subjects = collections.defaultdict(set)
    for path in rests:
        subject, session = path.split('/')[:2]
        session_subjects[session].add(subject)
    pre = session_subjects['ses-1']
    post = session_subjects['ses-2']
    paired = pre & post
    summary['ds008064'] = {
        'commit': tree['sha'], 'truncated': tree.get('truncated'),
        'rest_nifti_files': len(rests),
        'rest_by_session': dict(collections.Counter(p.split('/')[1] for p in rests)),
        'matched_rest_pairs_before_qc': len(paired),
        'matched_rest_subjects': sorted(paired),
        'pre_only_rest_subjects': sorted(pre - post),
        'post_only_rest_subjects': sorted(post - pre),
        'root_files': [p for p in paths if '/' not in p],
        'possible_allocation_files': [p for p in paths if any(k in p.lower() for k in ['participant', 'session', 'clinical', 'treatment', 'group', 'allocation'])],
        'sample_rest_files': rests[:4],
    }
    metadata_paths = [p for p in paths if p.endswith('.json') and ('task-rest' in p or p.endswith('sessions.json'))]
    with ThreadPoolExecutor(max_workers=4) as pool:
        meta = list(pool.map(lambda p: fetch('ds008064_' + p.replace('/', '_'), 'https://raw.githubusercontent.com/OpenNeuroDatasets/ds008064/1.0.0/' + p), metadata_paths))
    values = collections.defaultdict(set)
    relevant_fields = {}
    metadata_errors = {}
    for name, txt in meta:
        if not isinstance(txt, str):
            metadata_errors[name] = txt
            continue
        obj = json.loads(txt)
        for key in ['RepetitionTime', 'TaskName', 'Manufacturer', 'ManufacturersModelName', 'SeriesDescription', 'ProtocolName']:
            values[key].add(str(obj.get(key)))
        hits = {k: v for k, v in obj.items() if any(term in k.lower() for term in ['treatment', 'group', 'therapy', 'condition', 'randomiz', 'allocation'])}
        if hits:
            relevant_fields[name] = hits
    summary['ds008064']['rest_metadata_values'] = {k: sorted(v) for k, v in values.items()}
    summary['ds008064']['allocation_fields_in_rest_metadata'] = relevant_fields
    summary['ds008064']['rest_metadata_errors'] = metadata_errors
    participants = results.get('ds008064_participants.tsv')
    if isinstance(participants, str):
        rows = list(csv.DictReader(io.StringIO(participants), delimiter='\t'))
        summary['ds008064']['participant_columns'] = list(rows[0])
        summary['ds008064']['participant_rows'] = len(rows)
        summary['ds008064']['post_madrs_available'] = sum(r.get('MADRS Post') not in ['NA', 'n/a', ''] for r in rows)
        summary['ds008064']['participants_without_rest'] = sorted({r['participant_id'] for r in rows} - (pre | post))

for key in ['figshare_tbs_collection.json', 'figshare_tbs_articles.json', 'dmt_zenodo_record.json']:
    txt = results.get(key)
    if not isinstance(txt, str):
        continue
    obj = json.loads(txt)
    if isinstance(obj, list):
        summary[key] = [{k: a.get(k) for k in ['id', 'title', 'url', 'url_public_html']} for a in obj]
    elif key.startswith('figshare'):
        summary[key] = {k: obj.get(k) for k in ['id', 'title', 'description', 'doi', 'url_public_html']}
    else:
        files = obj.get('files', [])
        summary[key] = {
            'title': obj.get('metadata', {}).get('title'),
            'license': obj.get('metadata', {}).get('license'),
            'access_right': obj.get('metadata', {}).get('access_right'),
            'file_count': len(files), 'total_bytes': sum(f['size'] for f in files),
            'scale_file': [f for f in files if 'scale' in f['key'].lower()],
        }
pdc = results.get('pdc_completeness.csv')
if isinstance(pdc, str):
    rows = list(csv.DictReader(io.StringIO(pdc)))
    groups = {'e': 'ECT', 'k': 'ketamine', 'h': 'healthy_controls', 's': 'sleep_deprivation'}
    by_group = {key: collections.defaultdict(set) for key in groups}
    for row in rows:
        subject = row['src_subject_id']
        prefix = subject.split('_')[0]
        if prefix not in groups:
            continue  # The second CSV line describes the columns.
        count = row.get('RS-fMRI_Count', '')
        if count and float(count) > 0:
            by_group[prefix][row['visit']].add(subject)
    summary['pdc_completeness'] = {
        'note': 'At least one resting-state run present; no imaging QC or clinical complete-case filtering.',
        'rest_subjects_by_visit': {
            groups[key]: {visit: len(subjects) for visit, subjects in sorted(visits.items())}
            for key, visits in by_group.items()
        },
        'paired_rest_baseline_to_course_end_or_repeat': {
            groups[key]: {
                'baseline_visit': '1', 'post_visit': post_visit,
                'paired_subjects': len(by_group[key]['1'] & by_group[key][post_visit]),
            }
            for key, post_visit in [('e', '2'), ('k', '3'), ('h', '2')]
        },
    }
(OUT / 'summary.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
print(json.dumps(summary, indent=2))
