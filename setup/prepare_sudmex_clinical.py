"""Prepare published SUDMEX scores, preserving source values and MRI cohort rules."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'data/sudmex_tms/clinical/source/zenodo-10409461'
OUT = ROOT / 'outputs/tms_clinical_20261125'
SESSION_MAP = {'T0': 't0', 'T1': 't1', 'T1-4': 't14', 'T2': 't2', 'T3': 't3', 'T4': 't4'}
SPECS = [
    ('12_VAS Clinico.xlsx', 'Data', {'vas': 'craving_vas'}, None),
    ('5_CCQ-G.xlsx', 'CCQG', {'ccq_g': 'craving_ccq_general'}, set(range(1, 8))),
    ('5_CCQ-N.xlsx', 'CCQN', {'ccq_n': 'craving_ccq_now'}, set(range(1, 8))),
    ('7_BIS-11.xlsx', 'data', {'tot_score': 'impulsivity_bis_total', 'cog': 'impulsivity_bis_cognitive',
                            'mot': 'impulsivity_bis_motor', 'nonp': 'impulsivity_bis_nonplanning'}, {0, 1, 3, 4}),
]
METRICS = [metric for _, _, mapping, _ in SPECS for metric in mapping.values()]


def records(frame):
    return json.loads(frame.to_json(orient='records'))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((SOURCE / 'download_manifest.json').read_text(encoding='utf-8'))
    for entry in manifest['files']:
        if hashlib.sha256((SOURCE / entry['file']).read_bytes()).hexdigest() != entry['sha256']:
            raise ValueError('Source file changed: ' + entry['file'])
    demo = pd.read_excel(SOURCE / 'Demographic_data.xlsx', sheet_name='DEMOGRAPHIC')
    if demo.rid.duplicated().any() or not demo.group.isin([1, 2]).all():
        raise ValueError('Invalid demographic identities or randomized arms')
    demo = demo.set_index('rid')
    bids = pd.read_csv(ROOT / 'data/sudmex_tms/participants.tsv', sep='\t')
    bids['rid'] = bids.participant_id.str.removeprefix('sub-').astype(int)
    bids = bids.set_index('rid')
    for rid in bids.index:
        if any(demo.loc[rid, c] != bids.loc[rid, b] for c, b in [('group', 'group'), ('q1_age', 'age'), ('q1_sex', 'sex')]):
            raise ValueError(f'Demographic/BIDS conflict: {rid}')
    qc = pd.read_csv(ROOT / 'pipeline/results/qc/sudmex_tms/exclusions.csv', dtype={'subject': str, 'session': str})
    if qc.duplicated(['subject', 'session']).any():
        raise ValueError('Duplicate MRI QC records')
    qc_lookup = qc.set_index(['subject', 'session']).status.to_dict()
    # Motion QC was run on the configured MRI cohort; do not expand that cohort
    # merely because additional randomized arms are available in clinical files.
    cohort = {int(x) for x in qc.subject}
    if not cohort.issubset(bids.index) or not bids.loc[sorted(cohort), 'exclusion'].eq(0).all():
        raise ValueError('MRI QC cohort conflicts with BIDS participant exclusions')
    combined = None
    issues = []
    for filename, sheet, mapping, levels in SPECS:
        raw = pd.read_excel(SOURCE / filename, sheet_name=sheet)
        if raw.duplicated(['rid', 'stage']).any() or not raw.stage.isin(SESSION_MAP).all():
            raise ValueError(f'Duplicate/unknown clinical visits: {filename}')
        if not raw['group'].eq(raw.rid.map(demo.group)).all():
            raise ValueError(f'Clinical arm conflict: {filename}')
        frame = raw[['rid', 'stage', *mapping]].rename(columns=mapping).copy()
        qcols = [c for c in raw if c.startswith('q') and c[1:].isdigit()]
        flags = []
        for _, row in raw.iterrows():
            reasons = []
            for item in qcols:
                if pd.isna(row[item]):
                    reasons.append(f'missing_item:{item}')
                elif row[item] not in levels:
                    reasons.append(f'invalid_item:{item}={row[item]:g}')
            if 'tot_score' in mapping and row.tot_score != row[['cog', 'mot', 'nonp']].sum():
                reasons.append('total_subscale_mismatch')
            if 'vas' in mapping and not 0 <= row.vas <= 10:
                reasons.append('total_out_of_range')
            for col in mapping:
                if pd.isna(row[col]):
                    reasons.append('missing_score:' + col)
            flag = ';'.join(reasons) or 'ok'
            flags.append(flag)
            if reasons:
                issues.append({'file': filename, 'subject_id': f'sudmex_tms_{int(row.rid):03d}',
                               'stage': row.stage, 'issue': flag, 'policy': 'Preserve published score; exclude affected measure by default in loader.'})
        for metric in mapping.values():
            frame[metric + '_qc'] = flags
        frame = frame.set_index(['rid', 'stage'])
        combined = frame if combined is None else combined.join(frame, how='outer', validate='one_to_one')
    sessions = combined.reset_index()
    for metric in METRICS:
        sessions[metric + '_qc'] = sessions[metric + '_qc'].fillna('missing_assessment')
    sessions['dataset'] = 'sudmex_tms'
    sessions['subject_num'] = sessions.rid.astype(int)
    sessions['subject_id'] = sessions.rid.map(lambda n: f'sudmex_tms_{int(n):03d}')
    sessions['bids_subject_id'] = sessions.rid.map(lambda n: f'sub-{int(n):03d}')
    sessions['session'] = sessions.stage.map(SESSION_MAP)
    sessions['arm'] = sessions.rid.map(demo.group).map({1: 'sham', 2: 'active'})
    sessions['age'] = sessions.rid.map(demo.q1_age)
    sessions['sex'] = sessions.rid.map(demo.q1_sex).map({1: 'male', 2: 'female'})
    sessions['analysis_cohort'] = sessions.rid.isin(cohort).astype(int)
    sessions['mri_qc'] = [qc_lookup.get((f'{int(r):03d}', s), 'NOT_IN_MRI_COHORT') for r, s in zip(sessions.rid, sessions.session)]
    sessions['mri_qc_usable'] = sessions.mri_qc.isin(['INCLUDE', 'TRUNCATE']).astype(int)
    sessions['phase'] = sessions.session.map({'t0': 'baseline', 't1': 'double_blind_acute', 't14': 'open_label_crossover',
                                               't2': 'open_label_followup', 't3': 'open_label_followup', 't4': 'open_label_followup'})
    sessions['source_visit'] = sessions.stage
    sessions['_order'] = sessions.stage.map({s: i for i, s in enumerate(SESSION_MAP)})
    sessions = sessions.sort_values(['subject_num', '_order'])
    session_cols = ['dataset', 'subject_id', 'subject_num', 'bids_subject_id', 'arm', 'session', 'phase', 'age', 'sex',
                    'analysis_cohort', 'mri_qc', 'mri_qc_usable'] + METRICS + [m + '_qc' for m in METRICS] + ['source_visit']
    sessions = sessions[session_cols].reset_index(drop=True)
    by_visit = sessions.set_index(['subject_id', 'session'])
    effects = []
    for _, baseline in sessions[(sessions.session == 't0') & (sessions.analysis_cohort == 1)].iterrows():
        sid = baseline.subject_id
        follow = by_visit.loc[(sid, 't1')] if (sid, 't1') in by_visit.index else None
        item = {c: baseline[c] for c in ['dataset', 'subject_id', 'subject_num', 'bids_subject_id', 'arm', 'age', 'sex']}
        item.update(comparison='t1_minus_t0', baseline_session='t0', followup_session='t1',
                    mri_pair_usable=int(baseline.mri_qc_usable == 1 and follow is not None and follow.mri_qc_usable == 1),
                    baseline_mri_qc=baseline.mri_qc, followup_mri_qc=follow.mri_qc if follow is not None else 'MISSING', analysis_cohort=1)
        for metric in METRICS:
            b = baseline[metric]
            f = follow[metric] if follow is not None else np.nan
            item[metric + '_t0'] = b
            item[metric + '_t1'] = f
            item[metric + '_delta'] = f - b
            flags = []
            if baseline[metric + '_qc'] != 'ok':
                flags.append('t0:' + baseline[metric + '_qc'])
            if follow is None:
                flags.append('t1:missing_assessment')
            elif follow[metric + '_qc'] != 'ok':
                flags.append('t1:' + follow[metric + '_qc'])
            item[metric + '_qc'] = '|'.join(flags) or 'ok'
        effects.append(item)
    effects = pd.DataFrame(effects)
    first = ['dataset', 'subject_id', 'subject_num', 'bids_subject_id', 'arm', 'comparison', 'baseline_session', 'followup_session',
             'age', 'sex', 'analysis_cohort', 'mri_pair_usable', 'baseline_mri_qc', 'followup_mri_qc']
    effects = effects[first + [m + suffix for m in METRICS for suffix in ['_t0', '_t1', '_delta', '_qc']]]
    sessions.to_csv(ROOT / 'data/master_clinical_sessions.csv', index=False, na_rep='')
    effects.to_csv(ROOT / 'data/master_clinical_effects.csv', index=False, na_rep='')
    pd.DataFrame(issues).to_csv(OUT / 'source_issues.csv', index=False)
    dictionary = [
        ['Scope', 'TMS clinical outcomes', 'SUDMEX only. Aachen has no matching craving/impulsivity assessments in these sources.'],
        ['Source', 'Zenodo version 4.2', 'https://zenodo.org/records/10409461'],
        ['Identifiers', 'dataset + subject_id + comparison', 'GWT-style subject_id = sudmex_tms_001. BIDS sub-001 maps to source rid 1.'],
        ['Sessions', 'T0/T1/T1-4/T2/T3/T4', 'Mapped to t0/t1/t14/t2/t3/t4. T1-4 is the sham crossover visit, not T1.'],
        ['Arms', 'active / sham', 'Randomized assignment from source group 2 / 1, cross-checked against BIDS. Later visits do not retain a sham-controlled interpretation.'],
        ['Clinical changes', 't1_minus_t0', 'Two-week score minus baseline. Negative changes indicate lower craving or impulsivity. No drug-minus-placebo calculation or 0-1 normalization.'],
        ['Cohort', '44 configured participants', 'Existing MRI cohort is retained. All 53 clinical participants remain in Clinical sessions. New clinical arm information does not silently reinstate MRI exclusions.'],
        ['MRI', 'mri_pair_usable', '1 means both visits pass existing motion QC (INCLUDE or TRUNCATE). This is eligibility, not a claim that PhiID processing is complete.'],
        ['Missing', 'blank numeric cells', 'No imputation. Zero is a valid score and is never treated as missing. BIS assessment missing at one later visit.'],
        ['Quality', '*_qc', 'ok means no detected item-range/missingness/total-subscale issue; it is not full psychometric revalidation. Other flags preserve the reported score but mask the affected outcome in the default loader.'],
        ['Scores', 'Published totals', 'Totals and BIS subscales copied unchanged. We do not sum raw items or apply a second reverse-scoring operation. BIS dictionary has ambiguous response labels; item values 0,1,3,4 are used only for range checks.'],
        ['craving_vas', '0-10 cm', 'Published clinical visual analogue craving score; genuine zeros retained.'],
        ['craving_ccq_general', '45-item reported total', 'CCQ-General; higher indicates more craving. Items have source range 1-7.'],
        ['craving_ccq_now', '45-item reported total', 'CCQ-Now; higher indicates more craving. Items have source range 1-7.'],
        ['impulsivity_bis_total', 'Published BIS-11 total', 'Spanish source coding 0,1,3,4. Higher indicates more impulsivity; do not substitute a 1-4 scoring system.'],
        ['impulsivity_bis_cognitive', 'Published cognitive subscale', 'Exploratory outcome; preserved without rescoring.'],
        ['impulsivity_bis_motor', 'Published motor subscale', 'Exploratory outcome; preserved without rescoring.'],
        ['impulsivity_bis_nonplanning', 'Published nonplanning subscale', 'Exploratory outcome; preserved without rescoring.'],
        ['Analysis policy', 'Before outcome testing', 'Choose primary endpoints, multiplicity control and the handling of flagged records in the analysis plan. The sheet does not establish a primary endpoint.'],
        ['Reproduction', 'prepare_sudmex_clinical.py', 'Rebuild source CSVs with setup/prepare_sudmex_clinical.py; use pipeline/scripts/utils/clinical_data.py for validated CSV or XLSX loading.'],
    ]
    payload = {'metrics': METRICS, 'sessions': records(sessions), 'effects': records(effects),
               'session_columns': session_cols, 'effect_columns': list(effects), 'dictionary': dictionary,
               'issues': issues, 'source': manifest['record']}
    (OUT / 'workbook_data.json').write_text(json.dumps(payload, indent=2, allow_nan=False), encoding='utf-8')
    summary = {'created_utc': datetime.now(timezone.utc).isoformat(), 'source': manifest['record'],
               'clinical_participants': int(sessions.subject_id.nunique()), 'clinical_visits': len(sessions),
               'analysis_participants': len(effects), 'mri_eligible_pairs': int(effects.mri_pair_usable.sum()),
               'mri_eligible_by_arm': effects[effects.mri_pair_usable == 1].arm.value_counts().to_dict(),
               'flagged_source_assessments': len(issues), 'issues': issues,
               'available_clean_mri_pairs': {m: int(((effects.mri_pair_usable == 1) & (effects[m + '_qc'] == 'ok') & effects[m + '_delta'].notna()).sum()) for m in METRICS},
               'source_manifest_sha256': hashlib.sha256((SOURCE / 'download_manifest.json').read_bytes()).hexdigest()}
    (OUT / 'clinical_validation.json').write_text(json.dumps(summary, indent=2), encoding='utf-8')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
