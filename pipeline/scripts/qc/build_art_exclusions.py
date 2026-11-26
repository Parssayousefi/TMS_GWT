#!/usr/bin/env python3
"""
Build ART-based exclusion CSV from CONN ART outlier summary.

Reads ``results/timeseries_intermediate/{dataset}/art_outlier_summary.csv``
(written by the CONN ART denoise stage, e.g.
``step01b_conn_denoise_cambridge_art.m`` lines 438-451) and emits two files:

  results/qc/{dataset}/exclusions_art.csv
      QCExclusions-compatible CSV with one EXCLUDE row per scan whose
      fraction of ART-flagged outlier volumes meets or exceeds the
      threshold. Schema is a strict superset of the motion-screen
      exclusions.csv (same first columns) so QCExclusions can ingest it
      unchanged.

  results/qc/{dataset}/art_exclusion_report.csv
      Audit table with one row per scan in the input summary, joining the
      ART fraction, the pre-existing motion-screen status, and the
      worst-case final status. This is the supplementary-methods
      deliverable.

Threshold rationale
-------------------
The default 20 % threshold follows:
  * Power et al. 2014 (NeuroImage 84:320) — recommends excluding scans
    with more than ~20 % of volumes flagged for motion.
  * Parkes et al. 2018 (NeuroImage 171:415) — validates 20 % as the
    stringent ART/scrubbing threshold in resting-state.
  * Luppi et al. 2024 (eLife 12:e88173) and the broader Luppi DOC line
    use ART (FD>0.5 mm, GS-Z>3) and exclude scans with excessive
    flagging at 20-25 %.

Usage
-----
    python pipeline/scripts/qc/build_art_exclusions.py \
        --dataset cambridge_propofol --threshold 0.20
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_DIR / 'scripts' / 'utils'))

from dataset_handler import load_dataset_handler  # noqa: E402

SCAN_ID_RE = re.compile(
    r'^sub-(?P<subject>[^_]+)'
    r'(?:_ses-(?P<session>[^_]+))?'
    r'_task-(?P<task>[^_]+)'
    r'_run-(?P<run>\d+)$'
)


def parse_scan_id(scan_id: str) -> dict:
    """Split a scan_id like ``sub-04HD_task-restlight_run-1`` into fields."""
    m = SCAN_ID_RE.match(scan_id)
    if not m:
        raise ValueError(f"Could not parse scan_id '{scan_id}'")
    g = m.groupdict()
    return {
        'subject': g['subject'],
        'session': g.get('session'),
        'task': g['task'],
        'run': int(g['run']),
    }


def load_motion_exclusions(dataset: str, project_dir: Path) -> pd.DataFrame:
    """Load the existing motion-screen exclusions.csv (returns empty DF if absent)."""
    path = project_dir / 'results' / 'qc' / dataset / 'exclusions.csv'
    if not path.exists():
        return pd.DataFrame(columns=['subject', 'session', 'task', 'status', 'reason'])
    return pd.read_csv(path, dtype={'subject': str})


def lookup_motion_status(motion_df: pd.DataFrame, subject: str,
                         session: Optional[str], task: str) -> str:
    """Return motion-screen status for a (subject, session, task) tuple, or 'UNKNOWN'."""
    if motion_df.empty:
        return 'UNKNOWN'
    subj_norm = str(subject).lstrip('0') or '0'
    rows = motion_df.copy()
    rows['_subj_norm'] = rows['subject'].astype(str).str.lstrip('0').replace('', '0')
    mask = rows['_subj_norm'] == subj_norm
    if 'task' in rows.columns:
        mask &= rows['task'].fillna('').astype(str) == (task or '')
    if 'session' in rows.columns:
        sess_match = rows['session'].astype(str).where(rows['session'].notna(), '')
        target = '' if session is None else str(session)
        mask &= sess_match == target
    matched = rows.loc[mask, 'status']
    if matched.empty:
        return 'UNKNOWN'
    statuses = set(matched.str.upper().tolist())
    for s in ('EXCLUDE', 'TRUNCATE', 'INCLUDE'):
        if s in statuses:
            return s
    return 'UNKNOWN'


def worst_case(*statuses: str) -> str:
    """EXCLUDE > TRUNCATE > INCLUDE > UNKNOWN."""
    s = set(statuses)
    for tier in ('EXCLUDE', 'TRUNCATE', 'INCLUDE'):
        if tier in s:
            return tier
    return 'UNKNOWN'


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--dataset', required=True,
                   help='Dataset name registered in configs/dataset_configs.yaml')
    p.add_argument('--threshold', type=float, default=0.20,
                   help='Exclude scans with fraction_outliers >= this value '
                        '(default 0.20; Power 2014 / Parkes 2018)')
    p.add_argument('--art-summary', type=Path, default=None,
                   help='Override input art_outlier_summary.csv path')
    p.add_argument('--exclusions-out', type=Path, default=None,
                   help='Override output exclusions_art.csv path')
    p.add_argument('--report-out', type=Path, default=None,
                   help='Override output art_exclusion_report.csv path')
    args = p.parse_args()

    handler = load_dataset_handler(args.dataset)
    dataset_label = handler.config.get('name', args.dataset)
    project_dir = handler.project_dir

    art_path = args.art_summary or (
        project_dir / 'results' / 'timeseries_intermediate' / args.dataset
        / 'art_outlier_summary.csv'
    )
    if not art_path.exists():
        print(f"ERROR: ART summary not found at {art_path}")
        return 2

    out_dir = project_dir / 'results' / 'qc' / args.dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    exclusions_out = args.exclusions_out or (out_dir / 'exclusions_art.csv')
    report_out = args.report_out or (out_dir / 'art_exclusion_report.csv')

    art_df = pd.read_csv(art_path)
    motion_df = load_motion_exclusions(args.dataset, project_dir)

    # Build per-scan rows
    art_rows = []
    report_rows = []

    for _, row in art_df.iterrows():
        scan_id = row['scan_id']
        try:
            parsed = parse_scan_id(scan_id)
        except ValueError as e:
            print(f"  SKIP unparseable scan_id: {scan_id} ({e})")
            continue

        n_outliers = int(row['n_outliers'])
        fraction = float(row['fraction_outliers'])
        art_fail = fraction >= args.threshold
        art_status = 'EXCLUDE' if art_fail else 'INCLUDE'

        motion_status = lookup_motion_status(
            motion_df, parsed['subject'], parsed['session'], parsed['task']
        )

        final = worst_case(motion_status, art_status)

        # Only ART-EXCLUDE rows go in exclusions_art.csv (so worst-case
        # merge in QCExclusions only adds new EXCLUDEs; INCLUDE rows are
        # left to the motion-screen file).
        if art_fail:
            reason = (f"ART outliers {fraction:.3f} >= {args.threshold:.2f} "
                      f"(Power 2014; Parkes 2018)")
            art_rows.append({
                'dataset': dataset_label,
                'subject': parsed['subject'],
                'session': parsed['session'] if parsed['session'] else '',
                'task': parsed['task'],
                'status': 'EXCLUDE',
                'reason': reason,
                'truncate_start': '',
                'truncate_end': '',
                'art_fraction_outliers': fraction,
                'art_n_outliers': n_outliers,
            })

        report_rows.append({
            'scan_id': scan_id,
            'dataset': dataset_label,
            'subject': parsed['subject'],
            'session': parsed['session'] if parsed['session'] else '',
            'task': parsed['task'],
            'run': parsed['run'],
            'n_outliers': n_outliers,
            'fraction_outliers': fraction,
            'art_threshold': args.threshold,
            'art_status': art_status,
            'motion_status_pre_art': motion_status,
            'final_status': final,
            'included_in_analysis': final == 'INCLUDE'
                                    or (final == 'TRUNCATE'),
        })

    # Sort for stable output
    sort_keys = ['subject', 'task']
    art_rows_sorted = sorted(art_rows, key=lambda r: (r['subject'], r['task']))
    report_rows_sorted = sorted(report_rows, key=lambda r: (r['subject'], r['task']))

    # Write exclusions_art.csv with provenance header (as comments)
    header_lines = [
        f"# ART-based scan exclusions for dataset: {args.dataset}",
        f"# Threshold: fraction_outliers >= {args.threshold:.2f}",
        f"# Source: {art_path.relative_to(project_dir)}",
        f"# Citations: Power 2014 (NeuroImage 84:320); Parkes 2018 "
        f"(NeuroImage 171:415); Luppi 2024 (eLife 12:e88173).",
        f"# Schema is a superset of exclusions.csv so QCExclusions ingests it.",
    ]
    art_columns = [
        'dataset', 'subject', 'session', 'task', 'status', 'reason',
        'truncate_start', 'truncate_end',
        'art_fraction_outliers', 'art_n_outliers',
    ]
    with open(exclusions_out, 'w', newline='') as f:
        for line in header_lines:
            f.write(line + '\n')
        pd.DataFrame(art_rows_sorted, columns=art_columns).to_csv(f, index=False)

    # Write report
    report_columns = [
        'scan_id', 'dataset', 'subject', 'session', 'task', 'run',
        'n_outliers', 'fraction_outliers', 'art_threshold',
        'art_status', 'motion_status_pre_art', 'final_status',
        'included_in_analysis',
    ]
    pd.DataFrame(report_rows_sorted, columns=report_columns).to_csv(
        report_out, index=False
    )

    n_total = len(report_rows_sorted)
    n_excluded = sum(1 for r in report_rows_sorted if r['art_status'] == 'EXCLUDE')
    n_motion_already = sum(1 for r in report_rows_sorted
                            if r['motion_status_pre_art'] == 'EXCLUDE')
    n_included_final = sum(1 for r in report_rows_sorted
                            if r['final_status'] == 'INCLUDE')

    print("=" * 60)
    print(f"ART exclusion build — dataset: {args.dataset}")
    print(f"  Threshold: fraction_outliers >= {args.threshold:.2f}")
    print(f"  Scans in ART summary:        {n_total}")
    print(f"  Excluded by ART threshold:   {n_excluded}")
    print(f"  Already excluded by motion:  {n_motion_already}")
    print(f"  Final INCLUDE after merge:   {n_included_final}")
    print()
    print(f"  exclusions_art.csv: {exclusions_out.relative_to(project_dir)}")
    print(f"  report:             {report_out.relative_to(project_dir)}")
    print("=" * 60)

    if n_excluded == 0:
        print("\nNote: no scans crossed the threshold. exclusions_art.csv is empty "
              "(but kept for downstream config consistency).")

    print("\nART-failed scans:")
    for r in report_rows_sorted:
        if r['art_status'] == 'EXCLUDE':
            print(f"  sub-{r['subject']:<8} task-{r['task']:<14} "
                  f"frac={r['fraction_outliers']:.3f} (motion: {r['motion_status_pre_art']})")

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
