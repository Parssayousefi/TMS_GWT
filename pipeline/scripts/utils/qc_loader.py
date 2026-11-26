#!/usr/bin/env python3
"""
QC Loader - Utility for loading and parsing QC exclusion files
==============================================================
Provides functions to check if subjects/sessions should be excluded
based on motion QC criteria.

Usage:
    from qc_loader import QCExclusions

    qc = QCExclusions('berlin_ketamine')

    # Check if a session is excluded
    if qc.is_excluded(subject='03', session='ketamine'):
        print("Session excluded")

    # Get truncation range for salvageable session
    truncation = qc.get_truncation(subject='03', session='placebo')
    if truncation:
        print(f"Truncate to volumes {truncation['start']}-{truncation['end']}")

    # Get all included sessions
    included = qc.get_included_sessions()
"""

import os
import yaml
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass


@dataclass
class TruncationInfo:
    """Information about truncation for salvageable sessions."""
    start_volume: int
    end_volume: int
    n_volumes: int
    duration_minutes: float
    reason: str


class QCExclusions:
    """
    Manages QC exclusion data for a dataset.

    Loads exclusion decisions from the standardized exclusions.csv file
    and provides methods to check exclusion status and get truncation info.
    """

    def __init__(self, dataset_name: str, project_dir: Optional[Path] = None):
        """
        Initialize QC exclusions for a dataset.

        Parameters:
        -----------
        dataset_name : str
            Name of the dataset (e.g., 'berlin_ketamine')
        project_dir : Path, optional
            Project directory. If None, infers from script location.
        """
        self.dataset_name = dataset_name

        if project_dir is None:
            # Assume we're in scripts/utils/, go up two levels
            project_dir = Path(__file__).parent.parent.parent
        self.project_dir = Path(project_dir)

        # Load dataset config to get QC settings
        self._load_config()

        # Load exclusions file
        self._load_exclusions()

    def _load_config(self):
        """Load dataset configuration."""
        # Load registry
        registry_path = self.project_dir / 'configs' / 'dataset_configs.yaml'
        with open(registry_path, 'r') as f:
            registry = yaml.safe_load(f)

        # Get dataset config path
        if self.dataset_name not in registry['datasets']:
            available = ', '.join(registry['datasets'].keys())
            raise ValueError(f"Dataset '{self.dataset_name}' not found. Available: {available}")

        config_file = registry['datasets'][self.dataset_name]['config_file']
        config_path = self.project_dir / config_file

        with open(config_path, 'r') as f:
            self.config = yaml.safe_load(f)

        # Get QC settings
        self.qc_config = self.config.get('qc', {})
        self.motion_thresholds = self.qc_config.get('motion_thresholds', {
            'max_translation_mm': 3.0,
            'max_rotation_deg': 3.0
        })
        self.auto_truncate = self.qc_config.get('auto_truncate', False)

        # Get exclusions file path
        exclusions_file = self.qc_config.get('exclusions_file',
            f'results/qc/{self.dataset_name}/exclusions.csv')
        self.exclusions_path = self.project_dir / exclusions_file

        # Optional second exclusions file (e.g. ART-based scrubbing exclusions
        # for cambridge_propofol). When present, rows merge into the same
        # store with worst-case collapse (EXCLUDE > TRUNCATE > INCLUDE), so
        # an ART EXCLUDE can override a motion-screen INCLUDE for the same
        # (subject, session, task) without losing the audit trail (the
        # source file remains on disk).
        art_exclusions_file = self.qc_config.get('art_exclusions_file')
        self.art_exclusions_path = (
            self.project_dir / art_exclusions_file
            if art_exclusions_file else None
        )

    def _load_exclusions(self):
        """Load exclusions from the motion-screen CSV and (if configured) the
        ART-based CSV, merging both with worst-case collapse."""
        self.exclusions = {}
        self.truncations = {}

        if not self.exclusions_path.exists():
            print(f"QC exclusions file not found: {self.exclusions_path}")
            print("Run the QC script first to generate exclusions.")
        else:
            self._ingest_exclusion_csv(self.exclusions_path, source_label='motion')

        # Optional second source (ART). Worst-case collapse means an ART
        # EXCLUDE wins over a motion INCLUDE for the same key.
        if self.art_exclusions_path is not None:
            if self.art_exclusions_path.exists():
                self._ingest_exclusion_csv(self.art_exclusions_path,
                                            source_label='art')
            else:
                print(f"ART exclusions file configured but not found: "
                      f"{self.art_exclusions_path}")

        n_excluded = sum(1 for e in self.exclusions.values() if e['status'] == 'EXCLUDE')
        n_truncate = sum(1 for e in self.exclusions.values() if e['status'] == 'TRUNCATE')
        n_include = sum(1 for e in self.exclusions.values() if e['status'] == 'INCLUDE')

        print(f"Loaded QC exclusions: {n_include} INCLUDE, {n_excluded} EXCLUDE, {n_truncate} TRUNCATE")

    def _ingest_exclusion_csv(self, path: Path, source_label: str = 'motion'):
        """Read one exclusion CSV and merge rows into ``self.exclusions``.

        Subsequent calls merge with worst-case collapse so a stricter status
        from a later file (e.g. an ART EXCLUDE) overrides a more permissive
        status from an earlier file (e.g. motion INCLUDE) for the same key.

        ``comment='#'`` lets the ART file carry an inline provenance header
        without a separate sidecar; the motion file has no comments so this
        is a no-op for it.
        """
        df = pd.read_csv(path, dtype={'subject': str}, comment='#')

        for _, row in df.iterrows():
            subject = str(row['subject'])
            session = str(row['session']) if pd.notna(row.get('session')) else None
            task = str(row['task']) if pd.notna(row.get('task')) else None
            status = row['status']
            reason = row.get('reason', '')

            # Key includes task so datasets that encode conditions in the task
            # name (has_sessions=false) get one entry per row instead of
            # later-task-alphabetical overwriting earlier ones. See
            # cambridge_propofol for the motivating case.
            key = self._normalize_key(subject, session, task)

            existing = self.exclusions.get(key)
            if existing is None:
                self.exclusions[key] = {
                    'status': status,
                    'reason': reason,
                    'source': source_label,
                }
            else:
                merged_status = self._worst_status([existing['status'], status])
                # Only overwrite reason/source when the new row's status is
                # the one that wins; otherwise keep the prior decision's
                # provenance for audit clarity.
                if merged_status != existing['status']:
                    self.exclusions[key] = {
                        'status': merged_status,
                        'reason': reason,
                        'source': source_label,
                    }
                else:
                    self.exclusions[key]['status'] = merged_status

            # Store truncation info if this row is a TRUNCATE (regardless of
            # which file it came from).
            if status == 'TRUNCATE':
                start = row.get('truncate_start')
                end = row.get('truncate_end')
                if pd.notna(start) and pd.notna(end):
                    self.truncations[key] = TruncationInfo(
                        start_volume=int(start),
                        end_volume=int(end),
                        n_volumes=int(end) - int(start),
                        duration_minutes=(int(end) - int(start)) * 2.0 / 60,  # Assuming TR=2
                        reason=reason
                    )

    def _normalize_key(self, subject: str, session: Optional[str] = None,
                       task: Optional[str] = None) -> tuple:
        """Normalize (subject, session, task) key for consistent lookup.

        Subject IDs have leading zeros stripped so "001"/"01"/"1" match.
        """
        subj = str(subject).lstrip('0') or '0'
        sess = str(session) if session else None
        tsk = str(task) if task else None
        return (subj, sess, tsk)

    def _matching_keys(self, subject: str, session: Optional[str],
                       task: Optional[str]) -> list:
        """Return all stored keys matching the (subject, session, task) query.

        If task is None, returns every stored entry that matches (subject, session)
        regardless of stored task. Used for backward-compat callers that don't
        pass task. The public methods collapse multiple matches via a
        worst-case policy so old callers become conservative, not wrong.
        """
        subj = str(subject).lstrip('0') or '0'
        sess = str(session) if session else None
        if task is not None:
            key = self._normalize_key(subject, session, task)
            return [key] if key in self.exclusions else []
        return [k for k in self.exclusions
                if k[0] == subj and k[1] == sess]

    @staticmethod
    def _worst_status(statuses: list) -> str:
        """Conservative collapse: EXCLUDE > TRUNCATE > INCLUDE."""
        if 'EXCLUDE' in statuses:
            return 'EXCLUDE'
        if 'TRUNCATE' in statuses:
            return 'TRUNCATE'
        if 'INCLUDE' in statuses:
            return 'INCLUDE'
        return 'UNKNOWN'

    def is_excluded(self, subject: str, session: Optional[str] = None,
                    task: Optional[str] = None) -> bool:
        """
        Check if a (subject, session, task) should be excluded.

        If task is None, returns True if ANY matching row would be excluded
        (conservative fallback for old session-keyed callers).
        """
        status = self.get_status(subject, session, task)

        if status == 'EXCLUDE':
            return True
        elif status == 'TRUNCATE':
            return not self.auto_truncate
        else:  # INCLUDE or UNKNOWN
            return False

    def get_status(self, subject: str, session: Optional[str] = None,
                   task: Optional[str] = None) -> str:
        """
        Get QC status for a (subject, session, task).

        If task is None, returns the worst-case status across all matching
        rows (EXCLUDE dominates TRUNCATE dominates INCLUDE).

        Returns: 'INCLUDE', 'EXCLUDE', 'TRUNCATE', or 'UNKNOWN'
        """
        matches = self._matching_keys(subject, session, task)
        if not matches:
            return 'UNKNOWN'
        statuses = [self.exclusions[k]['status'] for k in matches]
        return self._worst_status(statuses)

    def get_truncation(self, subject: str, session: Optional[str] = None,
                       task: Optional[str] = None) -> Optional[TruncationInfo]:
        """
        Get truncation info for a salvageable session.

        If task is None and multiple matching rows exist, returns the first
        TRUNCATE row encountered (undefined order; task-aware callers should
        always pass task to get deterministic behaviour).
        """
        if task is not None:
            key = self._normalize_key(subject, session, task)
            return self.truncations.get(key)
        for k in self._matching_keys(subject, session, None):
            if k in self.truncations:
                return self.truncations[k]
        return None

    def get_included_sessions(self) -> List[Tuple[str, Optional[str], Optional[str]]]:
        """
        Get list of all included (subject, session, task) rows.

        Rows with status INCLUDE or TRUNCATE (if auto_truncate=True) are included.
        """
        included = []

        for key, info in self.exclusions.items():
            status = info['status']

            if status == 'INCLUDE':
                included.append(key)
            elif status == 'TRUNCATE' and self.auto_truncate:
                included.append(key)

        return sorted(included)

    def get_excluded_sessions(self) -> List[Tuple[str, Optional[str], Optional[str], str]]:
        """
        Get list of all excluded rows as (subject, session, task, reason) tuples.
        """
        excluded = []

        for key, info in self.exclusions.items():
            if self.is_excluded(key[0], key[1], key[2]):
                excluded.append((key[0], key[1], key[2], info['reason']))

        return sorted(excluded)

    def get_exclusion_reason(self, subject: str, session: Optional[str] = None,
                             task: Optional[str] = None) -> str:
        """Get the reason for exclusion. With task=None, joins reasons from all matching rows."""
        matches = self._matching_keys(subject, session, task)
        if not matches:
            return ''
        reasons = [self.exclusions[k].get('reason', '') for k in matches
                   if self.exclusions[k].get('reason')]
        return '; '.join(reasons)

    def summary(self) -> Dict:
        """Get summary statistics of QC exclusions."""
        stats = {
            'total': len(self.exclusions),
            'include': 0,
            'exclude': 0,
            'truncate': 0,
            'effective_include': 0,
            'effective_exclude': 0
        }

        for key, info in self.exclusions.items():
            status = info['status']
            stats[status.lower()] = stats.get(status.lower(), 0) + 1

            if not self.is_excluded(key[0], key[1], key[2]):
                stats['effective_include'] += 1
            else:
                stats['effective_exclude'] += 1

        return stats


def load_qc_exclusions(dataset_name: str, project_dir: Optional[Path] = None) -> QCExclusions:
    """
    Convenience function to load QC exclusions.

    Parameters:
    -----------
    dataset_name : str
        Name of the dataset
    project_dir : Path, optional
        Project directory

    Returns:
    --------
    QCExclusions
        Loaded QC exclusions object
    """
    return QCExclusions(dataset_name, project_dir)


# CLI for testing
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Test QC exclusions loader')
    parser.add_argument('--dataset', '-d', required=True, help='Dataset name')
    parser.add_argument('--subject', '-s', help='Check specific subject')
    parser.add_argument('--session', help='Check specific session')

    args = parser.parse_args()

    qc = QCExclusions(args.dataset)

    print(f"\nQC Exclusions Summary for {args.dataset}:")
    print("-" * 40)

    summary = qc.summary()
    print(f"Total sessions: {summary['total']}")
    print(f"  INCLUDE: {summary['include']}")
    print(f"  EXCLUDE: {summary['exclude']}")
    print(f"  TRUNCATE: {summary['truncate']}")
    print(f"  Auto-truncate: {qc.auto_truncate}")
    print(f"  Effective include: {summary['effective_include']}")
    print(f"  Effective exclude: {summary['effective_exclude']}")

    if args.subject:
        subject = args.subject
        session = args.session

        print(f"\nChecking sub-{subject}" + (f" ses-{session}" if session else "") + ":")
        print(f"  Status: {qc.get_status(subject, session)}")
        print(f"  Excluded: {qc.is_excluded(subject, session)}")
        print(f"  Reason: {qc.get_exclusion_reason(subject, session)}")

        truncation = qc.get_truncation(subject, session)
        if truncation:
            print(f"  Truncation: volumes {truncation.start_volume}-{truncation.end_volume}")
