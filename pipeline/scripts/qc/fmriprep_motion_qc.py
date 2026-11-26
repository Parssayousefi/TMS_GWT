#!/usr/bin/env python3
"""
Dataset-Agnostic fMRIPrep Motion QC Script
==========================================
Generates a quality control report based on Luppi et al. exclusion criteria:
- Max translation > 3mm (configurable)
- Max rotation > 3 degrees (configurable)

Reads thresholds from dataset configuration and outputs to standardized location.
Generates:
- motion_qc.csv: Raw motion metrics for all sessions
- exclusions.csv: Final exclusion decisions (INCLUDE/EXCLUDE/TRUNCATE)
- motion_qc_report.md: Human-readable report

Usage:
    python fmriprep_motion_qc.py --dataset berlin_ketamine
    python fmriprep_motion_qc.py --dataset london_lsd --env local

Author: Generated for GWT Pipeline
"""

import os
import glob
import platform
import pandas as pd
import numpy as np
import yaml
from pathlib import Path
from datetime import datetime
import argparse


def detect_clean_segments(confounds_file, tr=2.0, fd_threshold=0.5, trans_threshold=3.0, rot_threshold=3.0, min_duration_minutes=8, skip_first=0, skip_last=0):
    """
    Detect clean segments in a session that could potentially be salvaged.

    A volume is considered "clean" if:
    - FD < fd_threshold (or is NaN for first volume)
    - Cumulative translation from start < trans_threshold
    - Cumulative rotation from start < rot_threshold

    Parameters
    ----------
    skip_first : int
        Number of initial volumes to exclude from search (default: 0)
    skip_last : int
        Number of final volumes to exclude from search (default: 0)

    Returns list of clean segments with start/end volumes and duration.
    Volume indices are in **original coordinates** (not offset by skip_first).
    """
    try:
        df = pd.read_csv(confounds_file, sep='\t')
    except Exception as e:
        return []

    n_volumes = len(df)
    min_volumes = int((min_duration_minutes * 60) / tr)

    # Constrain search range to [skip_first, n_volumes - skip_last]
    search_start = skip_first
    search_end = n_volumes - skip_last if skip_last > 0 else n_volumes

    if search_end - search_start < min_volumes:
        return []

    # Get motion parameters
    fd = df['framewise_displacement'].fillna(0).values if 'framewise_displacement' in df.columns else np.zeros(n_volumes)

    # Calculate instantaneous translation and rotation at each timepoint
    trans_cols = ['trans_x', 'trans_y', 'trans_z']
    rot_cols = ['rot_x', 'rot_y', 'rot_z']

    # For each volume, check if motion from first volume exceeds threshold
    trans_max = np.zeros(n_volumes)
    rot_max = np.zeros(n_volumes)

    for col in trans_cols:
        if col in df.columns:
            # Relative to first volume
            trans_max = np.maximum(trans_max, np.abs(df[col].values - df[col].values[0]))

    for col in rot_cols:
        if col in df.columns:
            rot_max = np.maximum(rot_max, np.abs(np.degrees(df[col].values - df[col].values[0])))

    # Find clean segments by sliding window approach
    # A segment is clean if max motion within segment is below threshold
    segments = []

    # Try different starting points within the constrained range
    for start in range(search_start, search_end - min_volumes):
        # Check segment from this start point
        best_end = start

        for end in range(start + min_volumes, search_end + 1):
            segment_length = end - start

            # Get motion within this segment (relative to segment start)
            seg_trans_max = 0
            seg_rot_max = 0

            for col in trans_cols:
                if col in df.columns:
                    seg_vals = df[col].values[start:end]
                    seg_trans_max = max(seg_trans_max, np.abs(seg_vals - seg_vals[0]).max())

            for col in rot_cols:
                if col in df.columns:
                    seg_vals = df[col].values[start:end]
                    seg_rot_max = max(seg_rot_max, np.abs(np.degrees(seg_vals - seg_vals[0])).max())

            # Check if segment is clean
            if seg_trans_max <= trans_threshold and seg_rot_max <= rot_threshold:
                best_end = end
            else:
                break

        if best_end > start + min_volumes:
            duration_minutes = ((best_end - start) * tr) / 60

            # Calculate segment stats
            seg_fd = fd[start:best_end]
            seg_mean_fd = np.nanmean(seg_fd)

            # Volume indices remain in original coordinates
            segments.append({
                'start_volume': start,
                'end_volume': best_end,
                'n_volumes': best_end - start,
                'duration_minutes': duration_minutes,
                'mean_fd': seg_mean_fd,
                'start_time': f"{int(start * tr // 60):02d}:{int(start * tr % 60):02d}",
                'end_time': f"{int(best_end * tr // 60):02d}:{int(best_end * tr % 60):02d}"
            })

    # Remove overlapping segments, keep longest
    if not segments:
        return []

    # Sort by duration descending
    segments.sort(key=lambda x: x['duration_minutes'], reverse=True)

    # Keep only non-overlapping segments
    final_segments = []
    used_volumes = set()

    for seg in segments:
        seg_volumes = set(range(seg['start_volume'], seg['end_volume']))
        if not seg_volumes & used_volumes:  # No overlap
            final_segments.append(seg)
            used_volumes.update(seg_volumes)

    return final_segments[:3]  # Return top 3 segments


def get_motion_metrics(confounds_file, skip_first=0, skip_last=0, relative=False):
    """
    Extract motion metrics from fMRIPrep confounds file.

    Parameters
    ----------
    confounds_file : str
        Path to fMRIPrep confounds TSV file
    skip_first : int
        Number of initial volumes to skip in assessment (default: 0)
    skip_last : int
        Number of final volumes to skip in assessment (default: 0)
    relative : bool
        If True, compute trans/rot relative to the first volume of the
        assessed range (for within-segment metrics). If False, use absolute
        values from fMRIPrep (default, for Luppi criteria checking).

    Returns dict with:
    - max_trans_x/y/z: Maximum absolute translation in each axis (mm)
    - max_trans: Maximum translation across all axes (mm)
    - max_rot_x/y/z: Maximum absolute rotation in each axis (degrees)
    - max_rot: Maximum rotation across all axes (degrees)
    - mean_fd: Mean framewise displacement (mm)
    - max_fd: Maximum framewise displacement (mm)
    - n_motion_outliers: Number of fMRIPrep-flagged motion outliers
    - n_volumes_total: Total number of volumes in file
    - n_volumes_assessed: Number of volumes actually assessed
    - skip_first / skip_last: Boundary trim applied
    """
    try:
        df = pd.read_csv(confounds_file, sep='\t')
    except Exception as e:
        return {'error': str(e)}

    n_volumes_total = len(df)

    # Trim to assessment range
    end_idx = n_volumes_total - skip_last if skip_last > 0 else n_volumes_total
    if skip_first > 0 or skip_last > 0:
        df = df.iloc[skip_first:end_idx]

    metrics = {}
    metrics['n_volumes_total'] = n_volumes_total
    metrics['n_volumes_assessed'] = len(df)
    metrics['skip_first'] = skip_first
    metrics['skip_last'] = skip_last

    # Translation (mm)
    for axis in ['x', 'y', 'z']:
        col = f'trans_{axis}'
        if col in df.columns:
            if relative and len(df) > 0:
                metrics[f'max_trans_{axis}'] = (df[col] - df[col].iloc[0]).abs().max()
            else:
                metrics[f'max_trans_{axis}'] = df[col].abs().max()
        else:
            metrics[f'max_trans_{axis}'] = np.nan

    metrics['max_trans'] = max(
        metrics.get('max_trans_x', 0),
        metrics.get('max_trans_y', 0),
        metrics.get('max_trans_z', 0)
    )

    # Rotation (radians -> degrees)
    for axis in ['x', 'y', 'z']:
        col = f'rot_{axis}'
        if col in df.columns:
            if relative and len(df) > 0:
                metrics[f'max_rot_{axis}'] = np.degrees((df[col] - df[col].iloc[0]).abs().max())
            else:
                metrics[f'max_rot_{axis}'] = np.degrees(df[col].abs().max())
        else:
            metrics[f'max_rot_{axis}'] = np.nan

    metrics['max_rot'] = max(
        metrics.get('max_rot_x', 0),
        metrics.get('max_rot_y', 0),
        metrics.get('max_rot_z', 0)
    )

    # Framewise displacement
    if 'framewise_displacement' in df.columns:
        fd = df['framewise_displacement'].dropna()
        metrics['mean_fd'] = fd.mean()
        metrics['max_fd'] = fd.max()
        metrics['median_fd'] = fd.median()
    else:
        metrics['mean_fd'] = np.nan
        metrics['max_fd'] = np.nan
        metrics['median_fd'] = np.nan

    # Motion outliers (fMRIPrep flagged) — count from full file, not trimmed
    outlier_cols = [c for c in pd.read_csv(confounds_file, sep='\t', nrows=0).columns if c.startswith('motion_outlier')]
    metrics['n_motion_outliers'] = len(outlier_cols)

    # Legacy key for backward compatibility
    metrics['n_volumes'] = n_volumes_total

    # Percentage of outliers (relative to total)
    metrics['pct_outliers'] = (metrics['n_motion_outliers'] / n_volumes_total) * 100 if n_volumes_total > 0 else 0

    return metrics


def check_luppi_criteria(metrics, trans_threshold=3.0, rot_threshold=3.0):
    """
    Check if session should be excluded based on Luppi et al. criteria.

    Criteria:
    - Max translation > 3mm
    - Max rotation > 3 degrees
    - Sufficient volumes after boundary trimming
    """
    if 'error' in metrics:
        return 'ERROR', 'Could not read confounds file'

    # Check for insufficient data after boundary trimming
    n_assessed = metrics.get('n_volumes_assessed', metrics.get('n_volumes', 0))
    if n_assessed == 0 or np.isnan(metrics.get('max_trans', np.nan)) or np.isnan(metrics.get('max_rot', np.nan)):
        return 'EXCLUDE', f'insufficient_volumes ({metrics.get("n_volumes_total", 0)} total, {n_assessed} after trim)'

    reasons = []

    if metrics['max_trans'] > trans_threshold:
        reasons.append(f"max_trans={metrics['max_trans']:.2f}mm > {trans_threshold}mm")

    if metrics['max_rot'] > rot_threshold:
        reasons.append(f"max_rot={metrics['max_rot']:.2f}° > {rot_threshold}°")

    if reasons:
        return 'EXCLUDE', '; '.join(reasons)
    else:
        return 'INCLUDE', ''


def find_all_sessions(derivatives_dir, has_sessions=True, task_filter=None,
                      include_subjects=None, exclude_subjects=None,
                      include_sessions=None, run_selection='first'):
    """
    Find all confounds files in the derivatives directory.
    Returns list of dicts with subject, session, and file path.

    When multiple runs exist for the same subject/session/task, only one is
    kept based on run_selection strategy, ensuring the QC report matches
    what actually enters the analysis pipeline.

    Parameters:
    -----------
    derivatives_dir : str
        Path to fMRIPrep derivatives directory
    has_sessions : bool
        Whether dataset has session structure (default: True)
    task_filter : list of str, optional
        If provided, only include sessions whose task matches one of these names.
    include_subjects : list of str, optional
        If provided, only include these subjects.
    exclude_subjects : list of str, optional
        If provided, exclude these subjects.
    include_sessions : list of str, optional
        If provided, only include sessions matching these names.
    run_selection : str
        Run selection strategy: 'first', 'last', or 'all' (default: 'first').
        'first' keeps the lowest run number per subject/session/task.
        'last' keeps the highest run number.
        'all' keeps all runs (no deduplication).
    """
    sessions = []

    if has_sessions:
        # Pattern: sub-*/ses-*/func/*_desc-confounds_timeseries.tsv
        pattern = os.path.join(derivatives_dir, 'sub-*', 'ses-*', 'func', '*_desc-confounds_timeseries.tsv')
    else:
        # Pattern: sub-*/func/*_desc-confounds_timeseries.tsv
        pattern = os.path.join(derivatives_dir, 'sub-*', 'func', '*_desc-confounds_timeseries.tsv')

    for confounds_file in glob.glob(pattern):
        path = Path(confounds_file)

        # Parse subject, session from directory parts (exclude filename)
        dir_parts = path.parent.parts
        subject = None
        session = None
        task = None
        run = None

        for part in dir_parts:
            if part.startswith('sub-'):
                subject = part.replace('sub-', '')
            elif part.startswith('ses-'):
                session = part.replace('ses-', '')

        # Extract task and run from filename
        filename = path.name
        for part in filename.split('_'):
            if part.startswith('task-'):
                task = part.replace('task-', '')
            elif part.startswith('run-'):
                try:
                    run = int(part.replace('run-', ''))
                except ValueError:
                    run = part.replace('run-', '')

        # Apply task filter if specified
        if task_filter and task not in task_filter:
            continue

        # Apply subject filters
        if include_subjects and subject not in include_subjects:
            continue
        if exclude_subjects and subject in exclude_subjects:
            continue

        # Apply session filter
        if include_sessions and session not in include_sessions:
            continue

        if subject:
            sessions.append({
                'subject': subject,
                'session': session,  # Will be None if no sessions
                'task': task,
                'run': run,  # None if no run in filename
                'confounds_file': confounds_file,
                'filename': path.name
            })

    # Deduplicate by (subject, session, task) keeping one run per strategy
    if run_selection != 'all':
        from collections import defaultdict
        grouped = defaultdict(list)
        for sess in sessions:
            key = (sess['subject'], sess['session'], sess['task'])
            grouped[key].append(sess)

        deduped = []
        for key, entries in grouped.items():
            if len(entries) == 1:
                deduped.append(entries[0])
            else:
                # Sort by run number (None treated as 0)
                entries.sort(key=lambda x: x['run'] if x['run'] is not None else 0)
                if run_selection == 'first':
                    deduped.append(entries[0])
                elif run_selection == 'last':
                    deduped.append(entries[-1])
                else:
                    deduped.extend(entries)
        sessions = deduped

    # Sort by subject, then session (None sorts first)
    return sorted(sessions, key=lambda x: (x['subject'], x['session'] or '', x.get('task', '')))


def generate_qc_report(derivatives_dir, output_dir, dataset_name='unknown',
                       trans_threshold=3.0, rot_threshold=3.0,
                       min_salvage_minutes=8, tr=2.0, has_sessions=True,
                       skip_first=5, skip_last=5, task_filter=None,
                       include_subjects=None, exclude_subjects=None,
                       include_sessions=None, run_selection='first'):
    """
    Generate comprehensive QC report and exclusion files.

    For excluded sessions, also detects potentially salvageable clean segments
    that meet the Luppi criteria for at least min_salvage_minutes.

    Parameters:
    -----------
    derivatives_dir : str
        Path to fMRIPrep derivatives directory
    output_dir : str
        Output directory for QC files
    dataset_name : str
        Name of the dataset (for report header)
    trans_threshold : float
        Maximum translation threshold in mm
    rot_threshold : float
        Maximum rotation threshold in degrees
    min_salvage_minutes : float
        Minimum duration for salvageable segments
    tr : float
        Repetition time in seconds
    has_sessions : bool
        Whether dataset has session structure
    skip_first : int
        Number of initial volumes to skip in QC assessment (default: 5)
    skip_last : int
        Number of final volumes to skip in QC assessment (default: 5)

    Outputs:
    --------
    - motion_qc.csv: Raw motion metrics
    - exclusions.csv: Exclusion decisions (INCLUDE/EXCLUDE/TRUNCATE)
    - motion_qc_report.md: Human-readable report
    """
    sessions = find_all_sessions(
        derivatives_dir, has_sessions=has_sessions, task_filter=task_filter,
        include_subjects=include_subjects, exclude_subjects=exclude_subjects,
        include_sessions=include_sessions, run_selection=run_selection
    )

    if not sessions:
        print(f"No sessions found in {derivatives_dir}")
        return

    filters_applied = []
    if task_filter:
        filters_applied.append(f"tasks={task_filter}")
    if include_subjects:
        filters_applied.append(f"subjects={include_subjects}")
    if exclude_subjects:
        filters_applied.append(f"excluded subjects={exclude_subjects}")
    if include_sessions:
        filters_applied.append(f"sessions={include_sessions}")

    if filters_applied:
        print(f"Found {len(sessions)} sessions (filters: {', '.join(filters_applied)})")
    else:
        print(f"Found {len(sessions)} sessions")

    # Collect all metrics
    results = []
    salvage_analysis = []  # Track salvage potential for excluded sessions

    for sess in sessions:
        session_str = f"_ses-{sess['session']}" if sess['session'] else ""
        task_str = f"_task-{sess['task']}" if sess.get('task') else ""
        print(f"  Processing sub-{sess['subject']}{session_str}{task_str}...")

        metrics = get_motion_metrics(sess['confounds_file'], skip_first=skip_first, skip_last=skip_last)
        status, reason = check_luppi_criteria(metrics, trans_threshold, rot_threshold)

        result_row = {
            'Subject': f"sub-{sess['subject']}",
            'Session': sess['session'] if sess['session'] else '',
            'Task': sess.get('task', ''),
            'Max Trans (mm)': metrics.get('max_trans', np.nan),
            'Max Rot (deg)': metrics.get('max_rot', np.nan),
            'Mean FD (mm)': metrics.get('mean_fd', np.nan),
            'Max FD (mm)': metrics.get('max_fd', np.nan),
            'Motion Outliers': metrics.get('n_motion_outliers', 0),
            'N Volumes': metrics.get('n_volumes', 0),
            'N Volumes Assessed': metrics.get('n_volumes_assessed', metrics.get('n_volumes', 0)),
            'Status': status,
            'Exclusion Reason': reason,
            'Salvageable': '',
            'Best Segment Duration (min)': np.nan,
            'Best Segment Volumes': '',
            'Recommended Truncation': '',
            'Post-Trunc Max Trans (mm)': np.nan,
            'Post-Trunc Max Rot (deg)': np.nan,
            'Post-Trunc Mean FD (mm)': np.nan,
            'Post-Trunc N Volumes': np.nan,
            # "Usable" metrics: for INCLUDE sessions these equal the assessed metrics;
            # for TRUNCATE sessions they will be overwritten with post-truncation values
            'Usable Max Trans (mm)': metrics.get('max_trans', np.nan),
            'Usable Max Rot (deg)': metrics.get('max_rot', np.nan),
            'Usable Mean FD (mm)': metrics.get('mean_fd', np.nan),
            'Usable N Volumes': metrics.get('n_volumes_assessed', metrics.get('n_volumes', np.nan)),
        }

        # For excluded sessions, check for salvageable clean segments
        if status == 'EXCLUDE':
            print(f"    Checking for salvageable segments...")
            segments = detect_clean_segments(
                sess['confounds_file'],
                tr=tr,
                trans_threshold=trans_threshold,
                rot_threshold=rot_threshold,
                min_duration_minutes=min_salvage_minutes,
                skip_first=skip_first,
                skip_last=skip_last
            )

            if segments:
                best_seg = segments[0]  # Longest segment
                result_row['Salvageable'] = 'YES'
                result_row['Best Segment Duration (min)'] = best_seg['duration_minutes']
                result_row['Best Segment Volumes'] = f"{best_seg['start_volume']}-{best_seg['end_volume']}"
                result_row['Recommended Truncation'] = f"Keep volumes {best_seg['start_volume']}-{best_seg['end_volume']} ({best_seg['n_volumes']} vols, {best_seg['duration_minutes']:.1f} min)"

                # Recompute metrics on just the clean segment (relative to segment start)
                n_vols_total = metrics.get('n_volumes_total', metrics.get('n_volumes', 0))
                seg_skip_last = n_vols_total - best_seg['end_volume']
                post_trunc_metrics = get_motion_metrics(
                    sess['confounds_file'],
                    skip_first=best_seg['start_volume'],
                    skip_last=seg_skip_last if seg_skip_last > 0 else 0,
                    relative=True
                )
                result_row['Post-Trunc Max Trans (mm)'] = post_trunc_metrics.get('max_trans', np.nan)
                result_row['Post-Trunc Max Rot (deg)'] = post_trunc_metrics.get('max_rot', np.nan)
                result_row['Post-Trunc Mean FD (mm)'] = post_trunc_metrics.get('mean_fd', np.nan)
                result_row['Post-Trunc N Volumes'] = post_trunc_metrics.get('n_volumes_assessed', np.nan)
                # Usable metrics = post-truncation for salvaged sessions
                result_row['Usable Max Trans (mm)'] = post_trunc_metrics.get('max_trans', np.nan)
                result_row['Usable Max Rot (deg)'] = post_trunc_metrics.get('max_rot', np.nan)
                result_row['Usable Mean FD (mm)'] = post_trunc_metrics.get('mean_fd', np.nan)
                result_row['Usable N Volumes'] = post_trunc_metrics.get('n_volumes_assessed', np.nan)

                salvage_analysis.append({
                    'Subject': f"sub-{sess['subject']}",
                    'Session': sess['session'],
                    'Original Status': 'EXCLUDE',
                    'Original Reason': reason,
                    'N Segments Found': len(segments),
                    'Best Segment Start': best_seg['start_volume'],
                    'Best Segment End': best_seg['end_volume'],
                    'Best Segment Duration (min)': best_seg['duration_minutes'],
                    'Best Segment Mean FD': best_seg['mean_fd'],
                    'Start Time': best_seg['start_time'],
                    'End Time': best_seg['end_time'],
                    'All Segments': segments
                })
                print(f"      Found {len(segments)} clean segment(s), best: {best_seg['duration_minutes']:.1f} min")
            else:
                result_row['Salvageable'] = 'NO'
                print(f"      No salvageable segments found")

        results.append(result_row)

    df = pd.DataFrame(results)

    # Generate markdown report
    md_content = []

    md_content.append("# fMRIPrep Motion Quality Control Report")
    md_content.append("")
    md_content.append(f"**Dataset:** {dataset_name}")
    md_content.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    md_content.append(f"**Assessment Range:** Volumes [{skip_first} .. N-{skip_last}] (boundary trimmed)")
    md_content.append(f"**Exclusion Criteria (Luppi et al.):**")
    md_content.append(f"- Maximum translation > {trans_threshold} mm")
    md_content.append(f"- Maximum rotation > {rot_threshold} degrees")
    md_content.append(f"- First {skip_first} volumes skipped (steady-state magnetization)")
    md_content.append(f"- Last {skip_last} volumes skipped (end-of-scan artifact)")
    md_content.append("")

    # Summary statistics
    n_total = len(df)
    n_exclude = len(df[df['Status'] == 'EXCLUDE'])
    n_include = len(df[df['Status'] == 'INCLUDE'])

    md_content.append("## Summary")
    md_content.append("")
    md_content.append(f"| Metric | Value |")
    md_content.append(f"|--------|-------|")
    md_content.append(f"| Total sessions | {n_total} |")
    md_content.append(f"| Sessions to INCLUDE | {n_include} |")
    md_content.append(f"| Sessions to EXCLUDE | {n_exclude} |")
    md_content.append(f"| Exclusion rate | {100*n_exclude/n_total:.1f}% |")
    md_content.append("")

    # Group statistics — usable sessions (INCLUDE + TRUNCATE with post-trunc metrics)
    included = df[df['Status'] == 'INCLUDE']
    salvageable = df[df['Salvageable'] == 'YES']
    n_usable = len(included) + len(salvageable)

    md_content.append(f"### Motion Statistics (All Usable Sessions: {len(included)} INCLUDE + {len(salvageable)} TRUNCATE = {n_usable})")
    md_content.append("")
    md_content.append("*For TRUNCATE sessions, metrics are from the clean segment that enters analysis.*")
    md_content.append("")
    if n_usable > 0:
        md_content.append(f"| Metric | Mean | SD | Min | Max |")
        md_content.append(f"|--------|------|-----|-----|-----|")
        for col in ['Usable Max Trans (mm)', 'Usable Max Rot (deg)', 'Usable Mean FD (mm)']:
            # Combine INCLUDE and TRUNCATE usable metrics
            usable_vals = pd.concat([included[col], salvageable[col]]).dropna()
            if len(usable_vals) > 0:
                label = col.replace('Usable ', '')
                md_content.append(f"| {label} | {usable_vals.mean():.3f} | {usable_vals.std():.3f} | {usable_vals.min():.3f} | {usable_vals.max():.3f} |")
        # Usable volume counts (integer-valued; emphasizes the minimum)
        usable_vol_vals = pd.concat([included['Usable N Volumes'], salvageable['Usable N Volumes']]).dropna()
        if len(usable_vol_vals) > 0:
            md_content.append(
                f"| N Volumes (usable) | {usable_vol_vals.mean():.1f} | {usable_vol_vals.std():.1f} | "
                f"{int(usable_vol_vals.min())} | {int(usable_vol_vals.max())} |"
            )
    md_content.append("")

    # Excluded sessions detail
    if n_exclude > 0:
        md_content.append("## Excluded Sessions")
        md_content.append("")
        md_content.append("| Subject | Session | Max Trans (mm) | Max Rot (deg) | Mean FD (mm) | Reason | Salvageable |")
        md_content.append("|---------|---------|----------------|---------------|--------------|--------|-------------|")

        excluded = df[df['Status'] == 'EXCLUDE']
        for _, row in excluded.iterrows():
            salvage_status = row.get('Salvageable', 'N/A')
            md_content.append(
                f"| {row['Subject']} | {row['Session']} | "
                f"{row['Max Trans (mm)']:.2f} | {row['Max Rot (deg)']:.2f} | "
                f"{row['Mean FD (mm)']:.3f} | {row['Exclusion Reason']} | {salvage_status} |"
            )
        md_content.append("")

    # Salvage Analysis Section
    if salvage_analysis:
        md_content.append("## Salvage Analysis (Truncation Potential)")
        md_content.append("")
        md_content.append(f"The following excluded sessions contain clean segments of at least {min_salvage_minutes} minutes")
        md_content.append("that meet motion criteria and could potentially be salvaged via truncation.")
        md_content.append("")
        md_content.append("**Note:** Truncation preserves temporal continuity required for ΦID analysis,")
        md_content.append("unlike scrubbing which removes scattered volumes.")
        md_content.append("")

        n_salvageable = len(salvage_analysis)
        md_content.append(f"**Salvageable sessions:** {n_salvageable} of {n_exclude} excluded sessions")
        md_content.append("")

        md_content.append("| Subject | Session | Original Reason | Clean Segment | Duration (min) | Mean FD (mm) | Recommendation |")
        md_content.append("|---------|---------|-----------------|---------------|----------------|--------------|----------------|")

        for salvage in salvage_analysis:
            seg_range = f"vol {salvage['Best Segment Start']}-{salvage['Best Segment End']}"
            time_range = f"{salvage['Start Time']}-{salvage['End Time']}"
            md_content.append(
                f"| {salvage['Subject']} | {salvage['Session']} | "
                f"{salvage['Original Reason'][:40]}... | {seg_range} ({time_range}) | "
                f"{salvage['Best Segment Duration (min)']:.1f} | {salvage['Best Segment Mean FD']:.3f} | "
                f"Truncate to {salvage['Best Segment End'] - salvage['Best Segment Start']} volumes |"
            )
        md_content.append("")

        # Summary of salvage potential
        md_content.append("### Salvage Summary")
        md_content.append("")
        if n_salvageable > 0:
            total_salvage_minutes = sum(s['Best Segment Duration (min)'] for s in salvage_analysis)
            avg_salvage_minutes = total_salvage_minutes / n_salvageable
            md_content.append(f"- Sessions with salvageable segments: {n_salvageable}")
            md_content.append(f"- Average salvageable duration: {avg_salvage_minutes:.1f} minutes")
            md_content.append(f"- Minimum duration threshold: {min_salvage_minutes} minutes")
            md_content.append("")
            md_content.append("**If salvage is applied:**")
            md_content.append(f"- Original inclusion rate: {100*n_include/n_total:.1f}%")
            md_content.append(f"- Potential inclusion rate (with truncation): {100*(n_include + n_salvageable)/n_total:.1f}%")
        md_content.append("")

    # Full table for supplementary
    md_content.append("## Complete Motion Parameters (Supplementary Table)")
    md_content.append("")
    md_content.append(f"*Metrics computed on assessed volumes only (skipping first {skip_first}, last {skip_last})*")
    md_content.append("")
    md_content.append("| Subject | Session | Vols (assessed/total) | Max Trans (mm) | Max Rot (deg) | Mean FD (mm) | Max FD (mm) | Outliers | Status | Salvageable | Best Segment |")
    md_content.append("|---------|---------|----------------------|----------------|---------------|--------------|-------------|----------|--------|-------------|--------------|")

    for _, row in df.iterrows():
        status_emoji = "✓" if row['Status'] == 'INCLUDE' else "✗"
        salvage = row.get('Salvageable', '')
        best_seg = row.get('Best Segment Volumes', '') if pd.notna(row.get('Best Segment Duration (min)', np.nan)) else '-'
        if row['Status'] == 'INCLUDE':
            salvage = '-'
            best_seg = '-'
        n_assessed = row.get('N Volumes Assessed', row.get('N Volumes', 0))
        n_vols = row.get('N Volumes', 0)
        md_content.append(
            f"| {row['Subject']} | {row['Session']} | "
            f"{n_assessed}/{n_vols} | "
            f"{row['Max Trans (mm)']:.2f} | {row['Max Rot (deg)']:.2f} | "
            f"{row['Mean FD (mm)']:.3f} | {row['Max FD (mm)']:.3f} | "
            f"{row['Motion Outliers']} | {status_emoji} {row['Status']} | {salvage} | {best_seg} |"
        )
    md_content.append("")

    # QC checklist
    md_content.append("## Visual QC Checklist")
    md_content.append("")
    md_content.append("For each session, review the following fMRIPrep outputs:")
    md_content.append("")
    md_content.append("### Priority 1: Essential Checks")
    md_content.append("")
    md_content.append("| File | What to Check |")
    md_content.append("|------|---------------|")
    md_content.append("| `*_desc-carpetplot_bold.svg` | Look for vertical bands (motion), signal dropout, global signal shifts |")
    md_content.append("| `*_desc-coreg_bold.svg` | Verify BOLD-to-T1w alignment (red outlines should match brain edges) |")
    md_content.append("| `sub-XX.html` | Main report - check all sections for warnings/errors |")
    md_content.append("")
    md_content.append("### Priority 2: Additional Checks")
    md_content.append("")
    md_content.append("| File | What to Check |")
    md_content.append("|------|---------------|")
    md_content.append("| `*_desc-confoundcorr_bold.svg` | Correlations between confounds (high correlation = redundant regressors) |")
    md_content.append("| `*_desc-compcorvar_bold.svg` | Variance explained by CompCor components |")
    md_content.append("| `*_desc-rois_bold.svg` | ROI masks for CompCor (should cover WM/CSF properly) |")
    md_content.append("| `*_space-MNI152NLin2009cAsym_T1w.svg` | T1w to MNI normalization quality |")
    md_content.append("| `*_desc-reconall_T1w.svg` | FreeSurfer surface reconstruction quality |")
    md_content.append("")

    md_content.append("### Carpet Plot Interpretation")
    md_content.append("")
    md_content.append("- **Good**: Uniform horizontal bands, smooth transitions")
    md_content.append("- **Motion artifact**: Sharp vertical lines/bands")
    md_content.append("- **Signal dropout**: Dark patches or sudden intensity changes")
    md_content.append("- **Global signal shift**: All regions shift together (horizontal band)")
    md_content.append("")

    # Create output directory if needed
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Write markdown report
    report_file = output_path / 'motion_qc_report.md'
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write('\n'.join(md_content))
    print(f"\nReport saved to: {report_file}")

    # Save motion QC CSV with all metrics
    metrics_file = output_path / 'motion_qc.csv'
    df.to_csv(metrics_file, index=False)
    print(f"Motion metrics saved to: {metrics_file}")

    # Generate standardized exclusions.csv for pipeline integration
    exclusions_data = []
    for _, row in df.iterrows():
        # Subject column is now clean "sub-XX" format
        subject = str(row['Subject']).replace('sub-', '')
        session = row.get('Session') if has_sessions else None

        status = row['Status']
        reason = row.get('Exclusion Reason', '')

        # Handle truncation info
        truncate_start = None
        truncate_end = None

        if row.get('Salvageable') == 'YES' and status == 'EXCLUDE':
            # This session can be salvaged via truncation
            status = 'TRUNCATE'
            seg_volumes = row.get('Best Segment Volumes', '')
            if seg_volumes and '-' in str(seg_volumes):
                try:
                    parts = str(seg_volumes).split('-')
                    truncate_start = int(parts[0])
                    truncate_end = int(parts[1])
                except (ValueError, IndexError):
                    pass

        exclusions_data.append({
            'dataset': dataset_name,
            'subject': subject,
            'session': session if session else '',
            'task': row.get('Task', ''),
            'status': status,
            'reason': reason,
            'truncate_start': truncate_start if truncate_start else '',
            'truncate_end': truncate_end if truncate_end else '',
            'max_trans_mm': row.get('Max Trans (mm)', ''),
            'max_rot_deg': row.get('Max Rot (deg)', ''),
            'mean_fd_mm': row.get('Mean FD (mm)', ''),
            'usable_max_trans_mm': row.get('Usable Max Trans (mm)', ''),
            'usable_max_rot_deg': row.get('Usable Max Rot (deg)', ''),
            'usable_mean_fd_mm': row.get('Usable Mean FD (mm)', '')
        })

    exclusions_df = pd.DataFrame(exclusions_data)
    exclusions_file = output_path / 'exclusions.csv'
    exclusions_df.to_csv(exclusions_file, index=False)
    print(f"Exclusions file saved to: {exclusions_file}")

    # Print summary
    print(f"\n{'='*50}")
    print(f"SUMMARY")
    print(f"{'='*50}")
    print(f"Total sessions: {n_total}")
    print(f"To INCLUDE: {n_include}")
    print(f"To EXCLUDE: {n_exclude}")
    print(f"Exclusion rate: {100*n_exclude/n_total:.1f}%")

    if n_exclude > 0:
        print(f"\nExcluded sessions:")
        for _, row in excluded.iterrows():
            salvage_info = f" [SALVAGEABLE: {row['Best Segment Volumes']}]" if row.get('Salvageable') == 'YES' else ""
            print(f"  - {row['Subject']} {row['Session']}: {row['Exclusion Reason']}{salvage_info}")

    if salvage_analysis:
        print(f"\n{'='*50}")
        print(f"SALVAGE POTENTIAL")
        print(f"{'='*50}")
        print(f"Sessions with salvageable segments (>={min_salvage_minutes} min): {len(salvage_analysis)}")
        print(f"Potential new inclusion rate: {100*(n_include + len(salvage_analysis))/n_total:.1f}%")
        print(f"\nSalvageable sessions:")
        for salvage in salvage_analysis:
            print(f"  - {salvage['Subject']} {salvage['Session']}:")
            print(f"      Keep volumes {salvage['Best Segment Start']}-{salvage['Best Segment End']} "
                  f"({salvage['Best Segment Duration (min)']:.1f} min, mean FD={salvage['Best Segment Mean FD']:.3f})")

    return df, salvage_analysis


def generate_report_table(motion_qc_csv, output_path, dataset_config, tr):
    """
    Generate a clean, publication-ready report table from motion_qc.csv.

    Writes to results/qc/{dataset_id}/report_{dataset_id}.csv with
    human-friendly columns suitable for supplementary material.
    """
    df = pd.read_csv(motion_qc_csv)

    report_rows = []
    for _, row in df.iterrows():
        n_volumes = row.get('N Volumes', 0)
        scan_duration = (n_volumes * tr) / 60.0 if n_volumes > 0 else 0.0

        # Derive condition from session name
        # Preserve case for acronyms (LSD, PCB), title-case for words (ketamine -> Ketamine)
        session = str(row.get('Session', ''))
        if session:
            condition = session if session.isupper() else session.capitalize()
        else:
            condition = 'rest'

        n_assessed = row.get('N Volumes Assessed', n_volumes)

        report_rows.append({
            'Subject': row['Subject'],
            'Session': session,
            'Condition': condition,
            'Max Translation (mm)': round(row.get('Max Trans (mm)', np.nan), 3),
            'Max Rotation (deg)': round(row.get('Max Rot (deg)', np.nan), 3),
            'Mean FD (mm)': round(row.get('Mean FD (mm)', np.nan), 4),
            'Max FD (mm)': round(row.get('Max FD (mm)', np.nan), 3),
            'Outlier Volumes': int(row.get('Motion Outliers', 0)),
            'Total Volumes': int(n_volumes),
            'Assessed Volumes': int(n_assessed) if pd.notna(n_assessed) else int(n_volumes),
            'Scan Duration (min)': round(scan_duration, 1),
            'Status': row.get('Status', ''),
            'Exclusion Reason': row.get('Exclusion Reason', '')
        })

    report_df = pd.DataFrame(report_rows)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_df.to_csv(output_path, index=False)
    print(f"Report table saved to: {output_path}")
    return report_df


def generate_cross_dataset_summary(all_results, output_path):
    """
    Generate a cross-dataset summary table from results of multiple datasets.

    Parameters
    ----------
    all_results : list of dict
        Each dict has keys: dataset_id, dataset_name, motion_qc_csv, report_csv, tr
    output_path : str or Path
        Path to write report_summary.csv
    """
    summary_rows = []
    for res in all_results:
        try:
            df = pd.read_csv(res['motion_qc_csv'])
        except Exception:
            continue

        n_total = len(df)
        n_included = len(df[df['Status'] == 'INCLUDE'])
        n_excluded = len(df[df['Status'] == 'EXCLUDE'])
        n_truncated = len(df[df.get('Salvageable', pd.Series()) == 'YES']) if 'Salvageable' in df.columns else 0
        n_subjects = df['Subject'].nunique()
        n_sessions = n_total

        included = df[df['Status'] == 'INCLUDE']
        salvageable = df[df.get('Salvageable', pd.Series()) == 'YES'] if 'Salvageable' in df.columns else pd.DataFrame()

        # Usable metrics: INCLUDE assessed metrics + TRUNCATE post-trunc metrics
        if 'Usable Mean FD (mm)' in df.columns:
            usable = pd.concat([included, salvageable])
            mean_fd = round(usable['Usable Mean FD (mm)'].dropna().mean(), 4) if len(usable) > 0 else np.nan
            mean_max_trans = round(usable['Usable Max Trans (mm)'].dropna().mean(), 3) if len(usable) > 0 else np.nan
        else:
            # Fallback for older QC files without usable columns
            mean_fd = round(included['Mean FD (mm)'].mean(), 4) if len(included) > 0 else np.nan
            mean_max_trans = round(included['Max Trans (mm)'].mean(), 3) if len(included) > 0 else np.nan

        # Assessed volume counts: INCLUDE assessed volumes + TRUNCATE post-trunc volumes.
        # NOTE: this is the length of the assessed window (N Volumes - skip_first - skip_last,
        # or the post-truncation segment). It does NOT subtract motion-outlier volumes,
        # because the downstream pipeline does not scrub them (CLAUDE.md: "no ART, no despiking").
        # Motion-outlier counts are reported separately below.
        if 'Usable N Volumes' in df.columns:
            assessed_vols = pd.concat([included, salvageable])['Usable N Volumes'].dropna()
        elif 'N Volumes Assessed' in df.columns:
            assessed_vols = included['N Volumes Assessed'].dropna()
        else:
            assessed_vols = pd.Series(dtype=float)

        if len(assessed_vols) > 0:
            mean_assessed_vols = round(assessed_vols.mean(), 1)
            min_assessed_vols = int(assessed_vols.min())
            max_assessed_vols = int(assessed_vols.max())
            assessed_vols_range = f"{min_assessed_vols}-{max_assessed_vols}"
        else:
            mean_assessed_vols = np.nan
            min_assessed_vols = np.nan
            max_assessed_vols = np.nan
            assessed_vols_range = ''

        # Motion-flagged volumes (fMRIPrep motion_outlier_* columns) among usable sessions.
        # This surfaces per-scan residual motion variance that the assessed-window length hides
        # when a dataset has a uniform acquisition protocol (e.g. London LSD = 220 vols × run-3).
        if 'Motion Outliers' in df.columns:
            motion_flagged = pd.concat([included, salvageable])['Motion Outliers'].dropna()
        else:
            motion_flagged = pd.Series(dtype=float)

        if len(motion_flagged) > 0:
            mean_motion_flagged = round(motion_flagged.mean(), 1)
            min_motion_flagged = int(motion_flagged.min())
            max_motion_flagged = int(motion_flagged.max())
            motion_flagged_range = f"{min_motion_flagged}-{max_motion_flagged}"
        else:
            mean_motion_flagged = np.nan
            motion_flagged_range = ''

        n_usable = n_included + n_truncated

        summary_rows.append({
            'Dataset': res['dataset_name'],
            'N Subjects': n_subjects,
            'N Sessions': n_sessions,
            'Included': n_included,
            'Excluded': n_excluded,
            'Truncated': n_truncated,
            'Usable': n_usable,
            'Exclusion Rate (%)': round(100 * n_excluded / n_total, 1) if n_total > 0 else 0,
            'Usable Rate (%)': round(100 * n_usable / n_total, 1) if n_total > 0 else 0,
            'Mean FD Usable (mm)': mean_fd,
            'Mean Max Trans Usable (mm)': mean_max_trans,
            'Mean Assessed Volumes': mean_assessed_vols,
            'Min Assessed Volumes': min_assessed_vols,
            'Max Assessed Volumes': max_assessed_vols,
            'Assessed Volumes Range': assessed_vols_range,
            'Mean Motion-Flagged Volumes': mean_motion_flagged,
            'Motion-Flagged Volumes Range': motion_flagged_range
        })

    summary_df = pd.DataFrame(summary_rows)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(output_path, index=False)
    print(f"\nCross-dataset summary saved to: {output_path}")
    return summary_df


def auto_detect_environment():
    """Auto-detect environment based on OS."""
    system = platform.system()
    if system == 'Windows':
        return 'local'
    elif system == 'Linux':
        return 'hpc'
    else:
        return 'local'


def load_dataset_config(dataset_name, project_dir=None, env='auto'):
    """
    Load dataset configuration from YAML files.

    Parameters:
    -----------
    dataset_name : str
        Name of the dataset (e.g., 'berlin_ketamine')
    project_dir : Path, optional
        Project directory. If None, infers from script location.
    env : str
        Environment: 'local', 'hpc', or 'auto'

    Returns:
    --------
    dict
        Configuration dictionary with paths and QC settings
    """
    if project_dir is None:
        # Assume we're in scripts/qc/, go up two levels
        project_dir = Path(__file__).parent.parent.parent
    project_dir = Path(project_dir)

    # Determine environment
    if env == 'auto':
        env = auto_detect_environment()
    print(f"Using environment: {env}")

    # Load registry
    registry_path = project_dir / 'configs' / 'dataset_configs.yaml'
    with open(registry_path, 'r') as f:
        registry = yaml.safe_load(f)

    # Get dataset config path
    if dataset_name not in registry['datasets']:
        available = ', '.join(registry['datasets'].keys())
        raise ValueError(f"Dataset '{dataset_name}' not found. Available: {available}")

    config_file = registry['datasets'][dataset_name]['config_file']
    config_path = project_dir / config_file

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    # Resolve paths based on environment
    paths = config.get('paths', {})

    # Check for environment-specific paths
    if env in paths:
        fmriprep_dir = paths[env].get('fmriprep_dir', paths.get('fmriprep_dir'))
    else:
        fmriprep_dir = paths.get('fmriprep_dir')

    # Get QC settings with defaults
    qc_config = config.get('qc', {})
    motion_thresholds = qc_config.get('motion_thresholds', {})

    # Extract unique task names from config to filter QC sessions
    tasks_config = config.get('tasks', {})
    task_names = set()
    for task_info in tasks_config.values():
        bids_task = task_info.get('bids_task_name')
        if bids_task:
            task_names.add(bids_task)

    # Get subject and session filters from QC config
    include_sessions = qc_config.get('include_sessions')
    exclude_subjects = qc_config.get('exclude_subjects')
    # Also check subjects config for include/exclude
    subjects_config = config.get('subjects', {})
    include_subjects = subjects_config.get('include')
    if not exclude_subjects:
        exclude_subjects = subjects_config.get('exclude')

    # Get run selection strategy (matches what Script 01 uses)
    run_selection = config.get('run_selection', {}).get('strategy', 'first')

    return {
        'dataset_name': config.get('name', dataset_name),
        'dataset_id': config.get('dataset_id', dataset_name),
        'fmriprep_dir': fmriprep_dir,
        'has_sessions': config.get('bids_structure', {}).get('has_sessions', True),
        'trans_threshold': motion_thresholds.get('max_translation_mm', 3.0),
        'rot_threshold': motion_thresholds.get('max_rotation_deg', 3.0),
        'min_salvage_minutes': qc_config.get('min_salvage_duration_min', 8.0),
        'skip_first_volumes': config.get('preprocessing', {}).get('volumes_to_remove', 5),
        'skip_last_volumes': config.get('preprocessing', {}).get('volumes_to_remove_last', 0),
        'tr': config.get('preprocessing', {}).get('tr', 2.0),
        'project_dir': project_dir,
        'task_filter': list(task_names) if task_names else None,
        'include_sessions': include_sessions,
        'include_subjects': include_subjects,
        'exclude_subjects': exclude_subjects,
        'run_selection': run_selection
    }


def main():
    parser = argparse.ArgumentParser(
        description='Generate fMRIPrep motion QC report based on Luppi et al. criteria',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python fmriprep_motion_qc.py --dataset berlin_ketamine
  python fmriprep_motion_qc.py --dataset hcp --env local
  python fmriprep_motion_qc.py --all

Output files (in results/qc/{dataset}/):
  - motion_qc.csv: Raw motion metrics
  - exclusions.csv: Exclusion decisions for pipeline integration
  - motion_qc_report.md: Human-readable report
  - report_{dataset}.csv: Publication-ready report table

When using --all, also generates:
  - results/qc/report_summary.csv: Cross-dataset summary table
        """
    )
    parser.add_argument(
        '--dataset', '-d',
        type=str,
        default=None,
        help='Dataset name (e.g., berlin_ketamine, london_lsd, hcp)'
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help='Run QC on all configured datasets'
    )
    parser.add_argument(
        '--env', '-e',
        type=str,
        choices=['local', 'hpc', 'auto'],
        default='auto',
        help='Environment for paths (default: auto-detect based on OS)'
    )
    parser.add_argument(
        '--trans-threshold',
        type=float,
        default=None,
        help='Maximum translation threshold in mm (default: from config)'
    )
    parser.add_argument(
        '--rot-threshold',
        type=float,
        default=None,
        help='Maximum rotation threshold in degrees (default: from config)'
    )
    parser.add_argument(
        '--min-salvage-minutes',
        type=float,
        default=None,
        help='Minimum duration (minutes) for salvageable segment (default: from config)'
    )
    parser.add_argument(
        '--tr',
        type=float,
        default=None,
        help='Repetition time in seconds (default: from config)'
    )
    parser.add_argument(
        '--skip-first-volumes',
        type=int,
        default=None,
        help='Number of initial volumes to skip in QC assessment (default: from config, fallback 5)'
    )
    parser.add_argument(
        '--skip-last-volumes',
        type=int,
        default=None,
        help='Number of final volumes to skip in QC assessment (default: from config, fallback 5)'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Output directory (default: results/qc/{dataset}/)'
    )

    args = parser.parse_args()

    if not args.all and not args.dataset:
        parser.error("Either --dataset or --all is required")

    print("=" * 60)
    print("  Dataset-Agnostic Motion QC Script")
    print("=" * 60)
    print()

    # Determine which datasets to process
    if args.all:
        # Load registry to get all dataset names
        project_dir = Path(__file__).parent.parent.parent
        registry_path = project_dir / 'configs' / 'dataset_configs.yaml'
        with open(registry_path, 'r') as f:
            registry = yaml.safe_load(f)
        dataset_names = list(registry['datasets'].keys())
        print(f"Running QC on all {len(dataset_names)} datasets: {', '.join(dataset_names)}\n")
    else:
        dataset_names = [args.dataset]

    all_results = []

    for dataset_name in dataset_names:
        print()
        print("=" * 60)
        print(f"  Processing: {dataset_name}")
        print("=" * 60)
        print()

        # Load dataset configuration
        try:
            config = load_dataset_config(dataset_name, env=args.env)
        except Exception as e:
            print(f"ERROR: Failed to load dataset config for '{dataset_name}': {e}")
            if not args.all:
                return 1
            continue

        # Override with command line args if provided
        trans_threshold = args.trans_threshold if args.trans_threshold is not None else config['trans_threshold']
        rot_threshold = args.rot_threshold if args.rot_threshold is not None else config['rot_threshold']
        min_salvage_minutes = args.min_salvage_minutes if args.min_salvage_minutes is not None else config['min_salvage_minutes']
        tr = args.tr if args.tr is not None else config['tr']
        skip_first = args.skip_first_volumes if args.skip_first_volumes is not None else config['skip_first_volumes']
        skip_last = args.skip_last_volumes if args.skip_last_volumes is not None else config['skip_last_volumes']

        # Determine output directory
        if args.output_dir and not args.all:
            output_dir = Path(args.output_dir)
        else:
            output_dir = config['project_dir'] / 'results' / 'qc' / config['dataset_id']

        print(f"Dataset: {config['dataset_name']}")
        print(f"fMRIPrep directory: {config['fmriprep_dir']}")
        print(f"Has sessions: {config['has_sessions']}")
        print(f"Output directory: {output_dir}")
        print()
        print("QC Thresholds:")
        print(f"  Max translation: {trans_threshold} mm")
        print(f"  Max rotation: {rot_threshold} degrees")
        print(f"  Min salvage duration: {min_salvage_minutes} minutes")
        print(f"  TR: {tr} seconds")
        print(f"  Skip first volumes: {skip_first}")
        print(f"  Skip last volumes: {skip_last}")
        print()

        # Check if derivatives directory exists
        if not Path(config['fmriprep_dir']).exists():
            print(f"WARNING: fMRIPrep directory not found: {config['fmriprep_dir']}")
            print("Skipping this dataset.")
            if not args.all:
                return 1
            continue

        # Generate report with filters
        task_filter = config.get('task_filter')
        include_subjects = config.get('include_subjects')
        exclude_subjects = config.get('exclude_subjects')
        include_sessions = config.get('include_sessions')
        run_selection = config.get('run_selection', 'first')

        result = generate_qc_report(
            config['fmriprep_dir'],
            output_dir,
            dataset_name=config['dataset_name'],
            trans_threshold=trans_threshold,
            rot_threshold=rot_threshold,
            min_salvage_minutes=min_salvage_minutes,
            tr=tr,
            has_sessions=config['has_sessions'],
            skip_first=skip_first,
            skip_last=skip_last,
            task_filter=task_filter,
            include_subjects=include_subjects,
            exclude_subjects=exclude_subjects,
            include_sessions=include_sessions,
            run_selection=run_selection
        )

        if result is not None:
            # Generate clean report table
            motion_qc_csv = output_dir / 'motion_qc.csv'
            report_csv = output_dir / f"report_{config['dataset_id']}.csv"
            generate_report_table(motion_qc_csv, report_csv, config, tr)

            all_results.append({
                'dataset_id': config['dataset_id'],
                'dataset_name': config['dataset_name'],
                'motion_qc_csv': str(motion_qc_csv),
                'report_csv': str(report_csv),
                'tr': tr
            })

        print()
        print(f"  QC complete for {dataset_name}")
        print(f"  Output files in: {output_dir}")
        print("  - exclusions.csv: Use with pipeline scripts")
        print("  - motion_qc.csv: Full motion metrics")
        print("  - motion_qc_report.md: Human-readable report")
        print(f"  - report_{config['dataset_id']}.csv: Publication-ready report table")

    # Generate cross-dataset summary if multiple datasets processed
    if len(all_results) > 1:
        project_dir = Path(__file__).parent.parent.parent
        summary_path = project_dir / 'results' / 'qc' / 'report_summary.csv'
        generate_cross_dataset_summary(all_results, summary_path)

    print()
    print("=" * 60)
    print("  All QC Complete!")
    print("=" * 60)

    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main() or 0)
