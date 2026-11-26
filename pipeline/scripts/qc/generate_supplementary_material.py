#!/usr/bin/env python3
"""
Generate Supplementary Material Word Document for Motion QC
===========================================================
Dataset-agnostic script that generates a publication-quality Word document
covering all datasets with QC results, or a single specified dataset.

Always generates a combined document. Uses dataset configs for parameterization
(TR, thresholds, session structure, etc.).

Usage:
    python generate_supplementary_material.py                    # All datasets with QC results
    python generate_supplementary_material.py --dataset hcp      # Single dataset
    python generate_supplementary_material.py --all              # Explicit all

Output:
    results/qc/Supplementary_Motion_QC.docx

Author: Generated for GWT Pipeline
"""

import argparse
import yaml
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from docx import Document
from docx.shared import Pt, Cm
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.oxml import OxmlElement


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def set_cell_shading(cell, color):
    """Set background color of a table cell."""
    shading = OxmlElement('w:shd')
    shading.set(qn('w:fill'), color)
    cell._tc.get_or_add_tcPr().append(shading)


def set_table_font(table, font_name='Times New Roman', font_size=Pt(10)):
    """Set font for all cells in a table."""
    for row in table.rows:
        for cell in row.cells:
            for paragraph in cell.paragraphs:
                for run in paragraph.runs:
                    run.font.name = font_name
                    run.font.size = font_size


def add_formatted_table(doc, headers, rows, header_color='D9E2F3'):
    """
    Add a formatted table with headers and data rows.

    Parameters
    ----------
    doc : Document
    headers : list of str
    rows : list of list/tuple (one per data row, values as strings)
    header_color : str  hex color for header shading
    """
    table = doc.add_table(rows=1, cols=len(headers))
    table.style = 'Table Grid'

    # Header
    for i, h in enumerate(headers):
        cell = table.rows[0].cells[i]
        cell.text = h
        cell.paragraphs[0].runs[0].bold = True
        cell.paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.CENTER
        set_cell_shading(cell, header_color)

    # Data
    for row_vals in rows:
        row_cells = table.add_row().cells
        for i, val in enumerate(row_vals):
            row_cells[i].text = str(val) if val is not None else '-'
            row_cells[i].paragraphs[0].alignment = WD_ALIGN_PARAGRAPH.CENTER

    set_table_font(table)
    return table


def fmt(val, decimals=3):
    """Format a numeric value, returning '-' for NaN."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return '-'
    return f'{val:.{decimals}f}'


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------

def get_project_dir():
    return Path(__file__).parent.parent.parent


def load_registry():
    project_dir = get_project_dir()
    registry_path = project_dir / 'configs' / 'dataset_configs.yaml'
    with open(registry_path, 'r') as f:
        return yaml.safe_load(f)


def load_dataset_config(dataset_id):
    """Load a single dataset config YAML."""
    project_dir = get_project_dir()
    registry = load_registry()

    if dataset_id not in registry['datasets']:
        return None

    config_file = registry['datasets'][dataset_id]['config_file']
    config_path = project_dir / config_file
    if not config_path.exists():
        return None

    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def get_dataset_tr(config, global_config):
    """Get TR for a dataset, falling back to global default."""
    preprocessing = config.get('preprocessing', {})
    if isinstance(preprocessing, dict) and preprocessing.get('tr'):
        return preprocessing['tr']
    return global_config.get('global_preprocessing', {}).get('tr', 2.0)


def get_motion_thresholds(config):
    """Get motion thresholds from dataset config."""
    qc = config.get('qc', {})
    mt = qc.get('motion_thresholds', {})
    return {
        'trans_mm': mt.get('max_translation_mm', 3.0),
        'rot_deg': mt.get('max_rotation_deg', 3.0)
    }


def get_min_salvage_minutes(config):
    qc = config.get('qc', {})
    return qc.get('min_salvage_duration_min', 8.0)


# ---------------------------------------------------------------------------
# Discover datasets with QC results
# ---------------------------------------------------------------------------

def discover_datasets_with_qc(dataset_filter=None):
    """
    Find all datasets that have QC report tables.

    Returns list of dicts with dataset metadata + paths.
    """
    project_dir = get_project_dir()
    registry = load_registry()
    qc_base = project_dir / 'results' / 'qc'

    datasets = []

    dataset_ids = [dataset_filter] if dataset_filter else list(registry['datasets'].keys())

    for ds_id in dataset_ids:
        if ds_id not in registry['datasets']:
            continue

        # Check for report table (preferred) or motion_qc.csv (fallback)
        report_csv = qc_base / ds_id / f'report_{ds_id}.csv'
        motion_csv = qc_base / ds_id / 'motion_qc.csv'

        if not report_csv.exists() and not motion_csv.exists():
            continue

        config = load_dataset_config(ds_id)
        if config is None:
            continue

        tr = get_dataset_tr(config, registry)
        thresholds = get_motion_thresholds(config)
        min_salvage = get_min_salvage_minutes(config)

        datasets.append({
            'dataset_id': ds_id,
            'dataset_name': config.get('name', ds_id),
            'report_csv': str(report_csv) if report_csv.exists() else None,
            'motion_csv': str(motion_csv),
            'tr': tr,
            'thresholds': thresholds,
            'min_salvage_min': min_salvage,
            'has_sessions': config.get('bids_structure', {}).get('has_sessions', True)
        })

    return datasets


def load_qc_data(ds_info):
    """Load QC data for a dataset, preferring report CSV over raw motion CSV."""
    if ds_info.get('report_csv') and Path(ds_info['report_csv']).exists():
        return pd.read_csv(ds_info['report_csv']), 'report'
    elif Path(ds_info['motion_csv']).exists():
        return pd.read_csv(ds_info['motion_csv']), 'raw'
    return None, None


# ---------------------------------------------------------------------------
# Document generation
# ---------------------------------------------------------------------------

def generate_supplementary_document(datasets_info, output_path):
    """Generate the combined supplementary material Word document."""

    doc = Document()

    # -- Default font --
    style = doc.styles['Normal']
    font = style.font
    font.name = 'Times New Roman'
    font.size = Pt(11)

    # -- Title --
    title = doc.add_heading('Supplementary Material', level=0)
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER

    # ======================================================================
    # S1. fMRI Motion Quality Control
    # ======================================================================
    doc.add_heading('S1. fMRI Motion Quality Control', level=1)

    # ------------------------------------------------------------------
    # S1.1 Methods
    # ------------------------------------------------------------------
    doc.add_heading('S1.1 Methods', level=2)

    # Collect fmriprep version from global config
    registry = load_registry()
    fmriprep_version = registry.get('fmriprep', {}).get('version', '23.2.1')

    # Collect unique threshold sets
    threshold_sets = set()
    for ds in datasets_info:
        t = ds['thresholds']
        threshold_sets.add((t['trans_mm'], t['rot_deg']))

    if len(threshold_sets) == 1:
        t_mm, r_deg = threshold_sets.pop()
        threshold_text = (
            f"Sessions were excluded if they exceeded either of the following thresholds: "
            f"maximum translation > {t_mm} mm in any direction, or maximum rotation > {r_deg}\u00b0 "
            f"in any axis."
        )
    else:
        threshold_text = (
            "Dataset-specific motion thresholds were applied (see Table S0 below). "
            "Sessions were excluded if they exceeded the maximum translation or rotation thresholds."
        )

    methods_para = (
        f"Functional MRI data quality was assessed using motion parameters derived from "
        f"fMRIPrep (version {fmriprep_version}) preprocessing. Head motion was quantified "
        f"using framewise displacement (FD) and the six rigid-body motion parameters "
        f"(three translations: x, y, z; three rotations: pitch, roll, yaw). "
        f"Quality control criteria were based on Luppi et al. (2024) methodology for "
        f"integrated information decomposition analysis. {threshold_text} "
        f"These conservative thresholds ensure data quality suitable for information-theoretic "
        f"analyses, which are sensitive to motion-related artifacts that can artificially inflate "
        f"or distort connectivity estimates (Power et al., 2012)."
    )
    doc.add_paragraph(methods_para)

    # Salvage / truncation methods
    # Collect min salvage durations and TR values for the description
    salvage_descs = []
    for ds in datasets_info:
        min_min = ds['min_salvage_min']
        tr = ds['tr']
        min_vols = int((min_min * 60) / tr)
        salvage_descs.append((ds['dataset_name'], min_min, min_vols, tr))

    salvage_para = (
        "Additionally, excluded sessions were evaluated for potential salvage through "
        "temporal truncation. For each excluded session, continuous segments meeting "
        "the motion criteria were identified using a sliding window approach. "
        "Unlike volume censoring (\"scrubbing\"), which removes scattered timepoints and "
        "disrupts temporal continuity, truncation preserves the continuous timeseries "
        "required for integrated information decomposition (PhiID) analysis."
    )
    doc.add_paragraph(salvage_para)

    # TR / threshold parameter table
    if len(datasets_info) > 1:
        doc.add_paragraph()
        p = doc.add_paragraph('Table S0. ')
        p.runs[0].bold = True
        p.add_run('Dataset acquisition and QC parameters.')

        param_headers = [
            'Dataset', 'TR (s)', 'Max Trans (mm)', 'Max Rot (\u00b0)',
            'Min Salvage (min)', 'Min Salvage (vols)'
        ]
        param_rows = []
        for ds in datasets_info:
            t = ds['thresholds']
            min_vols = int((ds['min_salvage_min'] * 60) / ds['tr'])
            param_rows.append([
                ds['dataset_name'],
                fmt(ds['tr'], 2),
                fmt(t['trans_mm'], 1),
                fmt(t['rot_deg'], 2),
                fmt(ds['min_salvage_min'], 0),
                str(min_vols)
            ])
        add_formatted_table(doc, param_headers, param_rows)
        doc.add_paragraph()

    # ------------------------------------------------------------------
    # S1.2 Cross-Dataset Summary
    # ------------------------------------------------------------------
    doc.add_heading('S1.2 Cross-Dataset Summary', level=2)

    summary_headers = [
        'Dataset', 'N Subjects', 'N Sessions', 'Included', 'Truncated',
        'Excluded', 'Usable', 'Usable Rate (%)', 'Mean FD Usable (mm)'
    ]
    summary_rows = []

    for ds in datasets_info:
        df, src = load_qc_data(ds)
        if df is None:
            continue

        n_total = len(df)
        n_subjects = df['Subject'].nunique() if 'Subject' in df.columns else 0

        # Handle column name differences between report and raw CSVs
        status_col = 'Status'
        if status_col not in df.columns:
            continue

        n_include = len(df[df[status_col] == 'INCLUDE'])
        n_exclude = len(df[df[status_col] == 'EXCLUDE'])
        n_truncate = 0
        if 'Salvageable' in df.columns:
            n_truncate = len(df[df['Salvageable'] == 'YES'])

        n_usable = n_include + n_truncate

        # Mean FD for usable sessions (INCLUDE + TRUNCATE with post-trunc metrics)
        included = df[df[status_col] == 'INCLUDE']
        salvageable = df[df['Salvageable'] == 'YES'] if 'Salvageable' in df.columns else pd.DataFrame()
        if 'Usable Mean FD (mm)' in df.columns:
            usable = pd.concat([included, salvageable])
            mean_fd = usable['Usable Mean FD (mm)'].dropna().mean() if len(usable) > 0 else np.nan
        else:
            fd_col = 'Mean FD (mm)'
            mean_fd = included[fd_col].mean() if fd_col in included.columns and len(included) > 0 else np.nan

        usable_rate = 100 * n_usable / n_total if n_total > 0 else 0

        summary_rows.append([
            ds['dataset_name'],
            str(n_subjects),
            str(n_total),
            str(n_include),
            str(n_truncate),
            str(n_exclude - n_truncate),  # Truly excluded
            str(n_usable),
            fmt(usable_rate, 1),
            fmt(mean_fd, 4)
        ])

    if summary_rows:
        p = doc.add_paragraph('Table S1. ')
        p.runs[0].bold = True
        p.add_run('Cross-dataset motion quality control summary.')
        add_formatted_table(doc, summary_headers, summary_rows)
        doc.add_paragraph()

    # ------------------------------------------------------------------
    # S1.3+ Per-Dataset Results
    # ------------------------------------------------------------------
    table_counter = 2  # S1 was the summary table

    for idx, ds in enumerate(datasets_info):
        section_num = idx + 3  # S1.3, S1.4, ...
        doc.add_heading(f'S1.{section_num} {ds["dataset_name"]}', level=2)

        df, src = load_qc_data(ds)
        if df is None:
            doc.add_paragraph('No QC data available for this dataset.')
            continue

        status_col = 'Status'
        n_total = len(df)
        n_include = len(df[df[status_col] == 'INCLUDE'])
        n_exclude = len(df[df[status_col] == 'EXCLUDE'])
        n_subjects = df['Subject'].nunique()

        # Determine column names based on source type
        if src == 'report':
            trans_col = 'Max Translation (mm)'
            rot_col = 'Max Rotation (deg)'
            fd_col = 'Mean FD (mm)'
            maxfd_col = 'Max FD (mm)'
            outlier_col = 'Outlier Volumes'
            nvol_col = 'Total Volumes'
        else:
            trans_col = 'Max Trans (mm)'
            rot_col = 'Max Rot (deg)'
            fd_col = 'Mean FD (mm)'
            maxfd_col = 'Max FD (mm)'
            outlier_col = 'Motion Outliers'
            nvol_col = 'N Volumes'

        # Salvageable count
        n_salvageable = 0
        if 'Salvageable' in df.columns:
            n_salvageable = len(df[df['Salvageable'] == 'YES'])

        # Results text
        n_usable = n_include + n_salvageable
        n_truly_excluded = n_exclude - n_salvageable
        session_info = f" across {n_total} sessions" if ds['has_sessions'] else ""
        results_text = (
            f"A total of {n_subjects} subjects{session_info} were evaluated. "
            f"Based on the motion criteria, {n_include} sessions ({100*n_include/n_total:.1f}%) "
            f"met quality criteria directly (INCLUDE) and {n_exclude} ({100*n_exclude/n_total:.1f}%) "
            f"exceeded motion thresholds."
        )
        if n_salvageable > 0:
            results_text += (
                f" Of the {n_exclude} flagged sessions, {n_salvageable} contained salvageable "
                f"continuous segments meeting motion criteria (TRUNCATE), "
                f"yielding {n_usable} usable sessions ({100*n_usable/n_total:.1f}%) "
                f"and {n_truly_excluded} truly excluded."
            )
        doc.add_paragraph(results_text)

        # --- Table Sn.1: Motion statistics for usable sessions ---
        included = df[df[status_col] == 'INCLUDE']
        salvageable_df = df[df['Salvageable'] == 'YES'] if 'Salvageable' in df.columns else pd.DataFrame()
        usable_df = pd.concat([included, salvageable_df]) if len(salvageable_df) > 0 else included

        if len(usable_df) > 0:
            doc.add_paragraph()
            p = doc.add_paragraph(f'Table S{table_counter}. ')
            p.runs[0].bold = True
            caption = f'Motion statistics for usable sessions ({ds["dataset_name"]}). '
            if len(salvageable_df) > 0:
                caption += f'Includes {len(included)} INCLUDE + {len(salvageable_df)} TRUNCATE sessions. '
                caption += 'TRUNCATE metrics are from the clean segment entering analysis.'
            p.add_run(caption)
            table_counter += 1

            stat_headers = ['Metric', 'Mean', 'SD', 'Min', 'Max']
            stat_rows = []

            # Use "Usable" columns if available, fallback to original
            usable_cols = [
                ('Usable Max Trans (mm)', trans_col, 'Max Translation (mm)'),
                ('Usable Max Rot (deg)', rot_col, 'Max Rotation (deg)'),
                ('Usable Mean FD (mm)', fd_col, 'Mean FD (mm)')
            ]
            for usable_col, fallback_col, label in usable_cols:
                col_name = usable_col if usable_col in usable_df.columns else fallback_col
                if col_name in usable_df.columns:
                    vals = usable_df[col_name].dropna()
                    if len(vals) > 0:
                        stat_rows.append([
                            label,
                            fmt(vals.mean()),
                            fmt(vals.std()),
                            fmt(vals.min()),
                            fmt(vals.max())
                        ])
            add_formatted_table(doc, stat_headers, stat_rows)

        # --- Table Sn.2: Excluded sessions with reasons ---
        excluded = df[df[status_col] == 'EXCLUDE']
        if len(excluded) > 0:
            doc.add_paragraph()
            p = doc.add_paragraph(f'Table S{table_counter}. ')
            p.runs[0].bold = True
            p.add_run(f'Excluded sessions ({ds["dataset_name"]}).')
            table_counter += 1

            excl_headers = ['Subject', 'Session', 'Max Trans (mm)', 'Max Rot (deg)',
                            'Mean FD (mm)', 'Exclusion Reason']
            excl_rows = []
            for _, row in excluded.iterrows():
                session_val = str(row.get('Session', '')) if pd.notna(row.get('Session')) else ''
                reason_col = 'Exclusion Reason' if 'Exclusion Reason' in row.index else 'Reason'
                reason = str(row.get(reason_col, ''))
                excl_rows.append([
                    str(row['Subject']),
                    session_val,
                    fmt(row.get(trans_col, np.nan)),
                    fmt(row.get(rot_col, np.nan)),
                    fmt(row.get(fd_col, np.nan)),
                    reason
                ])
            add_formatted_table(doc, excl_headers, excl_rows)

        # --- Table Sn.3: Complete motion parameters ---
        doc.add_paragraph()
        p = doc.add_paragraph(f'Table S{table_counter}. ')
        p.runs[0].bold = True
        p.add_run(f'Complete motion parameters ({ds["dataset_name"]}).')
        table_counter += 1

        full_headers = ['Subject', 'Session', 'Max Trans (mm)', 'Max Rot (deg)',
                        'Mean FD (mm)', 'Max FD (mm)', 'Outlier Vols', 'Total Vols', 'Status']
        full_rows = []
        for _, row in df.iterrows():
            session_val = str(row.get('Session', '')) if pd.notna(row.get('Session')) else ''
            full_rows.append([
                str(row['Subject']),
                session_val,
                fmt(row.get(trans_col, np.nan)),
                fmt(row.get(rot_col, np.nan)),
                fmt(row.get(fd_col, np.nan), 4),
                fmt(row.get(maxfd_col, np.nan)),
                str(int(row.get(outlier_col, 0))) if pd.notna(row.get(outlier_col)) else '-',
                str(int(row.get(nvol_col, 0))) if pd.notna(row.get(nvol_col)) else '-',
                str(row.get(status_col, ''))
            ])
        add_formatted_table(doc, full_headers, full_rows)
        doc.add_paragraph()

    # ------------------------------------------------------------------
    # References
    # ------------------------------------------------------------------
    doc.add_heading('References', level=1)

    refs = [
        "Esteban, O., Markiewicz, C. J., Blair, R. W., et al. (2019). fMRIPrep: a robust "
        "preprocessing pipeline for functional MRI. Nature Methods, 16, 111\u2013116.",

        "Luppi, A. I., Mediano, P. A. M., Rosas, F. E., et al. (2024). A synergistic workspace "
        "for human consciousness revealed by Integrated Information Decomposition. eLife, 12, e88173.",

        "Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). "
        "Spurious but systematic correlations in functional connectivity MRI networks arise from "
        "subject motion. NeuroImage, 59(3), 2142\u20132154."
    ]
    for ref in refs:
        doc.add_paragraph(ref, style='List Bullet')

    # ------------------------------------------------------------------
    # Save
    # ------------------------------------------------------------------
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    doc.save(str(output_path))
    print(f"\nSupplementary material saved to: {output_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Generate publication-quality supplementary material for motion QC',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python generate_supplementary_material.py              # All datasets with QC results
  python generate_supplementary_material.py --dataset hcp
  python generate_supplementary_material.py --all

Output:
  results/qc/Supplementary_Motion_QC.docx
        """
    )
    parser.add_argument(
        '--dataset', '-d',
        type=str,
        default=None,
        help='Generate supplement for a single dataset'
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help='Generate combined supplement for all datasets with QC results'
    )
    parser.add_argument(
        '--output', '-o',
        type=str,
        default=None,
        help='Output path for the Word document'
    )

    args = parser.parse_args()

    # Default behaviour: generate for all datasets that have QC results
    dataset_filter = args.dataset if args.dataset and not args.all else None

    print("=" * 60)
    print("  Supplementary Material Generator")
    print("=" * 60)
    print()

    datasets_info = discover_datasets_with_qc(dataset_filter)

    if not datasets_info:
        if dataset_filter:
            print(f"ERROR: No QC results found for dataset '{dataset_filter}'")
            print("Run fmriprep_motion_qc.py first to generate QC data.")
        else:
            print("ERROR: No QC results found for any dataset.")
            print("Run fmriprep_motion_qc.py --all first.")
        return 1

    print(f"Found QC results for {len(datasets_info)} dataset(s):")
    for ds in datasets_info:
        print(f"  - {ds['dataset_name']} (TR={ds['tr']}s)")
    print()

    # Determine output path
    if args.output:
        output_path = Path(args.output)
    else:
        project_dir = get_project_dir()
        output_path = project_dir / 'results' / 'qc' / 'Supplementary_Motion_QC.docx'

    generate_supplementary_document(datasets_info, output_path)

    print()
    print("=" * 60)
    print("  Done!")
    print("=" * 60)

    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main() or 0)
