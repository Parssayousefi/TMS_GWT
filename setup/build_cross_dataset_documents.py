"""Build the two-dataset report using the bundled document runtime."""
from pathlib import Path
import csv
import hashlib
import json
import re
import sys
from docx.shared import Inches, Pt
from docx.enum.section import WD_SECTION_START, WD_ORIENT

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'setup'))
import build_corrected_manuscript_documents as helpers

OUT = ROOT / 'manuscript/tms_cross_dataset_report/20261126'
TITLE = 'Haemodynamic covariation of information integration across two prefrontal stimulation datasets'
MAIN = (OUT / 'manuscript_source.md').read_text(encoding='utf-8')
SUPP = (OUT / 'supplement_source.md').read_text(encoding='utf-8')
REFS = json.loads((OUT / 'reference_audit.json').read_text(encoding='utf-8'))['references']
KEYS = list(dict.fromkeys(re.findall(r'@([a-z]+\d+)', MAIN + SUPP)))
NUMBERS = {k: i + 1 for i, k in enumerate(KEYS)}
FIGURES = {'FIGURE1': 'Figure_1_hrf_width.png', 'FIGURE2': 'Figure_2_scope_associations.png', 'SUPPLEMENT_FIGURE1': 'Figure_S1_control_workspace_comparison.png'}


def read(name):
    with (OUT / 'tables' / name).open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))


def interval(row, scale=1):
    return f"{scale*float(row['estimate']):.5f} [{scale*float(row['ci_low']):.5f}, {scale*float(row['ci_high']):.5f}]"


def tables():
    names = {'aachen_itbs': 'Aachen', 'sudmex_tms': 'SUDMEX'}
    scopes = {'whole_brain': 'Whole brain', 'within_workspace': 'Within', 'outside_workspace': 'Outside'}
    width = []
    for r in read('hrf_association_results.csv'):
        if r['model'] != 'motion_arm_adjusted' or r['parameter'] != 'width_seconds' or r['outcome_scope'] != r['predictor_scope'] or r['outcome_scope'] not in scopes:
            continue
        if r['branch'] == 'normalized' and r['outcome_scope'] != 'outside_workspace':
            continue
        width.append([names[r['dataset']], r['branch'], scopes[r['outcome_scope']], interval(r, 100 if r['branch'] == 'normalized' else 1), f"{float(r['p_holm15']):.5f}"])
    spatial = []
    labels = {'synergy': 'Persistent synergy', 'redundancy': 'Persistent redundancy', 'synergy_total': 'Source synergy', 'phi': 'Φ', 'phi_r': 'ΦR'}
    for r in read('spatial_results.csv'):
        if r['branch'] != 'raw' or not ((r['dataset'] == 'aachen_itbs' and r['model'] == 'paired_change') or (r['dataset'] == 'sudmex_tms' and r['model'] == 'baseline_motion_adjusted')):
            continue
        p = r['p_maxT5'] if r['dataset'] == 'aachen_itbs' else r['p_holm5']
        spatial.append([names[r['dataset']], labels[r['measure']], interval(r), f'{float(p):.4f}'])
    return {'WIDTH_TABLE': (['Dataset', 'Scale', 'Scope', 'Estimate [95% CI]', 'Holm15 p'], width, [.72, .80, .88, 3.12, .75]), 'SPATIAL_TABLE': (['Dataset', 'Measure', 'Estimate [95% CI]', 'Adjusted p'], spatial, [.72, 1.38, 3.37, .80])}


def build(source, filename):
    doc = helpers.base.document()
    section = doc.sections[0]
    section.page_width = Inches(8.5); section.page_height = Inches(11)
    section.left_margin = section.right_margin = Inches(.75)
    section.top_margin = section.bottom_margin = Inches(.7)
    doc.core_properties.title = TITLE
    doc.core_properties.author = ''; doc.core_properties.last_modified_by = ''
    text = re.sub(r'\[@[^\]]+\]', lambda m: '(' + ', '.join(str(NUMBERS[k]) for k in re.findall(r'@([a-z]+\d+)', m[0])) + ')', source)
    used = list(dict.fromkeys(re.findall(r'@([a-z]+\d+)', source)))
    table_data = tables()
    for block in text.split('\n\n'):
        block = block.strip()
        if not block: continue
        tag = block.strip('{}')
        if tag in FIGURES:
            if tag == 'SUPPLEMENT_FIGURE1':
                sec = doc.add_section(WD_SECTION_START.NEW_PAGE)
                sec.orientation = WD_ORIENT.LANDSCAPE; sec.page_width = Inches(11); sec.page_height = Inches(8.5)
                sec.top_margin = sec.bottom_margin = Inches(.5)
                width = 9.5
            else:
                doc.add_page_break(); width = 6.95
            p = doc.add_paragraph(); p.paragraph_format.keep_with_next = True
            p.add_run().add_picture(str(OUT / 'figures' / FIGURES[tag]), width=Inches(width))
        elif tag in table_data:
            headers, rows, widths = table_data[tag]
            helpers.base.table(doc, headers, rows, widths)
        elif tag == 'REFERENCES':
            for k in sorted(used, key=lambda k: NUMBERS[k]):
                p = doc.add_paragraph(f'{NUMBERS[k]}. ' + helpers.base.reference_text(REFS[k]))
                p.paragraph_format.line_spacing = 1.; p.paragraph_format.space_after = Pt(5)
                for run in p.runs: run.font.size = Pt(10)
        elif block.startswith('# '): doc.add_paragraph(block[2:], style='Title')
        elif block.startswith('## '):
            heading = block[3:]
            if heading in ['Figures and tables', 'Supplementary figure']: continue
            if heading in ['References', 'Supplementary data']:
                if doc.sections[-1].orientation == WD_ORIENT.LANDSCAPE:
                    sec = doc.add_section(WD_SECTION_START.NEW_PAGE)
                    sec.orientation = WD_ORIENT.PORTRAIT; sec.page_width = Inches(8.5); sec.page_height = Inches(11)
                    sec.top_margin = sec.bottom_margin = Inches(.7)
                else: doc.add_page_break()
            doc.add_heading(heading, level=1)
        elif block.startswith('### '): doc.add_heading(block[4:], level=2)
        else:
            if block.startswith('Table '): doc.add_page_break()
            p = doc.add_paragraph(block)
            if block.startswith(('Figure ', 'Table ')):
                p.paragraph_format.line_spacing = 1.
                for run in p.runs: run.font.size = Pt(10.5)
                if block.startswith('Table '): p.paragraph_format.keep_with_next = True
    for table in doc.tables: helpers.style_table(table)
    assert '{{' not in '\n'.join(p.text for p in doc.paragraphs)
    doc.save(OUT / filename)
    resolved = text.replace('{{REFERENCES}}', '\n\n'.join(f'{NUMBERS[k]}. ' + helpers.base.reference_text(REFS[k]) for k in sorted(used, key=lambda k: NUMBERS[k])))
    for tag, name in FIGURES.items(): resolved = resolved.replace('{{' + tag + '}}', f'![{tag}](figures/{name})')
    for tag, (headers, rows, _) in table_data.items():
        md = '| ' + ' | '.join(headers) + ' |\n| ' + ' | '.join(['---'] * len(headers)) + ' |\n'
        md += '\n'.join('| ' + ' | '.join(row) + ' |' for row in rows)
        resolved = resolved.replace('{{' + tag + '}}', md)
    assert '{{' not in resolved
    (OUT / ('supplement.md' if 'Supplement' in filename else 'manuscript.md')).write_text(resolved, encoding='utf-8')
    return OUT / filename


if __name__ == '__main__':
    files = [build(MAIN, 'TMS_Cross_Dataset_Brief_Report.docx'), build(SUPP, 'TMS_Cross_Dataset_Supplement.docx')]
    receipt = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    (OUT / 'build_receipt.json').write_text(json.dumps(receipt, indent=2), encoding='utf-8')
    print(json.dumps(receipt))
