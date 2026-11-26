"""Build corrected DOCX deliverables with the bundled document runtime."""
from pathlib import Path
import json
import re
import sys
import csv
import hashlib
import ast
import types
from docx.shared import Inches,Pt,RGBColor
from docx.oxml.ns import qn
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.section import WD_SECTION_START, WD_ORIENT

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'setup'))
FOLDER=ROOT/'manuscript/aachen_brief_report/corrected_20261126'
base=types.ModuleType('document_helpers')
helper_source=ROOT/'setup/build_aachen_manuscript_documents.py'
tree=ast.parse(helper_source.read_text(encoding='utf-8'))
nodes=[n for n in tree.body if isinstance(n,(ast.Import,ast.ImportFrom)) or isinstance(n,ast.FunctionDef) and n.name in ['clean','reference_text','document','table','main_table']]
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(helper_source),'exec'),base.__dict__)
main_text=(FOLDER/'manuscript_source.md').read_text(encoding='utf-8')
supp_text=(FOLDER/'supplement_source.md').read_text(encoding='utf-8')
base.FOLDER=FOLDER
keys=list(dict.fromkeys(re.findall(r'@([a-z]+\d+)',main_text+supp_text)))
garza={'authors':['Garza-Villarreal EA','Alcala-Lozano R','Fernandez-Lozano S','Morelos-Santana E','Dávalos A','Villicaña V','Alcauter S','Castellanos FX','Gonzalez-Olvera JJ'],
    'bibauthors':['Garza-Villarreal, Eduardo A.','Alcala-Lozano, Ruth','Fernandez-Lozano, Sofia','Morelos-Santana, Erik','Dávalos, Alan','Villicaña, Viviana','Alcauter, Sarael','Castellanos, F. Xavier','Gonzalez-Olvera, Jorge J.'],
    'title':'Clinical and Functional Connectivity Outcomes of 5-Hz Repetitive Transcranial Magnetic Stimulation as an Add-on Treatment in Cocaine Use Disorder: A Double-Blind Randomized Controlled Trial',
    'journal':'Biological Psychiatry: Cognitive Neuroscience and Neuroimaging','year':2021,'volume':'6','pages':'745–757',
    'doi':'10.1016/j.bpsc.2021.01.003','url':'https://pubmed.ncbi.nlm.nih.gov/33508499/'}
base.refs=json.loads((FOLDER/'reference_audit.json').read_text(encoding='utf-8'))['references']
base.refs['garza2021']=garza
base.numbers={k:i+1 for i,k in enumerate(keys)}
def citations(text):
    return re.sub(r'\[@[^\]]+\]',lambda m:'('+', '.join(str(base.numbers[k]) for k in re.findall(r'@([a-z]+\d+)',m[0]))+')',text)
base.citations=citations


def style_table(t):
    for row in t.rows:
        for cell in row.cells:
            cell.vertical_alignment=WD_CELL_VERTICAL_ALIGNMENT.CENTER
            for border in cell._tc.iter(qn('w:tcBorders')):
                for edge in border:edge.set(qn('w:color'),'D9D9D9')


def build(source,name):
    doc=base.document();section=doc.sections[0]
    section.page_width=Inches(8.5);section.page_height=Inches(11)
    section.left_margin=section.right_margin=Inches(.75)
    section.top_margin=section.bottom_margin=Inches(.7)
    doc.core_properties.title='Information dynamics following prefrontal theta burst stimulation'
    figures={'FIGURE1':'Figure_1_participant_changes.png','FIGURE2':'Figure_2_gwt_expectations.png',
        'SUPPLEMENT_FIGURE1':'Figure_S1_control_workspace_comparison.png','SUPPLEMENT_FIGURE2':'Figure_S2_leave_one_out.png','SUPPLEMENT_FIGURE3':'Figure_S3_hrf_width_associations.png'}
    used=list(dict.fromkeys(re.findall(r'@([a-z]+\d+)',source)))
    text=base.citations(source)
    for block in text.split('\n\n'):
        block=block.strip()
        if not block:continue
        tag=block.strip('{}')
        if tag in figures:
            if tag=='SUPPLEMENT_FIGURE1':
                sec=doc.add_section(WD_SECTION_START.NEW_PAGE);sec.orientation=WD_ORIENT.LANDSCAPE
                sec.page_width=Inches(11);sec.page_height=Inches(8.5)
                sec.top_margin=sec.bottom_margin=Inches(.5)
                width=9.5
            elif tag=='SUPPLEMENT_FIGURE2':
                sec=doc.add_section(WD_SECTION_START.NEW_PAGE);sec.orientation=WD_ORIENT.PORTRAIT
                sec.page_width=Inches(8.5);sec.page_height=Inches(11)
                sec.top_margin=sec.bottom_margin=Inches(.7)
                width=6.95
            else:doc.add_page_break();width=6.95
            p=doc.add_paragraph();p.paragraph_format.keep_with_next=True
            p.add_run().add_picture(str(FOLDER/'figures'/figures[tag]),width=Inches(width))
        elif tag in ['TABLE1','TABLE2']:
            base.main_table(doc,'table_1_scope_results.csv' if tag=='TABLE1' else 'table_2_scope_results.csv',tag=='TABLE2')
        elif tag=='SUPPLEMENT_TABLE':
            mapping=json.loads((FOLDER/'table_inventory.json').read_text())
            base.table(doc,['Table','Contents'],[(a,b) for a,b,_ in mapping],[.85,6.0])
        elif tag=='REFERENCES':
            for k in sorted(used,key=lambda k:base.numbers[k]):
                p=doc.add_paragraph(f'{base.numbers[k]}. '+base.reference_text(base.refs[k]));p.paragraph_format.line_spacing=1.;p.paragraph_format.space_after=Pt(5)
                for r in p.runs:r.font.size=Pt(10)
        elif block.startswith('# '):doc.add_paragraph(block[2:],style='Title')
        elif block.startswith('## '):
            heading=block[3:]
            if heading in ['Figures and tables','Supplementary figures']:continue
            if heading=='References':doc.add_page_break()
            doc.add_heading(heading,level=1)
        elif block.startswith('### '):doc.add_heading(block[4:],level=2)
        else:
            if block.startswith(('Table 1.','Table 2.')):doc.add_page_break()
            p=doc.add_paragraph(block)
            if block.startswith(('Figure ','Table ')):
                p.paragraph_format.line_spacing=1.0
                for r in p.runs:r.font.size=Pt(10.5)
                if block.startswith('Table '):p.paragraph_format.keep_with_next=True
            if '[TODO:' in block:
                for r in p.runs:r.font.color.rgb=RGBColor.from_string('7F2A1D')
    for t in doc.tables:style_table(t)
    assert '{{' not in '\n'.join(p.text for p in doc.paragraphs)
    path=FOLDER/name;doc.save(path)
    resolved=text.replace('{{REFERENCES}}','\n\n'.join(f'{base.numbers[k]}. '+base.reference_text(base.refs[k]) for k in used))
    for tag,file in figures.items():resolved=resolved.replace('{{'+tag+'}}',f'![{tag}](figures/{file})')
    for tag,file in [('TABLE1','table_1_scope_results.csv'),('TABLE2','table_2_scope_results.csv')]:
        with (FOLDER/file).open(encoding='utf-8',newline='') as stream:rows=list(csv.DictReader(stream))
        cols=['scope','measure','estimate','ci_low','ci_high','p_maxT_family','p_maxT_whole_and_workspace10']
        md=['| '+' | '.join(cols)+' |','| '+' | '.join(['---']*len(cols))+' |']
        md+=['| '+' | '.join(r[c] for c in cols)+' |' for r in rows]
        resolved=resolved.replace('{{'+tag+'}}','\n'.join(md))
    mapping=json.loads((FOLDER/'table_inventory.json').read_text())
    resolved=resolved.replace('{{SUPPLEMENT_TABLE}}','\n'.join(f'- {a}: {b} ({c})' for a,b,c in mapping))
    assert '{{' not in resolved
    (FOLDER/('supplement.md' if 'Supplementary' in name else 'manuscript.md')).write_text(resolved,encoding='utf-8')
    return path


def main():
    files=[build(main_text,'Aachen_Brief_Research_Report_Corrected.docx'),build(supp_text,'Supplementary_Material_Corrected.docx')]
    (FOLDER/'reference_audit.json').write_text(json.dumps({'keys':keys,'references':base.refs,'garza_verified_source':'https://pubmed.ncbi.nlm.nih.gov/33508499/'},ensure_ascii=False,indent=2),encoding='utf-8')
    bib=[]
    for k,r in base.refs.items():bib.append('@article{'+k+',\n author = {'+' and '.join(r['bibauthors'])+'},\n title = {{'+r['title']+'}},\n journal = {'+r['journal']+'},\n year = {'+str(r['year'])+'},\n doi = {'+r['doi']+'}\n}')
    (FOLDER/'references.bib').write_text('\n\n'.join(bib)+'\n',encoding='utf-8')
    receipt={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    (FOLDER/'build_receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8');print(json.dumps(receipt))


if __name__=='__main__':main()
