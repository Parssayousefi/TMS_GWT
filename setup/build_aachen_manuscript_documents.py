"""Create manuscript/supplement DOCX using the bundled document runtime."""
from pathlib import Path
import argparse
import csv
import html
import json
import re
import shutil
import zipfile

from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.text import WD_ALIGN_PARAGRAPH

ROOT=Path(__file__).resolve().parents[1]
FOLDER=ROOT/'manuscript/aachen_brief_report'
SOURCES=FOLDER/'sources'
RUN=ROOT/'output/tms_workspace_scope'/(ROOT/'output/tms_workspace_scope/LATEST.txt').read_text().strip()
DATA=FOLDER/'submission_data'
DATA.mkdir(exist_ok=True)


def clean(text):
    return re.sub(r'\s+',' ',html.unescape(re.sub(r'<[^>]+>','',str(text)))).strip().replace('Conn :','Conn:')


def initials(given):
    return ''.join(s[0] for s in re.findall(r'[^\W\d_]+',given))


def reference(key):
    if key=='holm1979':
        return {'authors':['Holm S'],'bibauthors':['Holm, Sture'],'title':'A simple sequentially rejective multiple test procedure',
                'journal':'Scandinavian Journal of Statistics','year':1979,'volume':'6','pages':'65–70','doi':'','url':'https://www.jstor.org/stable/4615733'}
    if key in ['yousefi2026','blackburne2026']:
        rows=json.loads((SOURCES/(key+'_biorxiv.json')).read_text(encoding='utf-8'))['collection']
        row=max(rows,key=lambda x:int(x['version']))
        authors=[];bib=[]
        for author in row['authors'].split(';'):
            family,given=author.strip().split(',',1)
            authors.append(family+' '+initials(given));bib.append(family+', '+given.strip())
        return {'authors':authors,'bibauthors':bib,'title':clean(row['title']),'journal':'bioRxiv [Preprint]',
                'year':int(row['date'][:4]),'volume':'','pages':'version '+row['version'],
                'doi':row['doi']+'v'+row['version'],'url':'https://doi.org/'+row['doi']+'v'+row['version']}
    raw=json.loads((SOURCES/(key+'.json')).read_text(encoding='utf-8'));row=raw.get('message',raw)
    return {'authors':[a.get('family',a.get('name',''))+' '+initials(a.get('given','')) for a in row['author']],
            'bibauthors':[a.get('family',a.get('name',''))+', '+a.get('given','') for a in row['author']],
            'title':clean(row['title'][0]),'journal':clean(row['container-title'][0]),
            'year':row.get('published-print',row['published'])['date-parts'][0][0],
            'volume':row.get('volume',''),'pages':row.get('page',row.get('article-number','')),
            'doi':row.get('DOI',''),'url':'https://doi.org/'+row.get('DOI','')}


def reference_text(r):
    authors=', '.join(r['authors'][:6])+(', et al.' if len(r['authors'])>6 else '.')
    tail=f"{r['journal']}. {r['year']}"+(f";{r['volume']}" if r['volume'] else '')+(f":{r['pages']}" if r['pages'] else '')+'.'
    return f"{authors} {r['title']}. {tail}"+(f" doi:{r['doi']}." if r['doi'] else '')


main_text=(FOLDER/'manuscript_source.md').read_text(encoding='utf-8')
supp_text=(FOLDER/'supplement_source.md').read_text(encoding='utf-8')
keys=list(dict.fromkeys(re.findall(r'@([a-z]+\d+)',main_text+supp_text)))
refs={k:reference(k) for k in keys};numbers={k:i+1 for i,k in enumerate(keys)}


def citations(text):
    def sub(match):return '('+', '.join(str(numbers[k]) for k in re.findall(r'@([a-z]+\d+)',match[0]))+')'
    return re.sub(r'\[@[^\]]+\]',sub,text)


def prepare_data():
    mapping=[
      ('Table S1','Main endpoints','aachen_main_results.csv','TableS1_main.csv'),
      ('Table S2','Motion models for all spatial scopes','aachen_all_spatial_motion.csv','TableS2_motion.csv'),
      ('Table S3','Intensity and adjusted intensity comparisons','aachen_intensity_results.csv','TableS3_intensity.csv'),
      ('Table S4','Independent conditional-reference seeds','normalization_seed_stability.csv','TableS4_seeds.csv'),
      ('Table S5','Motion, bootstrap and signed-rank sensitivities','aachen_robustness.csv','TableS5_robustness.csv'),
      ('Table S6','Leave-one-participant-out results','aachen_leave_one_out.csv','TableS6_influence.csv'),
      ('Table S8a','Four spatial scopes and 20-test correction','aachen_all_spatial20_results.csv','TableS8a_spatial20.csv'),
      ('Table S8b','Direct outside-minus-within change contrasts','aachen_spatial_interactions.csv','TableS8b_interactions.csv'),
      ('Table S8c','Raw 40-endpoint sensitivity','raw_aachen_expanded40_sensitivity.csv','TableS8c_raw40.csv'),
      ('Table S8d','Normalized 40-endpoint sensitivity','normalized_aachen_expanded40_sensitivity.csv','TableS8d_normalized40.csv')]
    for _,_,source,target in mapping:shutil.copyfile(RUN/source,DATA/target)
    shutil.copyfile(FOLDER/'supplementary_data/Table_S7_null_model_calibration.csv',DATA/'TableS7_calibration.csv')
    mapping.insert(6,('Table S7','Finite-record null-model calibration','', 'TableS7_calibration.csv'))
    # Preserve the historical global bootstrap used in the main text, with its own seed/provenance.
    shutil.copyfile(FOLDER/'supplementary_data/Table_S5_additional_robustness.csv',DATA/'TableS5b_original_global.csv')
    mapping.append(('Table S5b','Original whole-brain bootstrap and checks','','TableS5b_original_global.csv'))
    for branch in ['raw','normalized']:
        for kind in ['scan_metrics','paired_changes']:
            source=RUN/(branch+'_'+kind+'.csv')
            with source.open(encoding='utf-8',newline='') as f:
                reader=csv.DictReader(f);fields=reader.fieldnames;rows=[r for r in reader if r['dataset']=='aachen_itbs']
            with (DATA/(branch+'_'+kind+'.csv')).open('w',encoding='utf-8',newline='') as f:
                writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
        shutil.copyfile(RUN/(branch+'_aachen_pooled_changes.csv'),DATA/(branch+'_participant_averages.csv'))
    for name in ['fixed_hcp_workspace.csv','independent_verification.json','analysis_specification.json','provenance.json','report_provenance.json']:
        shutil.copyfile(RUN/name,DATA/name)
    follow=ROOT/'output/aachen_spatial_followup'/(ROOT/'output/aachen_spatial_followup/LATEST.txt').read_text().strip()
    extra=['estimator_sensitivity_results','method_interactions','temporal_component_results','network_results',
           'visit_consistency_results','baseline_motion_sensitivity','nonparametric_checks','participant_influence',
           'network_pair_counts','network_node_counts']
    for i,name in enumerate(extra):shutil.copyfile(follow/(name+'.csv'),DATA/('TableS9'+chr(97+i)+'_'+name+'.csv'))
    mapping.append(('Table S9a–j','Signal representation, temporal components, network and visit diagnostics','','TableS9*.csv'))
    for i,name in enumerate(['hrf_parameter_tests','hrf_phi_associations']):
        shutil.copyfile(follow/(name+'.csv'),DATA/('TableS10'+chr(97+i)+'_'+name+'.csv'))
    mapping.append(('Table S10a–b','Estimated HRF parameters and their associations with Φ','','TableS10*.csv'))
    for name in ['signal_decomposition_scans','signal_participant_changes','network_profile_scans','network_participant_changes','hrf_parameter_scans']:
        shutil.copyfile(follow/(name+'.csv'),DATA/(name+'.csv'))
    for name in ['provenance.json','followup_verification.json','analysis_scope.json']:
        shutil.copyfile(follow/name,DATA/('spatial_followup_'+name))
    (DATA/'README.txt').write_text('Aachen manuscript supplementary data\n\n'
      'All manuscript participant outputs contain Aachen only. n=15 people, 60 scans.\n'
      'Raw quantities are nats; normalized CSV estimates are CDF units (multiply by 100 for percentile points).\n'
      'Main tables 1 and 2 already use the displayed manuscript units. Confidence intervals are marginal.\n'
      'Whole brain has 102831 unique pairs; within workspace 27495; outside workspace 23871; between 51465.\n'
      'Families: five information measures per scope; 10 across whole/workspace; 30 with original PC/roles;\n'
      '20 across four scopes; 40 with all PC/roles. No file controls the entire historical exploratory search.\n'
      'Table S7 is simulation calibration, not participant data. Table S5b uses the historical bootstrap seed.\n'
      'The provenance files document the full project audit, which also checked SUDMEX; those data are not included here.\n\n'
      +'\n'.join(f'{label}: {description} — {target}' for label,description,_,target in mapping)+'\n',encoding='utf-8')
    return mapping


def document():
    doc=Document();section=doc.sections[0]
    section.page_width=Inches(8.27);section.page_height=Inches(11.69)
    section.top_margin=Inches(.8);section.bottom_margin=Inches(.8)
    section.left_margin=Inches(.8);section.right_margin=Inches(.8)
    normal=doc.styles['Normal'];normal.font.name='Times New Roman';normal.font.size=Pt(11.5)
    normal.font.color.rgb=RGBColor(0,0,0)
    normal.paragraph_format.line_spacing=1.08;normal.paragraph_format.space_after=Pt(5)
    for name,size in [('Title',18),('Heading 1',14),('Heading 2',12)]:
        s=doc.styles[name];s.font.name='Times New Roman';s.font.size=Pt(size);s.font.color.rgb=RGBColor(0,0,0)
        s.paragraph_format.space_after=Pt(7);s.paragraph_format.keep_with_next=True
    # Remove inherited decorative title rules and theme font overrides.
    for s in doc.styles:
        for border in list(s.element.iter(qn('w:pBdr'))):border.getparent().remove(border)
        if s.name in ['Normal','Title','Heading 1','Heading 2']:
            fonts=s.element.get_or_add_rPr().get_or_add_rFonts()
            for key in ['asciiTheme','hAnsiTheme','eastAsiaTheme','cstheme','csTheme']:
                fonts.attrib.pop(qn('w:'+key),None)
            fonts.set(qn('w:ascii'),'Times New Roman');fonts.set(qn('w:hAnsi'),'Times New Roman')
    footer=section.footer.paragraphs[0];footer.alignment=WD_ALIGN_PARAGRAPH.CENTER
    run=footer.add_run();field=OxmlElement('w:fldSimple');field.set(qn('w:instr'),'PAGE');run._r.addnext(field)
    doc.core_properties.author='Zain Salamor; Parsa Yousefi; Hector Taylor; Abigail Calder; Gregor Hasler'
    doc.core_properties.title='Whole brain and workspace information dynamics following prefrontal theta burst stimulation'
    return doc


def table(doc,headers,rows,widths=None):
    t=doc.add_table(rows=1,cols=len(headers));t.autofit=False
    for i,h in enumerate(headers):t.rows[0].cells[i].text=str(h)
    for row in rows:
        cells=t.add_row().cells
        for cell,value in zip(cells,row):cell.text=str(value)
    repeat=OxmlElement('w:tblHeader');repeat.set(qn('w:val'),'true');t.rows[0]._tr.get_or_add_trPr().append(repeat)
    for row_index,row in enumerate(t.rows):
        cant=OxmlElement('w:cantSplit');row._tr.get_or_add_trPr().append(cant)
        for i,cell in enumerate(row.cells):
            if widths:cell.width=Inches(widths[i])
            props=cell._tc.get_or_add_tcPr()
            margin=OxmlElement('w:tcMar')
            for edge in ['top','left','bottom','right']:
                e=OxmlElement('w:'+edge);e.set(qn('w:w'),'65');e.set(qn('w:type'),'dxa');margin.append(e)
            props.append(margin)
            borders=OxmlElement('w:tcBorders')
            for edge in ['top','left','bottom','right']:
                e=OxmlElement('w:'+edge);e.set(qn('w:val'),'single');e.set(qn('w:sz'),'4');e.set(qn('w:color'),'CCCCCC');borders.append(e)
            props.append(borders)
            if row_index==0:
                shade=OxmlElement('w:shd');shade.set(qn('w:fill'),'EEEEEE');props.append(shade)
            for p in cell.paragraphs:
                p.paragraph_format.space_after=Pt(2);p.paragraph_format.line_spacing=1.0
                for r in p.runs:r.font.size=Pt(9.5);r.bold=(row_index==0)
    doc.add_paragraph().paragraph_format.space_after=Pt(0)
    return t


def main_table(doc,filename,normalized=False):
    with (FOLDER/filename).open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f))
    contents=[]
    for r in rows:
        nd=3 if normalized else 6
        contents.append([r['scope'],r['measure'],f"{float(r['estimate']):.{nd}f}",
          f"{float(r['ci_low']):.{nd}f} to {float(r['ci_high']):.{nd}f}",
          f"{float(r['p_maxT_family']):.4f}",f"{float(r['p_maxT_whole_and_workspace10']):.4f}"])
    table(doc,['Scope','Measure','Change','95% CI','Max-t p5','Max-t p10'],contents,[.95,1.25,.72,1.65,.72,.72])


def build(source,target,supplement=False,mapping=None):
    doc=document();used=list(dict.fromkeys(re.findall(r'@([a-z]+\d+)',source)))
    hold=FOLDER/'ANATOMICAL_REVIEW_PENDING.txt'
    if hold.exists():
        p=doc.add_paragraph('INTERNAL DRAFT — Aachen anatomical reprocessing is under review. Numerical results must be revalidated before submission.')
        for r in p.runs:r.font.color.rgb=RGBColor.from_string('7F2A1D');r.bold=True
    source=citations(source)
    for block in source.split('\n\n'):
        block=block.strip()
        if not block:continue
        if block=='{{FIGURE1}}':doc.add_picture(str(FOLDER/'figures/Figure_1_participant_changes.png'),width=Inches(6.65))
        elif block=='{{FIGURE2}}':
            doc.add_page_break();doc.add_picture(str(FOLDER/'figures/Figure_2_spatial_scope.png'),width=Inches(6.65))
        elif block=='{{TABLE1}}':main_table(doc,'table_1_scope_results.csv')
        elif block=='{{TABLE2}}':main_table(doc,'table_2_scope_results.csv',True)
        elif block=='{{SUPPLEMENT_TABLE}}':table(doc,['Table','Contents'],[(a,b) for a,b,_,_ in mapping],[1.,5.5])
        elif block=='{{SUPPLEMENT_FIGURE}}':doc.add_picture(str(FOLDER/'figures/Figure_S1_leave_one_out.png'),width=Inches(6.65))
        elif block=='{{REFERENCES}}':
            for k in sorted(used,key=lambda k:numbers[k]):
                p=doc.add_paragraph(f'{numbers[k]}. '+reference_text(refs[k]))
                p.paragraph_format.space_after=Pt(5);p.paragraph_format.line_spacing=1.0
                for r in p.runs:r.font.size=Pt(10.5)
        elif block.startswith('# '):doc.add_paragraph(block[2:],style='Title')
        elif block.startswith('## '):
            heading=block[3:]
            if heading in ['Introduction','Figures and tables','References','Supplementary figure','Supplementary tables and participant outputs']:doc.add_page_break()
            doc.add_heading(heading,level=1)
        elif block.startswith('### '):doc.add_heading(block[4:],level=2)
        else:
            if block.startswith(('Table 1.','Table 2.')):doc.add_page_break()
            p=doc.add_paragraph(block)
            if block.startswith(('Figure ','Table ')):
                for r in p.runs:r.font.size=Pt(10.5)
                if block.startswith('Table '):p.paragraph_format.keep_with_next=True
            if block.startswith('[TODO:'):
                for r in p.runs:r.font.color.rgb=RGBColor.from_string('7F2A1D')
    doc.save(target)
    assert '{{' not in '\n'.join(p.text for p in doc.paragraphs)
    return target


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--suffix',default='')
    args=parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9_-]*',args.suffix):raise ValueError('Invalid document suffix')
    mapping=prepare_data()
    main_path=build(main_text,FOLDER/('Aachen_Brief_Research_Report'+args.suffix+'.docx'))
    supp_path=build(supp_text,FOLDER/('Supplementary_Material'+args.suffix+'.docx'),True,mapping)
    main_keys=list(dict.fromkeys(re.findall(r'@([a-z]+\d+)',main_text)))
    resolved=citations(main_text).replace('{{REFERENCES}}','\n\n'.join(f'{numbers[k]}. '+reference_text(refs[k]) for k in main_keys))
    for number in [1,2]:
        figure='Figure_1_participant_changes.png' if number==1 else 'Figure_2_spatial_scope.png'
        resolved=resolved.replace('{{FIGURE'+str(number)+'}}',f'![Figure {number}](figures/{figure})')
        with (FOLDER/f'table_{number}_scope_results.csv').open(encoding='utf-8',newline='') as f:
            records=list(csv.DictReader(f))
        columns=['scope','measure','estimate','ci_low','ci_high','p_maxT_family','p_maxT_whole_and_workspace10']
        lines=['| '+' | '.join(columns)+' |','| '+' | '.join(['---']*len(columns))+' |']
        lines+=['| '+' | '.join(r[c] for c in columns)+' |' for r in records]
        resolved=resolved.replace('{{TABLE'+str(number)+'}}','\n'.join(lines))
    if (FOLDER/'ANATOMICAL_REVIEW_PENDING.txt').exists():
        resolved='**INTERNAL DRAFT: Aachen anatomical reprocessing is under review; results require revalidation before submission.**\n\n'+resolved
    assert '{{' not in resolved
    (FOLDER/'manuscript.md').write_text(resolved,encoding='utf-8')
    bib=[]
    for k,r in refs.items():
        bib.append('@article{'+k+',\n  author = {'+' and '.join(r['bibauthors'])+'},\n  title = {{'+r['title']+'}},\n  journal = {'+r['journal']+'},\n  year = {'+str(r['year'])+'},\n  volume = {'+r['volume']+'},\n  pages = {'+r['pages']+'},\n  doi = {'+r['doi']+'}\n}')
    (FOLDER/'references.bib').write_text('\n\n'.join(bib)+'\n',encoding='utf-8')
    (FOLDER/'reference_audit.json').write_text(json.dumps({'keys':keys,'references':refs},ensure_ascii=False,indent=2),encoding='utf-8')
    with zipfile.ZipFile(FOLDER/'Aachen_Supplementary_Data.zip','w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(DATA.iterdir()):
            if p.is_file():z.write(p,'Aachen_Supplementary_Data/'+p.name)
    count=lambda t:len(re.findall(r"\b[\w]+(?:[’'-][\w]+)*\b",t))
    abstract_words=count(main_text.split('## Abstract',1)[1].split('## Introduction',1)[0])
    body_words=count(main_text.split('## Introduction',1)[1].split('## Data availability',1)[0])
    hold_note='Aachen anatomical reprocessing is pending. Do not submit or treat the current numerical claims as final until the corrected pipeline has been verified. The apparent target-mask mismatch was traced to forced extraction of already skull-stripped anatomy; see PROVENANCE/AACHEN_ANATOMICAL_REPROCESSING_20261125.md.\n\n' if (FOLDER/'ANATOMICAL_REVIEW_PENDING.txt').exists() else ''
    notes='''# Submission notes

The scientific draft now includes the fixed-HCP-workspace reanalysis. The older narrative is retained in archive/ only. The current document is Aachen-only; the separate TMS project audit also reports SUDMEX.

The author order is Zain Salamor, Parsa Yousefi, Hector Taylor, Abigail Calder, Gregor Hasler. Zain and Parsa share first authorship; Parsa is corresponding author; Gregor is PI and last author.

Before submission, complete Zain's affiliation, all authors' CRediT contributions and approvals, funding/grant details, conflict declarations, and the persistent public code/materials repository URL. These fields are explicitly marked and have not been invented. Confirm the secondary-analysis ethics wording with the responsible institution if needed.

The manuscript has two figures and two tables, an ABSTRACT_COUNT-word abstract and approximately BODY_COUNT words from Introduction through Discussion. Figures are available in PNG and vector PDF. The supplement and data ZIP retain all previously reported numerical outputs, subject to the anatomical reprocessing review. The earlier supplementary_data/ folder contains historical preparation inputs; the current ZIP is mirrored in submission_data/.

Suggested article type: Brief Research Report. Intended Research Topic: Psychedelic Medicine and Neuromodulation From Bench to Bedside. Check the live submission page and requirements before submission. The files have not been submitted or publicly uploaded.

The prior candidate claim is a whole-brain pairwise Φ change with limited workspace evidence. Signal-representation diagnostics are summarized briefly in the main text and documented fully in Supplementary Results. Do not relabel the result increased synergy, increased consciousness, therapeutic efficacy, or a direct replication of Luppi/Yousefi. Spatial follow-up evidence outside the workspace is reported alongside uncorrected and five-measure-corrected direct contrasts. The component identity suggests a relative shift in temporal predictability, not a demonstrated increase in total joint information. Target-based Aachen analyses are separate and suspended pending anatomical validation; SUDMEX is documented separately.
'''
    notes=notes.replace('# Submission notes\n\n','# Submission notes\n\n'+hold_note,1).replace('ABSTRACT_COUNT',str(abstract_words)).replace('BODY_COUNT',f'{body_words:,}')
    (FOLDER/'submission_notes.md').write_text(notes,encoding='utf-8')
    print('Created:',main_path,supp_path,sep='\n')
    print('References:',len(keys))


if __name__=='__main__':main()
