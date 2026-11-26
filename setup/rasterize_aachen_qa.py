"""Rasterize the Word-exported PDFs with the bundled Poppler runtime."""
from pathlib import Path
import argparse
import json
import re
from pdf2image import convert_from_path
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'manuscript/aachen_brief_report/qa'
POPPLER=Path(r'C:\Program Files\WindowsApps\OpenAI.CodexPrimaryRuntime.v26-1123-915-0_26.1123.915.0_x64__3k8sg7r9htsxt\dependencies\native\poppler\Library\bin')
parser=argparse.ArgumentParser()
parser.add_argument('--kind',choices=['all','main','supplement'],default='all')
args=parser.parse_args()
receipt_path=BASE/'render_receipt.json'
receipt=json.loads(receipt_path.read_text(encoding='utf-8')) if receipt_path.exists() else {}
receipt.pop('renderer',None)
for kind in (['main','supplement'] if args.kind=='all' else [args.kind]):
    folder=BASE/kind
    reader=PdfReader(folder/'document.pdf')
    pages=convert_from_path(str(folder/'document.pdf'),dpi=140,poppler_path=str(POPPLER))
    texts=[]
    for i,(page,pdfpage) in enumerate(zip(pages,reader.pages),1):
        page.save(folder/f'page-{i}.png')
        text=pdfpage.extract_text()
        texts.append({'page':i,'characters':len(text),'first':text[:110],'last':text[-110:]})
        assert '\ufffd' not in text and '{{' not in text and '[@' not in text
    receipt[kind]={'pages':len(pages),'text':texts}
(BASE/'render_receipt.json').write_text(json.dumps({'renderer':'Word native PDF export plus bundled Poppler; packaged render_docx failed because soffice is absent',**receipt},indent=2,ensure_ascii=False),encoding='utf-8')
for k,v in receipt.items():
    print(k,v['pages'],'pages')
    for p in v['text']:print(p['page'],p['characters'],re.sub(r'\s+',' ',p['first']), '| END:',re.sub(r'\s+',' ',p['last']))
