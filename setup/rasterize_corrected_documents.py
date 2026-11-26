from pathlib import Path
import json
import hashlib
from pdf2image import convert_from_path
from pypdf import PdfReader
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'manuscript/aachen_brief_report/corrected_20261126/qa'
POPPLER=Path(r'C:\Program Files\WindowsApps\OpenAI.CodexPrimaryRuntime.v26-1123-915-0_26.1123.915.0_x64__3k8sg7r9htsxt\dependencies\native\poppler\Library\bin')
receipt={'renderer':'Native Word PDF export plus bundled Poppler; packaged render_docx failed because bundled LibreOffice is unavailable on Windows.'}
for kind in ['main','supplement']:
    folder=BASE/kind;reader=PdfReader(folder/'document.pdf')
    pages=convert_from_path(str(folder/'document.pdf'),dpi=120,poppler_path=str(POPPLER))
    checks=[]
    for i,(page,pdfpage) in enumerate(zip(pages,reader.pages),1):
        page.save(folder/f'page-{i}.png');text=pdfpage.extract_text()
        assert '\ufffd' not in text and '{{' not in text and '[@' not in text
        checks.append({'page':i,'characters':len(text),'first':text[:75],'last':text[-75:]})
    receipt[kind]={'pages':len(pages),'page_checks':checks};print(kind,len(pages),'pages')
    for check in checks:print(check)
(BASE/'render_receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False),encoding='utf-8')
