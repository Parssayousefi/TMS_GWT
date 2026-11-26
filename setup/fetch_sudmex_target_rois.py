"""Fetch the published SUDMEX localization workbook and mask archive."""
from pathlib import Path, PurePosixPath
from datetime import datetime, timezone
import hashlib
import json
from urllib.parse import quote
from urllib.request import urlopen
import zipfile

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/sudmex_target_rois'
FILES={'14_DLPFC localization information.xlsx':'6edba5c2655187515c8472a6ae84c61b',
       'SUDMEX_TMS_masks.zip':'53544e0a33eaecf7208f150b0478ceb8'}


def main():
    OUT.mkdir(exist_ok=True);rows=[]
    for name,digest in FILES.items():
        path=OUT/name;url='https://zenodo.org/records/10409461/files/'+quote(name)+'?download=1'
        if not path.exists():
            with urlopen(url,timeout=45) as r: content=r.read()
            assert hashlib.md5(content).hexdigest()==digest
            path.write_bytes(content)
        assert hashlib.md5(path.read_bytes()).hexdigest()==digest
        rows.append({'file':name,'url':url,'md5':digest,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    with zipfile.ZipFile(OUT/'SUDMEX_TMS_masks.zip') as z:
        names=z.namelist();print('Archive entries:',len(names))
        for name in names[:14]:print(name)
        for info in z.infolist():
            rel=PurePosixPath(info.filename)
            if rel.is_absolute() or '..' in rel.parts:raise ValueError('Unsafe archive path')
            dest=(OUT/'extracted'/Path(*rel.parts)).resolve()
            assert str(dest).startswith(str(OUT.resolve())+'\\')
            if info.is_dir():dest.mkdir(parents=True,exist_ok=True)
            else:
                dest.parent.mkdir(parents=True,exist_ok=True)
                dest.write_bytes(z.read(info.filename))
    (OUT/'download_manifest.json').write_text(json.dumps({'created_utc':datetime.now(timezone.utc).isoformat(),'files':rows},indent=2),encoding='utf-8')


if __name__=='__main__':main()
