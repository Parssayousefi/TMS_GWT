"""Check the staged publication snapshot without printing secret candidates."""
from pathlib import Path
import subprocess
import re
import json
import zipfile
ROOT=Path(__file__).resolve().parents[1]
paths=subprocess.check_output(['git','ls-files','--cached','-z'],cwd=ROOT).decode().split('\0')
paths=[p for p in paths if p]
assert paths, 'No staged files; refusing to report a successful empty snapshot audit.'
entries=subprocess.check_output(['git','ls-files','--stage','-z'],cwd=ROOT).decode().split('\0')
index={entry.split('\t',1)[1]:entry.split('\t',1)[0].split() for entry in entries if entry}
for relative in paths:
    mode,object_id,stage=index[relative]
    assert mode in ('100644','100755') and stage=='0', 'Unsupported staged entry: '+relative
    # Check the committed bytes, not only the working files; this also detects
    # unintentional line-ending conversion of hash-bound release sources.
    actual=subprocess.check_output(['git','hash-object','--no-filters',relative],cwd=ROOT).decode().strip()
    if relative!='PROVENANCE/REPOSITORY_SNAPSHOT_AUDIT_20261126.json':
        assert actual==object_id, 'Staged and working bytes differ: '+relative
forbidden=[];oversized=[];secret_files=[];binary=[]
patterns=[re.compile(rb'gh[pousr]_[A-Za-z0-9]{30,}'),re.compile(rb'github_pat_[A-Za-z0-9_]{30,}'),re.compile(rb'AKIA[0-9A-Z]{16}'),re.compile(rb'-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----')]
safe_numpy={'pipeline/results/group_analysis/hcp/rest/hcp_rest_workspace_mask.npy','pipeline/results/group_analysis/hcp/rest/hcp_rest_synergy_redundancy_gradient.npy'}
for relative in paths:
    p=ROOT/relative
    if relative.startswith(('data/','logs/','reruns/','tmp/','outputs/','tier1_robustness/')) or relative=='license.txt':forbidden.append(relative)
    if p.suffix.lower() in ['.mat','.nii','.gz','.npy','.npz','.h5','.pkl'] and relative not in safe_numpy:forbidden.append(relative)
    if p.stat().st_size>25*1024**2:oversized.append(relative)
    if p.suffix.lower() in ['.zip','.docx']:
        with zipfile.ZipFile(p) as z:
            for name in z.namelist():
                if p.suffix.lower()=='.zip' and Path(name).suffix.lower() not in ['.csv','.json','.txt','.md']:forbidden.append(relative+'::'+name)
                if Path(name).suffix.lower() in ['.xml','.json','.txt','.csv','.md'] and any(q.search(z.read(name)) for q in patterns):secret_files.append(relative+'::'+name)
    elif p.suffix.lower() not in ['.png','.jpg','.jpeg','.npy']:
        if any(q.search(p.read_bytes()) for q in patterns):secret_files.append(relative)
assert not forbidden, 'Unexpected data or excluded files: '+repr(forbidden)
assert not oversized, 'Oversized files: '+repr(oversized)
assert not secret_files, 'Credential-like patterns in files (content not printed): '+repr(secret_files)
record={'status':'passed','staged_files':len(paths),'total_bytes':sum((ROOT/p).stat().st_size for p in paths),
        'raw_imaging_and_credentials_excluded':True,'small_frozen_group_reference_files':sorted(safe_numpy),
        'public_dataset_scalar_zip_reviewed':True}
print(json.dumps(record,indent=2))
(ROOT/'PROVENANCE/REPOSITORY_SNAPSHOT_AUDIT_20261126.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
