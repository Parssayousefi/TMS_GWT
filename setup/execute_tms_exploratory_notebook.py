"""Execute only the lightweight TMS exploratory-results notebook."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys

import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager

ROOT=Path(__file__).resolve().parents[1]
runtime=ROOT/'tmp/jupyter-notebook/exploratory-runtime'
runtime.mkdir(parents=True,exist_ok=True)
os.environ['JUPYTER_RUNTIME_DIR']=str(runtime)
path=ROOT/'output/jupyter-notebook/TMS_Exploratory_Findings.ipynb'
notebook=nbformat.read(path,as_version=4)
name='gwt' if 'gwt' in KernelSpecManager().find_kernel_specs() else 'python3'
manager=KernelManager(kernel_name=name)
manager.kernel_spec.argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}']
client=NotebookClient(notebook,timeout=120,startup_timeout=60,kernel_name=name,km=manager,
                      resources={'metadata':{'path':str(ROOT)}},allow_errors=False)
try:
    print('Executing the exploratory results notebook with the gwt kernel.',flush=True)
    client.execute()
finally:
    nbformat.write(notebook,path)
    if manager.has_kernel: manager.shutdown_kernel(now=True)
code=[cell for cell in notebook.cells if cell.cell_type=='code']
assert all(cell.execution_count is not None for cell in code)
errors=[o for cell in code for o in cell.get('outputs',[]) if o.output_type=='error']
assert not errors
nbformat.validate(notebook)
record={'status':'executed_without_errors','created_utc':datetime.now(timezone.utc).isoformat(),
        'code_cells':len(code),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'notebook':str(path),'kernel_python':sys.executable}
(ROOT/'output/tms_exploratory_report/20261125/notebook_execution.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2),flush=True)
