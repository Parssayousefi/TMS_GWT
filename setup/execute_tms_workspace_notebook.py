"""Execute only the TMS workspace review notebook using the gwt environment."""
from pathlib import Path
import os
import sys
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager

root=Path(__file__).resolve().parents[1]
runtime=root/'tmp/jupyter-notebook/runtime'
runtime.mkdir(parents=True,exist_ok=True)
os.environ['JUPYTER_RUNTIME_DIR']=str(runtime)
path=root/'output/jupyter-notebook/TMS_Fixed_Workspace_Reanalysis.ipynb'
nb=nbformat.read(path,as_version=4)
kernel='gwt' if 'gwt' in KernelSpecManager().find_kernel_specs() else 'python3'
manager=KernelManager(kernel_name=kernel)
manager.kernel_spec.argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}']
client=NotebookClient(nb,timeout=120,startup_timeout=60,kernel_name=kernel,km=manager,
                      resources={'metadata':{'path':str(root)}},allow_errors=False)
try:client.execute()
finally:
    nbformat.write(nb,path)
    if manager.has_kernel:manager.shutdown_kernel(now=True)
assert all(o.output_type!='error' for c in nb.cells for o in c.get('outputs',[]))
run=root/'output/tms_workspace_scope'/(root/'output/tms_workspace_scope/LATEST.txt').read_text().strip()
nbformat.write(nb,run/path.name)
print('Executed',sum(c.cell_type=='code' for c in nb.cells),'cells:',path)
