"""Execute and save the TMS notebook using the current Python environment."""
from pathlib import Path
import os
import sys
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
from jupyter_client.kernelspec import KernelSpecManager

ROOT = Path(__file__).resolve().parents[1]
runtime = ROOT / 'tmp/jupyter-notebook/runtime'
runtime.mkdir(parents=True, exist_ok=True)
os.environ['JUPYTER_RUNTIME_DIR'] = str(runtime)
path = ROOT / 'output/jupyter-notebook/TMS_Global_Workspace_Publication.ipynb'
notebook = nbformat.read(path, as_version=4)
kernel_name = 'gwt' if 'gwt' in KernelSpecManager().find_kernel_specs() else 'python3'
manager = KernelManager(kernel_name=kernel_name)
manager.kernel_spec.argv = [sys.executable, '-m', 'ipykernel_launcher', '-f', '{connection_file}']
client = NotebookClient(notebook, timeout=180, startup_timeout=60, kernel_name=kernel_name, km=manager,
                        resources={'metadata': {'path': str(ROOT)}}, allow_errors=False)
try:
    print('Starting gwt kernel and executing notebook cells.', flush=True)
    client.execute()
finally:
    nbformat.write(notebook, path)
    if manager.has_kernel:
        manager.shutdown_kernel(now=True)
errors = [output for cell in notebook.cells for output in cell.get('outputs', []) if output.output_type == 'error']
if errors:
    raise RuntimeError('Notebook contains execution errors')
print(f'Executed {sum(c.cell_type == "code" for c in notebook.cells)} code cells: {path}')
