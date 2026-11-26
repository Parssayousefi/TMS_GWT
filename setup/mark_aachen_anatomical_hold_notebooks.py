"""Label current notebooks without modifying their historical numerical outputs."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import shutil

import nbformat


ROOT = Path(__file__).resolve().parents[1]
MARKER = "aachen_anatomical_review_hold_20261125"
NOTICE = """> **Aachen results on hold — anatomical QC update, 25 November 2026.**
>
> The Aachen estimates below were calculated from the original preprocessing, which repeated skull stripping on already stripped T1 images and removed substantial outer cortex. These are historical results and must not be used as current evidence for a TMS effect. Corrected processing of all 60 scans has completed through PhiID in an isolated branch; the estimates below require statistical recomputation before the hold can be lifted. The SUDMEX results are not affected by this particular error, although their separate audit identifies regional functional-coverage limitations.
>
> See `PROVENANCE/AACHEN_ANATOMICAL_REPROCESSING_20261125.md` and `docs/TMS_CONTROL_WORKSPACE_QC_20261125.md` for the audit and control-workspace QC. This notice leaves the original cells and numerical outputs intact.
"""


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    archive = ROOT / "output/aachen_alignment_review/notebook_hold_archive"
    archive.mkdir(parents=True, exist_ok=True)
    records = []
    for name in ["TMS_Global_Workspace_Publication.ipynb",
                 "TMS_Fixed_Workspace_Reanalysis.ipynb", "TMS_Exploratory_Findings.ipynb"]:
        path = ROOT / "output/jupyter-notebook" / name
        notebook = nbformat.read(path, as_version=4)
        if any(cell.metadata.get(MARKER) for cell in notebook.cells):
            continue
        backup = archive / name
        assert not backup.exists(), "Inspect the existing archive before overwriting"
        before_hash = sha256(path)
        before_cells = json.dumps(notebook.cells, sort_keys=True)
        shutil.copy2(path, backup)
        assert sha256(backup) == before_hash
        notice = nbformat.v4.new_markdown_cell(NOTICE, metadata={MARKER: True})
        notebook.cells.insert(0, notice)
        nbformat.validate(notebook)
        nbformat.write(notebook, path)
        saved = nbformat.read(path, as_version=4)
        assert json.dumps(saved.cells[1:], sort_keys=True) == before_cells
        records.append({"notebook": str(path), "archive": str(backup),
                        "before_sha256": before_hash, "after_sha256": sha256(path),
                        "historical_cells_and_outputs_unchanged": True})
    receipt = {"updated_utc": datetime.now(timezone.utc).isoformat(), "records": records,
               "validation": "Notebook schema valid; original cells and outputs unchanged. No historical analysis rerun."}
    if records:
        (archive / "receipt.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
