"""Copy the frozen GWT psychedelic pipeline into the independent TMS_GWT project.

What it does
------------
1. Copies the production code, configs, atlases, tools, frozen psychedelic results and the
   reference analysis notebooks from D:/projects/GWT into D:/projects/TMS_GWT, keeping the
   same relative layout (scripts resolve their project root relative to __file__, so the
   layout is what makes the copy run unchanged).
2. Hashes every source and destination file (SHA-256) and writes PROVENANCE/MANIFEST.tsv.
3. Rewrites absolute path literals that point into the GWT tree (D:/projects/GWT,
   D:\\projects\\GWT, /d/projects/GWT, /mnt/d/projects/GWT) to the TMS_GWT root, in code and
   config files only (.py .m .sh .ps1 .yaml .yml .json). Notebooks, Markdown and every
   results/output file stay byte-identical: they are the record of the psychedelic runs.
   Every rewritten line is logged in PROVENANCE/PATCHES_paths.tsv.

The source tree is opened read-only; nothing is ever written under D:/projects/GWT.
Re-running is safe: identical files are skipped, a conflicting destination file stops the run.

Usage:  python setup/copy_from_gwt.py [--dry-run]
"""
from __future__ import annotations

import argparse
import csv
import fnmatch
import hashlib
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path

SRC = Path(r"D:\projects\GWT")
DST = Path(r"D:\projects\TMS_GWT")
PROV = DST / "PROVENANCE"

# (source path relative to SRC, include globs or None for everything, exclude globs).
# Globs are matched against the path relative to the item root, with forward slashes.
COMMON_EXCLUDE = ["*__pycache__*", "*.pyc", "*.ipynb_checkpoints*"]
ITEMS = [
    # --- production code, configs, atlases, external tools shipped in the repo
    ("pipeline/scripts", None, []),
    ("pipeline/configs", None, []),
    ("pipeline/atlases", None, []),
    ("pipeline/tools/integrated-info-decomp", None, []),
    ("pipeline/tools/npy-matlab", None, []),
    ("pipeline/tools/rsHRF", None, []),
    ("pipeline/tools/spm12", None, []),
    # Luppi/Mediano MATLAB PhiIDFull.m, required by step04_compute_phiid_matlab.m
    ("pipeline/docs/luppi_supps/software/MMI_PhiID", None, []),
    ("pipeline/docs/Luppi_2023.md", None, []),
    ("pipeline/docs/luppi_2022.md", None, []),
    # --- reference analysis notebook and its documentation
    ("pipeline/Global_Workspace_Complete_Pipeline_publication.ipynb", None, []),
    ("pipeline/Global_Workspace_Complete_Pipeline_publication.html", None, []),
    ("pipeline/README_publication_notebook.md", None, []),
    ("pipeline/requirements_publication.txt", None, []),
    ("pipeline/CLAUDE.md", None, []),
    # --- frozen psychedelic and HCP results (inputs to the TMS comparison; never recomputed)
    ("pipeline/results/phiid_matrices", None, []),
    ("pipeline/results/group_analysis", None, []),
    ("pipeline/results/qc", None, []),
    ("output/publication-notebook-2026-09-26", None, []),
    # --- tier-1 robustness: code, notebook and final tables (not its 4 GB of copied inputs)
    ("tier1_robustness/README.md", None, []),
    ("tier1_robustness/RESULTS_SUMMARY.md", None, []),
    ("tier1_robustness/INPUT_MANIFEST.csv", None, []),
    ("tier1_robustness/requirements.txt", None, []),
    ("tier1_robustness/notebook_outputs", None, []),
    ("tier1_robustness/pipeline/Global_Workspace_Complete_Pipeline_tier1.ipynb", None, []),
    ("tier1_robustness/pipeline/Global_Workspace_Complete_Pipeline_tier1.html", None, []),
    ("tier1_robustness/pipeline/scripts", None, []),
    ("tier1_robustness/pipeline/atlases", None, []),
    ("tier1_robustness/pipeline/plot_assets", None, []),
    ("tier1_robustness/tier1/code", None, []),
    ("tier1_robustness/tier1/source_notebook", None, []),
    # --- NuMIT-PhiID workstream: code, notes and result tables only (no npy/npz arrays)
    ("numit-phiid", ["*.py", "*.m", "*.md", "*.csv", "*.txt", "*.json", "*.yaml", "*.yml"], []),
    # --- FreeSurfer licence used by fMRIPrep
    ("license.txt", None, []),
]

EMPTY_DIRS = [
    "tmp/templateflow",          # TemplateFlow cache mounted by run_fmriprep.ps1
    "tmp/fmriprep_work",         # fMRIPrep work dirs (delete after each run)
    "pipeline/logs",
    "pipeline/results/timeseries",
]

# Path rewriting: only these suffixes, only under these copied roots.
PATCH_SUFFIXES = {".py", ".m", ".sh", ".ps1", ".yaml", ".yml", ".json"}
PATCH_ROOTS = ["pipeline/scripts", "pipeline/configs", "tier1_robustness/tier1/code",
               "tier1_robustness/pipeline/scripts", "numit-phiid"]
WIN_RE = re.compile(r"(?i)\b(D:)([\\/]+)(projects)([\\/]+)(GWT)(?![A-Za-z0-9_])")
POSIX_RE = re.compile(r"(?i)((?:/mnt)?/d)(/projects/)(GWT)(?![A-Za-z0-9_])")


def rewrite_paths(text: str) -> str:
    text = WIN_RE.sub(lambda m: f"{m.group(1)}{m.group(2)}{m.group(3)}{m.group(4)}TMS_GWT", text)
    return POSIX_RE.sub(lambda m: f"{m.group(1)}{m.group(2)}TMS_GWT", text)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def matches(rel: str, globs) -> bool:
    return any(fnmatch.fnmatch(rel, g) or fnmatch.fnmatch(Path(rel).name, g) for g in globs)


def iter_files(src_root: Path, include, exclude):
    if src_root.is_file():
        yield src_root, src_root.name
        return
    for p in sorted(src_root.rglob("*")):
        if not p.is_file():
            continue
        rel = p.relative_to(src_root).as_posix()
        if matches(rel, COMMON_EXCLUDE + exclude):
            continue
        if include is not None and not matches(rel, include):
            continue
        yield p, rel


def in_patch_scope(dst_rel: str) -> bool:
    return Path(dst_rel).suffix.lower() in PATCH_SUFFIXES and any(
        dst_rel == r or dst_rel.startswith(r + "/") for r in PATCH_ROOTS)


def decode(raw: bytes):
    try:
        return raw.decode("utf-8"), "utf-8"
    except UnicodeDecodeError:
        return raw.decode("latin-1"), "latin-1"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not SRC.is_dir():
        sys.exit(f"Source not found: {SRC}")
    PROV.mkdir(parents=True, exist_ok=True)
    manifest_rows, patch_rows = [], []
    n_copied = n_skipped = 0
    total_bytes = 0

    for item, include, exclude in ITEMS:
        src_root = SRC / item
        if not src_root.exists():
            sys.exit(f"Missing source item: {src_root}")
        for src_file, rel in iter_files(src_root, include, exclude):
            dst_rel = item if src_root.is_file() else f"{item}/{rel}"
            dst_file = DST / dst_rel
            size = src_file.stat().st_size
            total_bytes += size
            if args.dry_run:
                continue
            src_hash = sha256(src_file)
            raw = src_file.read_bytes() if in_patch_scope(dst_rel) else None
            new_raw = None
            if raw is not None:
                text, enc = decode(raw)
                new_text = rewrite_paths(text)
                if new_text != text:
                    new_raw = new_text.encode(enc)
                    old_lines = text.splitlines()
                    new_lines = new_text.splitlines()
                    for i, (a, b) in enumerate(zip(old_lines, new_lines), start=1):
                        if a != b:
                            patch_rows.append([dst_rel, i, a.strip(), b.strip()])
            if dst_file.exists():
                dst_hash = sha256(dst_file)
                expected = hashlib.sha256(new_raw).hexdigest() if new_raw is not None else src_hash
                if dst_hash != expected:
                    sys.exit(f"Conflict: {dst_file} exists with different content; not overwriting.")
                n_skipped += 1
            else:
                dst_file.parent.mkdir(parents=True, exist_ok=True)
                if new_raw is not None:
                    dst_file.write_bytes(new_raw)
                    shutil.copystat(src_file, dst_file)
                else:
                    shutil.copy2(src_file, dst_file)
                dst_hash = sha256(dst_file)
                n_copied += 1
            manifest_rows.append([dst_rel, f"{item}" if src_root.is_file() else f"{item}/{rel}",
                                  size, src_hash, dst_hash,
                                  "path-rewritten" if new_raw is not None else "identical"])

    if args.dry_run:
        print(f"DRY RUN: {total_bytes / 1e9:.2f} GB would be copied")
        return 0

    for d in EMPTY_DIRS:
        (DST / d).mkdir(parents=True, exist_ok=True)

    with open(PROV / "MANIFEST.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["path_in_TMS_GWT", "source_path_in_GWT", "bytes", "sha256_source",
                    "sha256_copy", "status"])
        w.writerows(manifest_rows)
    with open(PROV / "PATCHES_paths.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["file", "line", "before", "after"])
        w.writerows(patch_rows)

    n_patched = sum(1 for r in manifest_rows if r[5] == "path-rewritten")
    summary = (f"{datetime.now().isoformat(timespec='seconds')}  files: {len(manifest_rows)} "
               f"(copied {n_copied}, already present {n_skipped}); bytes: {total_bytes / 1e9:.2f} GB; "
               f"path-rewritten files: {n_patched}; rewritten lines: {len(patch_rows)}")
    (PROV / "COPY_SUMMARY.txt").write_text(summary + "\n", encoding="utf-8")
    print(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
