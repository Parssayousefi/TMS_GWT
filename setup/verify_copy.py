"""Verify that TMS_GWT holds exact copies of the frozen GWT pipeline.

Checks, for every row of PROVENANCE/MANIFEST.tsv:
  1. the copy exists and still has the recorded SHA-256;
  2. the GWT source still has its recorded SHA-256 (the original was not touched);
  3. "identical" rows: copy hash == source hash;
     "path-rewritten" rows: applying the documented path rule to the source reproduces the
     copy byte for byte, so the only differences are GWT path literals.
Then scans every code/config file in the rewrite scope for leftover GWT path literals.
Writes PROVENANCE/VERIFY_REPORT.txt and exits non-zero on any failure.

Usage:  python setup/verify_copy.py
"""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import sys
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("copy_from_gwt", HERE / "copy_from_gwt.py")
cfg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cfg)

SRC, DST, PROV = cfg.SRC, cfg.DST, cfg.PROV


def load_post_copy_edits() -> dict:
    """Deliberate edits recorded by setup/record_post_copy_edits.py (path -> row)."""
    p = PROV / "POST_COPY_EDITS.tsv"
    if not p.is_file():
        return {}
    return {r["path_in_TMS_GWT"]: r for r in csv.DictReader(open(p, encoding="utf-8"), delimiter="\t")}


def main() -> int:
    rows = list(csv.DictReader(open(PROV / "MANIFEST.tsv", encoding="utf-8"), delimiter="\t"))
    edits = load_post_copy_edits()
    failures = []
    n_ident = n_rewritten = n_edited = 0
    for r in rows:
        dst = DST / r["path_in_TMS_GWT"]
        src = SRC / r["source_path_in_GWT"]
        if not dst.is_file():
            failures.append(f"missing copy: {dst}")
            continue
        edit = edits.get(r["path_in_TMS_GWT"])
        if edit is not None:
            # Recorded post-copy edit: the as-copied hash must match the manifest and the
            # file must still be exactly the recorded edited version.
            n_edited += 1
            if edit["sha256_as_copied"] != r["sha256_copy"]:
                failures.append(f"post-copy edit record does not match the manifest: {dst}")
            if cfg.sha256(dst) != edit["sha256_now"]:
                failures.append(f"edited file changed after its edit was recorded: {dst}")
        elif cfg.sha256(dst) != r["sha256_copy"]:
            failures.append(f"copy changed since copying: {dst}")
        if not src.is_file() or cfg.sha256(src) != r["sha256_source"]:
            failures.append(f"GWT source missing or changed: {src}")
        if r["status"] == "identical":
            n_ident += 1
            if r["sha256_copy"] != r["sha256_source"]:
                failures.append(f"identical row with different hashes: {dst}")
        else:
            n_rewritten += 1
            # The as-copied file must equal the GWT source with only the path rule applied.
            # Compared by hash against the manifest, so recorded post-copy edits do not mask it.
            text, enc = cfg.decode(src.read_bytes())
            if hashlib.sha256(cfg.rewrite_paths(text).encode(enc)).hexdigest() != r["sha256_copy"]:
                failures.append(f"rewritten file differs beyond the path rule: {dst}")

    leftovers = []
    for root in cfg.PATCH_ROOTS:
        base = DST / root
        files = [base] if base.is_file() else base.rglob("*")
        for p in files:
            if p.is_file() and p.suffix.lower() in cfg.PATCH_SUFFIXES:
                text, _ = cfg.decode(p.read_bytes())
                for i, line in enumerate(text.splitlines(), start=1):
                    if cfg.WIN_RE.search(line) or cfg.POSIX_RE.search(line):
                        leftovers.append(f"{p.relative_to(DST).as_posix()}:{i}: {line.strip()}")
    failures += [f"GWT path left in code/config: {x}" for x in leftovers]

    lines = [f"Verified {datetime.now().isoformat(timespec='seconds')}",
             f"Manifest rows: {len(rows)} ({n_ident} identical, {n_rewritten} path-rewritten; "
             f"{n_edited} with a recorded post-copy edit, see POST_COPY_EDITS.tsv)",
             f"Failures: {len(failures)}"] + failures
    (PROV / "VERIFY_REPORT.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines[:3] + failures[:20]))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
