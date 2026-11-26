"""Record deliberate edits made to copied GWT files after the copy (dataset registration).

For each file in EDITS it rebuilds the copy-time content (GWT source + the documented path
rule), writes a unified diff against the current file to PROVENANCE/post_copy_edits.diff and
the hashes to PROVENANCE/POST_COPY_EDITS.tsv. setup/verify_copy.py accepts exactly these
recorded edits and nothing else.

Usage:  python setup/record_post_copy_edits.py
"""
from __future__ import annotations

import csv
import difflib
import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("copy_from_gwt", HERE / "copy_from_gwt.py")
cfg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cfg)

# path in TMS_GWT -> (date, reason)
EDITS = {
    "pipeline/configs/dataset_configs.yaml":
        ("2026-11-05", "Registered the TMS datasets sudmex_tms and aachen_itbs (2 registry entries added)"),
    "pipeline/scripts/preprocessing/step02_hrf_deconvolution_matlab.m":
        ("2026-11-05", "Added TR 2.0 s for sudmex_tms and aachen_itbs to the dataset_tr map, as done for each new GWT cohort"),
}


def main() -> int:
    manifest = {r["path_in_TMS_GWT"]: r for r in
                csv.DictReader(open(cfg.PROV / "MANIFEST.tsv", encoding="utf-8"), delimiter="\t")}
    rows, diffs = [], []
    for rel, (date, reason) in EDITS.items():
        m = manifest[rel]
        src_text, enc = cfg.decode((cfg.SRC / m["source_path_in_GWT"]).read_bytes())
        at_copy = cfg.rewrite_paths(src_text) if m["status"] == "path-rewritten" else src_text
        now_path = cfg.DST / rel
        now_text, _ = cfg.decode(now_path.read_bytes())
        diffs += difflib.unified_diff(at_copy.splitlines(keepends=True), now_text.splitlines(keepends=True),
                                      fromfile=f"{rel} (as copied)", tofile=f"{rel} (now)")
        rows.append([rel, date, reason, m["sha256_copy"], cfg.sha256(now_path)])
    with open(cfg.PROV / "POST_COPY_EDITS.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["path_in_TMS_GWT", "date", "reason", "sha256_as_copied", "sha256_now"])
        w.writerows(rows)
    (cfg.PROV / "post_copy_edits.diff").write_text("".join(diffs), encoding="utf-8")
    print(f"Recorded {len(rows)} edited files; diff lines: {len(diffs)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
