"""Check that the TMS datasets are registered correctly in the copied pipeline.

Read-only. For each TMS dataset it
  1. loads the config through the pipeline's own fmriprep_config_loader (as run_fmriprep.ps1 does),
  2. lists the subjects fMRIPrep would process and the extra fMRIPrep arguments,
  3. queries the raw BIDS tree with pybids the way fMRIPrep does, applying bids_filter.json when
     the config uses one, and reports T1w and rest-BOLD inputs per subject,
  4. checks that the sessions on disk match the conditions in the config.

Usage:  python setup/check_tms_registration.py
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import yaml
from bids import BIDSLayout

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "pipeline" / "scripts" / "preprocessing"))
import fmriprep_config_loader as fl  # noqa: E402

DATASETS = ["sudmex_tms", "aachen_itbs"]


def check(name: str) -> list[str]:
    problems = []
    configs = fl.load_configs(name, ROOT / "pipeline")
    paths = fl.get_paths(configs, "local")
    params = fl.get_fmriprep_params(configs)
    subjects = fl.discover_subjects(configs, paths["bids_root"])
    cfg = configs["dataset"]
    print(f"\n=== {name}: {cfg['name']}")
    print(f"bids_root {paths['bids_root']}  ->  fmriprep_dir {paths['fmriprep_dir']}")
    print(f"license {paths['license_file']} exists={Path(paths['license_file']).is_file()}")
    print(f"fMRIPrep {params.get('version')}  spaces {params.get('output_spaces')}  "
          f"extra args: {params.get('extra_args') or '(none)'}")
    print(f"subjects to process: {len(subjects)}")

    queries = {"t1w": {"datatype": "anat", "suffix": "T1w"},
               "bold": {"datatype": "func", "suffix": "bold"}}
    if "--bids-filter-file" in (cfg.get("fmriprep_extra_args") or ""):
        filt = json.loads((Path(paths["bids_root"]) / "bids_filter.json").read_text())
        for key in queries:
            queries[key].update(filt.get(key, {}))
        print(f"bids_filter.json applied: {filt}")

    layout = BIDSLayout(paths["bids_root"], validate=False)
    conditions = set(cfg["tasks"])
    n_t1 = Counter()
    bold_sessions = Counter()
    for s in subjects:
        t1 = layout.get(subject=s, extension=[".nii", ".nii.gz"], return_type="filename", **queries["t1w"])
        bold = layout.get(subject=s, extension=[".nii", ".nii.gz"], **queries["bold"])
        n_t1[len(t1)] += 1
        if len(t1) != 1:
            problems.append(f"{name} sub-{s}: {len(t1)} T1w inputs (expected 1)")
        if not bold:
            problems.append(f"{name} sub-{s}: no BOLD inputs")
        for b in bold:
            ses = b.entities.get("session")
            bold_sessions[ses] += 1
            if ses not in conditions:
                problems.append(f"{name} sub-{s}: session '{ses}' is not a condition in the config")
    print(f"T1w inputs per subject: {dict(n_t1)}")
    print(f"rest BOLD runs per session: {dict(sorted(bold_sessions.items()))}  "
          f"(total {sum(bold_sessions.values())})")
    qc_sessions = set(cfg.get("qc", {}).get("include_sessions") or [])
    if qc_sessions and qc_sessions != conditions:
        problems.append(f"{name}: qc.include_sessions {sorted(qc_sessions)} != conditions {sorted(conditions)}")
    return problems


def main() -> int:
    registry = yaml.safe_load((ROOT / "pipeline" / "configs" / "dataset_configs.yaml").read_text())
    problems = [f"{d} not in dataset_configs.yaml" for d in DATASETS if d not in registry["datasets"]]
    for d in DATASETS:
        problems += check(d)
    print(f"\nProblems: {len(problems)}")
    for p in problems[:40]:
        print("  -", p)
    if len(problems) > 40:
        print(f"  ... {len(problems) - 40} more")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
