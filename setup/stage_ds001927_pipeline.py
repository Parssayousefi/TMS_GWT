"""Prepare a ds001927-only copy; leave production and frozen cohorts untouched.

Only run paths, the CONN dataset setting, and the rsHRF TR map change in copied
processing code. No stage is executed by this preparation script.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import shutil
import yaml

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "pipeline"
DEST = ROOT / "reruns/ds001927_20261126/pipeline"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    DEST.resolve().relative_to(ROOT.resolve())
    if DEST.exists():
        raise RuntimeError("Isolated pipeline exists; verify it instead of overwriting.")
    DEST.mkdir(parents=True)
    for name in ("scripts", "atlases", "tools/spm12", "tools/rsHRF", "tools/npy-matlab",
                 "docs/luppi_supps/software/MMI_PhiID"):
        shutil.copytree(SOURCE / name, DEST / name, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    for name in ("configs/datasets", "results/timeseries", "results/phiid_matrices",
                 "results/qc/ds001927", "data/conn_input/ds001927", "results/timeseries_intermediate"):
        (DEST / name).mkdir(parents=True, exist_ok=True)
    edits = []
    names = ("scripts/preprocessing/step01b_conn_denoise.m",
             "scripts/preprocessing/step02_hrf_deconvolution_matlab.m",
             "scripts/analysis/step04_compute_phiid_matlab.m")
    for name in names:
        before = (SOURCE / name).read_bytes()
        old, new = str(SOURCE).encode(), str(DEST).encode()
        assert old in before
        after = before.replace(old, new)
        replacements = []
        if name.endswith("step01b_conn_denoise.m"):
            replacements = [(b"DATASET   = 'missouri_psilo'", b"DATASET   = 'ds001927'")]
        elif name.endswith("step02_hrf_deconvolution_matlab.m"):
            replacements = [(b"'missouri_psilo', 'sudmex_tms', 'aachen_itbs'}, ...",
                             b"'missouri_psilo', 'sudmex_tms', 'aachen_itbs', 'ds001927'}, ..."),
                            (b"0.91, 1.761, 2.0, 2.0});", b"0.91, 1.761, 2.0, 2.0, 1.25});")]
        for a, b in replacements:
            assert after.count(a) == 1, (name, a)
            after = after.replace(a, b)
        reconstructed = after
        for a, b in reversed(replacements):
            reconstructed = reconstructed.replace(b, a)
        assert reconstructed.replace(new, old) == before
        (DEST / name).write_bytes(after)
        edits.append(dict(file=name, source_sha256=digest(before), staged_sha256=digest(after),
                          change="Run directory; ds001927 CONN selection / TR1.25 lookup only where applicable"))

    config = dict(dataset_id="ds001927", name="Castrillon three-site rTMS (ds001927)",
                  paths={"local": {"bids_root": (ROOT / "data/ds001927").as_posix(),
                                   "fmriprep_dir": (ROOT / "data/ds001927/derivatives/fmriprep_23_2_1").as_posix()}},
                  bids_structure={"has_sessions": True},
                  subjects={"auto_discover": True, "include": [], "exclude": []},
                  tasks={"restpre": {"bids_task_name": "restpre", "is_control": True,
                                     "description": "Pre-stimulation; select stimulation site separately"},
                         "restpost": {"bids_task_name": "restpost", "is_control": False,
                                      "description": "Post-stimulation; select stimulation site separately"}},
                  run_selection={"strategy": "first"},
                  file_patterns={"bold": "sub-{subject}_ses-{session}_task-{task}_space-MNI152NLin2009cAsym_res-2_desc-preproc_bold.nii.gz",
                                 "mask": "sub-{subject}_ses-{session}_task-{task}_space-MNI152NLin2009cAsym_res-2_desc-brain_mask.nii.gz",
                                 "confounds": "sub-{subject}_ses-{session}_task-{task}_desc-confounds_timeseries.tsv"},
                  qc={"motion_thresholds": {"max_translation_mm": 3.0, "max_rotation_deg": 3.15},
                      "auto_truncate": False, "include_sessions": ["ctr", "fro", "occ"],
                      "exclusions_file": "results/qc/ds001927/exclusions.csv",
                      "required_external_gate": "setup/qc_ds001927_derivatives.py plus anatomical/coverage release"},
                  preprocessing={"tr": 1.25, "volumes_to_remove": 8, "volumes_to_remove_last": 0,
                                 "bandpass_low": 0.008, "bandpass_high": 0.09,
                                 "compcor_wm": 5, "compcor_csf": 5, "global_signal_regression": False},
                  fmriprep_extra_args="--skull-strip-t1w force --longitudinal --ignore fieldmaps slicetiming",
                  notes="Use reviewed subject-specific launcher/filter. sub-s022 ses-ctr T1 excluded only. Six scans per participant. Never pool sites using generic group task loader. Nonsteady-state extensions and actual motion decisions require explicit QC adapter before CONN preparation.")
    config["paths"].update(config["paths"]["local"])
    (DEST / "configs/datasets/ds001927.yaml").write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    registry = yaml.safe_load((SOURCE / "configs/dataset_configs.yaml").read_text())
    registry["datasets"] = {"ds001927": {"name": config["name"], "config_file": "configs/datasets/ds001927.yaml",
                                         "description": "Three stimulation sites by pre/post; isolated processing only"}}
    registry["default_dataset"] = "ds001927"
    registry["global_paths"] = {"local": {**registry["global_paths"]["local"],
                                         "project_dir": DEST.as_posix(),
                                         "work_dir": (ROOT / "tmp/fmriprep_ds001927_20261126").as_posix(),
                                         "log_dir": (DEST.parent / "logs").as_posix()}}
    registry["global_preprocessing"].update(config["preprocessing"])
    registry["fmriprep"]["skip_bids_validation"] = False
    (DEST / "configs/dataset_configs.yaml").write_text(yaml.safe_dump(registry, sort_keys=False), encoding="utf-8")
    unchanged = []
    for path in (DEST / "scripts").rglob("*"):
        if path.is_file() and path.suffix in (".py", ".m", ".ps1", ".sh"):
            rel = path.relative_to(DEST).as_posix()
            if rel not in names:
                assert path.read_bytes() == (SOURCE / rel).read_bytes(), rel
                unchanged.append(rel)
    receipt = dict(created_utc=datetime.now(timezone.utc).isoformat(), source=str(SOURCE), stage=str(DEST),
                   script_edits=edits, unchanged_scripts=unchanged,
                   status="prepared only; anatomical and explicit motion release required before downstream processing",
                   frozen_cohort_data_copied=False, processing_started=False,
                   registry_datasets=["ds001927"], qc_policy="PROVENANCE/ds001927_preprocessing_20261126/qc_policy.json")
    (DEST.parent / "stage_manifest.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    (ROOT / "PROVENANCE/ds001927_preprocessing_20261126/isolated_pipeline.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    print(f"Prepared {DEST}; {len(unchanged)} scripts unchanged; {len(edits)} settings-only script copies")


if __name__ == "__main__":
    main()
