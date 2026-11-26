# Provenance of the copied GWT pipeline

Copied from `D:\projects\GWT` on 5 November 2026 by `setup/copy_from_gwt.py`; checked by `setup/verify_copy.py`.

| File | Contents |
| --- | --- |
| `MANIFEST.tsv` | One row per copied file: path here, path in GWT, bytes, SHA-256 of source and copy, status (`identical` or `path-rewritten`) |
| `PATCHES_paths.tsv` | Every line changed by the path rule: file, line number, before, after |
| `COPY_SUMMARY.txt` | File, byte and rewrite counts of the copy run |
| `VERIFY_REPORT.txt` | Result of `setup/verify_copy.py` |
| `POST_COPY_EDITS.tsv`, `post_copy_edits.diff` | Deliberate edits to copied files after the copy (TMS dataset registration), with hashes before and after; written by `setup/record_post_copy_edits.py` and accepted by `verify_copy.py` |
| `environment_gwt_pip_freeze.txt`, `environment_gwt_conda_export.txt` | Package versions of the `gwt` conda environment used for the psychedelic runs |

## The one permitted change

Absolute path literals pointing into the GWT tree (`D:/projects/GWT`, `D:\projects\GWT`, `/d/projects/GWT`, `/mnt/d/projects/GWT`, any case) were replaced by the same path into `D:\projects\TMS_GWT`. The rule ran only on code and config files (`.py .m .sh .ps1 .yaml .yml .json`) under `pipeline/scripts`, `pipeline/configs`, `tier1_robustness/tier1/code`, `tier1_robustness/pipeline/scripts` and `numit-phiid`. `verify_copy.py` re-applies the rule to each GWT source and requires the result to equal the copy byte for byte.

Everything else is identical to the source, including notebooks, Markdown and all result files. Result metadata (for example the YAML sidecars next to ΦID matrices) still names GWT paths, because it records where those results were produced.

## Machine-wide tools, not copied

- Docker image `nipreps/fmriprep:23.2.1`
- MATLAB R2025a (`C:\Program Files\MATLAB\R2025a`) and CONN 25b (`C:\Users\YousefiP\software\conn25b\conn`)
- conda environment `gwt` (`C:\Miniforge3\envs\gwt`, Python 3.9.23)
