"""Copy verified pilot derivatives to a dedicated, larger-volume output tree.

Never removes the source or overwrites different destination bytes. This is an
operational copy, not a preprocessing change. The destination is deliberately
fixed and the copy is restartable only when existing files match exactly.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import shutil

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/ds001927/derivatives/fmriprep_23_2_1"
STORAGE = Path("E:/data/TMS_GWT/ds001927_20261126")
DEST = STORAGE / "derivatives/fmriprep_23_2_1"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for data in iter(lambda: stream.read(8 * 1024**2), b""):
            h.update(data)
    return h.hexdigest()


def main():
    if not Path("E:/").is_dir() or shutil.disk_usage("E:/").free < 400 * 1024**3:
        raise RuntimeError("Dedicated E: storage needs at least 400 GiB free")
    files = sorted(path for path in SOURCE.rglob("*") if path.is_file())
    if not files or any(p.is_symlink() for p in SOURCE.rglob("*")):
        raise RuntimeError("Expected ordinary, nonempty derivative files")
    for subject in ("sub-s001", "sub-s022"):
        if len(list((SOURCE / subject).glob("ses-*/func/*desc-preproc_bold.nii.gz"))) != 6:
            raise RuntimeError("Both six-run pilots must be complete before staging")
    DEST.mkdir(parents=True, exist_ok=True)
    manifest = []
    for i, source in enumerate(files):
        relative = source.relative_to(SOURCE)
        target = DEST / relative
        target.resolve().relative_to(STORAGE.resolve())
        target.parent.mkdir(parents=True, exist_ok=True)
        digest = sha(source)
        if target.exists():
            if target.stat().st_size != source.stat().st_size or sha(target) != digest:
                raise RuntimeError(f"Different destination file preserved: {relative}")
        else:
            # An interrupted copy is recognizable; do not silently overwrite it.
            temporary = target.with_name(target.name + ".copying")
            if temporary.exists():
                raise RuntimeError(f"Incomplete copy requires inspection: {temporary}")
            shutil.copy2(source, temporary)
            if sha(temporary) != digest:
                raise RuntimeError(f"Copy checksum mismatch: {relative}")
            temporary.rename(target)
        manifest.append(dict(relative_path=relative.as_posix(), bytes=source.stat().st_size, sha256=digest))
        if i % 50 == 0:
            print(f"Verified {i+1}/{len(files)} pilot derivative files", flush=True)
    out = ROOT / "output/ds001927_storage_20261126"
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = out / "copied_pilot_file_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    record = dict(created_utc=datetime.now(timezone.utc).isoformat(), status="copied_and_all_hashes_verified",
                  source=str(SOURCE), derivatives=str(DEST), storage_root=str(STORAGE),
                  work=str(STORAGE / "work"), file_count=len(manifest),
                  bytes=sum(x["bytes"] for x in manifest), source_deleted=False,
                  input_image_bytes_changed=False, manifest_sha256=sha(manifest_path),
                  script_sha256=sha(Path(__file__)))
    (ROOT / "PROVENANCE/ds001927_preprocessing_20261126/storage.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
    print(json.dumps(record, indent=2), flush=True)


if __name__ == "__main__":
    main()
