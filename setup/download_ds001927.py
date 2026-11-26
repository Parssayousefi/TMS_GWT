"""Download the pinned public ds001927 release and verify every file.

Regular metadata are checked against Git blob IDs; annex files against the
release's embedded content digest and size. Supports interrupted transfers.
Downloaded study code is retained as data and is never executed.
"""
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import hashlib
import json
import re
import shutil
import threading
import time

import nibabel as nib
import numpy as np
import requests

ROOT = Path(__file__).resolve().parents[1]
DATASET = "ds001927"
VERSION = "2.1.0"
COMMIT = "3239a60164f7a230e02b6a42eb36a13ff4f97f57"
TARGET = ROOT / "data" / DATASET
TARGET_RESOLVED = TARGET.resolve()
AUDIT = ROOT / "output/additional_dataset_audit_20261126"
OUT = ROOT / "output/ds001927_download_20261126"
LOCK = threading.Lock()
STATE = {"completed_files": 0, "verified_bytes": 0, "received_bytes_this_run": 0,
         "active": [], "failed": [], "started_utc": datetime.now(timezone.utc).isoformat()}


def destination(relative):
    assert not Path(relative).is_absolute() and ".." not in Path(relative).parts, relative
    # Windows can transiently fail canonicalization while another worker
    # creates a shared parent. Recheck; never bypass the containment guard.
    for attempt in range(4):
        path = (TARGET / relative).resolve()
        if path.is_relative_to(TARGET_RESOLVED):
            return path
        time.sleep(.05)
    raise AssertionError(f"Path containment failed: {path!s} under {TARGET_RESOLVED!s}")


def git_hash(content):
    return hashlib.sha1(b"blob " + str(len(content)).encode() + b"\0" + content).hexdigest()


def write_status():
    tmp = OUT / "status.json.tmp"
    tmp.write_text(json.dumps(STATE, indent=2), encoding="utf-8")
    tmp.replace(OUT / "status.json")


def bounded_get(url, limit=5_000_000):
    for attempt in range(4):
        try:
            with requests.get(url, timeout=(15, 45), stream=True) as response:
                response.raise_for_status()
                content = response.raw.read(limit + 1, decode_content=True)
            assert len(content) <= limit, url
            return content
        except requests.RequestException:
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def prepare_one(entry):
    relative = entry["path"]
    destination(relative)  # Reject unsafe paths before any write.
    cache = AUDIT / (DATASET + "_" + relative.replace("/", "_"))
    content = cache.read_bytes() if cache.exists() else None
    if content is None or git_hash(content) != entry["sha"]:
        content = bounded_get(f"https://raw.githubusercontent.com/OpenNeuroDatasets/{DATASET}/{COMMIT}/{relative}")
    assert git_hash(content) == entry["sha"], (relative, "Git blob hash mismatch")
    item = {"path": relative, "git_blob": entry["sha"], "git_mode": entry["mode"]}
    if entry["mode"] == "120000":
        pointer = content.decode("utf-8").strip()
        match = re.search(r"/(MD5E|SHA256E)-s(\d+)--([a-fA-F0-9]+)(?:\.[^/]*)?$", pointer)
        assert match, (relative, "Unsupported annex digest")
        item.update(kind="annex", digest_algorithm="md5" if match[1] == "MD5E" else "sha256",
                    bytes=int(match[2]), expected_digest=match[3].lower(),
                    url=f"https://s3.amazonaws.com/openneuro.org/{DATASET}/{relative}")
    else:
        assert entry["mode"] in ["100644", "100755"]
        item.update(kind="metadata", bytes=len(content), expected_digest=hashlib.sha256(content).hexdigest())
        metadata = OUT / "metadata" / relative
        metadata.parent.mkdir(parents=True, exist_ok=True)
        metadata.write_bytes(content)
    return item


def digest_file(path, algorithm):
    digest, sha = hashlib.new(algorithm), hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024**2), b""):
            digest.update(chunk); sha.update(chunk)
    return digest.hexdigest(), sha.hexdigest()


def download_one(item):
    relative, size = item["path"], item["bytes"]
    path = destination(relative)
    path.parent.mkdir(parents=True, exist_ok=True)
    algorithm = item.get("digest_algorithm", "sha256")
    if path.exists():
        digest, sha = digest_file(path, algorithm)
        assert path.stat().st_size == size and digest == item["expected_digest"], (relative, "Existing file differs; preserved for review")
    elif item["kind"] == "metadata":
        content = (OUT / "metadata" / relative).read_bytes()
        assert git_hash(content) == item["git_blob"]
        path.write_bytes(content)
        sha = hashlib.sha256(content).hexdigest()
    else:
        partial = path.with_name(path.name + ".part")
        with LOCK:
            STATE["active"].append(relative)
        try:
            for attempt in range(5):
                have = partial.stat().st_size if partial.exists() else 0
                assert have <= size, (relative, "Oversized partial file")
                if have == size:
                    break
                try:
                    headers = {"Range": f"bytes={have}-"} if have else {}
                    with requests.get(item["url"], headers=headers, timeout=(15, 60), stream=True) as response:
                        response.raise_for_status()
                        if have:
                            assert response.status_code == 206 and response.headers.get("Content-Range", "").startswith(f"bytes {have}-"), (relative, "Resume range not honored")
                        else:
                            assert response.status_code in [200, 206]
                        with partial.open("ab" if have else "wb") as handle:
                            for chunk in response.iter_content(chunk_size=4 * 1024**2):
                                if not chunk:
                                    continue
                                have += len(chunk)
                                assert have <= size, (relative, "Download exceeds pinned release size")
                                handle.write(chunk)
                                with LOCK:
                                    STATE["received_bytes_this_run"] += len(chunk)
                    assert partial.stat().st_size == size, (relative, "Incomplete object")
                    break
                except (requests.RequestException, AssertionError):
                    if attempt == 4:
                        raise
                    time.sleep(min(2 ** attempt, 8))
            digest, sha = digest_file(partial, algorithm)
            assert partial.stat().st_size == size and digest == item["expected_digest"], (relative, "Annex content checksum mismatch; partial preserved")
            partial.replace(path)
        finally:
            with LOCK:
                STATE["active"].remove(relative)
    with LOCK:
        STATE["completed_files"] += 1
        STATE["verified_bytes"] += size
    return {**item, "sha256": sha, "verified": True}


def monitor(stop):
    while not stop.wait(15):
        with LOCK:
            write_status()
            print(f"Verified {STATE['completed_files']}/{STATE['total_files']} files; "
                  f"{STATE['verified_bytes']/1e9:.2f}/{STATE['total_bytes']/1e9:.2f} GB; "
                  f"received {STATE['received_bytes_this_run']/1e9:.2f} GB this run", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    assert 1 <= args.workers <= 128
    OUT.mkdir(parents=True, exist_ok=True)
    TARGET.mkdir(parents=True, exist_ok=True)
    tree = json.loads((AUDIT / f"{DATASET}_tree.json").read_text())
    assert tree["sha"] == COMMIT and not tree.get("truncated")
    entries = [x for x in tree["tree"] if x["type"] == "blob"]
    manifest = OUT / "download_manifest.json"
    if manifest.exists():
        record = json.loads(manifest.read_text())
        assert record["commit"] == COMMIT
        items = record["files"]
        assert {v["path"]: v["git_blob"] for v in items} == {v["path"]: v["sha"] for v in entries}
    else:
        print(f"Preparing pinned manifest for {len(entries)} files", flush=True)
        with ThreadPoolExecutor(max_workers=4) as executor:
            items = list(executor.map(prepare_one, entries))
        manifest.write_text(json.dumps({"dataset": DATASET, "version": VERSION, "commit": COMMIT, "files": items}, indent=2), encoding="utf-8")
    total = sum(v["bytes"] for v in items)
    remaining = sum(v["bytes"] for v in items if not destination(v["path"]).exists())
    free = shutil.disk_usage(TARGET).free
    assert free > remaining + 5 * 1024**3, "Insufficient free disk space for download and working margin"
    STATE.update(total_files=len(items), total_bytes=total, dataset=DATASET, version=VERSION,
                 commit=COMMIT, target=str(TARGET), free_bytes_at_start=free, state="downloading")
    write_status()
    print(f"Downloading {len(items)} files, {total/1e9:.3f} GB to {TARGET}", flush=True)
    stop = threading.Event()
    monitor_thread = threading.Thread(target=monitor, args=(stop,), daemon=True)
    monitor_thread.start()
    results, failures = [], []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        jobs = {executor.submit(download_one, item): item for item in items}
        for future in as_completed(jobs):
            item = jobs[future]
            try:
                results.append(future.result())
            except Exception as exc:
                failures.append({"path": item["path"], "error": repr(exc)})
                with LOCK:
                    STATE["failed"] = list(failures)
                print("FAILED:", item["path"], repr(exc), flush=True)
    stop.set(); monitor_thread.join()
    (OUT / "verified_files.json").write_text(json.dumps(sorted(results, key=lambda v: v["path"]), indent=2), encoding="utf-8")
    STATE["failed"] = failures
    if failures:
        STATE["state"] = "incomplete"
        write_status()
        raise RuntimeError(f"{len(failures)} downloads failed; valid completed files and partials retained for resume")
    headers = []
    for item in items:
        if item["path"].endswith((".nii.gz", ".nii")):
            image = nib.load(destination(item["path"]))
            assert np.isfinite(image.affine).all() and abs(np.linalg.det(image.affine[:3, :3])) > 0
            if "_bold.nii" in item["path"]:
                assert image.shape == (64, 64, 40, 600), (item["path"], image.shape)
                assert abs(image.header.get_zooms()[3] - 1.25) < 1e-6
            headers.append({"path": item["path"], "shape": list(image.shape), "zooms": [float(v) for v in image.header.get_zooms()]})
    assert len(headers) == 207
    assert STATE["completed_files"] == len(items) and STATE["verified_bytes"] == total
    assert not list(TARGET.rglob("*.part"))
    (OUT / "nifti_header_checks.json").write_text(json.dumps(headers, indent=2), encoding="utf-8")
    STATE.update(state="complete", completed_utc=datetime.now(timezone.utc).isoformat(),
                 nifti_headers_verified=len(headers), rest_runs=138, t1w_images=69,
                 participants=23, sessions_per_participant=3,
                 content_checks="Pinned Git blobs and annex digests; SHA-256 recorded for every local file",
                 quality_limit="Header/content integrity checks only; no motion, registration or signal-coverage QC and no preprocessing")
    write_status()
    print(json.dumps(STATE, indent=2), flush=True)


if __name__ == "__main__":
    main()
