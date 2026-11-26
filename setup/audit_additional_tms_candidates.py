"""Public metadata and 64-KiB access probes; no full imaging download."""
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
import collections
import csv
import hashlib
import io
import json
import re
import struct
import xml.etree.ElementTree as ET
import zlib
import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/additional_dataset_audit_20261126"
OUT.mkdir(parents=True, exist_ok=True)
SOURCES = {}


def fetch(name, url):
    with requests.get(url, timeout=40, stream=True) as response:
        response.raise_for_status()
        content = response.raw.read(20_000_001, decode_content=True)
    assert len(content) <= 20_000_000, "Metadata exceeds 20 MB"
    (OUT / name).write_bytes(content)
    SOURCES[name] = {"url": url, "sha256": hashlib.sha256(content).hexdigest(), "bytes": len(content)}
    return content


def probe(dataset, path):
    url = f"https://s3.amazonaws.com/openneuro.org/{dataset}/{path}"
    with requests.get(url, headers={"Range": "bytes=0-65535"}, timeout=40, stream=True) as response:
        response.raise_for_status()
        header = response.raw.read(65536)
        record = {"url": url, "http_status": response.status_code,
                  "content_range": response.headers.get("Content-Range"),
                  "bytes_read": len(header), "sha256_partial": hashlib.sha256(header).hexdigest()}
    assert response.status_code == 206 and record["content_range"], "Server did not honor range request"
    (OUT / f"{dataset}_sample_first64KiB.bin").write_bytes(header)
    nifti = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(header, 540) if path.endswith(".gz") else header
    endian = "<" if struct.unpack("<i", nifti[:4])[0] == 348 else ">"
    assert struct.unpack(endian + "i", nifti[:4])[0] == 348, "Expected NIfTI-1 header"
    record["dimensions"] = struct.unpack(endian + "8h", nifti[40:56])
    record["pixdim"] = struct.unpack(endian + "8f", nifti[76:108])
    record["header_only_not_full_integrity_check"] = True
    return record


def audit(dataset, ref):
    tree = json.loads(fetch(f"{dataset}_tree.json", f"https://api.github.com/repos/OpenNeuroDatasets/{dataset}/git/trees/{ref}?recursive=1"))
    assert not tree.get("truncated")
    commit = tree["sha"]
    paths = [x["path"] for x in tree["tree"] if x["type"] == "blob"]
    base = f"https://raw.githubusercontent.com/OpenNeuroDatasets/{dataset}/{commit}/"
    root = [p for p in paths if "/" not in p and not p.startswith(".")]
    bodies = {}
    for p in root:
        if p.endswith((".json", ".tsv", ".md")) or p in ["README", "CHANGES"]:
            bodies[p] = fetch(f"{dataset}_{p}", base + p).decode("utf-8-sig")
    rest = [p for p in paths if "rest" in p.lower() and p.endswith(("_bold.nii.gz", "_bold.nii"))]
    anat = [p for p in paths if p.endswith(("_T1w.nii.gz", "_T1w.nii"))]
    metadata = [p for p in paths if p.endswith(".json") and "rest" in p.lower()]
    def sidecar(p):
        return p, json.loads(fetch(dataset + "_" + p.replace("/", "_"), base + p))
    with ThreadPoolExecutor(max_workers=4) as executor:
        sidecars = dict(executor.map(sidecar, metadata))
    values = collections.defaultdict(set)
    for p, entry in sidecars.items():
        for k in ["RepetitionTime", "EchoTime", "TaskName", "ProtocolName", "SeriesDescription"]:
            values[k].add(str(entry.get(k)))
    sessions = collections.defaultdict(list)
    for p in rest:
        sessions["/".join(p.split("/")[:2])].append(p)
    participants = list(csv.DictReader(io.StringIO(bodies.get("participants.tsv", "")), delimiter="\t"))
    sizes = {}
    # S3 lists actual imaging object lengths; git-annex tree 'size' is only a
    # small pointer length and must not be reported as imaging download size.
    try:
        marker = ""
        for page in range(20):
            suffix = "&continuation-token=" + requests.utils.quote(marker, safe="") if marker else ""
            payload = fetch(f"{dataset}_s3_listing_{page:02d}.xml", f"https://s3.amazonaws.com/openneuro.org?list-type=2&prefix={dataset}/&max-keys=1000" + suffix)
            xml = ET.fromstring(payload)
            ns = {"s": "http://s3.amazonaws.com/doc/2006-03-01/"}
            for item in xml.findall("s:Contents", ns):
                key = item.find("s:Key", ns).text
                sizes[key[len(dataset)+1:]] = int(item.find("s:Size", ns).text)
            if xml.find("s:IsTruncated", ns).text != "true":
                break
            marker = xml.find("s:NextContinuationToken", ns).text
        else:
            raise RuntimeError("S3 listing exceeded metadata page cap")
    except Exception as exc:
        sizes = {}
        size_error = str(exc)
    summary = {"dataset": dataset, "ref": ref, "commit": commit, "root_files": root,
               "participant_rows": len(participants), "participant_columns": list(participants[0]) if participants else [],
               "rest_files": len(rest), "rest_subjects": len({p.split('/')[0] for p in rest}),
               "rest_by_session": dict(collections.Counter(p.split('/')[1] for p in rest)),
               "rest_session_scan_counts": {s: len(v) for s, v in sessions.items()},
               "t1w_files": len(anat), "sidecar_values": {k: sorted(v) for k,v in values.items()},
               "possible_condition_metadata": [p for p in paths if any(w in p.lower() for w in ['participants', 'sessions', 'events', 'condition', 'target', 'stim', 'readme']) and not p.endswith(('.nii.gz', '.nii'))],
               "license": json.loads(bodies.get('dataset_description.json', '{}')).get('License'),
               "sample_access": probe(dataset, rest[0]) if rest else None,
               "rest_plus_t1_bytes": sum(sizes[p] for p in rest+anat) if sizes and all(p in sizes for p in rest+anat) else None,
               "size_error": locals().get('size_error'), "full_recordings_downloaded": False}
    if dataset == "ds001927":
        summary["rest_filenames"] = [Path(p).name for p in rest[:12]]
    manifest = []
    for p in rest + anat:
        manifest.append({"dataset": dataset, "path": p, "bytes": sizes.get(p), "public_url": f"https://s3.amazonaws.com/openneuro.org/{dataset}/{p}", "selected_for_imaging_download": False})
    with (OUT / f"{dataset}_rest_t1_manifest.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(manifest[0]))
        writer.writeheader(); writer.writerows(manifest)
    (OUT / f"{dataset}_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


if __name__ == "__main__":
    results = []
    for dataset, ref in [("ds001927", "2.1.0"), ("ds008893", "main"), ("ds008064", "main")]:
        try:
            result = audit(dataset, ref)
            results.append(result)
            print(json.dumps({k: result[k] for k in ['dataset', 'commit', 'participant_rows', 'participant_columns', 'rest_files', 'rest_subjects', 'rest_by_session', 'rest_plus_t1_bytes', 'sample_access']}), flush=True)
        except Exception as exc:
            results.append({"dataset": dataset, "error": str(exc)})
            print(dataset, "FAILED", str(exc), flush=True)
    (OUT / "source_manifest.json").write_text(json.dumps({"checked_utc": datetime.now(timezone.utc).isoformat(), "sources": SOURCES}, indent=2), encoding="utf-8")
    (OUT / "summary.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
