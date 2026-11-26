"""Download and verify the published SUDMEX-TMS clinical source workbooks."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "data/sudmex_tms/clinical/source/zenodo-10409461"
FILES = {
    "5_CCQ-G.xlsx": "e0bac99afe01379a36eacc958a568f1b",
    "5_CCQ-N.xlsx": "aaf95692874edc5399dd8b5aa0ae4b25",
    "7_BIS-11.xlsx": "745cada170678b664b3c4839d045674e",
    "12_VAS Clinico.xlsx": "ec685444ae5f6f9fde46d0c88fbe0daa",
    "Demographic_data.xlsx": "10ff721407d9bb97f1b2c5f8532d3f2b",
}


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    records = []
    for name, expected_md5 in FILES.items():
        path = DEST / name
        url = "https://zenodo.org/records/10409461/files/" + quote(name) + "?download=1"
        if path.exists():
            content = path.read_bytes()
        else:
            print(f"Downloading {name}", flush=True)
            request = Request(url, headers={"User-Agent": "TMS-GWT-research/1.0"})
            with urlopen(request, timeout=45) as response:
                content = response.read()
        digest = hashlib.md5(content).hexdigest()
        if digest != expected_md5 or not content.startswith(b"PK"):
            raise RuntimeError(f"Source checksum or XLSX format mismatch: {name}: {digest}")
        if not path.exists():
            path.write_bytes(content)
        records.append({"file": name, "url": url, "bytes": len(content),
                        "md5": digest, "sha256": hashlib.sha256(content).hexdigest()})
        print(f"Verified {name}: {len(content)} bytes", flush=True)
    (DEST / "download_manifest.json").write_text(json.dumps({
        "record": "https://zenodo.org/records/10409461", "version": "4.2",
        "verified_utc": datetime.now(timezone.utc).isoformat(), "files": records,
    }, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
