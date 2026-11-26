"""Retrieve public left/right DLPFC masks without altering the BIDS or pipeline trees."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import urllib.request
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'data/aachen_target_rois'
NS = {'s': 'http://s3.amazonaws.com/doc/2006-03-01/'}
LIST_URL = 'https://s3.amazonaws.com/openneuro.org?list-type=2&prefix=ds001832/derivatives/'


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(LIST_URL, timeout=30) as response:
        content = response.read()
    tree = ET.fromstring(content)
    assert tree.findtext('s:IsTruncated', namespaces=NS) == 'false'
    (DEST / 'public_listing.xml').write_bytes(content)
    manifest = []
    for item in tree.findall('s:Contents', NS):
        key = item.findtext('s:Key', namespaces=NS)
        if not re.search(r'rFC[lr].*\.nii$', key):
            continue
        local = DEST / Path(key).name
        url = 'https://s3.amazonaws.com/openneuro.org/' + key
        if not local.exists():
            with urllib.request.urlopen(url, timeout=30) as response:
                local.write_bytes(response.read())
        assert local.stat().st_size == int(item.findtext('s:Size', namespaces=NS))
        manifest.append({'key': key, 'url': url, 'file': local.name,
                         'etag': item.findtext('s:ETag', namespaces=NS),
                         'sha256': hashlib.sha256(local.read_bytes()).hexdigest()})
    report = {'retrieved_utc': datetime.now(timezone.utc).isoformat(),
              'source': LIST_URL, 'dataset': 'ds001832',
              'note': 'Public S3 release; snapshot identified by object ETags and SHA256. Session numbers are original BIDS numbering.',
              'files': manifest}
    (DEST / 'download_manifest.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print('Downloaded/verified', len(manifest), 'public DLPFC masks.')
    print('Subjects:', sorted({re.search(r'sub-\d+', r['file'])[0] for r in manifest}))


if __name__ == '__main__':
    main()
