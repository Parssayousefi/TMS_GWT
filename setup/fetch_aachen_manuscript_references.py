"""Fetch public DOI metadata and the cited Yousefi preprint; no participant data."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import json
import time
import requests

OUT=Path(__file__).resolve().parents[1]/'manuscript/aachen_brief_report/sources'
OUT.mkdir(parents=True,exist_ok=True)
DOIS={
 'blumberger2018':'10.1016/S0140-6736(18)30295-2',
 'alkhasli2019':'10.3389/fnhum.2019.00190',
 'huang2005':'10.1016/j.neuron.2004.12.033',
 'barrett2011':'10.1371/journal.pcbi.1001052',
 'mediano2025':'10.1073/pnas.2423297122',
 'luppi2024':'10.7554/eLife.88173',
 'esteban2019':'10.1038/s41592-018-0235-4',
 'whitfield2012':'10.1089/brain.2012.0073',
 'behzadi2007':'10.1016/j.neuroimage.2007.04.042',
 'schaefer2018':'10.1093/cercor/bhx179',
 'tian2020':'10.1038/s41593-020-00711-6',
 'wu2013':'10.1016/j.media.2013.01.003',
 'liardi2025':'10.1371/journal.pcbi.1013629',
 'winkler2014':'10.1016/j.neuroimage.2014.01.060',
 'power2012':'10.1016/j.neuroimage.2011.10.018',
 'ben2025':'10.1016/j.neurot.2025.e00556',
}

def fetch(item):
 key,doi=item
 target=OUT/(key+'.json')
 try:
  if target.exists():
   data=json.loads(target.read_text(encoding='utf-8'))
  else:
   response=requests.get('https://api.crossref.org/works/'+doi,timeout=30)
   if response.status_code==429:
    time.sleep(1)
    response=requests.get('https://api.crossref.org/works/'+doi,timeout=30)
   response.raise_for_status()
   data=response.json()['message']
   target.write_text(json.dumps(data,indent=2),encoding='utf-8')
  return key, {'doi':doi,'title':data.get('title'), 'year':data.get('published',data.get('issued'))}
 except Exception as exc:
  return key,{'doi':doi,'error':str(exc)}

with ThreadPoolExecutor(max_workers=1) as pool:
 results=dict(pool.map(fetch,DOIS.items()))
for key,doi in [('yousefi2026','10.64898/2026.09.28.755172'),('blackburne2026','10.64898/2026.10.01.755915')]:
 try:
  response=requests.get('https://api.biorxiv.org/details/biorxiv/'+doi,timeout=30)
  response.raise_for_status()
  data=response.json()
  (OUT/(key+'_biorxiv.json')).write_text(json.dumps(data,indent=2),encoding='utf-8')
  # User's Yousefi citation is v1; use the latest version of the Mediano-group preprint.
  versions=data['collection']
  chosen=versions[0] if key=='yousefi2026' else versions[-1]
  results[key]={'doi':doi,'title':chosen['title'],'date':chosen['date'],'version':chosen['version']}
  xml=requests.get(chosen['jatsxml'],timeout=30)
  xml.raise_for_status()
  (OUT/(key+'.xml')).write_bytes(xml.content)
 except Exception as exc:
  results[key]={'doi':doi,'error':str(exc)}
(OUT/'reference_fetch_receipt.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
print(json.dumps(results,indent=2))
