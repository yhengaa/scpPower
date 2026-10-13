"""Fetch only documentary sources for the four registered SCP datasets."""
import urllib.request, urllib.parse, json, hashlib, time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
OUT=Path('.local/sources'); OUT.mkdir(parents=True,exist_ok=True)
DOIS={'specht':'10.1186/s13059-021-02267-5','leduc':'10.1186/s13059-022-02817-5','derks':'10.1038/s41587-022-01389-w','brunner':'10.15252/msb.202110798'}
def fetch(item):
 name,doi=item;start=time.monotonic()
 search='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':'DOI:'+doi,'format':'json'})
 data=json.load(urllib.request.urlopen(search,timeout=90)); pmcid=next(x['pmcid'] for x in data['resultList']['result'] if x.get('pmcid'))
 url=f'https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML'

 try: content=urllib.request.urlopen(url,timeout=90).read()
 except Exception:
  url='https://pmc.ncbi.nlm.nih.gov/articles/'+pmcid+'/?report=xml'
  content=urllib.request.urlopen(url,timeout=90).read()
 path=OUT/(name+'-fulltext.xml');path.write_bytes(content)
 return dict(name=name,doi=doi,pmcid=pmcid,url=url,sha256=hashlib.sha256(content).hexdigest(),bytes=len(content),seed=20261013,wall_seconds=time.monotonic()-start)
with ThreadPoolExecutor(max_workers=4) as pool:
 records=list(pool.map(fetch,DOIS.items()))
(OUT/'resumed-source-downloads.json').write_text(json.dumps(records,indent=2)+'\n')
print('Downloaded and hashed all four source papers.')

url='https://raw.githubusercontent.com/cvanderaa/SCoPE2/master/code/SCoPE2_analysis.R'
start=time.monotonic();content=urllib.request.urlopen(url,timeout=90).read()
(OUT/'SCoPE2_analysis.R').write_bytes(content)
(OUT/'specht-analysis-source.json').write_text(json.dumps(dict(url=url,sha256=hashlib.sha256(content).hexdigest(),seed=20261013,wall_seconds=time.monotonic()-start)))

# Quantities and documentary inputs for the already-selected Walzik cohort.
inputs={'PXD058573-output.tsv':'https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD058573/output.tsv',
'exercise-fulltext.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12775541/fullTextXML',
'PXD058573-setup.txt':'https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD058573/setup.txt'}
extra=[]
for name,url in inputs.items():
 start=time.monotonic();path=OUT/name
 if not path.exists(): path.write_bytes(urllib.request.urlopen(url,timeout=180).read())
 content=path.read_bytes()
 extra.append(dict(name=name,url=url,sha256=hashlib.sha256(content).hexdigest(),bytes=len(content),seed=20261013,wall_seconds=time.monotonic()-start))
(OUT/'resumed-bulk-downloads.json').write_text(json.dumps(extra,indent=2)+'\n')
base='https://raw.githubusercontent.com/UCLouvain-CBIO/scpdata/master/'
extra=[]
for dataset in ['specht2019v3','leduc2022_pSCoPE','derks2022','brunner2022']:
 for rel in ['man/'+dataset+'.Rd','inst/scripts/make-data_'+dataset+'.R']:
  start=time.monotonic();url=base+rel;content=urllib.request.urlopen(url,timeout=90).read()
  name=rel.replace('/','__');(OUT/name).write_bytes(content)
  extra.append(dict(name=name,url=url,sha256=hashlib.sha256(content).hexdigest(),seed=20261013,wall_seconds=time.monotonic()-start))
(OUT/'resumed-assay-doc-downloads.json').write_text(json.dumps(extra,indent=2)+'\n')
