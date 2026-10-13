"""Parse M3 durations from the same downloaded primary source papers."""
from pathlib import Path
from lxml import etree,html
import re,json,hashlib,time
start=time.monotonic();out=Path('reports/M3');out.mkdir(exist_ok=True)
records={}
rules={
 'specht':('SCoPE2 samples were separated','Loading and separating each sample took',r'took\s+(\d+)\s*[  ]*min', 'https://doi.org/10.1186/s13059-021-02267-5','Reported total loading and separation method, not active gradient alone'),
 'leduc':('All samples except for plexDIA single cells','95-min total',r'(\d+)-min total length method','https://doi.org/10.1186/s13059-022-02817-5','Reported total TMT method, not the plexDIA method in the same paper'),
 'derks_qe':('plexDIA single-cell sets and 100-cell standards were injected','minute 53',r'until minute\s+(\d+)','https://doi.org/10.1038/s41587-022-01389-w','Q Exactive single-cell full reported gradient, including re-equilibration'),
 'derks_tims':('The single-cell plexDIA sets were separated','separated over 30 minutes',r'separated over\s+(\d+)\s+minutes','https://doi.org/10.1038/s41587-022-01389-w','timsTOF single-cell reported full gradient'),
 'brunner':('For all other proteome analyses','35',r'novel\s+(\d+)\s*min','https://doi.org/10.15252/msb.202110798','Reported single-cell EvoSep gradient; unreported loading/turnover overhead excluded')}
for key,(anchor,check,pattern,url,meaning) in rules.items():
 name='derks' if key.startswith('derks') else key
 p=Path('.local/sources')/(name+'-fulltext.xml')
 try:root=etree.fromstring(p.read_bytes())
 except etree.XMLSyntaxError:root=html.fromstring(p.read_bytes())
 paragraphs=[' '.join(' '.join(x.itertext()).split()) for x in root.xpath('//p')]
 candidates=[x for x in paragraphs if anchor in x]
 if len(set(candidates))!=1:raise RuntimeError(f'Ambiguous gradient section: {key}')
 text=candidates[0];match=re.search(pattern,text)
 if not match:raise RuntimeError(f'Gradient duration cannot be parsed: {key}')
 records[key]=dict(minutes_per_run=int(match.group(1)),source=url,definition=meaning,
   paper_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),section_anchor=anchor,
   section_sha256=hashlib.sha256(text.encode()).hexdigest())
vs=json.loads(Path('reports/M2-revised/run-manifest.json').read_text())['versions']
m=dict(seed=20261013,versions=vs,sources=records,wall_seconds=time.monotonic()-start,
   input_hashes={k:r['paper_sha256'] for k,r in records.items()},reproduction='Run fetch_resumed_sources.py before this parser; source URLs and hashes are preserved.')
(out/'instrument-times.json').write_text(json.dumps(m,indent=2)+'\n')
print({k:r['minutes_per_run'] for k,r in records.items()})
