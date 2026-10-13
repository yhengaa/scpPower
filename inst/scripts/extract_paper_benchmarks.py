"""Extract comparison numbers from downloaded primary papers; no fitted tuning."""
from pathlib import Path
import json,re,time,hashlib,xml.etree.ElementTree as E
from html.parser import HTMLParser
start=time.monotonic()
class Text(HTMLParser):
 def __init__(self):super().__init__();self.parts=[]
 def handle_data(self,d):self.parts.append(d)
def text(name):
 p=Path('.local/sources/'+name+'-fulltext.xml')
 try:t=' '.join(E.parse(p).getroot().itertext())
 except E.ParseError:
  h=Text();h.feed(p.read_text());t=' '.join(h.parts)
 return re.sub(r'\s+',' ',t)
def number(t,pattern):
 m=re.search(pattern,t,re.I)
 if not m:raise RuntimeError('Cannot reproduce paper benchmark: '+pattern)
 return float(m.group(1).replace(',',''))
s,l,d,b=map(text,['specht','leduc','derks','brunner'])
bench=dict(specht=dict(proteins=number(s,r'diverse set of (3042) proteins'),cells=number(s,r'(1490) single cells have lower'),mean_cells_per_protein=number(s,r'per protein is about (305)'),within_identified_peptide_missing=number(s,r'only (10)% of the data'),source='https://doi.org/10.1186/s13059-021-02267-5'),leduc=dict(mean_proteins_per_cell=number(l,r'on average (997) proteins per cell'),proteins=number(l,r'with (2844) proteins quantified'),cells=number(l,r'across the (1556) single cells'),cv_threshold=number(l,r'threshold of (0\.41)'),source='https://doi.org/10.1186/s13059-022-02817-5'),derks=dict(within_plex_completeness=number(d,r'(98)% data completeness'),proteins_per_cell=number(d,r'quantified\s*[~∼]?\s*([1,]+000) proteins per cell'),source='https://doi.org/10.1038/s41587-022-01389-w'),brunner=dict(min_proteins_per_cell=number(b,r'at least (600) protein'),median_cv=number(b,r'median coefficient of variation \(CV\) of (0\.3)'),starting_cells=number(b,r'The (434) single'),source='https://doi.org/10.15252/msb.202110798'))
# All data-derived paper numbers in the note come from these parsed fields.
bench['specht']['global_missing_fraction']=1-bench['specht']['mean_cells_per_protein']/bench['specht']['cells']
bench['specht']['implied_mean_proteins_per_cell']=bench['specht']['proteins']*bench['specht']['mean_cells_per_protein']/bench['specht']['cells']
bench['leduc']['implied_missing_fraction']=1-bench['leduc']['mean_proteins_per_cell']/bench['leduc']['proteins']
manifest=json.loads(Path('reports/M1-resumed/run-manifest.json').read_text())
Path('reports/M1-resumed/paper-benchmarks.json').write_text(json.dumps(dict(benchmarks=bench,seed=20261013,versions=manifest['versions'],input_hashes={n:hashlib.sha256(Path('.local/sources/'+n+'-fulltext.xml').read_bytes()).hexdigest() for n in ['specht','leduc','derks','brunner']},wall_seconds=time.monotonic()-start),indent=2)+'\n')
print('Reproduced all paper comparison numbers.')
