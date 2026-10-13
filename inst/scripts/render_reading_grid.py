"""Export the requested Reading the grid section as one reproducible page."""
from pathlib import Path
import json,hashlib,time,subprocess,tempfile
start=time.monotonic()
source=Path('reports/M3/power-report.md')
text=source.read_text().split('## Reading the grid',1)[1]
text='# Reading the grid\n\n13 October 2026\n'+text
with tempfile.TemporaryDirectory() as tmp:
 md=Path(tmp)/'reading.md';md.write_text(text)
 subprocess.run(['pandoc',str(md),'-o','reports/M3/reading-grid.pdf','--pdf-engine=pdflatex','-V','fontsize=10pt','-V','geometry:margin=0.6in'],check=True)
 info=subprocess.check_output(['pdfinfo','reports/M3/reading-grid.pdf'],text=True)
 assert any(line.split()==['Pages:','1'] for line in info.splitlines()),info
versions=json.loads(Path('reports/M3/run-manifest.json').read_text())['versions']
m={'seed':20261013,'versions':versions,'input_hashes':{'report':hashlib.sha256(source.read_bytes()).hexdigest(),'renderer':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},'output_sha256':hashlib.sha256(Path('reports/M3/reading-grid.pdf').read_bytes()).hexdigest(),'pages':1,'pandoc_version':subprocess.check_output(['pandoc','--version'],text=True).splitlines()[0],'wall_seconds':time.monotonic()-start}
Path('reports/M3/reading-grid-render-manifest.json').write_text(json.dumps(m,indent=2)+'\n')
print('Reading the grid: one-page PDF verified')
