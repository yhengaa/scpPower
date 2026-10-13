"""Package the final M1 deliverables for the authorised GitHub release."""
from pathlib import Path
import zipfile
root=Path(__file__).resolve().parents[2]
with zipfile.ZipFile(root/'scpPower_M1_v0.1.0_evidence.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted((root/'reports/M1-final').iterdir()):
  if p.is_file():z.write(p,str(p.relative_to(root)))
 for name in ['scpPower_mandate_v1.1.md','NEWS.md','renv.lock','README.md']:
  z.write(root/name,name)
 for p in sorted((root/'inst/scripts').glob('*')):
  if p.is_file():z.write(p,str(p.relative_to(root)))
