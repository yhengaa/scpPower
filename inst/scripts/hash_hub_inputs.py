"""Record the actual downloaded ExperimentHub files, not just R objects.

Run after inspect_scp.R. The hub downloads its public RDA files; this script
records the input file hash, public URL, and software versions in the manifest.
"""
import hashlib
import json
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CACHE = Path.home()/".cache"/"R"/"ExperimentHub"
connection = sqlite3.connect("file:"+str(CACHE/"BiocFileCache.sqlite")+"?mode=ro", uri=True)
connection.row_factory = sqlite3.Row
rows = connection.execute("SELECT rname, rpath, fpath, last_modified_time, etag FROM resource").fetchall()
# These resource IDs are returned by the four registered scpdata loaders.
datasets = {"3943": "specht2019v3", "7762": "leduc2022_pSCoPE",
            "7763": "derks2022", "8360": "brunner2022"}
manifest_path = ROOT/"reports"/"M1"/"run-manifest.json"
manifest = json.loads(manifest_path.read_text())
sources = []
for row in rows:
    item = dict(row)
    identifier = item["rname"].removeprefix("EH")
    if identifier not in datasets:
        identifier = str(item["fpath"]).rstrip("/").rsplit("/", 1)[-1]
    if identifier not in datasets:
        continue
    path = Path(item["rpath"])
    if not path.is_absolute():
        path = CACHE/path
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024*1024):
            hasher.update(block)
    sources.append(dict(dataset=datasets[identifier], hub_resource=item["rname"].split(" : ")[0],
        fetch_resource=identifier,
        public_url=item["fpath"], sha256=hasher.hexdigest(), bytes=path.stat().st_size,
        etag=item["etag"], last_modified=item["last_modified_time"],
        seed=20261013, versions=manifest["versions"]))
if len(sources) != 4:
    raise RuntimeError("Could not identify all four downloaded public hub inputs")
manifest["downloaded_hub_inputs"] = sources
manifest_path.write_text(json.dumps(manifest, indent=2)+"\n")
(ROOT/"reports"/"M1"/"public-input-hashes.json").write_text(json.dumps(sources, indent=2)+"\n")
print("Hashed all four downloaded ExperimentHub inputs.")
