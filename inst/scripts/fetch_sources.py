"""Download public documentary inputs and record SHA-256 and wall times.

Run from the package root. Only the registered SCP datasets, the authorised
bulk-cohort search candidates, and the named calibration paper are fetched.
"""
import hashlib
import json
import time
import urllib.request
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / ".local" / "sources"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261013


def download(name, url):
    started = time.monotonic()
    path = OUT / name
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "scpPower-M1"})
        with urllib.request.urlopen(req, timeout=90) as response:
            content = response.read()
        path.write_bytes(content)
        if name.endswith((".xlsx", ".zip")) and not content.startswith(b"PK"):
            return dict(name=name, url=url,
                        sha256=hashlib.sha256(content).hexdigest(), bytes=len(content),
                        wall_seconds=time.monotonic()-started, status="unavailable",
                        error="Server returned a browser challenge or other non-ZIP content",
                        seed=SEED)
        return dict(name=name, url=url, sha256=hashlib.sha256(content).hexdigest(),
                    bytes=len(content), wall_seconds=time.monotonic()-started,
                    status="downloaded", seed=SEED)
    except Exception as exc:
        return dict(name=name, url=url, wall_seconds=time.monotonic()-started,
                    status="unavailable", error=str(exc), seed=SEED)


if __name__ == "__main__":
    base = "https://raw.githubusercontent.com/UCLouvain-CBIO/scpdata/master/"
    sources = []
    for dataset in ["specht2019v3", "leduc2022_pSCoPE", "derks2022", "brunner2022"]:
        for path in [f"R/{dataset}.R", f"man/{dataset}.Rd",
                     f"inst/scripts/make-data_{dataset}.R"]:
            sources.append((path.replace("/", "__"), base + path))
    for accession in ["PXD058573", "PXD004352"]:
        api = "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/" + accession
        sources.append((accession+"-project.json", api))
        # PRIDE caps file-list pages at 100; request all pages, not one page.
        for page in range(15):
            sources.append((f"{accession}-files-{page}.json",
                            api+f"/files?page={page}&pageSize=100"))
    sources += [
        ("cheung-pubmed.html", "https://pubmed.ncbi.nlm.nih.gov/33288958/"),
        ("cheung-publisher.html", "https://www.nature.com/articles/s41592-020-01002-5"),
        ("exercise-publisher.html", "https://www.nature.com/articles/s41467-025-68101-9"),
        ("rieckmann-publisher.html", "https://www.nature.com/articles/ni.3693"),
        ("monocyte-12donors.html", "https://pmc.ncbi.nlm.nih.gov/articles/PMC5874767/"),
        ("orcid-person.json", "https://pub.orcid.org/v3.0/0000-0002-2758-171X/person"),
        ("PXD058573-output.tsv", "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD058573/output.tsv"),
        ("PXD058573-sdrf.tsv", "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD058573/PXD058573_community_annotated.sdrf.tsv"),
        ("PXD058573-setup.txt", "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/01/PXD058573/setup.txt"),
        ("exercise-supplements.zip", "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12775541/supplementaryFiles"),
        ("exercise-fulltext.xml", "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12775541/fullTextXML"),
    ]
    for n in [4, 5, 6]:
        sources.append((f"cheung-dataset{n-2}.xlsx",
                        "https://media.springernature.com/original/springer-static/esm/"
                        "art%3A10.1038%2Fs41592-020-01002-5/MediaObjects/"
                        f"41592_2020_1002_MOESM{n}_ESM.xlsx"))
    with ThreadPoolExecutor(max_workers=6) as executor:
        results = list(executor.map(lambda x: download(*x), sources))
    (OUT/"download_manifest.json").write_text(json.dumps(results, indent=2)+"\n")
    for result in results:
        if result["status"] == "unavailable":
            print(result["name"], result["error"])
    print("Public-source manifest written.")
