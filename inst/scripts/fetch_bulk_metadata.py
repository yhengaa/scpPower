"""Fetch the authorised cohort's donor map without unrelated inline figures.

Europe PMC streams ZIP entries. Stop once the complete metadata entry is read
and verify its MD5 against the article XML, then record its SHA-256. This avoids
requiring the much larger, unrelated figure and analysis supplements.
"""
import hashlib
import json
import re
import struct
import time
import urllib.request
import zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / ".local" / "sources"
OUT.mkdir(parents=True, exist_ok=True)
XML_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12775541/fullTextXML"
ZIP_URL = ("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12775541/"
           "supplementaryFiles?includeInlineImage=false")
NAME = "41467_2025_68101_MOESM3_ESM.xlsx"


def extract_member(buffer, target):
    pos = 0
    while buffer[pos:pos+4] == b"PK\x03\x04":
        if len(buffer) < pos+30:
            return None
        fields = struct.unpack_from("<4s5H3I2H", buffer, pos)
        flag, method, name_len, extra_len = fields[2], fields[3], fields[-2], fields[-1]
        if len(buffer) < pos+30+name_len+extra_len:
            return None
        name = buffer[pos+30:pos+30+name_len].decode()
        start = pos+30+name_len+extra_len
        if method != 8:
            raise RuntimeError("Unexpected ZIP method; cannot recover verified metadata")
        decoder = zlib.decompressobj(-15)
        content = decoder.decompress(buffer[start:])
        if not decoder.eof:
            return None
        used = len(buffer[start:])-len(decoder.unused_data)
        if name == target:
            return content
        pos = start+used
        if flag & 8:
            pos += 16 if buffer[pos:pos+4] == b"PK\x07\x08" else 12
    return None


if __name__ == "__main__":
    start = time.monotonic()
    with urllib.request.urlopen(XML_URL, timeout=90) as response:
        xml = response.read()
    (OUT/"exercise-fulltext.xml").write_bytes(xml)
    block = re.search(rb'<supplementary-material[^>]*id="MOESM3".*?</supplementary-material>', xml, re.S)
    expected = re.search(rb'<\?suppdata-md5 ([a-f0-9]+)\?>', block.group()).group(1).decode()
    buffer = bytearray()
    result = None
    with urllib.request.urlopen(ZIP_URL, timeout=120) as response:
        while len(buffer) < 32*1024*1024:
            chunk = response.read(65536)
            if not chunk:
                break
            buffer.extend(chunk)
            result = extract_member(bytes(buffer), NAME)
            if result is not None:
                break
    if result is None:
        raise RuntimeError("Donor metadata could not be recovered from the public supplement")
    if hashlib.md5(result).hexdigest() != expected:
        raise RuntimeError("Donor metadata differs from the checksum published in the article")
    target = OUT/"exercise-supplements"/NAME
    target.parent.mkdir(exist_ok=True)
    target.write_bytes(result)
    manifest = dict(seed=20261013, url=ZIP_URL, member=NAME, article_xml_url=XML_URL,
        article_xml_sha256=hashlib.sha256(xml).hexdigest(),
        md5_published=expected, md5_verified=True,
        sha256=hashlib.sha256(result).hexdigest(), bytes=len(result),
        wall_seconds=time.monotonic()-start, status="downloaded and checksum verified")
    (OUT/"bulk-metadata-download-manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print("Author donor metadata downloaded and checksum verified.")
