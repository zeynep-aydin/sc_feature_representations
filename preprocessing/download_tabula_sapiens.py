#!/usr/bin/env python3
"""Download Tabula Sapiens deposited h5ad from CellxGene Data Portal.

Direct file download — all obs/var columns included as deposited by authors.
Supports resume if a partial file already exists.

Output: raw_data/tabula_sapiens/counts.h5ad (~45 GB)
"""

import os
import sys
from pathlib import Path

URL = "https://datasets.cellxgene.cziscience.com/5a495302-b7cd-4bf9-853e-95627b00bb03.h5ad"
EXPECTED_BYTES = 45013946701

PROJ_ROOT = Path(os.getenv("PROJ_ROOT", Path(__file__).resolve().parents[1]))
OUT_PATH = PROJ_ROOT / "raw_data" / "tabula_sapiens" / "counts.h5ad"


def download():
    import urllib.request

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    existing = OUT_PATH.stat().st_size if OUT_PATH.exists() else 0
    if existing == EXPECTED_BYTES:
        print(f"Already complete: {OUT_PATH} ({existing:,} bytes)")
        return
    if existing > 0:
        print(f"Resuming from byte {existing:,} / {EXPECTED_BYTES:,}")

    req = urllib.request.Request(URL)
    if existing > 0:
        req.add_header("Range", f"bytes={existing}-")

    mode = "ab" if existing > 0 else "wb"
    downloaded = existing

    with urllib.request.urlopen(req) as resp, open(OUT_PATH, mode) as fh:
        total = int(resp.headers.get("Content-Length", EXPECTED_BYTES - existing))
        chunk = 1024 * 1024  # 1 MB
        while True:
            buf = resp.read(chunk)
            if not buf:
                break
            fh.write(buf)
            downloaded += len(buf)
            pct = 100 * downloaded / EXPECTED_BYTES
            print(f"\r  {downloaded/1e9:.2f} / {EXPECTED_BYTES/1e9:.2f} GB  ({pct:.1f}%)",
                  end="", flush=True)

    print()
    final = OUT_PATH.stat().st_size
    if final != EXPECTED_BYTES:
        print(f"WARNING: expected {EXPECTED_BYTES:,} bytes, got {final:,}", file=sys.stderr)
    else:
        print(f"Download complete: {OUT_PATH}")


if __name__ == "__main__":
    download()
