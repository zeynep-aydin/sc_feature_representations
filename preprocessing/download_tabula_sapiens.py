#!/usr/bin/env python3
"""Download Tabula Sapiens from CellxGene Census, filtered to 10x assays only."""

import sys
from pathlib import Path

import cellxgene_census
import anndata

COLLECTION_ID = "e5f58829-1a66-40b5-a624-9046778e74f5"
TENX_ASSAYS = ["10x 3' v3", "10x 5' v2"]

PROJ_ROOT = Path(__file__).resolve().parents[1]
OUT_PATH = PROJ_ROOT / "data" / "tabula_sapiens" / "counts.h5ad"

ALL_CELLS_DATASET_ID = "53d208b0-2cfd-4366-9866-c3c6114081bc"

OBS_COLS = [
    "cell_type",
    "tissue",
    "donor_id",
    "assay",
]

def main():
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    assay_str = ", ".join(f'"{a}"' for a in TENX_ASSAYS)
    obs_filter = (
        f'dataset_id == "{ALL_CELLS_DATASET_ID}"'
        f" and assay in [{assay_str}]"
    )

    with cellxgene_census.open_soma(census_version="2025-11-08") as census:
        print(f"Filter: {obs_filter}")
        print("Streaming anndata from Census...")

        adata = cellxgene_census.get_anndata(
            census,
            organism="Homo sapiens",
            obs_value_filter=obs_filter,
            obs_column_names=OBS_COLS,
            X_name="raw",
        )

    print(f"\nDownloaded {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"Assay breakdown:\n{adata.obs['assay'].value_counts().to_string()}")
    print(f"Donors: {adata.obs['donor_id'].nunique()}")

    print(f"\nWriting to {OUT_PATH} ...")
    adata.write_h5ad(OUT_PATH)
    print("Done.")

if __name__ == "__main__":
    main()
