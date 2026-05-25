#!/usr/bin/env python3
#
# Filter Tabula Sapiens deposited h5ad into pipeline format.
#
# Filters: 10x assays only, published_2022 == False (exclude v1 training data),
#          manually_annotated == True, bad labels, min_cells=50.
# Output: data/tabula_sapiens/counts.h5ad  (filtered; var_names = ENSEMBL IDs)
# All obs columns preserved in the output h5ad for downstream use.

import os
import re
import argparse
import numpy as np
import h5py
import pandas as pd
import scipy.sparse as sp

DOUBLET_PATTERN = re.compile(
    r'doublet|contaminant|contamination|low[._]quality|'
    r'unassigned|unknown|debris|multiplet',
    re.IGNORECASE
)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--proj_root',
        default=os.getenv('PROJ_ROOT',
                          os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    args = parser.parse_args()

    import scanpy as sc

    proj_root = args.proj_root
    raw_h5ad = os.path.join(proj_root, 'raw_data', 'tabula_sapiens', 'counts.h5ad')
    out_dir = os.path.join(proj_root, 'data', 'tabula_sapiens')
    os.makedirs(out_dir, exist_ok=True)

    print(f"Loading obs from {raw_h5ad} (backed)", flush=True)
    adata = sc.read_h5ad(raw_h5ad, backed='r')
    print(f"Raw: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes", flush=True)

    # 10x assays only
    is_10x = adata.obs['assay'].astype(str).str.contains('10x', case=False)
    print(f"Keeping {is_10x.sum():,} 10x cells (dropping {(~is_10x).sum():,})", flush=True)
    print(f"  Assay breakdown:\n{adata.obs.loc[is_10x.values, 'assay'].value_counts().to_string()}", flush=True)

    # exclude v1 cells (published_2022 == True were in SCimilarity / Census training)
    is_new = adata.obs['published_2022'].astype(str) == 'False'
    print(f"Keeping {is_new.sum():,} non-v1 cells (dropping {(~is_new).sum():,} published_2022==True)", flush=True)

    # manually annotated only
    is_manual = adata.obs['manually_annotated'].astype(str) == 'True'
    print(f"Keeping {is_manual.sum():,} manually annotated cells "
          f"(dropping {(~is_manual).sum():,})", flush=True)

    ct_raw = adata.obs['cell_type']
    ct = ct_raw.astype(str)

    # bad labels
    bad_label = ct_raw.isna() | ct.str.contains(DOUBLET_PATTERN, na=False)
    if bad_label.sum():
        print(f"Dropping {bad_label.sum():,} cells by label pattern", flush=True)

    # min_cells=50 (evaluated on the cells passing all prior filters)
    keep_pre = is_10x & is_new & is_manual & ~bad_label
    ct_pre = ct[keep_pre]
    type_counts = ct_pre.value_counts()
    drop_types = set(type_counts[type_counts < 50].index.tolist())
    if drop_types:
        n_drop_cells = ct_pre.isin(drop_types).sum()
        print(f"Dropping {len(drop_types)} type(s) with <50 cells "
              f"({n_drop_cells:,} cells): {sorted(drop_types)}", flush=True)

    keep = keep_pre & ~ct.isin(drop_types)
    keep_idx = np.where(keep.values)[0]
    obs_sub = adata.obs.iloc[keep_idx].copy()

    # Close the backed file and immediately free the full 1.1M-cell obs from RAM
    # before starting the heavy h5py read.
    adata.file.close()
    del adata
    import gc; gc.collect()

    print(f"Loading raw counts for {len(keep_idx):,} cells via h5py "
          f"(skipping normalized X)...", flush=True)
    with h5py.File(raw_h5ad, 'r') as f:
        grp = f['raw/X']
        indptr = grp['indptr'][:].astype(np.int64)
        n_cells_full = len(indptr) - 1

        var_names = [
            v.decode() if isinstance(v, bytes) else v
            for v in f['raw/var/_index'][:]
        ]
        n_genes = len(var_names)

        var_cols = {}
        for col in f.get('raw/var', {}).keys():
            if col == '_index':
                continue
            try:
                arr = f[f'raw/var/{col}'][:]
                var_cols[col] = [
                    v.decode() if isinstance(v, bytes) else v for v in arr
                ]
            except Exception:
                pass

        # Pass 1: count total NNZ for selected rows from the indptr alone (no data read).
        total_nnz = int((indptr[keep_idx + 1] - indptr[keep_idx]).sum())
        print(f"  Total NNZ for filtered cells: {total_nnz:,}", flush=True)

        # Pre-allocate output arrays — one copy of the filtered matrix, never two.
        # The old approach (accumulate chunks + vstack) held two copies at peak,
        # which exceeded 128 GB for this dataset.
        out_data = np.empty(total_nnz, dtype=grp['data'].dtype)
        out_indices = np.empty(total_nnz, dtype=grp['indices'].dtype)
        out_indptr = np.zeros(len(keep_idx) + 1, dtype=np.int64)

        CHUNK = 10_000
        out_row = 0
        out_pos = 0
        for block_start in range(0, n_cells_full, CHUNK):
            block_end = min(block_start + CHUNK, n_cells_full)
            in_block = (keep_idx >= block_start) & (keep_idx < block_end)
            local_keep = keep_idx[in_block] - block_start
            if len(local_keep) == 0:
                continue
            ptr_s = int(indptr[block_start])
            ptr_e = int(indptr[block_end])
            block_data = grp['data'][ptr_s:ptr_e]
            block_cols = grp['indices'][ptr_s:ptr_e]
            block_ptr = indptr[block_start:block_end + 1] - ptr_s
            block_csr = sp.csr_matrix(
                (block_data, block_cols, block_ptr),
                shape=(block_end - block_start, n_genes),
            )
            sub = block_csr[local_keep]
            del block_csr, block_data, block_cols, block_ptr

            n_rows = len(local_keep)
            n_nnz = sub.nnz
            out_data[out_pos:out_pos + n_nnz] = sub.data
            out_indices[out_pos:out_pos + n_nnz] = sub.indices
            # sub.indptr[1:] is cumulative within this sub-block; shift by out_pos
            # to get absolute positions in the output arrays.
            out_indptr[out_row + 1:out_row + 1 + n_rows] = sub.indptr[1:].astype(np.int64) + out_pos
            out_row += n_rows
            out_pos += n_nnz
            del sub

    X_raw = sp.csr_matrix((out_data, out_indices, out_indptr), shape=(len(keep_idx), n_genes))
    del out_data, out_indices, out_indptr
    X_raw.check_format(full_check=True)  # validates indptr, index bounds, sorted indices

    var_df = pd.DataFrame(var_cols, index=var_names)
    adata = sc.AnnData(X=X_raw, obs=obs_sub, var=var_df)

    n_types = adata.obs['cell_type'].nunique()
    n_donors = adata.obs['donor_id'].nunique()
    n_tissues = adata.obs['tissue_in_publication'].nunique()
    print(f"Final: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes "
          f"| {n_types} types | {n_donors} donors | {n_tissues} tissues", flush=True)

    out_h5ad = os.path.join(out_dir, 'counts.h5ad')
    print(f"Writing {out_h5ad}...", flush=True)
    adata.write_h5ad(out_h5ad)
    print("Done.", flush=True)

if __name__ == '__main__':
    main()
