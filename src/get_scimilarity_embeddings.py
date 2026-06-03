import scanpy as sc
import resource
import torch
from scimilarity import CellAnnotation
from scimilarity.utils import align_dataset, lognorm_counts
import h5py
import numpy as np
import pandas as pd
import requests
import argparse
import os
import random
from datetime import datetime

BIOMART_URL = "https://may2025.archive.ensembl.org/biomart/martservice"


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate scimilarity embeddings")
    parser.add_argument("--data_dir", required=True, help="Path to dataset directory (contains counts.h5ad)")
    parser.add_argument("--indices_file", required=True, help="H5 file containing train/test indices")
    parser.add_argument("--output_file", required=True, help="Output H5 file for embeddings and metadata")
    parser.add_argument("--model_dir", default=os.environ.get("SCIMILARITY_MODEL_DIR"))
    parser.add_argument("--seed", type=int)
    return parser.parse_args()


def load_indices(filepath):
    with h5py.File(filepath, 'r') as f:
        return f['train_indices'][:], f['test_indices'][:]


def build_ensembl_to_symbol_map(ensembl_ids, model_dir):
    cache = os.path.join(model_dir, "ensembl_to_symbol.tsv")
    if os.path.exists(cache):
        print(f"Loading ENSEMBL->symbol map from cache: {cache}")
        df = pd.read_csv(cache, sep="\t", dtype=str)
        return dict(zip(df["ensembl_id"], df["symbol"]))

    print("Querying Ensembl BioMart for ENSEMBL->symbol (batches of 500)...")
    eid_to_sym = {}
    ids = list(ensembl_ids)
    for i in range(0, len(ids), 500):
        batch = ids[i:i + 500]
        query = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Filter name="ensembl_gene_id" value="{",".join(batch)}"/>
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="external_gene_name"/>
  </Dataset>
</Query>"""
        r = requests.post(BIOMART_URL, data={"query": query}, timeout=300)
        r.raise_for_status()
        for line in r.text.strip().split("\n")[1:]:
            parts = line.split("\t")
            if len(parts) >= 2 and parts[0].strip() and parts[1].strip():
                eid_to_sym[parts[0].strip()] = parts[1].strip()
        if i % 5000 == 0:
            print(f"  batch {i // 500 + 1}/{-(-len(ids) // 500)}: {len(eid_to_sym)} mapped")

    with open(cache, "w") as f:
        f.write("ensembl_id\tsymbol\n")
        for eid, sym in sorted(eid_to_sym.items()):
            f.write(f"{eid}\t{sym}\n")
    print(f"Cached ENSEMBL->symbol map -> {cache} ({len(eid_to_sym)} mappings)")
    return eid_to_sym


def remap_to_symbols(adata, model_dir):
    if not all(g.startswith("ENSG") for g in list(adata.var_names[:20])):
        print(f"var_names don't look like ENSEMBL IDs — skipping remap ({adata.n_vars} genes)")
        return adata
    eid_to_sym = build_ensembl_to_symbol_map(list(adata.var_names), model_dir)
    new_names = [eid_to_sym.get(eid, eid) for eid in adata.var_names]
    n_mapped = sum(1 for n in new_names if not n.startswith("ENSG"))
    adata.var_names = new_names
    adata.var_names_make_unique()
    print(f"ENSEMBL->symbol: {n_mapped}/{len(new_names)} genes mapped; {len(new_names) - n_mapped} kept as ENSEMBL IDs")
    return adata


def main():
    args = parse_arguments()

    np.random.seed(args.seed)
    random.seed(args.seed)
    torch.manual_seed(args.seed)

    train_indices, test_indices = load_indices(args.indices_file)

    data_file = os.path.join(args.data_dir, "counts.h5ad")
    preprocess_start = datetime.now()
    print("Loading expression matrix...")
    mtx = sc.read(data_file)
    mtx.obs_names_make_unique()
    mtx_train = mtx[train_indices].copy()
    mtx_test = mtx[test_indices].copy()
    del mtx
    preprocess_time = (datetime.now() - preprocess_start).total_seconds()

    print("Loading scimilarity model, remapping genes, aligning and log-normalizing...")
    if not os.path.isdir(args.model_dir):
        # Download before running (replace <dir> with parent of model_v1.1/):
        # curl -L -o <dir>/model_v1.1.tar.gz \
        #   https://zenodo.org/records/10685499/files/model_v1.1.tar.gz?download=1
        # tar -xzvf <dir>/model_v1.1.tar.gz -C <dir>/
        raise FileNotFoundError(f"SCimilarity model directory not found: {args.model_dir}")
    embedding_start = datetime.now()
    use_gpu = torch.cuda.is_available()
    ca = CellAnnotation(model_path=args.model_dir, use_gpu=use_gpu)
    print(f"Model device: {next(ca.model.parameters()).device}")

    mtx_train = remap_to_symbols(mtx_train, args.model_dir)
    mtx_test = remap_to_symbols(mtx_test, args.model_dir)

    mtx_train = align_dataset(mtx_train, ca.gene_order)
    mtx_test = align_dataset(mtx_test, ca.gene_order)

    lognorm_start = datetime.now()
    mtx_train.layers["counts"] = mtx_train.X.copy()
    mtx_train = lognorm_counts(mtx_train)
    mtx_test.layers["counts"] = mtx_test.X.copy()
    mtx_test = lognorm_counts(mtx_test)
    lognorm_time = (datetime.now() - lognorm_start).total_seconds()

    print("Generating embeddings...")
    train_embeddings = ca.get_embeddings(mtx_train.X)
    test_embeddings = ca.get_embeddings(mtx_test.X)
    embedding_time = (datetime.now() - embedding_start).total_seconds()

    peak_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)

    with h5py.File(args.output_file, 'w') as f:
        f.create_dataset('train_embeddings', data=train_embeddings)
        f.create_dataset('test_embeddings', data=test_embeddings)
        f.attrs['train_n_cells'] = train_embeddings.shape[0]
        f.attrs['train_n_features'] = train_embeddings.shape[1]
        f.attrs['test_n_cells'] = test_embeddings.shape[0]
        f.attrs['test_n_features'] = test_embeddings.shape[1]
        f.attrs['data_file'] = data_file
        f.attrs['embedding_time'] = embedding_time
        f.attrs['lognorm_time'] = lognorm_time
        f.attrs['preprocess_time'] = preprocess_time
        f.attrs['peak_memory_gb'] = peak_gb


if __name__ == "__main__":
    main()
