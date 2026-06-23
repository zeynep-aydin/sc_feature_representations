import scanpy as sc
import resource
import torch
from scimilarity import CellAnnotation
from scimilarity.utils import align_dataset, lognorm_counts
import h5py
import numpy as np
import pandas as pd
import argparse
import os
import random
from datetime import datetime

DEFAULT_MODEL_DIR = os.environ.get("SCIMILARITY_MODEL_DIR")


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate scimilarity embeddings")
    parser.add_argument("--data_dir", required=True, help="Path to dataset directory (contains counts.h5ad)")
    parser.add_argument("--indices_file", required=True, help="H5 file containing train/test indices")
    parser.add_argument("--output_file", required=True, help="Output H5 file for embeddings and metadata")
    parser.add_argument("--model_dir", default=DEFAULT_MODEL_DIR, help="Local directory containing the SCimilarity model")
    parser.add_argument("--seed", type=int, default=4853)
    return parser.parse_args()


def load_indices(filepath):
    with h5py.File(filepath, 'r') as f:
        return f['train_indices'][:], f['test_indices'][:]


def load_eid_to_sym():
    gene_map = os.environ.get("GENE_MAP_TSV")
    if not gene_map:
        raise RuntimeError("GENE_MAP_TSV is not set — point it at the gene_id/gene_name TSV (see preprocessing/build_gene_map.py).")
    df = pd.read_csv(gene_map, sep="\t", dtype=str)
    return dict(zip(df["gene_id"], df["gene_name"]))


def remap_to_symbols(adata, model_gene_order, eid_to_sym):
    """Match the model's symbols against the dataset.

    Model gene_order is gene symbols. If the dataset is already symbols we
    hand them straight to align_dataset. If the dataset is ENSEMBL we translate
    each column to a symbol and keep only the ones the model asks for; columns
    whose symbol the model never uses are dropped.
    """
    model_set = set(model_gene_order)
    head = list(adata.var_names[:50])

    if not all(g.startswith("ENSG") for g in head):
        n_overlap = sum(1 for g in adata.var_names if g in model_set)
        print(f"dataset is symbols ({adata.n_vars} genes); "
              f"{n_overlap}/{len(model_gene_order)} "
              f"({100 * n_overlap / len(model_gene_order):.1f}%) of model genes present")
        if n_overlap == 0:
            raise ValueError("0 model genes present in dataset — check gene namespace in counts.h5ad")
        return adata

    sym_to_col = {}
    for i, eid in enumerate(adata.var_names):
        sym = eid_to_sym.get(eid.split(".")[0])
        if sym is not None and sym in model_set:
            sym_to_col.setdefault(sym, i)  # first ENSEMBL column wins per symbol

    if not sym_to_col:
        raise ValueError("0 model genes matched after ENSEMBL->symbol — check gene namespace / gene map")

    cols = list(sym_to_col.values())
    names = list(sym_to_col.keys())
    adata = adata[:, cols].copy()
    adata.var_names = names
    print(f"ENSEMBL->symbol: {len(names)}/{len(model_gene_order)} "
          f"({100 * len(names) / len(model_gene_order):.1f}%) of model genes matched")
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

    eid_to_sym = load_eid_to_sym()
    mtx = remap_to_symbols(mtx, ca.gene_order, eid_to_sym)
    mtx = align_dataset(mtx, ca.gene_order)

    mtx_train = mtx[train_indices].copy()
    mtx_test = mtx[test_indices].copy()
    del mtx

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
        f.attrs['gene_map'] = os.environ.get("GENE_MAP_TSV", "")
        f.attrs['embedding_time'] = embedding_time
        f.attrs['lognorm_time'] = lognorm_time
        f.attrs['preprocess_time'] = preprocess_time
        f.attrs['peak_memory_gb'] = peak_gb


if __name__ == "__main__":
    main()
