import scvi
import scanpy as sc
import h5py
import numpy as np
import pandas as pd
import torch
import argparse
import os
import resource
from datetime import datetime

DEFAULT_MODEL_DIR = os.environ.get("SCVI_MODEL_DIR")


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate pretrained Census scVI embeddings (zero-shot)")
    parser.add_argument("--data_dir", required=True)
    parser.add_argument("--indices_file", required=True)
    parser.add_argument("--output_file", required=True)
    parser.add_argument("--model_dir", default=DEFAULT_MODEL_DIR, help="Local directory containing model.pt")
    return parser.parse_args()


def load_indices(filepath):
    with h5py.File(filepath, 'r') as f:
        return f['train_indices'][:], f['test_indices'][:]


def load_model_var_names(model_dir):
    ckpt = torch.load(os.path.join(model_dir, "model.pt"), map_location="cpu", weights_only=False)
    return list(ckpt["var_names"])


def load_eid_to_sym():
    df = pd.read_csv(os.environ["GENE_MAP_TSV"], sep="\t", dtype=str)
    return dict(zip(df["gene_id"], df["gene_name"]))


def remap_to_ensembl(adata, model_var_names, eid_to_sym):
    """Map the dataset's genes into the model's ENSEMBL namespace.

    If the dataset is already ENSEMBL, strip versions and return it as-is; the
    intersection with the model is left to prepare_query_anndata. If the dataset
    is symbols, map each model gene back to its symbol and keep the dataset
    columns that match, renamed to the model's ENSEMBL IDs.
    """
    model_set = set(model_var_names)
    head = list(adata.var_names[:50])

    if all(g.startswith("ENSG") for g in head):
        adata.var_names = [g.split(".")[0] for g in adata.var_names]
        adata.var_names_make_unique()
        n_overlap = sum(1 for g in adata.var_names if g in model_set)
        print(f"Input already ENSEMBL ({adata.n_vars} genes); "
              f"{n_overlap}/{len(model_var_names)} "
              f"({100 * n_overlap / len(model_var_names):.1f}%) of model genes present")
        return adata

    sym_to_col = {}
    for i, sym in enumerate(adata.var_names):
        sym_to_col.setdefault(sym, i)  # first dataset column wins per symbol

    cols = []
    names = []
    used_syms = set()
    for eid in model_var_names:
        sym = eid_to_sym.get(eid.split(".")[0])
        if sym is None or sym in used_syms:
            continue
        col = sym_to_col.get(sym)
        if col is None:
            continue
        cols.append(col)
        names.append(eid)
        used_syms.add(sym)

    adata = adata[:, cols].copy()
    adata.var_names = names
    print(f"Symbol->ENSEMBL (model-first): {len(names)}/{len(model_var_names)} "
          f"({100 * len(names) / len(model_var_names):.1f}%) of model genes matched")
    return adata


def main():
    args = parse_arguments()

    if args.model_dir is None or not os.path.exists(os.path.join(args.model_dir, "model.pt")):
        # Download manually before running:
        # wget -O <model_dir>/model.pt \
        #   https://cellxgene-contrib-public.s3.us-west-2.amazonaws.com/models/scvi/2024-07-01/homo_sapiens/model.pt
        raise FileNotFoundError(f"model.pt not found in {args.model_dir}. Set --model_dir or export SCVI_MODEL_DIR.")
    print(f"Using model: {os.path.join(args.model_dir, 'model.pt')}")

    train_indices, test_indices = load_indices(args.indices_file)

    data_file = os.path.join(args.data_dir, "counts.h5ad")
    preprocess_start = datetime.now()
    print(f"Loading {data_file} ...")
    adata = sc.read(data_file)
    adata.obs_names_make_unique()
    preprocess_time = (datetime.now() - preprocess_start).total_seconds()

    print("Remapping genes, aligning to model gene set, and loading scVI model...")
    embedding_start = datetime.now()
    model_var_names = load_model_var_names(args.model_dir)
    eid_to_sym = load_eid_to_sym()
    adata = remap_to_ensembl(adata, model_var_names, eid_to_sym)

    scvi.model.SCVI.prepare_query_anndata(adata, args.model_dir)
    print(f"After prepare_query_anndata: {adata.n_vars} genes")
    if adata.n_vars == 0:
        raise ValueError("prepare_query_anndata produced 0 genes — check ENSEMBL ID format in counts.h5ad")

    adata.obs["batch"] = "query"
    model = scvi.model.SCVI.load_query_data(adata, args.model_dir, freeze_dropout=True)
    model.is_trained_ = True

    print("Generating embeddings...")
    all_embeddings = model.get_latent_representation()
    train_embeddings = all_embeddings[train_indices]
    test_embeddings = all_embeddings[test_indices]
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
        f.attrs['model_dir'] = args.model_dir
        f.attrs['gene_map'] = os.environ.get("GENE_MAP_TSV", "")
        f.attrs['n_genes_after_align'] = adata.n_vars
        f.attrs['embedding_time'] = embedding_time
        f.attrs['preprocess_time'] = preprocess_time
        f.attrs['peak_memory_gb'] = peak_gb


if __name__ == "__main__":
    main()
