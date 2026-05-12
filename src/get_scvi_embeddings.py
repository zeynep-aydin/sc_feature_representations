import scvi
import scanpy as sc
import h5py
import numpy as np
import argparse
import os
import resource
import random
import torch
from datetime import datetime


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate scVI embeddings")
    parser.add_argument("--data_dir", required=True, help="Path to dataset directory (contains counts.h5ad)")
    parser.add_argument("--indices_file", required=True, help="H5 file containing train/test indices")
    parser.add_argument("--output_file", required=True, help="Output H5 file for embeddings and metadata")
    parser.add_argument("--n_latent", type=int, default=50, help="scVI latent dimensions")
    parser.add_argument("--max_epochs", type=int, default=400, help="Max training epochs")
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def load_indices(filepath):
    with h5py.File(filepath, 'r') as f:
        return f['train_indices'][:], f['test_indices'][:]


def main():
    args = parse_arguments()

    np.random.seed(args.seed)
    random.seed(args.seed)
    torch.manual_seed(args.seed)
    scvi.settings.seed = args.seed

    train_indices, test_indices = load_indices(args.indices_file)

    data_file = os.path.join(args.data_dir, "counts.h5ad")
    preprocess_start = datetime.now()

    print("Loading expression matrix...")
    mtx = sc.read(data_file)
    mtx.obs_names_make_unique()
    mtx.var_names_make_unique()

    adata_train = mtx[train_indices].copy()
    adata_test = mtx[test_indices].copy()
    del mtx

    scvi.model.SCVI.setup_anndata(adata_train)
    preprocess_time = (datetime.now() - preprocess_start).total_seconds()

    print(f"Training scVI model (n_latent={args.n_latent}, max_epochs={args.max_epochs})...")
    embedding_start = datetime.now()
    model = scvi.model.SCVI(adata_train, n_latent=args.n_latent)
    
    if torch.cuda.is_available():
        accelerator = "cuda"
    else:
        accelerator = "cpu"
    print(f"Using accelerator: {accelerator}")
    
    model.train(
        max_epochs=args.max_epochs,
        enable_progress_bar=True,
        accelerator=accelerator,
    )

    train_embeddings = model.get_latent_representation(adata_train)
    test_embeddings = model.get_latent_representation(adata_test)
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
        f.attrs['preprocess_time'] = preprocess_time
        f.attrs['peak_memory_gb'] = peak_gb


if __name__ == "__main__":
    main()
