import scanpy as sc
import resource
from scimilarity import CellAnnotation
from scimilarity.utils import align_dataset, lognorm_counts
import h5py
import numpy as np
import argparse
import os
import random
from datetime import datetime


def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate scimilarity embeddings")
    parser.add_argument("--data_dir", required=True, help="Path to dataset directory (contains counts.h5ad)")
    parser.add_argument("--indices_file", required=True, help="H5 file containing train/test indices")
    parser.add_argument("--output_file", required=True, help="Output H5 file for embeddings and metadata")
    parser.add_argument("--model_dir", default="/scratch/zeynepaydin21/scimilarity/models/model_v1.1")
    parser.add_argument("--seed", type=int)
    return parser.parse_args()


def load_indices(filepath):
    with h5py.File(filepath, 'r') as f:
        return f['train_indices'][:], f['test_indices'][:]


def main():
    args = parse_arguments()

    np.random.seed(args.seed)
    random.seed(args.seed)

    train_indices, test_indices = load_indices(args.indices_file)

    data_file = os.path.join(args.data_dir, "counts.h5ad")
    preprocess_start = datetime.now()

    print("Loading expression matrix...")
    mtx = sc.read(data_file)
    mtx.obs_names_make_unique()
    mtx.var_names_make_unique()

    mtx_train = mtx[train_indices].copy()
    mtx_test = mtx[test_indices].copy()

    print("Loading scimilarity model...")
    ca = CellAnnotation(model_path=args.model_dir)

    print("Aligning and log-normalizing...")
    mtx_train = align_dataset(mtx_train, ca.gene_order)
    mtx_test = align_dataset(mtx_test, ca.gene_order)

    mtx_train.layers["counts"] = mtx_train.X.copy()
    mtx_train = lognorm_counts(mtx_train)
    mtx_test.layers["counts"] = mtx_test.X.copy()
    mtx_test = lognorm_counts(mtx_test)

    preprocess_time = (datetime.now() - preprocess_start).total_seconds()

    print("Generating embeddings...")
    embedding_start = datetime.now()
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
        f.attrs['preprocess_time'] = preprocess_time
        f.attrs['peak_memory_gb'] = peak_gb


if __name__ == "__main__":
    main()
