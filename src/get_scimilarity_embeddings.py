import scanpy as sc
import resource
from scimilarity import CellAnnotation
from scimilarity.utils import align_dataset, lognorm_counts
import h5py
import numpy as np
import pandas as pd
import argparse
import os
import sys
import random
from datetime import datetime

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate scimilarity embeddings")
    parser.add_argument("--data_dir", required=True, help="Path to data directory")
    parser.add_argument("--feature_mode", required=True, help="Feature mode (e.g., thresh08)")
    parser.add_argument("--indices_file", required=True, help="Input H5 file containing train/test indices")
    parser.add_argument("--output_file", required=True, help="Output H5 file for embeddings and metadata")
    parser.add_argument("--seed", type=int, help="Random seed for reproducibility")
    return parser.parse_args()

def load_indices(filepath):
   with h5py.File(filepath, 'r') as f:
        train_indices = f['train_indices'][:]
        test_indices = f['test_indices'][:]
        return train_indices, test_indices

def get_data_filename(data_dir, feature_mode):
  dataset = os.path.basename(data_dir.rstrip('/'))
  if dataset == "scea":
    filename = f"scea_8tissues_{feature_mode}_expression_matrix_symbol.h5ad"
  else:
    filename = f"{dataset}_intersection_expression_matrix_symbol.h5ad"

  return os.path.join(data_dir, filename), dataset

def main():
  args = parse_arguments()

  model_dir = "path/to/your/model"

  np.random.seed(args.seed)
  random.seed(args.seed)

  train_indices, test_indices = load_indices(args.indices_file)

  data_file, dataset = get_data_filename(args.data_dir, args.feature_mode)

  preprocess_start_time = datetime.now()

  print("Loading and splitting expression matrix...")
  mtx = sc.read(data_file)
  mtx.obs_names_make_unique()
  mtx.var_names_make_unique()

  print("Splitting data into train/test sets...")
  mtx_train = mtx[train_indices].copy()
  mtx_test = mtx[test_indices].copy()

  print("Loading scimilarity model...")
  ca = CellAnnotation(model_path = model_dir)

  print("Aligning datasets with model gene order and log-normalizing...")
  mtx_train = align_dataset(mtx_train, ca.gene_order)
  mtx_test = align_dataset(mtx_test, ca.gene_order)

  mtx_train.layers["counts"] = mtx_train.X.copy()
  mtx_train = lognorm_counts(mtx_train)

  mtx_test.layers["counts"] = mtx_test.X.copy()
  mtx_test = lognorm_counts(mtx_test)

  preprocess_time = (datetime.now() - preprocess_start_time).total_seconds()

  embedding_start_time = datetime.now()

  print("Generating embeddings...")
  train_embeddings = ca.get_embeddings(mtx_train.X)
  test_embeddings = ca.get_embeddings(mtx_test.X)

  embedding_time = (datetime.now() - embedding_start_time).total_seconds()
  peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
  peak_gb = peak_kb / (1024 * 1024)

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
