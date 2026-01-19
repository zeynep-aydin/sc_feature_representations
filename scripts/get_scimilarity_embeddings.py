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
    parser.add_argument("--model_dir", default="/scratch/zeynepaydin21/scimilarity/models/model_v1.1", 
                        help="Path to scimilarity model directory")
    parser.add_argument("--mode", choices=['train_test', 'clustering'], default='train_test',
                        help="Splitting mode: 'train_test' for train/test split or 'clustering' for single indices")
    parser.add_argument("--seed", type=int, help="Random seed for reproducibility")
    return parser.parse_args()
  
def load_indices(filepath, mode):
   with h5py.File(filepath, 'r') as f:
        if mode == 'train_test':
            train_indices = f['train_indices'][:]
            test_indices = f['test_indices'][:]
            return train_indices, test_indices
        else:  # clustering mode
            indices = f['indices'][:]
            return indices
          
def get_data_filename(data_dir, feature_mode):
  dataset = os.path.basename(data_dir.rstrip('/'))
  if dataset == "scea":
    filename = f"scea_8tissues_{feature_mode}_expression_matrix_symbol.h5ad"
  else:
    filename = f"{dataset}_intersection_expression_matrix_symbol.h5ad"
        
  return os.path.join(data_dir, filename), dataset
  
def main():
  args = parse_arguments()
  
  np.random.seed(args.seed)
  random.seed(args.seed)
  
  if args.mode == 'train_test':
    train_indices, test_indices = load_indices(args.indices_file, mode = 'train_test')
  else:  
    clustering_indices = load_indices(args.indices_file, mode = 'clustering')
  
  data_file, dataset = get_data_filename(args.data_dir, args.feature_mode)
  
  preprocess_start_time = datetime.now()
   
  print("Loading and splitting expression matrix...")
  mtx = sc.read(data_file)
  mtx.obs_names_make_unique()
  mtx.var_names_make_unique()
  
  if args.mode == 'train_test':
    print("Splitting data into train/test sets...")
    mtx_train = mtx[train_indices].copy()
    mtx_test = mtx[test_indices].copy()
  else:  
    print("Extracting data for clustering...")
    mtx_clustering = mtx[clustering_indices].copy()
  
  print("Loading scimilarity model...")
  ca = CellAnnotation(model_path = args.model_dir)
  
  print("Aligning datasets with model gene order and log-normalizing...")
  if args.mode == 'train_test':
    mtx_train = align_dataset(mtx_train, ca.gene_order)
    mtx_test = align_dataset(mtx_test, ca.gene_order)
    
    mtx_train.layers["counts"] = mtx_train.X.copy()
    mtx_train = lognorm_counts(mtx_train)
    
    mtx_test.layers["counts"] = mtx_test.X.copy()
    mtx_test = lognorm_counts(mtx_test)
  else: 
    mtx_clustering = align_dataset(mtx_clustering, ca.gene_order)
    
    mtx_clustering.layers["counts"] = mtx_clustering.X.copy()
    mtx_clustering = lognorm_counts(mtx_clustering)
  
  preprocess_time = (datetime.now() - preprocess_start_time).total_seconds()
  
  embedding_start_time = datetime.now()
  
  print("Generating embeddings...")
  if args.mode == 'train_test':
    train_embeddings = ca.get_embeddings(mtx_train.X)
    test_embeddings = ca.get_embeddings(mtx_test.X)
  else: 
    clustering_embeddings = ca.get_embeddings(mtx_clustering.X)
  
  embedding_time = (datetime.now() - embedding_start_time).total_seconds()
  peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
  peak_gb = peak_kb / (1024 * 1024)
  
  with h5py.File(args.output_file, 'w') as f:
    if args.mode == 'train_test':
      f.create_dataset('train_embeddings', data=train_embeddings)
      f.create_dataset('test_embeddings', data=test_embeddings)
      
      f.attrs['train_n_cells'] = train_embeddings.shape[0]
      f.attrs['train_n_features'] = train_embeddings.shape[1]
      f.attrs['test_n_cells'] = test_embeddings.shape[0]
      f.attrs['test_n_features'] = test_embeddings.shape[1]
    else:  
      f.create_dataset('clustering_embeddings', data=clustering_embeddings)
      
      f.attrs['clustering_n_cells'] = clustering_embeddings.shape[0]
      f.attrs['clustering_n_features'] = clustering_embeddings.shape[1]
    
    f.attrs['data_file'] = data_file
    f.attrs['embedding_time'] = embedding_time
    f.attrs['preprocess_time'] = preprocess_time
    f.attrs['peak_memory_gb'] = peak_gb

    
if __name__ == "__main__":
    main()
