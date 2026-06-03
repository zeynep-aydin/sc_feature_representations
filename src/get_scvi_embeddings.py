import scvi
import scanpy as sc
import h5py
import numpy as np
import pandas as pd
import requests
import torch
import argparse
import os
import resource
from datetime import datetime

DEFAULT_MODEL_DIR = os.environ.get("SCVI_MODEL_DIR")
BIOMART_URL = "https://may2025.archive.ensembl.org/biomart/martservice"


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


def build_symbol_map(model_dir):
    cache = os.path.join(model_dir, "ensembl_symbol_to_ensembl_id.tsv")
    if os.path.exists(cache):
        print(f"Loading gene map from cache: {cache}")
        df = pd.read_csv(cache, sep="\t", dtype=str)
        return dict(zip(df["symbol"], df["ensembl_id"]))

    print("Querying Ensembl BioMart for gene symbols (batches of 500)...")
    ckpt = torch.load(os.path.join(model_dir, "model.pt"), map_location="cpu", weights_only=False)
    var_names = list(ckpt["var_names"])

    eid_to_sym = {}
    for i in range(0, len(var_names), 500):
        batch = var_names[i:i + 500]
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
        print(f"  batch {i // 500 + 1}/{-(-len(var_names) // 500)}: {len(eid_to_sym)} symbols so far")

    sym_to_eid = {sym: eid for eid, sym in eid_to_sym.items()} # invert mapping; some symbols may map to multiple ENSEMBL IDs
    with open(cache, "w") as f:
        f.write("symbol\tensembl_id\n")
        for sym, eid in sorted(sym_to_eid.items()):
            f.write(f"{sym}\t{eid}\n")
    print(f"Cached gene map -> {cache} ({len(sym_to_eid)} mappings)")
    return sym_to_eid


def remap_to_ensembl(adata, model_dir):
    if all(g.startswith("ENSG") for g in list(adata.var_names[:20])):
        print(f"var_names already ENSEMBL IDs — skipping remap ({adata.n_vars} genes)")
        return adata
    sym_to_eid = build_symbol_map(model_dir)
    mapped = [sym_to_eid.get(g) for g in adata.var_names]
    valid = [i for i, e in enumerate(mapped) if e is not None]
    adata = adata[:, valid].copy()
    adata.var_names = [mapped[i] for i in valid]
    print(f"Symbol->ENSEMBL: {len(valid)} / {len(mapped)} genes mapped")
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
    adata = remap_to_ensembl(adata, args.model_dir)

    scvi.model.SCVI.prepare_query_anndata(adata, args.model_dir)
    print(f"After prepare_query_anndata: {adata.n_vars} genes")
    if adata.n_vars == 0:
        raise ValueError("prepare_query_anndata produced 0 genes — check ENSEMBL ID format in counts.h5ad")

    adata.obs["batch"] = "query"  # required by model registry (batch_key="batch"), treated as new batch
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
        f.attrs['n_genes_after_align'] = adata.n_vars
        f.attrs['embedding_time'] = embedding_time
        f.attrs['preprocess_time'] = preprocess_time
        f.attrs['peak_memory_gb'] = peak_gb


if __name__ == "__main__":
    main()
