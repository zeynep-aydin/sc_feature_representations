import scanpy as sc
import scipy.io
import argparse


def main():
    parser = argparse.ArgumentParser(description="Save sparse count matrix as h5ad")
    parser.add_argument("--mtx", required=True, help="Matrix Market file (cells x genes)")
    parser.add_argument("--genes", required=True, help="Gene symbols, one per line")
    parser.add_argument("--out", required=True, help="Output .h5ad path")
    args = parser.parse_args()

    X = scipy.io.mmread(args.mtx).tocsr()
    with open(args.genes) as f:
        genes = [line.strip() for line in f]

    adata = sc.AnnData(X)
    adata.var_names = genes
    adata.write_h5ad(args.out)
    print(f"Saved {adata.shape[0]} cells x {adata.shape[1]} genes -> {args.out}")


if __name__ == "__main__":
    main()
