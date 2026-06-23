# sc_feature_representations

Benchmarking **feature representation strategies** for single-cell RNA-seq classification.

We compare a full gene set **reference** against linear (PCA), nonlinear kernel approximation (Random Fourier Features — Laplacian + Gaussian), and pretrained deep models (SCimilarity, zero-shot Census scVI) on cell-type and tissue classification across several datasets. The comparison axis is the **representation**, not the classifier, the classifier (glmnet / SVM) is just the readout.

## Methods compared

| Method | Pretrained | Supervised | Source |
|--------|-----------|------------|--------|
| Reference | No | No | fit on train split |
| PCA | No | No | fit on train split |
| RFF (Laplacian / Gaussian) | No | No | fit on train split |
| scVI | Yes | No | CZI S3 Census 2024-07-01, unsupervised VAE (50-dim, zero-shot) |
| SCimilarity | Yes | Yes (cell-type labels) | HuggingFace |

## Datasets

All confirmed out-of-distribution for the pretrained models (not in SCimilarity training or CellxGene Census).

| Dataset | Cells | Tasks |
|---------|-------|-------|
| `zilionis_lung` | ~51K | celltype (22) |
| `tabula_sapiens` | ~421K | celltype (broad 31 / fine 121), tissue (26) |
| `crc_icb` | ~540K | celltype (6), compound (12) |
| `lodi_pancancer` | ~611K | celltype = 9-way cancer-type-of-origin |

## Repository layout

```
src/             active pipeline scripts
preprocessing/   data-prep scripts
data/            prepared per-dataset inputs (gitignored)
experiments/     raw results (gitignored)
```

### Pipeline stages (`src/`)

1. **Load + split** (`io.R`) — `load_dataset()`; `make_split()` reserves a donor-pure test set, then takes nested cell-level train prefixes. Classes with <2 training cells in the smallest (40%) prefix are dropped so `test ⊆ train` at every fraction.
2. **Preprocess** (`utils.R`) — log-normalize (library size scale factor 10000) → `log1p` → per-gene min-max scale; optional HVG selection (deep learning excluded).
3. **Reduction** (`reduction.R`) — `rff_reduce()`, `pca_reduce()`, `scimilarity_embed()`, `scvi_embed()` (the last two call the Python scripts on their respective `scimilarity` / `scvi` conda envs).
4. **Train + predict** (`models.R`) — `train_glmnet()` / `train_svm()`.
5. **Metrics** (`metrics.R`) — (base R + pROC): macro precision/recall/F1 over test-present classes, micro accuracy, AUROC.

## Usage

Run from the repo root so `data/` and `experiments/` resolve (`PROJ_ROOT` defaults to `.`):

```bash
Rscript src/classification.R \
  -r <run_id>            # replicate; seeds the split + model training
  -d <dataset>           # zilionis_lung | tabula_sapiens | crc_icb | lodi_pancancer
  -t <celltype|tissue>   # task (tissue only where tissue.qs exists)
  -a <glmnet|svm>        # classifier (default glmnet)
  -s <40|60|80>          # train percentage
  -m <method>            # rff_lapl | rff_gauss | pca | scimilarity | scvi; OMIT for the reference (full gene set)
  -n <int>               # output dimensions (RFF uses n_dim/2 internal projections)
```

Additional flags:
- `--n_hvg <int>` — optional HVG selection (e.g. 2000); applied to reference/pca/rff, not scimilarity/scvi
- `--label_level <fine|compound>` — alternate label set
- `--no_gpu` — disable GPU for RFF projection

## Dependencies

**R:** `Matrix`, `MatrixExtra`, `sparseMatrixStats`, `qs2`, `rhdf5`, `RSpectra`, `glmnet`, `e1071`, `pROC`, `argparse`

**Python:** `scanpy`, `torch`, `h5py`, `numpy`, plus `scimilarity` and `scvi-tools` — managed via two separate conda envs (`scimilarity` and `scvi`); each embedding script invokes its own env's interpreter.
