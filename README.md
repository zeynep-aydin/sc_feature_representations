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
| SCimilarity | Yes | Yes (cell-type labels) | Zenodo (`model_v1.1`) |

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
3. **Reduction** (`reduction.R`) — `rff_reduce()`, `pca_reduce()`, `scimilarity_embed()`, `scvi_embed()` (the last two shell out to the Python embedding scripts via the interpreters named by `SCIMILARITY_PYTHON` / `SCVI_PYTHON`).
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

## Environment variables

Only the pretrained-model methods (`scimilarity`, `scvi`) need these — reference/PCA/RFF run with no setup. `PROJ_ROOT` is optional everywhere (defaults to the current directory).

| Variable | Used by | Purpose |
|----------|---------|---------|
| `PROJ_ROOT` | all | repo root for resolving `data/` and `experiments/` (default `.`) |
| `SCIMILARITY_MODEL_DIR` | scimilarity | directory containing the SCimilarity model |
| `SCVI_MODEL_DIR` | scvi | directory containing the pretrained scVI model (`model.pt`) |
| `GENE_MAP_TSV` | scimilarity, scvi | `gene_id`/`gene_name` TSV for ENSEMBL↔symbol remap (build with `preprocessing/build_gene_map.py`) |
| `SCIMILARITY_PYTHON` | scimilarity | Python interpreter for the embedding script (default: `python` on `PATH`) |
| `SCVI_PYTHON` | scvi | Python interpreter for the embedding script (default: `python` on `PATH`) |

If `SCIMILARITY_PYTHON` / `SCVI_PYTHON` are unset, the scripts run via whatever `python` is on your `PATH`. Point them at a specific interpreter (e.g. a conda env's `bin/python`) if the required packages live in a dedicated environment.

## Pretrained models (scimilarity / scvi only)

These two methods depend on externally-downloaded model weights (not included in the repo). Download them once, then point the corresponding `*_MODEL_DIR` variable at the directory.

**SCimilarity** — `model_v1.1` from Zenodo. Set `SCIMILARITY_MODEL_DIR` to the extracted `model_v1.1/` directory.

```bash
curl -L -o <dir>/model_v1.1.tar.gz \
  https://zenodo.org/records/10685499/files/model_v1.1.tar.gz?download=1
tar -xzvf <dir>/model_v1.1.tar.gz -C <dir>/
```

**scVI** — zero-shot CellxGene Census model, release **2024-07-01**, from CZI's public S3. Set `SCVI_MODEL_DIR` to the directory containing the downloaded `model.pt`.

```bash
wget -O <model_dir>/model.pt \
  https://cellxgene-contrib-public.s3.us-west-2.amazonaws.com/models/scvi/2024-07-01/homo_sapiens/model.pt
```

SCimilarity and scVI have conflicting dependencies, so they must live in **separate** environments. Specs are in `envs/`:

```bash
conda env create -f envs/scimilarity.yml   # then set SCIMILARITY_PYTHON=$(conda run -n scimilarity which python)
conda env create -f envs/scvi.yml          # then set SCVI_PYTHON=$(conda run -n scvi which python)
```

## Dependencies

**R:** `Matrix`, `MatrixExtra`, `sparseMatrixStats`, `qs2`, `rhdf5`, `RSpectra`, `glmnet`, `e1071`, `pROC`, `argparse`

**Python:** `scanpy`, `torch`, `h5py`, `numpy`, plus `scimilarity` and `scvi-tools`).

## Contact

For questions about the code or data:
- Zeynep Aydın: zeynepaydin21@ku.edu.tr
- Mehmet Gönen: mehmetgonen@ku.edu.tr
