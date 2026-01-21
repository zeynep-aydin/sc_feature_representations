# A comparative benchmark of linear, nonlinear and deep learning methods for single-cell transcriptomics

This repository provides the scripts for our benchmarking study comparing three feature representation strategies for single-cell RNA sequencing data:

1. **Linear dimensionality reduction** - Principal Component Analysis (PCA)
2. **Nonlinear kernel approximation** - Random Fourier Features (RFF) with Gaussian and Laplacian kernels
3. **Deep learning foundation models** - SCimilarity

We evaluate these methods across four diverse scRNA-seq datasets on cell type and tissue classification tasks, assessing classification accuracy, computational runtime, and memory efficiency.

## Repository Contents

```
├── classification.R              # Main classification pipeline
├── calculate_sigma.R             # Laplacian kernel bandwidth estimation
├── calculate_sigma_gaussian.R    # Gaussian kernel bandwidth estimation
├── get_scimilarity_embeddings.py # SCimilarity embedding generation
└── utils.R                       # Preprocessing and utility functions
```

**Data Availability:** Preprocessed data matrices are available upon request.


## Running Experiments

The main classification pipeline is executed via `classification.R`:

```bash
Rscript classification.R \
  --replicate <REP> \
  --dataset <DATASET> \
  --method <METHOD> \
  --algorithm <ALGORITHM> \
  --task <TASK> \
  --train_split <SPLIT>
```

**All Parameters:**
- `--dataset`: `scea`, `he_organs`, `zhao_immune`, `zilionis_lung`
- `--method`: `baseline`, `pca`, `rff`, `rff_gaussian`, `scimilarity`
- `--algorithm`: `glmnet`, `svm`
- `--feature_mode`: Feature mode selection for SCEA, either `intersection` or `thresh08`
- `--task`: `celltype` or `tissue`
- `--train_split`: Training proportion (0.5, 0.7, 0.8) or use 
- `--n_cells`: for fixed training sizes
- `--D`: RFF dimensions (100, 500, 1000)
- `--n_components`: PCA components (200, 1000, 2000)
- `--replicate`: Replicate number (1-30 in our study)


## Output Structure

Results are saved in:
```
experiments/{dataset}/classification/{method}/{train_config}/rep{X}/
```

Each replicate generates two files:
- `*_metadata.RData`: Performance metrics, timing, memory usage, parameters
- `*_model.RData`: Trained model object


## Contact

For questions about the code or data:
- Zeynep Aydın: zeynepaydin21@ku.edu.tr
- Mehmet Gönen: mehmetgonen@ku.edu.tr
