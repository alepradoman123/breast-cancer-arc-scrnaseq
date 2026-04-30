# Step-by-step: VS Code, environment, GitHub

Follow these in order. Commands are for macOS (your setup).

## 1. Install prerequisites

If you don't have them:

```bash
# Homebrew
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Miniconda (smaller than Anaconda, recommended)
brew install --cask miniconda
conda init zsh   # or: conda init bash, depending on your shell
# Restart your terminal

# Git
brew install git

# VS Code
brew install --cask visual-studio-code
```

In VS Code, install these extensions (Extensions pane, `Cmd+Shift+X`):

- **Python** (Microsoft)
- **Jupyter** (Microsoft)
- **Pylance** (Microsoft, for type hints)
- **GitLens** (optional, nice git UI)

## 2. Get the project

```bash
# Option A: if you received the project as a folder, just cd into it
cd ~/path/to/breast-cancer-arc-scrnaseq

# Option B: clone from your own GitHub (after step 6 below)
# git clone https://github.com/YOUR_USERNAME/breast-cancer-arc-scrnaseq.git
# cd breast-cancer-arc-scrnaseq
```

## 3. Create the conda environment

```bash
conda env create -f environment.yml
conda activate arc-breast
python -m ipykernel install --user --name arc-breast --display-name "Python (arc-breast)"
```

This takes 5-10 minutes. If it fails on any pip package, try:

```bash
conda create -n arc-breast python=3.11 -y
conda activate arc-breast
pip install scanpy harmonypy leidenalg python-igraph scrublet decoupler pydeseq2 adjustText GEOparse openpyxl jupyter ipykernel
python -m ipykernel install --user --name arc-breast --display-name "Python (arc-breast)"
```

## 4. Open in VS Code

```bash
code .
```

When VS Code opens:

1. `Cmd+Shift+P` → "Python: Select Interpreter" → pick `arc-breast`
2. Open any `.ipynb` file in `notebooks/` → click "Select Kernel" top-right → choose "Python (arc-breast)"

## 5. Run the pipeline

Run notebooks in order:

1. `01_data_download.ipynb` — downloads GSE176078 (~500 MB). Takes 5-10 min.
2. `02_qc_filtering.ipynb` — QC, doublet removal
3. `03_normalization_hvg.ipynb` — normalization + subsampling (in `DEV_MODE`)
4. `04_integration_clustering.ipynb` — Harmony + UMAP + Leiden
5. `05_annotation.ipynb` — cell type annotation
6. `06_arc_expression.ipynb` — the biological analysis
7. `07_figures.ipynb` (optional) — polish figures for a report

### If you hit memory errors on the 8 GB Mac

Edit `src/config.py` and lower `CELLS_PER_PATIENT_DEV` (e.g., 1000 or 800).

### For the full dataset run

Set `DEV_MODE = False` in `src/config.py` and run on Google Colab or a bigger machine.

## 6. Initialize the git repository and push to GitHub

```bash
cd ~/path/to/breast-cancer-arc-scrnaseq

# First-time git setup (if you haven't already)
git config --global user.name "Your Name"
git config --global user.email "you@example.com"

# Initialize the repo
git init
git add .
git commit -m "Initial commit: scRNA-seq pipeline for ARC in breast cancer"

# Create a new empty repo on GitHub via the web UI
# (do NOT initialize it with README/gitignore - you already have both)
# Then link and push:

git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/breast-cancer-arc-scrnaseq.git
git push -u origin main
```

### What gets pushed vs. excluded

Pushed to GitHub:
- All `.py` modules in `src/`
- All notebooks in `notebooks/`
- `environment.yml`, `README.md`, `.gitignore`, `SETUP.md`
- Empty `.gitkeep` placeholders for `figures/`, `results/`, `data/`

Excluded (via `.gitignore`):
- `data/raw/*` and `data/processed/*` (dataset files, too large)
- `*.h5ad` files (can be hundreds of MB)
- `figures/*.png`, `results/*.csv` by default (can be force-added with `git add -f` if you want to commit small ones)
- `__pycache__`, `.ipynb_checkpoints`, `.DS_Store`

## 7. Subsequent commits

As you work:

```bash
git status                                  # see what changed
git add src/some_file.py notebooks/06*.ipynb   # stage specific files
git commit -m "Refine pseudobulk test for cell-type-level stats"
git push
```

### Tip: commit notebooks without output

Nb outputs make diffs ugly. Strip outputs before committing:

```bash
pip install nbstripout
nbstripout --install             # sets it up as a git hook for this repo
```

From then on, `git add notebook.ipynb` will automatically strip outputs.

## 8. Running on Google Colab (for the full dataset)

1. Push the repo to GitHub (step 6).
2. In Colab: `File → Open notebook → GitHub → paste your repo URL`.
3. First cell of each notebook, add:

   ```python
   !git clone https://github.com/YOUR_USERNAME/breast-cancer-arc-scrnaseq.git
   %cd breast-cancer-arc-scrnaseq
   !pip install scanpy harmonypy leidenalg scrublet decoupler pydeseq2 adjustText
   ```

4. Set `DEV_MODE = False` in `src/config.py` via the Colab file browser, or inline:

   ```python
   import src.config as cfg
   cfg.DEV_MODE = False
   ```

5. Run as normal. Mount Google Drive if you want to persist `data/` across sessions:

   ```python
   from google.colab import drive
   drive.mount('/content/drive')
   ```
