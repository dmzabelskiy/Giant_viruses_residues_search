# gv_project

Pipeline for analysis of viral gene families against cellular references. The workflow extracts viral members of a curated family, builds a combined viral-plus-cellular alignment, calls family-discriminating residues, ranks representative viral sequences, and can generate PyMOL highlighting scripts for shortlisted sites.

## Scope

The current codebase is organized around generalized family-level scripts rather than a single protein family workflow. The main pipeline supports any family present in the curated metadata tables, with optional short acronyms such as `acon`, `gapdh`, or `sod`.

## Repository Layout

- `scripts/`: pipeline scripts and helpers
- `README.md`: project overview and usage
- `requirements.txt`: Python package requirements
- `initial_raw/`: source data inputs, alignments, and raw sequence collections
- `data/results/`: per-family outputs
- `figures/`, `logs/`: generated artifacts

The repository is intended to track code and documentation. Large raw inputs and generated outputs are ignored via `.gitignore`.

## Requirements

- Python 3.10+
- MAFFT available on `PATH` for alignment steps
- Python packages from `requirements.txt`
- Optional: PyMOL for visualization

Install the Python dependencies with:

```bash
python3 -m pip install -r requirements.txt
```

If the project is checked out outside the original working location, set:

```bash
export GV_PROJECT_ROOT=/path/to/gv_project
```

## Main Family Pipeline

The easiest entry point is:

```bash
python3 scripts/run_family_pipeline.py \
  --family <exact_gene_family_name> \
  --results-subdir <short_slug> \
  --alignment-source <path/to/source_alignment.fa>
```

Example:

```bash
python3 scripts/run_family_pipeline.py \
  --family TCA_cycle_Aconitase \
  --results-subdir acon \
  --alignment-source initial_raw/marker_gene_phylogenies/Phylogenetic_tree_alignments_and_tree_files/Aconitase_align.fa
```

This wrapper runs steps `03` through `08`:

1. extract viral family sequences
2. parse family membership from the source alignment
3. filter viral sequences by family-specific length thresholds
4. build the combined ungapped dataset
5. run MAFFT
6. call viral-versus-cellular residues
7. rank and choose a representative viral sequence

If you want to run the family pipeline manually, the family-based scripts are:

1. `scripts/03_extract_family_sequences.py`
2. `scripts/04_extract_family_from_alignment.py`
3. `scripts/05_filter_family_lengths.py`
4. `scripts/06_prepare_family_main_dataset.py`
5. `scripts/07_call_family_residues.py`
6. `scripts/08_choose_representative_viral_family.py`

There are also upstream data-preparation utilities:

1. `scripts/01_build_metadata.py`
2. `scripts/02_index_curated_gene_trees.py`

## Output Layout

Each family gets its own result directory:

```text
data/results/<family_slug>/
```

Top-level files are reserved for final user-facing outputs, including:

- `<family_slug>_main_ungapped.faa`
- `<family_slug>_main_mafft.trimmed.faa`
- `<family_slug>_viral_vs_cell_residue_calls.tsv`
- `<family_slug>_viral_vs_cell_residue_shortlist.tsv`
- `viral_<family_slug>_full_length.tsv`
- `viral_<family_slug>_representative_ranking.tsv`
- `viral_<family_slug>_representative_for_colabfold.faa`

Intermediate and QC files are written under:

```text
data/results/<family_slug>/intermediate/
```

## Residue Shortlist Logic

Residue shortlisting is based on the existing evidence axes:

- BH-adjusted Fisher q-value
- sequence-harmony-like score
- viral-versus-cellular residue-frequency difference
- viral top-residue frequency
- gap filtering

The canonical shortlist logic now uses rank percentiles across eligible positions instead of fixed top-15 or top-30 cutoffs. The main output columns include:

- `rank_q_pct`
- `rank_harmony_pct`
- `rank_freq_diff_pct`
- `rank_percentile_mean`
- `category`
- `ranking_score`
- `site_weight_score`

Canonical categories are:

- `rankpct_high_confidence`
- `rankpct_moderate_confidence`

This keeps the previous biological framing while scaling more smoothly across proteins of different lengths.

## PyMOL Helper

Use `scripts/highlight_top_family_sites_pymol.py` to map shortlisted family-specific sites onto a representative structure model and write a PyMOL script.

Example:

```bash
python3 scripts/highlight_top_family_sites_pymol.py \
  --shortlist data/results/acon/acon_viral_vs_cell_residue_shortlist.tsv \
  --alignment data/results/acon/acon_main_mafft.trimmed.faa \
  --representative-id "<sequence_id>" \
  --pdb model.pdb \
  --split-by-confidence \
  --init ~/Work/init.py
```

When `--split-by-confidence` is used, high- and moderate-confidence sites are colored separately.

## Preparing for GitHub

Before publishing:

- review whether any small example data should be added separately
- choose a license for the repository
- initialize git and create the remote

Typical setup commands:

```bash
git init
git add README.md requirements.txt .gitignore scripts
git commit -m "Initial commit"
git branch -M main
git remote add origin <your-github-repo-url>
git push -u origin main
```
