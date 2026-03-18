# alignment-free-hgt-detection

Alignment-free detection of horizontal gene transfer (HGT) candidates via protein similarity graphs, species-pair robust normalization, and graph anomaly scoring.

Paper title:
**Alignment-Free Detection of Horizontal Gene Transfer Candidates via Protein Similarity Graphs and Anomaly Scoring**

## Repository layout

- `src/hgt_pipeline/`: core pipeline package
- `src/graph_construction/`: pre-pipeline graph-construction steps (FASTA parsing, k-mer candidates, pruning)
- `src/hgt_pipeline/stages/`: graph-analysis stages used by `pipeline.py`
- `tools/reporting/`: reporting/explanation scripts
- `golden/`: regression baselines (canonical reference artifacts)
- `artifacts/`: historical and transient run outputs
- `tests/`: regression test suite
- `docs/`: manuscript and supporting report material

## Quickstart

From repository root (PowerShell):

```powershell
pushd src
```

### 1) Graph Construction (runs first)

This stage produces the pruned edge list consumed by `hgt_pipeline.pipeline`.

```powershell
python -m graph_construction.orchestrator construct-edges --manifest ..\golden\reference_inputs\manifest_tiny_set.tsv --downloads_dir ..\golden\reference_inputs\downloads_tiny --out_candidates ..\tmp_candidates.tsv --out_edges ..\tmp_edges.tsv
```

For the full canonical run, the repository already includes the pruned input:
- `golden/reference_inputs/edges_PRUNED_JACCARD_92790.tsv`

### 2) HGT Pipeline (graph input -> features/candidates)

With betweenness:

```powershell
python -m hgt_pipeline.pipeline --in_edges ..\golden\reference_inputs\edges_PRUNED_JACCARD_92790.tsv --out_dir ..\tmp_run_bw
```

Without betweenness (faster):

```powershell
python -m hgt_pipeline.pipeline --in_edges ..\golden\reference_inputs\edges_PRUNED_JACCARD_92790.tsv --out_dir ..\tmp_run_no_bw --no_betweenness
popd
```

Reviewer-facing, module-only E2E command lines are documented in [`REPRODUCE.md`](REPRODUCE.md).

```powershell
pushd src
python -m graph_construction.refseq_fetch_proteins --help
python -m graph_construction.orchestrator --help
python -m hgt_pipeline.pipeline --help
popd
```

## Reproduce Component Plots (5, 8, 32)

From repository root (PowerShell):

```powershell
python tools/reporting/plot_components.py `
  --edges artifacts/legacy_runs/out/edge_features.tsv `
  --protein_features artifacts/legacy_runs/out/protein_features.tsv `
  --hgt_candidates artifacts/legacy_runs/out/hgt_candidates.tsv `
  --component_ids 5,8,32 `
  --node_size_mode constant `
  --out_dir artifacts/updated_plots/simplified

python tools/reporting/plot_components.py `
  --edges artifacts/legacy_runs/out/edge_features.tsv `
  --protein_features artifacts/legacy_runs/out/protein_features.tsv `
  --hgt_candidates artifacts/legacy_runs/out/hgt_candidates.tsv `
  --component_ids 5,8,32 `
  --node_size_mode score `
  --out_dir artifacts/updated_plots/detailed
```

## Regression tests

Default fast gate (recommended while refactoring):

```powershell
python tests/test_regression_baselines.py
```

Modes:

- `graph` (default): no-betweenness pipeline regression + fast checks
- `graph_bw`: includes betweenness-on regression
- `full`: includes heavy full k-mer regression (only if full local data exists)

Examples:

```powershell
python tests/test_regression_baselines.py --mode graph
python tests/test_regression_baselines.py --mode graph_bw
```

## Data policy

Large raw data is not tracked in git (`data/` ignored).  
Regression relies on `golden/` fixtures, including small deterministic test inputs.

`data/assembly_summary_refseq.txt` is intentionally not tracked (very large; ~216 MB).
Fetch it directly from NCBI when needed (see `REPRODUCE.md`).

## Citation

Citation metadata is in `CITATION.cff`.
