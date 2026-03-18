# alignment-free-hgt-detection

Alignment-free detection of horizontal gene transfer (HGT) candidates via protein similarity graphs, species-pair robust normalization, and graph anomaly scoring.

Paper title:
**Alignment-Free Detection of Horizontal Gene Transfer Candidates via Protein Similarity Graphs and Anomaly Scoring**

## Repository layout

- `src/hgt_pipeline/`: core pipeline package
- `src/graph_construction/`: pre-pipeline graph-construction steps (FASTA parsing, k-mer candidates, pruning)
- `src/hgt_pipeline/stages/`: graph-analysis stages used by `pipeline.py`
- `src/reporting/`: reporting/explanation scripts
- `src/golden/`: regression baselines (canonical reference artifacts)
- `tests/`: regression test suite
- `docs/`: manuscript and supporting report material

## Quickstart

From repository root (PowerShell):

```powershell
$env:PYTHONPATH = "src"
python -m hgt_pipeline.pipeline --in_edges src/golden/reference_inputs/edges_PRUNED_JACCARD_92790.tsv --out_dir tmp_run --no_betweenness
```

Graph-construction orchestration (manifest/downloads -> candidates -> pruned edges):

```powershell
$env:PYTHONPATH = "src"
python -m graph_construction.orchestrator construct-edges --manifest src/golden/reference_inputs/manifest_tiny_set.tsv --downloads_dir src/golden/reference_inputs/downloads_tiny --out_candidates tmp_candidates.tsv --out_edges tmp_edges.tsv
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

Large raw data is not tracked in git (`src/data/` ignored).  
Regression relies on `src/golden/` fixtures, including small deterministic test inputs.

## Citation

Citation metadata is in `CITATION.cff`.
