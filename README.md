# Alignment-Free Detection of Horizontal Gene Transfer

Alignment-free detection of horizontal gene transfer (HGT) candidates via protein similarity graphs, species-pair robust normalization, and graph anomaly scoring.

## Overview

The method treats proteins as nodes in a cross-species similarity graph, then scores proteins/components for HGT-like behavior using graph structure and species-pair-normalized edge surprise.

Technically, the project has two parts:
- `graph_construction`: builds candidate protein similarity edges and prunes them into a graph input.
- `hgt_pipeline`: consumes a pruned edge list and produces edge/node/component features plus ranked HGT candidates.

There are two practical entry paths:
- Minimal-input path: start from `data/assembly_summary_refseq.txt` + `config/species.txt`, build manifest/downloads/candidates/pruned edges, then run the pipeline.
- Shortcut path: use the preincluded canonical pruned graph `golden/reference_inputs/edges_PRUNED_JACCARD_92790.tsv` and run `hgt_pipeline` directly.

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

From repository root (PowerShell), install once:

```powershell
python -m pip install -e .
```

Run the canonical pipeline input (preincluded in repo):

With betweenness:

```powershell
python -m hgt_pipeline.pipeline --in_edges golden\reference_inputs\edges_PRUNED_JACCARD_92790.tsv --out_dir tmp_run_bw
```

Without betweenness (faster):

```powershell
python -m hgt_pipeline.pipeline --in_edges golden\reference_inputs\edges_PRUNED_JACCARD_92790.tsv --out_dir tmp_run_no_bw --no_betweenness
```

## Full E2E Recipe (Reviewers)

From repository root (PowerShell):

```powershell
python -m pip install -e .
```

1. Prepare `assembly_summary_refseq.txt` (if missing):

```powershell
curl.exe -L -o data/assembly_summary_refseq.txt https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt
```

2. Build manifest from `config/species.txt` + assembly summary:

```powershell
python -m graph_construction.refseq_fetch_proteins `
  --assembly_summary data/assembly_summary_refseq.txt `
  --species_list config/species.txt `
  --out_dir tmp/e2e/graph_construction `
  --max_assemblies_per_species 2 `
  --require_latest
```

3. Construct candidates and pruned edges:

```powershell
python -m graph_construction.orchestrator construct-edges `
  --manifest tmp/e2e/graph_construction/manifest.tsv `
  --downloads_dir data/out_refseq/downloads `
  --out_candidates tmp/e2e/candidates.tsv `
  --out_edges tmp/e2e/edges_pruned.tsv `
  --k 6 --min_len 50 --max_postings 100 --min_shared 6 --top_m 10 --q 0.9 --top_x 20
```

4. Run HGT pipeline (with and without betweenness):

```powershell
python -m hgt_pipeline.pipeline `
  --in_edges tmp/e2e/edges_pruned.tsv `
  --out_dir tmp/e2e/pipeline_bw

python -m hgt_pipeline.pipeline `
  --in_edges tmp/e2e/edges_pruned.tsv `
  --out_dir tmp/e2e/pipeline_no_bw `
  --no_betweenness
```

5. Reporting outputs:

```powershell
python tools/reporting/top_anomaly_edges.py `
  --edges tmp/e2e/pipeline_bw/edge_features.tsv `
  --top_n 25 `
  --out_dir tmp/e2e/pipeline_bw/results

python tools/reporting/summarize_global_stats.py `
  --component_features tmp/e2e/pipeline_bw/component_features.tsv `
  --protein_features tmp/e2e/pipeline_bw/protein_features.tsv `
  --edge_features tmp/e2e/pipeline_bw/edge_features.tsv `
  --hgt_candidates tmp/e2e/pipeline_bw/hgt_candidates.tsv `
  --out_prefix tmp/e2e/pipeline_bw/results/global_stats
```

6. Explanations:

```powershell
python tools/reporting/explain_component.py --component_id 5 --edges tmp/e2e/pipeline_bw/edge_features.tsv --protein_features tmp/e2e/pipeline_bw/protein_features.tsv --component_features tmp/e2e/pipeline_bw/component_features.tsv --hgt_candidates tmp/e2e/pipeline_bw/hgt_candidates.tsv --top_nodes 20 --top_edges 25
python tools/reporting/explain_component.py --component_id 8 --edges tmp/e2e/pipeline_bw/edge_features.tsv --protein_features tmp/e2e/pipeline_bw/protein_features.tsv --component_features tmp/e2e/pipeline_bw/component_features.tsv --hgt_candidates tmp/e2e/pipeline_bw/hgt_candidates.tsv --top_nodes 20 --top_edges 25
python tools/reporting/explain_component.py --component_id 32 --edges tmp/e2e/pipeline_bw/edge_features.tsv --protein_features tmp/e2e/pipeline_bw/protein_features.tsv --component_features tmp/e2e/pipeline_bw/component_features.tsv --hgt_candidates tmp/e2e/pipeline_bw/hgt_candidates.tsv --top_nodes 20 --top_edges 25
python tools/reporting/explain_top_candidates.py --edges tmp/e2e/pipeline_bw/edge_features.tsv --protein_features tmp/e2e/pipeline_bw/protein_features.tsv --component_features tmp/e2e/pipeline_bw/component_features.tsv --hgt_candidates tmp/e2e/pipeline_bw/hgt_candidates.tsv --top_n 20 --top_k_neighbors 12
```

## Reproduce Component Plots (5, 8, 32)

From repository root (PowerShell):

```powershell
python tools/reporting/plot_components.py `
  --edges golden/bw_pipeline/rerun_pruned/edge_features.tsv `
  --protein_features golden/bw_pipeline/rerun_pruned/protein_features.tsv `
  --hgt_candidates golden/bw_pipeline/rerun_pruned/hgt_candidates.tsv `
  --component_ids 5,8,32 `
  --z_min_highlight 3 `
  --node_size_mode score `
  --out_dir artifacts/updated_plots/detailed
```

## Data policy

Large raw data is not tracked in git (`data/` ignored).  
Regression relies on `golden/` fixtures, including small deterministic test inputs.

`data/assembly_summary_refseq.txt` is intentionally not tracked (very large; ~216 MB).
Fetch it directly from NCBI when needed:
`https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt`

## Citation

Citation metadata is in `CITATION.cff`.
