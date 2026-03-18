#!/usr/bin/env python3
"""
plot_components.py (fixed)

Readable component network plots.

Encodings:
  - Node color: top-N species (distinct colors), all others = light gray ("Other")
  - Node size: HGT score (from hgt_candidates.tsv), percentile-scaled
  - Node outline: top-K scored nodes in the component
  - Edge style:
      * faint gray edges for context
      * highlighted edges for z >= z_min_highlight with width ~ z
  - Legend included (species + edge types)

Outputs:
  out_dir/component_<id>.png
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


def parse_float(x: str) -> Optional[float]:
    x = (x or "").strip()
    if x == "":
        return None
    return float(x)


@dataclass
class NodeRow:
    u: str
    species: str
    component_id: int


@dataclass
class EdgeRow:
    u: str
    v: str
    species_u: str
    species_v: str
    jac: float
    shared: int
    z: Optional[float]


def load_nodes(protein_features: Path) -> Dict[str, NodeRow]:
    out: Dict[str, NodeRow] = {}
    with open(protein_features, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            u = row["u"]
            out[u] = NodeRow(u=u, species=row["species"], component_id=int(row["component_id"]))
    return out


def load_scores(hgt_candidates: Optional[Path]) -> Dict[str, float]:
    if hgt_candidates is None:
        return {}
    out: Dict[str, float] = {}
    with open(hgt_candidates, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            out[row["u"]] = float(row["score"])
    return out


def load_edges_for_component(edges_path: Path, nodes_in_comp: set[str]) -> List[EdgeRow]:
    edges: List[EdgeRow] = []
    with open(edges_path, "r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        fields = set(r.fieldnames or [])
        if "z_robust" not in fields:
            raise ValueError("edges file must include z_robust (use out/edge_features.tsv).")

        for row in r:
            u = row["u"]
            v = row["v"]
            if u not in nodes_in_comp or v not in nodes_in_comp:
                continue
            edges.append(
                EdgeRow(
                    u=u,
                    v=v,
                    species_u=row.get("species_u", ""),
                    species_v=row.get("species_v", ""),
                    jac=float(row["jaccard"]),
                    shared=int(row["shared_kmers"]),
                    z=parse_float(row.get("z_robust", "")),
                )
            )
    return edges


def percentile_rank(xs: List[float], x: float) -> float:
    """Return empirical percentile of x in xs in [0,1]."""
    if not xs:
        return 0.0
    # count <= x
    le = 0
    for v in xs:
        if v <= x:
            le += 1
    return le / len(xs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--edges", type=Path, required=True)
    ap.add_argument("--protein_features", type=Path, required=True)
    ap.add_argument("--hgt_candidates", type=Path, default=None)
    ap.add_argument("--component_ids", type=str, required=True, help="comma-separated, e.g. 5,32")
    ap.add_argument("--out_dir", type=Path, default=Path("figs"))

    # readability knobs
    ap.add_argument("--top_species_colors", type=int, default=8,
                    help="Only top-N species in component get distinct colors; others are gray.")
    ap.add_argument("--top_score_outline", type=int, default=6,
                    help="Outline top-K scored nodes (per component).")
    ap.add_argument("--z_min_highlight", type=float, default=5.0,
                    help="Highlight edges with z_robust >= this threshold.")
    ap.add_argument("--max_highlight_edges", type=int, default=60,
                    help="Cap number of highlighted edges (keeps plot clean).")
    ap.add_argument("--color_mode", type=str, default="highlight",
                    choices=["highlight", "frequency"],
                    help="Species coloring strategy: 'highlight' colors species that appear in highlighted edges; "
                         "'frequency' colors most frequent species.")
    ap.add_argument("--max_colored_species", type=int, default=12,
                    help="Max distinct species colors to show (rest -> Other/gray).")

    ap.add_argument("--fade_non_highlight", action="store_true",
                    help="Fade nodes/edges that are not incident to highlighted (high-z) edges.")
    ap.add_argument("--fade_alpha", type=float, default=0.20,
                    help="Alpha for faded nodes (if --fade_non_highlight).")
    ap.add_argument("--fade_scale", type=float, default=0.55,
                    help="Size multiplier for faded nodes (if --fade_non_highlight).")
    ap.add_argument("--context_edge_alpha", type=float, default=0.06,
                    help="Alpha for background (context) edges.")

    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--figsize", type=str, default="12,10")
    ap.add_argument("--backbone_only", action="store_true",
                    help="Plot only highlighted edges and their incident nodes (high-z backbone).")

    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    nodes = load_nodes(args.protein_features)
    scores = load_scores(args.hgt_candidates)
    have_scores = bool(scores)

    comp_ids = [int(x.strip()) for x in args.component_ids.split(",") if x.strip()]
    comp_to_nodes: Dict[int, List[str]] = defaultdict(list)
    for u, n in nodes.items():
        comp_to_nodes[n.component_id].append(u)

    fig_w, fig_h = (float(x) for x in args.figsize.split(","))

    for cid in comp_ids:
        comp_nodes = comp_to_nodes.get(cid, [])
        if not comp_nodes:
            print(f"[WARN] component {cid} not found.")
            continue
        comp_set = set(comp_nodes)
        edges = load_edges_for_component(args.edges, comp_set)

        # Build graph
        G = nx.Graph()
        for u in comp_nodes:
            sp = nodes[u].species
            sc = scores.get(u, 0.0)
            G.add_node(u, species=sp, score=sc)

        for e in edges:
            G.add_edge(e.u, e.v, z=e.z, jac=e.jac, shared=e.shared)

        # ---- Species coloring: top-N colored, rest gray ----
        sp_list = [nodes[u].species for u in comp_nodes]
        sp_counts = Counter(sp_list)

        # Determine which edges are highlighted first (we need this to choose species)
        edges_with_z = []
        for u, v, d in G.edges(data=True):
            z = d.get("z")
            if z is not None:
                edges_with_z.append((u, v, z))
        edges_with_z.sort(key=lambda t: t[2], reverse=True)

        # Choose highlighted edges either by threshold or top-K (keep your current behavior)
        highlight_edges = []
        for (u, v, z) in edges_with_z:
            if z >= args.z_min_highlight:
                highlight_edges.append((u, v, z))
        # cap
        highlight_edges = highlight_edges[:args.max_highlight_edges]

        if args.color_mode == "highlight":
            # Species that appear in highlighted edges (most relevant for the story)
            sp_in_hl = Counter()
            for u, v, z in highlight_edges:
                sp_in_hl[nodes[u].species] += 1
                sp_in_hl[nodes[v].species] += 1

            # Also include species of outlined nodes (top scored) so they aren't gray
            if have_scores:
                topk = sorted(comp_nodes, key=lambda u: scores.get(u, 0.0), reverse=True)[: args.top_score_outline]
                for u in topk:
                    sp_in_hl[nodes[u].species] += 1

            top_species = [sp for sp, _ in sp_in_hl.most_common(args.max_colored_species)]
        else:
            # Old behavior: most frequent species
            top_species = [sp for sp, _ in sp_counts.most_common(args.max_colored_species)]

        palette = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        sp_to_color: Dict[str, str] = {}
        for i, sp in enumerate(top_species):
            sp_to_color[sp] = palette[i % len(palette)]
        other_color = "#D0D0D0"  # light gray

        node_colors = []
        for u in G.nodes():
            sp = nodes[u].species
            node_colors.append(sp_to_color.get(sp, other_color))

        # ---- Node size: score percentile-scaled ----
        node_sizes = []
        comp_scores = [scores.get(u, 0.0) for u in comp_nodes] if have_scores else []
        for u in G.nodes():
            if have_scores:
                s = scores.get(u, 0.0)
                p = percentile_rank(comp_scores, s)  # 0..1
                # map percentile to size range
                size = 80 + 900 * (p ** 1.7)
            else:
                size = 80 + 20 * G.degree(u)
            node_sizes.append(size)

        # Top-K nodes by score (outline)
        outlined = set()
        if have_scores:
            topk = sorted(comp_nodes, key=lambda u: scores.get(u, 0.0), reverse=True)[: args.top_score_outline]
            outlined = set(topk)

        # ---- Edge styling ----
        # base edges: light gray
        base_edges = list(G.edges())
        base_width = 0.6
        base_alpha = 0.10

        # highlighted edges: z>=threshold, width proportional to z
        hl = [(u, v) for (u, v, z) in highlight_edges]
        hl_widths = [min(6.0, 1.2 + 0.25 * z) for (u, v, z) in highlight_edges]
        hl_nodes = set()
        for (u, v, z) in highlight_edges:
            hl_nodes.add(u)
            hl_nodes.add(v)

        # cap highlighted edges (keep top by z)
        if len(hl) > args.max_highlight_edges:
            hl_sorted = sorted(hl, key=lambda e: (G[e[0]][e[1]].get("z") or -1.0), reverse=True)
            hl = hl_sorted[: args.max_highlight_edges]
            hl_widths = [min(6.0, 1.2 + 0.25 * (G[u][v].get("z") or 0.0)) for (u, v) in hl]

        # ---- Layout ----
        pos = nx.spring_layout(G, seed=args.seed)

        # ---- Draw ----
        plt.figure(figsize=(fig_w, fig_h))

        # edges: base
        nx.draw_networkx_edges(G, pos, edgelist=base_edges, width=base_width, alpha=base_alpha, edge_color="#808080")

        # edges: highlighted
        if hl:
            nx.draw_networkx_edges(G, pos, edgelist=hl, width=hl_widths, alpha=0.9, edge_color="black")

        # nodes: draw in two passes so outlines are clear
        # non-outlined nodes
        def node_size(u: str) -> float:
            idx = list(G.nodes()).index(u)
            return node_sizes[idx]

        def node_col(u: str) -> str:
            return sp_to_color.get(nodes[u].species, other_color)

        non_outline_nodes = [u for u in G.nodes() if u not in outlined]
        outline_nodes = [u for u in G.nodes() if u in outlined]

        if args.fade_non_highlight:
            # keep highlighted-incident nodes normal; fade the rest
            non_outline_hi = [u for u in non_outline_nodes if u in hl_nodes]
            non_outline_lo = [u for u in non_outline_nodes if u not in hl_nodes]
        else:
            non_outline_hi = non_outline_nodes
            non_outline_lo = []

        # draw faded non-highlight nodes
        if non_outline_lo:
            nx.draw_networkx_nodes(
                G, pos,
                nodelist=non_outline_lo,
                node_color=[node_col(u) for u in non_outline_lo],
                node_size=[node_size(u) * args.fade_scale for u in non_outline_lo],
                alpha=args.fade_alpha,
                linewidths=0.0,
            )

        # draw normal non-outlined nodes
        if non_outline_hi:
            nx.draw_networkx_nodes(
                G, pos,
                nodelist=non_outline_hi,
                node_color=[node_col(u) for u in non_outline_hi],
                node_size=[node_size(u) for u in non_outline_hi],
                alpha=1.0,
                linewidths=0.0,
            )

        # outlined nodes always drawn last, always visible
        if outline_nodes:
            nx.draw_networkx_nodes(
                G, pos,
                nodelist=outline_nodes,
                node_color=[node_col(u) for u in outline_nodes],
                node_size=[node_size(u) for u in outline_nodes],
                alpha=1.0,
                linewidths=1.5,
                edgecolors="black",
            )

        # outlined nodes
        outline_nodes = [u for u in G.nodes() if u in outlined]
        if outline_nodes:
            nx.draw_networkx_edges(G, pos, edgelist=base_edges, width=base_width,
                                   alpha=args.context_edge_alpha, edge_color="#808080")

        # Title
        n = G.number_of_nodes()
        m = G.number_of_edges()
        dens = (2.0 * m) / (n * (n - 1.0)) if n > 1 else 0.0
        title = f"Component {cid} | nodes={n} edges={m} density={dens:.3f}"
        if have_scores:
            title += f" | node size = score (percentile), outlined = top-{args.top_score_outline}"
        title += f" | bold edges: z ≥ {args.z_min_highlight}"
        plt.title(title)
        plt.axis("off")

        # ---- Legend ----
        handles: List = []
        # species patches for top species only + other
        for sp in top_species:
            handles.append(Patch(facecolor=sp_to_color[sp], edgecolor="none", label=sp))
        handles.append(Patch(facecolor=other_color, edgecolor="none", label="Other species"))

        # edge legend
        handles.append(Line2D([0], [0], color="#808080", lw=2, alpha=0.25, label="All edges (context)"))
        handles.append(Line2D([0], [0], color="black", lw=3, alpha=0.9, label=f"High-surprise edges (z ≥ {args.z_min_highlight})"))

        plt.legend(handles=handles, loc="upper left", frameon=True, fontsize=9)

        out_path = args.out_dir / f"component_{cid}.png"
        plt.tight_layout()
        plt.savefig(out_path, dpi=240)
        plt.close()
        print(f"[OK] wrote {out_path}")


if __name__ == "__main__":
    main()
