#!/usr/bin/env python3
"""
Plot a plasmid similarity network from a NUCmer all_coords.tsv file.

Key features:
- Reads pairs of contigs/sequences from all_coords.tsv and builds an undirected graph.
- Caches parsed nodes/edges to accelerate subsequent runs (pickle/TSV outputs).
- Colors nodes by: (a) groups from group_report.tsv, (b) genes present in blocks_annotated.tsv,
  or (c) selected Block_IDs presence in blocks_annotated.tsv.
- Supports multi-seed layout runs (-n) with deterministic seeding, saving one image per seed.
- Node coloring supports multiple memberships per node via pie-slice rendering.
- Option to filter and hide small isolated groups (<= threshold size) that have no external connections.

Author: ChatGPT (GPT-5 Thinking)
Date: 2025-08-16
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import pickle
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
import numpy as np
import pandas as pd

# ------------------------------
# Data containers
# ------------------------------

@dataclass
class Edge:
    u: str
    v: str
    block_count: int
    total_len1: int
    total_len2: int
    total_ref_span: int
    total_qry_span: int
    mean_pid: float


# ------------------------------
# Utilities
# ------------------------------

def _hash_for_cache(path: str, extra: str = "") -> str:
    h = hashlib.sha256()
    h.update(os.path.abspath(path).encode())
    if os.path.exists(path):
        h.update(str(os.path.getmtime(path)).encode())
        h.update(str(os.path.getsize(path)).encode())
    if extra:
        h.update(extra.encode())
    return h.hexdigest()[:16]


def ensure_dir(path: str) -> None:
    if path and not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)


# ------------------------------
# Parsing & caching
# ------------------------------

def parse_all_coords(
    tsv_path: str,
    min_len: int = 1,
    min_pid: float = 0.0,
) -> Tuple[Set[str], Dict[Tuple[str, str], Edge]]:
    df = pd.read_csv(tsv_path, sep="\t")
    required_cols = {"REF_TAG", "QRY_TAG", "LEN_1", "LEN_2", "%_IDY", "LEN_R", "LEN_Q"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in all_coords.tsv: {sorted(missing)}")

    df = df[(df["LEN_1"] >= min_len) & (df["%_IDY"] >= min_pid)].copy()

    nodes: Set[str] = set(df["REF_TAG"]).union(df["QRY_TAG"])
    edges: Dict[Tuple[str, str], Edge] = {}

    for (ref, qry), sub in df.groupby(["REF_TAG", "QRY_TAG"], sort=False):
        u, v = sorted((str(ref), str(qry)))
        key = (u, v)
        block_count = len(sub)
        total_len1 = int(sub["LEN_1"].sum())
        total_len2 = int(sub["LEN_2"].sum())
        total_ref_span = int(sub["LEN_R"].sum())
        total_qry_span = int(sub["LEN_Q"].sum())
        mean_pid = float(sub["%_IDY"].mean()) if block_count else 0.0
        if key in edges:
            prev = edges[key]
            edges[key] = Edge(
                u=u,
                v=v,
                block_count=prev.block_count + block_count,
                total_len1=prev.total_len1 + total_len1,
                total_len2=prev.total_len2 + total_len2,
                total_ref_span=prev.total_ref_span + total_ref_span,
                total_qry_span=prev.total_qry_span + total_qry_span,
                mean_pid=(prev.mean_pid + mean_pid) / 2.0,
            )
        else:
            edges[key] = Edge(
                u=u,
                v=v,
                block_count=block_count,
                total_len1=total_len1,
                total_len2=total_len2,
                total_ref_span=total_ref_span,
                total_qry_span=total_qry_span,
                mean_pid=mean_pid,
            )

    return nodes, edges


def load_or_build_cache(
    all_coords: str,
    cache_file: str,
    min_len: int,
    min_pid: float,
    force: bool = False,
) -> Tuple[Set[str], Dict[Tuple[str, str], Edge]]:
    ensure_dir(os.path.dirname(cache_file) or ".")
    meta_key = _hash_for_cache(all_coords, extra=f"minlen={min_len}|minpid={min_pid}")

    if (not force) and os.path.exists(cache_file):
        try:
            with open(cache_file, "rb") as fh:
                payload = pickle.load(fh)
            if payload.get("meta_key") == meta_key:
                return payload["nodes"], payload["edges"]
        except Exception:
            pass

    nodes, edges = parse_all_coords(all_coords, min_len=min_len, min_pid=min_pid)
    with open(cache_file, "wb") as fh:
        pickle.dump({"meta_key": meta_key, "nodes": nodes, "edges": edges}, fh)

    return nodes, edges


# ------------------------------
# Graph build & filtering
# ------------------------------

def build_graph(nodes: Set[str], edges: Dict[Tuple[str, str], Edge], min_blocks: int = 1) -> nx.Graph:
    G = nx.Graph()
    for n in nodes:
        G.add_node(n)
    for (u, v), e in edges.items():
        if e.block_count >= min_blocks:
            G.add_edge(u, v, weight=float(e.total_len1), blocks=int(e.block_count), mean_pid=float(e.mean_pid))
    return G


def prune_small_components(G: nx.Graph, max_size: int) -> Tuple[nx.Graph, List[List[str]]]:
    """Remove connected components with size <= max_size and return pruned graph + list of components removed."""
    if max_size <= 0:
        return G, []
    removed: List[List[str]] = []
    nodes_to_remove: Set[str] = set()
    for comp in nx.connected_components(G):
        comp_nodes = list(comp)
        if len(comp_nodes) <= max_size:
            removed.append(sorted(comp_nodes))
            nodes_to_remove.update(comp_nodes)
    if not nodes_to_remove:
        return G, []
    G2 = G.copy()
    G2.remove_nodes_from(nodes_to_remove)
    return G2, removed


# ------------------------------
# Labels
# ------------------------------

def parse_group_report(path: str) -> Dict[str, Set[str]]:
    if not path:
        return {}
    df = pd.read_csv(path, sep="\t")
    if not set(["group", "members"]).issubset(df.columns):
        raise ValueError("group_report.tsv must have columns: group, members")
    node_to_groups: Dict[str, Set[str]] = defaultdict(set)
    for _, row in df.iterrows():
        grp = str(row["group"])  # keep as-is
        members = str(row["members"]).split(",") if not pd.isna(row["members"]) else []
        for m in members:
            m = m.strip()
            if m:
                node_to_groups[m].add(grp)
    return node_to_groups


def _split_multi_field(val: str) -> List[str]:
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return []
    parts = [p.strip() for p in str(val).split("<>")]
    return [p for p in parts if p]


def parse_blocks_for_labels(
    blocks_path: str,
    target_block_ids: Optional[Set[str]] = None,
    target_genes: Optional[Set[str]] = None,
) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
    if not blocks_path:
        return {}, {}
    df = pd.read_csv(blocks_path, sep="\t")
    required = {"Block_ID", "sequence_name"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"blocks_annotated.tsv missing columns: {sorted(missing)}")

    block_map: Dict[str, Set[str]] = defaultdict(set)
    gene_map: Dict[str, Set[str]] = defaultdict(set)

    target_genes_norm = {g.lower() for g in (target_genes or set())}

    for _, row in df.iterrows():
        node = str(row["sequence_name"]).strip()
        bid = str(row["Block_ID"]).strip()
        if target_block_ids is None or bid in target_block_ids:
            block_map[node].add(bid)
        if "gene" in df.columns:
            genes = _split_multi_field(row.get("gene"))
            if genes:
                if target_genes_norm:
                    for g in genes:
                        if g.lower() in target_genes_norm:
                            gene_map[node].add(g)
                else:
                    for g in genes:
                        gene_map[node].add(g)
    return block_map, gene_map


def generate_color_palette(keys: Sequence[str]) -> Dict[str, tuple]:
    base = plt.get_cmap("tab20").colors
    colors: Dict[str, tuple] = {}
    for i, k in enumerate(sorted(keys)):
        colors[k] = base[i % len(base)]
    return colors


def assemble_node_labels(
    G: nx.Graph,
    color_by: str,
    node_to_groups: Dict[str, Set[str]],
    node_to_blocks: Dict[str, Set[str]],
    node_to_genes: Dict[str, Set[str]],
) -> Tuple[Dict[str, Set[str]], Dict[str, tuple]]:
    node2labels: Dict[str, Set[str]] = defaultdict(set)
    if color_by == "group":
        for n in G.nodes:
            for g in node_to_groups.get(n, set()):
                node2labels[n].add(g)
        label_keys = sorted({g for s in node2labels.values() for g in s})
    elif color_by == "block":
        for n in G.nodes:
            for b in node_to_blocks.get(n, set()):
                node2labels[n].add(b)
        label_keys = sorted({b for s in node2labels.values() for b in s})
    elif color_by == "gene":
        for n in G.nodes:
            for g in node_to_genes.get(n, set()):
                node2labels[n].add(g)
        label_keys = sorted({g for s in node2labels.values() for g in s})
    else:
        label_keys = []
    color_map = generate_color_palette(label_keys)
    return node2labels, color_map


# ------------------------------
# Drawing
# ------------------------------

def draw_pie_nodes(
    ax,
    pos: Dict[str, Tuple[float, float]],
    node2labels: Dict[str, Set[str]],
    label_colors: Dict[str, tuple],
    node_size: int = 300,
    default_color=(0.7, 0.7, 0.7),
    node_edge_color="k",
    linewidths: float = 0.5,
):
    radius = math.sqrt(max(node_size, 10)) / 200.0
    legend_handles: Dict[str, mpatches.Patch] = {}
    for n, (x, y) in pos.items():
        labels = list(node2labels.get(n, []))
        if not labels:
            circ = plt.Circle((x, y), radius=radius, facecolor=default_color, edgecolor=node_edge_color, linewidth=linewidths)
            ax.add_patch(circ)
            continue
        seen = set()
        labels = [l for l in labels if not (l in seen or seen.add(l))][:8]
        total = len(labels)
        start = 0.0
        for lab in labels:
            color = label_colors.get(lab, default_color)
            theta1 = start * 360.0 / total
            theta2 = (start + 1) * 360.0 / total
            wedge = mpatches.Wedge((x, y), r=radius, theta1=theta1, theta2=theta2, facecolor=color, edgecolor=node_edge_color, linewidth=linewidths)
            ax.add_patch(wedge)
            start += 1
            if lab not in legend_handles:
                legend_handles[lab] = mpatches.Patch(color=color, label=lab)
    return list(legend_handles.values())


def draw_graph(
    G: nx.Graph,
    pos: Dict[str, Tuple[float, float]],
    out_png: str,
    node2labels: Dict[str, Set[str]],
    label_colors: Dict[str, tuple],
    with_labels: bool = False,
    node_size: int = 300,
    edge_alpha: float = 0.25,
    title: Optional[str] = None,
):
    fig, ax = plt.subplots(figsize=(12, 10), dpi=200)

    weights = np.array([d.get("weight", 1.0) for _, _, d in G.edges(data=True)], dtype=float)
    if weights.size:
        wmin, wmax = float(weights.min()), float(weights.max())
        if wmax > 0:
            widths = 0.5 + 2.5 * (weights - wmin) / (wmax - wmin + 1e-12)
        else:
            widths = [1.0 for _ in weights]
    else:
        widths = []

    nx.draw_networkx_edges(G, pos, ax=ax, width=widths, alpha=edge_alpha, edge_color="0.3")

    legend_handles = draw_pie_nodes(
        ax,
        pos,
        node2labels=node2labels,
        label_colors=label_colors,
        node_size=node_size,
    )

    if with_labels:
        for n, (x, y) in pos.items():
            ax.text(x, y, s=n, fontsize=6, ha="center", va="center")

    if legend_handles:
        ax.legend(handles=legend_handles, bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.0, frameon=False, title="Legend")

    if title:
        ax.set_title(title)

    ax.set_axis_off()
    ensure_dir(os.path.dirname(out_png) or ".")
    plt.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote plot: {out_png}")


# ------------------------------
# Main
# ------------------------------

def main():
    p = argparse.ArgumentParser(description="Plot plasmid network from all_coords.tsv with group/gene/block coloring and caching.")
    p.add_argument("all_coords", help="Path to all_coords.tsv")
    p.add_argument("outdir", help="Output directory for figures and caches")

    # Caching and parsing
    p.add_argument("--cache-file", default=None, help="Path to cache pickle. Default: <outdir>/allcoords_cache.pkl")
    p.add_argument("--force-recompute", action="store_true", help="Ignore cache and recompute nodes/edges")
    p.add_argument("--min-len", type=int, default=1, help="Minimum LEN_1 to include a block (default: 1)")
    p.add_argument("--min-pid", type=float, default=0.0, help="Minimum %_IDY to include a block (default: 0.0)")
    p.add_argument("--min-blocks", type=int, default=1, help="Minimum aggregated block count to keep an edge (default: 1)")

    # Component pruning
    p.add_argument("--hide-small-groups", type=int, default=0, help="Hide connected components with size <= N (default: 0, keep all)")

    # Coloring sources
    p.add_argument("--group-report", default=None, help="group_report.tsv with columns: group, size, members")
    p.add_argument("--blocks-annotated", default=None, help="blocks_annotated.tsv with Block_ID and sequence_name (and optionally gene)")
    p.add_argument("--block-id", action="append", default=None, help="Block_ID to highlight (can be given multiple times)")
    p.add_argument("--genes", default=None, help="Comma-separated list of gene names to color by (requires --blocks-annotated)")

    # Plotting controls
    p.add_argument("--color-by", choices=["group", "gene", "block", "none"], default="group", help="What to color nodes by")
    p.add_argument("-n", "--num-seeds", type=int, default=1, help="Number of seeded layouts to try; saves one PNG per seed")
    p.add_argument("--layout", choices=["kk", "spring", "fr"], default="kk", help="Graph layout algorithm (default: kk -> spring init + Kamada-Kawai refine)")
    p.add_argument("--node-size", type=int, default=300, help="Node size (default: 300)")
    p.add_argument("--with-labels", action="store_true", help="Draw node labels")
    p.add_argument("--prefix", default="network", help="Output filename prefix (default: network)")

    args = p.parse_args()

    ensure_dir(args.outdir)
    cache_file = args.cache_file or os.path.join(args.outdir, "allcoords_cache.pkl")

    print("[INFO] Loading or computing cache…")
    nodes, edges = load_or_build_cache(
        args.all_coords,
        cache_file=cache_file,
        min_len=args.min_len,
        min_pid=args.min_pid,
        force=args.force_recompute,
    )

    # Persist nodes/edges as TSVs
    nodes_tsv = os.path.join(args.outdir, f"{args.prefix}_nodes.tsv")
    edges_tsv = os.path.join(args.outdir, f"{args.prefix}_edges.tsv")
    pd.DataFrame(sorted(list(nodes)), columns=["node"]).to_csv(nodes_tsv, sep="\t", index=False)
    edges_records = []
    for (u, v), e in edges.items():
        edges_records.append({
            "u": u,
            "v": v,
            "block_count": e.block_count,
            "total_len1": e.total_len1,
            "total_len2": e.total_len2,
            "total_ref_span": e.total_ref_span,
            "total_qry_span": e.total_qry_span,
            "mean_pid": round(e.mean_pid, 3),
        })
    pd.DataFrame(edges_records).to_csv(edges_tsv, sep="\t", index=False)
    print(f"[INFO] Wrote nodes to {nodes_tsv}")
    print(f"[INFO] Wrote edges to {edges_tsv}")

    # Build graph
    G = build_graph(nodes, edges, min_blocks=args.min_blocks)
    print(f"[INFO] Graph built: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges before pruning")

    # Prune small isolated components
    G, removed_components = prune_small_components(G, args.hide_small_groups)
    if removed_components:
        pruned_tsv = os.path.join(args.outdir, f"{args.prefix}_pruned_components.tsv")
        with open(pruned_tsv, "w") as fh:
            fh.write("component_id\tsize\tmembers\n")
            for i, comp in enumerate(removed_components, 1):
                fh.write(f"C{i}\t{len(comp)}\t{','.join(comp)}\n")
        print(f"[INFO] Pruned {len(removed_components)} components (<= {args.hide_small_groups}). Details: {pruned_tsv}")
    print(f"[INFO] Graph after pruning: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    if G.number_of_nodes() == 0:
        print("[WARN] No nodes to plot after pruning/filters.")
        return

    # Labels
    node_to_groups = parse_group_report(args.group_report) if args.group_report else {}
    target_blocks = set(args.block_id) if args.block_id else None
    target_genes = set([g.strip() for g in args.genes.split(",") if g.strip()]) if args.genes else None
    node_to_blocks, node_to_genes = parse_blocks_for_labels(args.blocks_annotated, target_block_ids=target_blocks, target_genes=target_genes)

    node2labels, label_colors = assemble_node_labels(
        G,
        color_by=args.color_by,
        node_to_groups=node_to_groups,
        node_to_blocks=node_to_blocks,
        node_to_genes=node_to_genes,
    )

    # Layouts & plotting
    seeds = list(range(args.num_seeds))
    for i, seed in enumerate(seeds):
        if args.layout == "spring":
            pos = nx.spring_layout(G, seed=seed, weight="weight")
        elif args.layout == "fr":
            pos = nx.fruchterman_reingold_layout(G, seed=seed, weight="weight")
        else:
            init = nx.spring_layout(G, seed=seed, weight="weight")
            pos = nx.kamada_kawai_layout(G, pos=init, weight="weight")

        out_png = os.path.join(args.outdir, f"{args.prefix}_{args.color_by}_seed{i}.png")
        title = f"{args.color_by} | seed={seed} | nodes={G.number_of_nodes()} edges={G.number_of_edges()} (hide<= {args.hide_small_groups})"
        draw_graph(
            G,
            pos,
            out_png=out_png,
            node2labels=node2labels,
            label_colors=label_colors,
            with_labels=args.with_labels,
            node_size=args.node_size,
            title=title,
        )

    # Save the graph object and run settings
    gpickle_path = os.path.join(args.outdir, f"{args.prefix}.gpickle")
    nx.write_gpickle(G, gpickle_path)
    settings = {
        "all_coords": os.path.abspath(args.all_coords),
        "outdir": os.path.abspath(args.outdir),
        "cache_file": os.path.abspath(cache_file),
        "min_len": args.min_len,
        "min_pid": args.min_pid,
        "min_blocks": args.min_blocks,
        "group_report": os.path.abspath(args.group_report) if args.group_report else None,
        "blocks_annotated": os.path.abspath(args.blocks_annotated) if args.blocks_annotated else None,
        "block_id": sorted(list(target_blocks)) if target_blocks else None,
        "genes": sorted(list(target_genes)) if target_genes else None,
        "color_by": args.color_by,
        "num_seeds": args.num_seeds,
        "layout": args.layout,
        "node_size": args.node_size,
        "with_labels": args.with_labels,
        "prefix": args.prefix,
        "hide_small_groups": args.hide_small_groups,
    }
    run_json = os.path.join(args.outdir, f"{args.prefix}_run.json")
    with open(run_json, "w") as fh:
        json.dump(settings, fh, indent=2)
    print(f"[INFO] Saved graph to {gpickle_path}")
    print(f"[INFO] Saved run settings to {run_json}")


if __name__ == "__main__":
    main()
