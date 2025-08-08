#!/usr/bin/env python3
"""Identify plasmid cliques from ``all_vs_all_nucmer.py`` output.

The script reads the ``all_coords.tsv`` table produced by ``all_vs_all_nucmer.py``
and groups plasmids that share alignment blocks of at least a user supplied
length and percent identity.  Each group represents a clique in which every pair
of plasmids shares a block meeting the criteria.

For every clique a folder named ``group_XXkb_YY`` is created containing copies of
all member FASTA files.  Two reports are generated in the current directory:

``groups_report.tsv``
    One row per group with its size and member plasmids.
``multi_membership_report.tsv``
    Plasmids that belong to more than one group and their groups.

Optionally a network graph of plasmid relationships can be saved as HTML using
``mpld3`` or as a PNG if ``mpld3`` is unavailable.
"""
from __future__ import annotations

import argparse
import csv
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import networkx as nx

try:  # Optional imports for plotting
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover - matplotlib may be unavailable
    plt = None


def parse_coords(coords_file: Path, min_len: int, min_identity: float) -> Dict[Tuple[str, str], int]:
    """Return edges (pairs of plasmids) meeting the criteria.

    The maximum qualifying alignment length for each plasmid pair is retained."""
    edges: Dict[Tuple[str, str], int] = {}
    with coords_file.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                idy = float(row["%_IDY"])
                len1 = int(row["LEN_1"])
                len2 = int(row["LEN_2"])
            except KeyError as e:
                raise ValueError(f"Missing expected column in coordinates file: {e}")
            aln_len = min(len1, len2)
            if idy < min_identity or aln_len < min_len:
                continue
            pair = tuple(sorted((row["ref_file"], row["query_file"])))
            edges[pair] = max(edges.get(pair, 0), aln_len)
    return edges


def build_graph(edges: Dict[Tuple[str, str], int]) -> nx.Graph:
    """Create an undirected graph from qualifying edges."""
    g = nx.Graph()
    for (a, b), length in edges.items():
        g.add_edge(a, b, length=length)
    return g


def write_reports(groups: List[Tuple[str, List[str]]]) -> None:
    """Write summary reports for groups and multi-membership plasmids."""
    with open("groups_report.tsv", "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["group", "size", "members"])
        for name, members in groups:
            writer.writerow([name, len(members), ",".join(members)])

    membership: Dict[str, List[str]] = defaultdict(list)
    for name, members in groups:
        for m in members:
            membership[m].append(name)

    with open("multi_membership_report.tsv", "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["plasmid", "groups"])
        for plasmid, grps in membership.items():
            if len(grps) > 1:
                writer.writerow([plasmid, ",".join(grps)])


def copy_group_fastas(groups: List[Tuple[str, List[str]]], fasta_dir: Path) -> None:
    """Copy member FASTA files into group specific folders."""
    for name, members in groups:
        outdir = Path(name)
        outdir.mkdir(exist_ok=True)
        for member in members:
            src = fasta_dir / member
            if src.exists():
                shutil.copy(src, outdir / member)
            else:  # pragma: no cover - depends on user input
                raise FileNotFoundError(f"FASTA file not found for {member}: {src}")


def plot_graph(g: nx.Graph, k: float, iterations: int) -> None:
    """Plot the plasmid relationship graph using matplotlib."""
    if plt is None:  # pragma: no cover - plotting optional
        raise RuntimeError("matplotlib is required for plotting")

    fig, ax = plt.subplots()
    pos = nx.spring_layout(g, k=k, iterations=iterations)
    nx.draw_networkx(g, pos, ax=ax, node_size=500, font_size=8)
    ax.set_axis_off()
    try:  # HTML output if mpld3 is available
        import mpld3

        mpld3.save_html(fig, "network_graph.html")
    except Exception:  # pragma: no cover - mpld3 may be missing
        fig.savefig("network_graph.png")
    finally:
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coords_file", type=Path, help="Path to all_coords.tsv")
    parser.add_argument("fasta_dir", type=Path, help="Directory containing FASTA files")
    parser.add_argument("--min-length", type=int, required=True, help="Minimum shared block length (bp)")
    parser.add_argument(
        "--min-identity", type=float, default=90.0, help="Minimum percent identity to consider"
    )
    parser.add_argument("--plot", action="store_true", help="Plot plasmid relationship graph")
    parser.add_argument("--spring-k", type=float, default=None, help="Repulsion parameter for spring layout")
    parser.add_argument("--spring-iterations", type=int, default=50, help="Iterations for spring layout")
    args = parser.parse_args()

    edges = parse_coords(args.coords_file, args.min_length, args.min_identity)
    graph = build_graph(edges)

    cliques = [clq for clq in nx.find_cliques(graph) if len(clq) > 1]

    length_kb = int(round(args.min_length / 1000))
    groups: List[Tuple[str, List[str]]] = []
    for idx, members in enumerate(cliques, start=1):
        group_name = f"group_{length_kb:02d}kb_{idx:02d}"
        groups.append((group_name, sorted(members)))

    copy_group_fastas(groups, args.fasta_dir)
    write_reports(groups)

    if args.plot and graph.number_of_edges() > 0:
        plot_graph(graph, args.spring_k, args.spring_iterations)


if __name__ == "__main__":
    main()
