#!/usr/bin/env python3
"""Identify plasmid groups sharing long high-identity blocks.

Consumes ``all_coords.tsv`` (from all_vs_all_nucmer.py) and identifies cliques
of plasmids that share blocks of at least a user-specified length and percent
identity (default 90%).

For each clique a directory named ``group_XXkb_<n>`` is created. Groups that are
near-duplicates of a larger group (<= --max-missing members missing AND >= --overlap
coverage of the larger group) are recorded as alternatives and named
``group_XXkb_<n>_alt<m>``.

Outputs:
- group_report.tsv: members per group and alt groups
- multi_group_members.tsv: plasmids present in >1 group
- group_gene_counts.tsv: counts of genes (taken after "bla" in FASTA filename)
- network_graph.tsv: node/edge table for the full graph (including singletons)

Optionally: a network plot where nodes are plasmids and edges connect plasmids
sharing a qualifying block. Nodes can be coloured by primary gene or by main group.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import networkx as nx


# ------------------------- Helpers -------------------------

def find_fasta_files(directory: Path) -> List[Path]:
    """Return a sorted list of FASTA files within *directory*."""
    exts = {".fa", ".fasta", ".fna", ".ffn"}
    return sorted(p for p in directory.iterdir() if p.suffix.lower() in exts and p.is_file())


def parse_genes(filename: str) -> List[str]:
    """Extract gene names from a FASTA filename as tokens after 'bla', split by 'and'."""
    stem = Path(filename).stem
    if "bla" not in stem:
        return []
    gene_part = stem.split("bla", 1)[1]
    gene_part = gene_part.strip("_-. ")
    return [tok for tok in gene_part.split("and") if tok]


def _union_length(intervals: List[Tuple[int, int]]) -> int:
    """Return total length of the union of 1-based closed intervals [(s,e), ...]."""
    if not intervals:
        return 0
    ivals = sorted((min(s, e), max(s, e)) for s, e in intervals)
    total = 0
    cur_s, cur_e = ivals[0]
    for s, e in ivals[1:]:
        if s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            total += (cur_e - cur_s + 1)
            cur_s, cur_e = s, e
    total += (cur_e - cur_s + 1)
    return total


# ------------------------- Core logic -------------------------

def build_graph(coords_file: Path, min_len: int, min_id: float) -> nx.Graph:
    """
    Build graph using the UNION of covered intervals per pair (ref,qry).
    Each block must have ident >= min_id to be counted.

    Edge attrs:
      - shared_bp: min(union_ref_bp, union_qry_bp)
      - blocks: number of qualifying blocks
      - mean_id: length-weighted mean identity over qualifying blocks
    """
    ref_intervals = defaultdict(list)   # (ref,qry) -> [(S1,E1), ...]
    qry_intervals = defaultdict(list)   # (ref,qry) -> [(S2,E2), ...]
    id_len_sum = defaultdict(float)     # (ref,qry) -> sum(id * seg_len)
    len_sum = defaultdict(int)          # (ref,qry) -> sum(seg_len)
    block_count = defaultdict(int)

    with coords_file.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ident = float(row["%_IDY"])
            if ident < min_id:
                continue
            ref = row["ref_file"]
            qry = row["query_file"]

            S1, E1 = int(row["S1"]), int(row["E1"])
            S2, E2 = int(row["S2"]), int(row["E2"])
            len1 = int(row["LEN_1"])
            len2 = int(row["LEN_2"])
            seg_len = min(len1, len2)

            key = (ref, qry)
            ref_intervals[key].append((S1, E1))
            qry_intervals[key].append((S2, E2))
            id_len_sum[key] += ident * seg_len
            len_sum[key] += seg_len
            block_count[key] += 1

    G = nx.Graph()
    for key in ref_intervals.keys() | qry_intervals.keys():
        union_ref = _union_length(ref_intervals.get(key, []))
        union_qry = _union_length(qry_intervals.get(key, []))
        shared_bp = min(union_ref, union_qry)
        if shared_bp >= min_len:
            mean_id = (id_len_sum[key] / len_sum[key]) if len_sum[key] else 0.0
            ref, qry = key
            G.add_edge(ref, qry, shared_bp=shared_bp, blocks=block_count[key], mean_id=mean_id)
    return G


def find_groups(
    G: nx.Graph, max_missing: int, overlap: float
) -> List[Dict[str, Iterable[Set[str]]]]:
    """Find maximal cliques and annotate alternative groups."""
    cliques = [set(c) for c in nx.find_cliques(G) if len(c) > 1]
    cliques.sort(key=len, reverse=True)

    groups: List[Dict[str, Iterable[Set[str]]]] = []
    for clique in cliques:
        is_alt = False
        for main in groups:
            main_members: Set[str] = main["members"]
            missing = len(main_members - clique)
            coverage = len(main_members & clique) / len(main_members)
            if missing <= max_missing and coverage >= overlap:
                main.setdefault("alts", []).append(clique)
                is_alt = True
                break
        if not is_alt:
            groups.append({"members": clique, "alts": []})
    return groups


def copy_files(members: Set[str], fasta_dir: Path, dest: Path) -> None:
    """Copy FASTAs for members into dest; skip if dest exists and is non-empty."""
    if dest.exists() and any(dest.iterdir()):
        print(f"[INFO] Skipping copy to {dest} — already exists and not empty.")
        return
    dest.mkdir(parents=True, exist_ok=True)
    for name in members:
        src = fasta_dir / name
        if src.exists():
            shutil.copy(src, dest / name)


def write_reports(
    outdir: Path,
    group_names: Dict[frozenset, str],
    groups: List[Dict[str, Iterable[Set[str]]]],
    genes: Dict[str, List[str]],
    min_length: int,
) -> None:
    length_kb = min_length // 1000
    if (outdir / "group_report.tsv").exists():
        print(f"[INFO] Skipping report generation for {length_kb}kb — outputs already exist.")
        return

    group_report = outdir / "group_report.tsv"
    multi_report = outdir / "multi_group_members.tsv"
    gene_report = outdir / "group_gene_counts.tsv"

    plasmid_to_groups: Dict[str, List[str]] = defaultdict(list)

    with group_report.open("w") as gr, gene_report.open("w") as gg:
        gr.write("group\tsize\tmembers\n")
        gg.write("group\tgene\tcount\n")
        for grp in groups:
            main_name = group_names[frozenset(grp["members"])]
            members = sorted(grp["members"])
            gr.write(f"{main_name}\t{len(members)}\t{','.join(members)}\n")
            for m in members:
                plasmid_to_groups[m].append(main_name)
            counts = Counter(g for m in members for g in genes.get(m, []))
            for gene, count in counts.items():
                gg.write(f"{main_name}\t{gene}\t{count}\n")

            for idx, alt in enumerate(grp["alts"], 1):
                alt_name = f"{main_name}_alt{idx}"
                alt_members = sorted(alt)
                gr.write(f"{alt_name}\t{len(alt_members)}\t{','.join(alt_members)}\n")
                for m in alt_members:
                    plasmid_to_groups[m].append(alt_name)
                counts = Counter(g for m in alt_members for g in genes.get(m, []))
                for gene, count in counts.items():
                    gg.write(f"{alt_name}\t{gene}\t{count}\n")

    with multi_report.open("w") as mr:
        mr.write("plasmid\tgroups\n")
        for plasmid, grps in plasmid_to_groups.items():
            if len(grps) > 1:
                mr.write(f"{plasmid}\t{','.join(sorted(grps))}\n")


# ------------------------- Plotting -------------------------

def plot_graph(
    G: nx.Graph,
    genes: Dict[str, List[str]],
    group_color: Dict[str, str],      # node → color
    group_to_color: Dict[str, str],   # group → color
    group_members: Dict[str, Set[str]],
    args: argparse.Namespace,
) -> None:
    """Plot network; hide singletons; legend only for visible groups."""
    # Keep only nodes with at least one edge
    H = G.subgraph([n for n, deg in G.degree() if deg > 0]).copy()
    if H.number_of_nodes() == 0:
        print("[INFO] No edges to plot (all nodes are singletons).")
        return

    plt.figure(figsize=(12, 10), dpi=300)

    # Kamada–Kawai layout
    pos = nx.kamada_kawai_layout(H)

    # Node coloring
    if args.color_by == "gene":
        palette = plt.get_cmap("tab20")
        gene_names = sorted({(genes.get(n, ["unknown"])[0]) for n in H.nodes})
        denom = max(1, len(gene_names) - 1)
        color_map = {g: palette(i / denom) for i, g in enumerate(gene_names)}
        node_colors = [color_map.get(genes.get(n, ["unknown"])[0], "lightgrey") for n in H.nodes]
        legend_labels = gene_names
        legend_colors = [color_map[g] for g in legend_labels]
    else:
        node_colors = [group_color.get(n, "lightgrey") for n in H.nodes]
        # Legend only for groups that have at least one visible member
        visible_groups = []
        for gname, members in group_members.items():
            if any(m in H for m in members):
                visible_groups.append(gname)
        legend_labels = [g for g in group_to_color.keys() if g in visible_groups]
        legend_colors = [group_to_color[g] for g in legend_labels]

    nx.draw_networkx_edges(H, pos, alpha=0.5)
    nx.draw_networkx_nodes(
        H, pos, node_color=node_colors, node_size=args.node_size,
        edgecolors="black", linewidths=0.8
    )
    if args.show_labels:
        nx.draw_networkx_labels(H, pos, font_size=args.label_size)

    legend_handles = [mpatches.Patch(color=c, label=l) for l, c in zip(legend_labels, legend_colors)]
    if legend_handles:
        plt.legend(
            handles=legend_handles,
            title=args.color_by.capitalize(),
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            borderaxespad=0
        )

    plt.axis("off")

    if args.output_figure:
        plt.savefig(args.output_figure, bbox_inches="tight")
    else:
        plt.show()


# ------------------------- CLI -------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coords", type=Path, help="Path to all_coords.tsv")
    parser.add_argument("fasta_dir", type=Path, help="Directory containing plasmid FASTA files")
    parser.add_argument("--min-length", type=int, required=True, help="Minimum shared block length (bp)")
    parser.add_argument("--min-identity", type=float, default=90.0, help="Minimum percent identity")
    parser.add_argument("--output-dir", type=Path, default=Path("."), help="Directory for results")
    parser.add_argument("--max-missing", type=int, default=20, help="Maximum members missing to call alt group")
    parser.add_argument("--overlap", type=float, default=0.90, help="Fractional overlap to call alt group")
    parser.add_argument("--plot", action="store_true", help="Generate network plot")
    parser.add_argument("--color-by", choices=["gene", "group"], default="gene", help="Colour nodes by gene or group")
    parser.add_argument("--node-size", type=int, default=300, help="Node size for plot")
    parser.add_argument("--label-size", type=int, default=8, help="Label font size")
    parser.add_argument("--show-labels", action="store_true", help="Display node labels in plot")
    parser.add_argument("--output-figure", type=Path, help="Path to save network figure instead of showing")
    args = parser.parse_args()

    # Build graph
    G = build_graph(args.coords, args.min_length, args.min_identity)

    # Ensure all plasmids appear as nodes
    fasta_files = find_fasta_files(args.fasta_dir)
    for f in fasta_files:
        G.add_node(f.name)

    # Parse gene labels
    genes = {f.name: parse_genes(f.name) for f in fasta_files}

    # Find cliques and derive main/alt groups
    groups = find_groups(G, args.max_missing, args.overlap)

    # Assign names + copy files + color by largest membership group
    length_kb = args.min_length // 1000
    group_names: Dict[frozenset, str] = {}
    group_to_color: Dict[str, str] = {}
    group_members: Dict[str, Set[str]] = {}   # group name -> set of members (main only)
    group_color: Dict[str, str] = {}          # node -> chosen color
    node_best_size = defaultdict(int)         # node -> best group size seen

    palette = plt.get_cmap("tab20")

    for idx, grp in enumerate(groups, 1):
        name = f"group_{length_kb}kb_{idx}"
        members = set(grp["members"])
        group_names[frozenset(members)] = name
        group_members[name] = members
        copy_files(members, args.fasta_dir, args.output_dir / name)

        color = palette((idx - 1) % 20)
        group_to_color[name] = color

        size = len(members)
        for m in members:
            if size > node_best_size[m]:
                group_color[m] = color
                node_best_size[m] = size

        # Alternatives: copy and consider for coloring (same color as main group)
        for alt_idx, alt in enumerate(grp["alts"], 1):
            alt_name = f"{name}_alt{alt_idx}"
            alt_members = set(alt)
            copy_files(alt_members, args.fasta_dir, args.output_dir / alt_name)
            alt_size = len(alt_members)
            for m in alt_members:
                if alt_size > node_best_size[m]:
                    group_color[m] = color
                    node_best_size[m] = alt_size

    # Reports
    write_reports(args.output_dir, group_names, groups, genes, args.min_length)

    # === Export full graph with singletons included ===
    network_graph_file = args.output_dir / "network_graph.tsv"
    with network_graph_file.open("w") as out:
        out.write("node1\tnode2\tshared_bp\tblocks\tmean_id\tgroup1\tgroup2\tgene1\tgene2\n")
        for u in G.nodes:
            neighbors = list(G.neighbors(u))
            if neighbors:
                for v in neighbors:
                    if u < v:  # avoid duplicate undirected edges
                        d = G[u][v]
                        group1 = next((g for g, members in group_names.items() if u in members), "")
                        group2 = next((g for g, members in group_names.items() if v in members), "")
                        gene1 = genes.get(u, [""])[0] if genes.get(u) else ""
                        gene2 = genes.get(v, [""])[0] if genes.get(v) else ""
                        out.write(
                            f"{u}\t{v}\t{d.get('shared_bp','')}\t{d.get('blocks','')}\t{d.get('mean_id','')}"
                            f"\t{group1}\t{group2}\t{gene1}\t{gene2}\n"
                        )
            else:
                # singleton — no edges
                group1 = next((g for g, members in group_names.items() if u in members), "")
                gene1 = genes.get(u, [""])[0] if genes.get(u) else ""
                out.write(f"{u}\t\t\t\t\t{group1}\t\t{gene1}\t\n")

    print(f"Full graph with singletons saved to {network_graph_file}")
    # ==========================================

    if args.plot:
        plot_graph(G, genes, group_color, group_to_color, group_members, args)


if __name__ == "__main__":
    main()
