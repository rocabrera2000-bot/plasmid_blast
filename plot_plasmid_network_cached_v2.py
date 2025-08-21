#!/usr/bin/env python3
from __future__ import annotations
"""
Plot a plasmid similarity network from nucmer all_coords.tsv (same unbiased edge logic),
with caching and flexible coloring (group/gene/block/meta TSV). Pie-slice nodes show
multiple labels. Meta TSV mode keeps top-N categories; the rest collapse to WHITE 'minor'.
"""

import argparse, csv, hashlib, json, logging
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict, List, Set, Tuple, Callable, Optional

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
import re

# -------------------- logging --------------------

def setup_logging(verbose: bool, quiet: bool):
    level = logging.INFO
    if verbose: level = logging.DEBUG
    if quiet: level = logging.WARNING
    logging.basicConfig(level=level, format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")

# -------------------- name normalization --------------------

def norm_exact(name: str) -> str: return name
def norm_basename(name: str) -> str: return Path(name).name
def norm_stem(name: str) -> str: return Path(name).stem

def make_normalizer(mode: str) -> Callable[[str], str]:
    return {"exact": norm_exact, "basename": norm_basename, "stem": norm_stem}[mode]

# -------------------- helpers --------------------

def canon_pair(a: str, b: str) -> Tuple[str,str]: return (a, b) if a <= b else (b, a)

def parse_genes_from_name(filename: str) -> List[str]:
    stem = Path(filename).stem
    if "bla" not in stem: return []
    gene_part = stem.split("bla", 1)[1].strip("_-. ")
    return [tok for tok in gene_part.split("and") if tok]

def merge_union_weighted(intervals_with_id: List[Tuple[int,int,float]]) -> Tuple[int, float, int]:
    if not intervals_with_id: return 0, 0.0, 0
    ivals = [(min(s,e), max(s,e), float(iden)) for s,e,iden in intervals_with_id]
    ivals.sort(key=lambda x: (-x[2], x[0], x[1]))
    covered: List[Tuple[int,int,float]] = []
    for s, e, ident in ivals:
        a, b = s, e
        new_parts = []
        for cs, ce, _ in covered:
            if ce < a or cs > b: continue
            if a < cs: new_parts.append((a, cs - 1))
            a = max(a, ce + 1)
            if a > b: break
        if a <= b: new_parts.append((a, b))
        for x, y in new_parts: covered.append((x, y, ident))
    covered.sort(key=lambda x: (x[0], x[1]))
    union_len = 0; id_len_sum = 0.0; pieces = 0
    for s, e, ident in covered:
        seglen = e - s + 1
        union_len += seglen
        id_len_sum += ident * seglen
        pieces += 1
    return union_len, id_len_sum, pieces

# -------------------- read coords and build edges --------------------

def read_coords_build_edges(coords_path: Path, min_identity: float, edge_min_len: int,
                            log_interval: int, name_norm: Callable[[str], str]):
    per_side_with_id: Dict[Tuple[str,str], Dict[str, List[Tuple[int,int,float]]]] = defaultdict(lambda: defaultdict(list))
    all_nodes: Set[str] = set()
    total = kept = 0
    with coords_path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            total += 1
            if total % log_interval == 0:
                logging.info(f"[coords] read {total:,} rows… kept {kept:,} (≥{min_identity}% id)")
            ident = float(row["%_IDY"])
            if ident < min_identity: continue
            kept += 1
            ref = name_norm(row["ref_file"]); qry = name_norm(row["query_file"])
            s1, e1 = int(row["S1"]), int(row["E1"]); s2, e2 = int(row["S2"]), int(row["E2"])
            u, v = canon_pair(ref, qry)
            per_side_with_id[(u, v)][ref].append((s1, e1, ident))
            per_side_with_id[(u, v)][qry].append((s2, e2, ident))
            all_nodes.add(ref); all_nodes.add(qry)
    edges = []
    for (u, v), sides in per_side_with_id.items():
        u_len, u_idsum, u_pieces = merge_union_weighted(sides.get(u, []))
        v_len, v_idsum, v_pieces = merge_union_weighted(sides.get(v, []))
        shared_bp = min(u_len, v_len)
        if shared_bp < edge_min_len: continue
        if u_len == 0 and v_len == 0:
            mean_id, blocks = 0.0, 0
        elif u_len <= v_len:
            mean_id, blocks = (u_idsum / max(1, u_len)), u_pieces
        else:
            mean_id, blocks = (v_idsum / max(1, v_len)), v_pieces
        edges.append({"u": u, "v": v, "shared_bp": int(shared_bp), "mean_id": float(mean_id), "blocks": int(blocks)})
    logging.info(f"[graph] nodes: {len(all_nodes):,}; edges kept (≥{edge_min_len} bp): {len(edges):,}")
    return sorted(all_nodes), edges

# -------------------- cache --------------------

def cache_key(coords_path: Path, min_identity: float, edge_min_len: int, name_mode: str) -> str:
    st = coords_path.stat()
    payload = f"{coords_path.resolve()}|{st.st_size}|{int(st.st_mtime)}|{min_identity}|{edge_min_len}|{name_mode}"
    return hashlib.sha1(payload.encode()).hexdigest()

def load_or_build_cache(coords_path: Path, min_identity: float, edge_min_len: int,
                        log_interval: int, cache_path: Path, recompute: bool,
                        name_norm: Callable[[str], str], name_mode: str):
    ck = cache_key(coords_path, min_identity, edge_min_len, name_mode)
    if cache_path.exists() and not recompute:
        with cache_path.open() as fh: blob = json.load(fh)
        if blob.get("key") == ck:
            logging.info(f"[cache] using cached edges: {cache_path}")
            return blob["nodes"], blob["edges"]
        else:
            logging.info("[cache] cache mismatch; rebuilding…")
    nodes, edges = read_coords_build_edges(coords_path, min_identity, edge_min_len, log_interval, name_norm)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with cache_path.open("w") as fh: json.dump({"key": ck, "nodes": nodes, "edges": edges}, fh)
    logging.info(f"[cache] wrote: {cache_path}")
    return nodes, edges

# -------------------- groups / blocks / genes --------------------

ALT_RE = re.compile(r"(.+?)_alt\d+$")
def canonical_group(g: str) -> str:
    m = ALT_RE.match(g); return m.group(1) if m else g

def read_groups(groups_path: Path, name_norm):
    node2groups: Dict[str, List[str]] = defaultdict(list)
    group_to_canon: Dict[str, str] = {}
    canon_group_sizes: Dict[str, int] = defaultdict(int)
    with groups_path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            raw = row["group"].strip()
            size = int(row["size"])
            canon = canonical_group(raw)
            group_to_canon[raw] = canon
            if size > canon_group_sizes[canon]: canon_group_sizes[canon] = size
            members = [name_norm(m.strip()) for m in row["members"].split(",") if m.strip()]
            for m in members: node2groups[m].append(raw)
    return node2groups, group_to_canon, dict(canon_group_sizes)

def read_blocks_have_nodes(blocks_path: Path, block_ids: List[str], name_norm: Callable[[str], str]) -> Dict[str, Set[str]]:
    wanted = set(block_ids); out: Dict[str, Set[str]] = {bid: set() for bid in wanted}
    with blocks_path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            bid = row["Block_ID"]
            if bid in wanted:
                node = name_norm(row["sequence_name"])
                out[bid].add(node)
    return out

def list_fasta_files(directory: Path) -> List[Path]:
    exts = {".fa", ".fasta", ".fna", ".ffn"}
    return sorted(p for p in directory.iterdir() if p.suffix.lower() in exts and p.is_file())

def build_node_genes(fasta_dir: Path, name_norm: Callable[[str], str]) -> Dict[str, List[str]]:
    genes = {}
    for f in list_fasta_files(fasta_dir):
        genes[name_norm(f.name)] = parse_genes_from_name(f.name)
    return genes

# -------------------- NEW: metadata TSV coloring --------------------

def read_meta_tsv(meta_path: Path, field: str, key_col: str, name_norm: Callable[[str], str]) -> Dict[str, str]:
    """
    Returns node -> category string for the chosen field.
    """
    node2meta: Dict[str, str] = {}
    with meta_path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        if key_col not in rdr.fieldnames:
            raise ValueError(f"[meta] key column '{key_col}' not found in {meta_path}. Columns: {rdr.fieldnames}")
        if field not in rdr.fieldnames:
            raise ValueError(f"[meta] field '{field}' not found in {meta_path}. Columns: {rdr.fieldnames}")
        for row in rdr:
            k = name_norm(row[key_col])
            v = row[field].strip()
            if k: node2meta[k] = v if v else ""
    return node2meta

def top_n_categories_for_nodes(node2meta: Dict[str,str], nodes: List[str], top_n: int) -> Set[str]:
    """
    Count categories among 'nodes' and return the top-N category names.
    """
    cnt = Counter(node2meta.get(n, "") for n in nodes if node2meta.get(n, "") != "")
    return set([lab for lab, _ in cnt.most_common(top_n)]) if top_n > 0 else set(cnt.keys())

# -------------------- graph assembly & filtering --------------------

def build_graph_from_cache(nodes: List[str], edges: List[dict]) -> nx.Graph:
    G = nx.Graph()
    for n in nodes: G.add_node(n)
    for e in edges:
        G.add_edge(e["u"], e["v"], shared_bp=e["shared_bp"], mean_id=e["mean_id"], blocks=e["blocks"])
    return G

def filter_small_isolated_components(G: nx.Graph, max_size: int) -> nx.Graph:
    if max_size <= 0: return G.copy()
    H = G.copy()
    for comp in list(nx.connected_components(H)):
        if len(comp) <= max_size: H.remove_nodes_from(comp)
    return H

# -------------------- coloring logic --------------------

def choose_color_palette():
    cmap = plt.get_cmap("tab20")
    return [cmap(i / 19.0) for i in range(20)]

def color_map_for_labels(all_labels: List[str]) -> Dict[str, tuple]:
    palette = choose_color_palette()
    color_map = {}
    for i, lab in enumerate(sorted(all_labels)):
        color_map[lab] = palette[i % len(palette)]
    return color_map

def node_label_bundles(
    nodes: List[str],
    color_sources: List[str],
    node2groups: Dict[str, List[str]],
    node2genes: Dict[str, List[str]],
    block2nodes: Dict[str, Set[str]],
    group_to_canon: Dict[str, str] = None,
    canon_group_sizes: Dict[str, int] = None,
    group_major_min: int = 10,
    gene_mode: str = "first",
    # meta
    node2meta: Optional[Dict[str,str]] = None,
    meta_top_labels: Optional[Set[str]] = None,
) -> Dict[str, List[str]]:
    group_to_canon = group_to_canon or {}
    canon_group_sizes = canon_group_sizes or {}
    out: Dict[str, List[str]] = {}

    blocks_only = (set(color_sources) == {"block"})
    meta_only   = (set(color_sources) == {"meta"})

    for n in nodes:
        labs: List[str] = []

        # groups
        if "group" in color_sources:
            raw_grps = node2groups.get(n, [])
            canon_grps = list(dict.fromkeys([group_to_canon.get(g, g) for g in raw_grps]))
            majors = [g for g in canon_grps if canon_group_sizes.get(g, 0) >= group_major_min]
            if majors: labs.extend(majors)
            elif canon_grps: labs.append("__MINOR_GROUP__")

        # genes
        if "gene" in color_sources:
            genes = node2genes.get(n, [])
            if genes:
                if gene_mode == "all": labs.extend(genes)
                else: labs.append(genes[0])

        # blocks
        if "block" in color_sources and block2nodes:
            block_hits = [bid for bid, members in block2nodes.items() if n in members]
            if block_hits: labs.extend(block_hits)
            elif blocks_only: labs.append("__NO_BLOCK__")

        # meta
        if "meta" in color_sources and node2meta is not None:
            cat = node2meta.get(n, "")
            if cat and (meta_top_labels is None or cat in meta_top_labels):
                labs.append(f"meta:{cat}")
            else:
                # collapse to minor if not in top-N or missing
                labs.append("__MINOR_META__")

        labs = list(dict.fromkeys(labs))
        out[n] = labs if labs else ["unlabeled"]
    return out

# -------------------- drawing --------------------

def draw_pie_nodes(ax, pos: Dict[str, Tuple[float,float]],
                   node2labels: Dict[str, List[str]],
                   label2color: Dict[str, tuple],
                   node_size: int, edgecolor="black", linewidth=0.8):
    max_slices = 6
    for n, (x, y) in pos.items():
        labs = node2labels.get(n, ["unlabeled"])
        if len(labs) > max_slices:
            labs = labs[:max_slices-1] + ["…"]
            if "…" not in label2color:
                label2color["…"] = (0.7,0.7,0.7,1.0)
        r = (node_size ** 0.5) / 200.0
        theta0 = 90.0
        frac = 1.0 / len(labs)
        for lab in labs:
            col = label2color.get(lab, (0.85,0.85,0.85,1.0))
            wedge = mpatches.Wedge(center=(x, y), r=r, theta1=theta0, theta2=theta0 + 360.0*frac,
                                   facecolor=col, edgecolor=edgecolor, linewidth=linewidth)
            ax.add_patch(wedge)
            theta0 += 360.0 * frac

# -------------------- plot --------------------

def plot_network(
    G: nx.Graph,
    node2labels: Dict[str, List[str]],
    label2color: Dict[str, tuple],
    node_size: int,
    label_font: int,
    show_labels: bool,
    layout_seed: int,
    output_path: Path,
    title: str,
    legend_max: int = 0,
    legend_ncol: int = 1,
):
    if G.number_of_edges() == 0 and G.number_of_nodes() == 0:
        logging.warning("[plot] nothing to plot"); return

    H = G.subgraph([n for n, deg in G.degree() if deg > 0]).copy()
    if H.number_of_nodes() == 0: H = G.copy()

    # Seeded, reproducible layout: spring init -> KK refine
    #pos0 = nx.spring_layout(H, seed=layout_seed, dim=2)
    pos = nx.kamada_kawai_layout(H)

    fig, ax = plt.subplots(figsize=(12, 10), dpi=300)
    nx.draw_networkx_edges(H, pos, alpha=0.5, width=1.0, ax=ax)
    draw_pie_nodes(ax, pos, node2labels, label2color, node_size=node_size)
    if show_labels: nx.draw_networkx_labels(H, pos, font_size=label_font, ax=ax)

    present_labels = sorted({lab for labs in node2labels.values() for lab in labs if lab in label2color})
    show = present_labels if legend_max == 0 else present_labels[:legend_max]
    handles = [plt.Line2D([0],[0], marker='o', linestyle='', markersize=8,
                          markerfacecolor=label2color[l], markeredgecolor='black', label=l)
               for l in show]
    if handles:
        ax.legend(handles=handles, title="Labels", loc="center left", bbox_to_anchor=(1.02,0.5),
                  borderaxespad=0, fontsize=8, ncol=legend_ncol)

    ax.set_title(title); ax.axis("off"); fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight"); plt.close(fig)
    logging.info(f"[plot] saved: {output_path}")

# -------------------- debug dump --------------------

def dump_name_debug(outdir: Path, nodes_from_coords: List[str],
                    node2groups: Dict[str, List[str]], node2labels: Dict[str, List[str]]):
    outdir.mkdir(parents=True, exist_ok=True)
    coords_set = set(nodes_from_coords); groups_set = set(node2groups.keys())
    inter = coords_set & groups_set
    unmapped = sorted(n for n in nodes_from_coords if n not in node2groups)
    (outdir / "debug_nodes_from_coords.txt").write_text("\n".join(sorted(coords_set)) + "\n")
    (outdir / "debug_nodes_from_groups.txt").write_text("\n".join(sorted(groups_set)) + "\n")
    (outdir / "debug_nodes_intersection.txt").write_text("\n".join(sorted(inter)) + "\n")
    (outdir / "debug_nodes_unmapped_in_groups.txt").write_text("\n".join(unmapped) + "\n")
    with (outdir / "debug_groups_per_node.tsv").open("w") as fh:
        fh.write("node\tgroups\n")
        for n in sorted(coords_set):
            grps = ",".join(sorted(node2groups.get(n, [])))
            fh.write(f"{n}\t{grps}\n")
    logging.info("[debug] coords nodes: %d | grouped nodes: %d | intersection: %d | unmapped: %d",
                 len(coords_set), len(groups_set), len(inter), len(unmapped))
    for n in unmapped[:10]:
        logging.info("[debug] unmapped example: %s (labels=%s)", n, node2labels.get(n, []))

# -------------------- CLI --------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("coords", type=Path, help="Path to nucmer all_coords.tsv")
    ap.add_argument("fasta_dir", type=Path, help="Directory with plasmid FASTA files (for gene parsing)")
    ap.add_argument("--output-dir", type=Path, default=Path("plots"), help="Directory for outputs")
    # edge logic
    ap.add_argument("--min-identity", type=float, default=90.0, help="Minimum %% identity for nucmer rows")
    ap.add_argument("--edge-min-length", type=int, required=True, help="Minimum shared length (bp) to keep an edge")
    ap.add_argument("--log-interval", type=int, default=100_000)
    # caching
    ap.add_argument("--cache", type=Path, default=None, help="Path to cache JSON (defaults to output-dir/edge_cache.json)")
    ap.add_argument("--recompute", action="store_true", help="Ignore cache and recompute edges")
    # color sources
    ap.add_argument("--groups-file", type=Path, help="group_report.tsv with columns: group size members")
    ap.add_argument("--blocks-file", type=Path, help="blocks_annotated.tsv")
    ap.add_argument("--block-id", action="append", default=[], help="Block_ID to highlight (repeatable)")
    ap.add_argument("--color-sources", default="group", help="Comma-separated subset of: group,gene,block,meta")
    ap.add_argument("--gene-mode", choices=["first","all"], default="first",
                    help="When coloring by gene, use only the first gene or all genes (pie-slices)")
    ap.add_argument("--group-major-min", type=int, default=10,
                    help="Groups with >= this many members are 'major' (colored); others become MINOR (white).")
    # meta TSV
    ap.add_argument("--meta-file", type=Path, help="TSV with per-node metadata")
    ap.add_argument("--meta-field", type=str, help="Which column in --meta-file to color by")
    ap.add_argument("--meta-key-col", type=str, default="sequence_name", help="Column with node IDs in --meta-file")
    ap.add_argument("--meta-top", type=int, default=10, help="Keep only top-N categories for meta; others become MINOR (white)")
    # plotting
    ap.add_argument("-n", "--num-layouts", type=int, default=1, help="How many seeds to try")
    ap.add_argument("--seed-base", type=int, default=1, help="Base seed; seeds are seed_base..seed_base+N-1")
    ap.add_argument("--node-size", type=int, default=300)
    ap.add_argument("--label-size", type=int, default=8)
    ap.add_argument("--show-labels", action="store_true")
    ap.add_argument("--title", type=str, default="", help="Custom plot title prefix")
    # filtering
    ap.add_argument("--hide-small-components", type=int, default=0,
                    help="Hide connected components of size <= this number")
    # name matching & debug
    ap.add_argument("--name-match-mode", choices=["stem","basename","exact"], default="stem",
                    help="How to normalize names before matching between coords, groups, blocks, and meta/gene files")
    ap.add_argument("--debug-dump-names", action="store_true", help="Dump name-mapping diagnostics to output-dir")
    # legend controls
    ap.add_argument("--legend-max", type=int, default=0, help="Max labels to show in legend (0 = show all)")
    ap.add_argument("--legend-ncol", type=int, default=1, help="Legend columns")
    # logging verbosity
    ap.add_argument("--verbose", action="store_true", help="Enable debug logging")
    ap.add_argument("--quiet", action="store_true", help="Minimal logging")

    args = ap.parse_args()
    setup_logging(args.verbose, args.quiet)

    name_norm = make_normalizer(args.name_match_mode)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cache_path = args.cache if args.cache else (args.output_dir / "edge_cache.json")

    # 1) Load or compute edges/nodes
    nodes, edges = load_or_build_cache(
        coords_path=args.coords, min_identity=args.min_identity, edge_min_len=args.edge_min_length,
        log_interval=args.log_interval, cache_path=cache_path, recompute=args.recompute,
        name_norm=name_norm, name_mode=args.name_match_mode
    )

    # 2) Build G and apply component filter
    G = build_graph_from_cache(nodes, edges)
    if args.hide_small_components > 0:
        G = filter_small_isolated_components(G, args.hide_small_components)
        nodes = sorted(G.nodes())

    # 3) Load metadata sources
    node2groups, group_to_canon, canon_group_sizes = (defaultdict(list), {}, {})
    if args.groups_file and args.groups_file.exists():
        node2groups, group_to_canon, canon_group_sizes = read_groups(args.groups_file, name_norm)
    node2genes = build_node_genes(args.fasta_dir, name_norm)
    block2nodes: Dict[str, Set[str]] = {}
    if args.blocks_file and args.block_id:
        block2nodes = read_blocks_have_nodes(args.blocks_file, args.block_id, name_norm)

    node2meta: Optional[Dict[str,str]] = None
    meta_top_labels: Optional[Set[str]] = None
    if args.meta_file:
        if not args.meta_field:
            raise SystemExit("--meta-field is required when --meta-file is provided")
        node2meta = read_meta_tsv(args.meta_file, args.meta_field, args.meta_key_col, name_norm)
        meta_top_labels = top_n_categories_for_nodes(node2meta, nodes, args.meta_top)

    # 4) Labels per node
    color_sources = [t.strip() for t in args.color_sources.split(",") if t.strip()]
    node2labels = node_label_bundles(
        nodes, color_sources, node2groups, node2genes, block2nodes,
        group_to_canon=group_to_canon, canon_group_sizes=canon_group_sizes,
        group_major_min=args.group_major_min, gene_mode=args.gene_mode,
        node2meta=node2meta, meta_top_labels=meta_top_labels,
    )

    if args.debug_dump_names:
        dump_name_debug(outdir=args.output_dir, nodes_from_coords=nodes,
                        node2groups=node2groups, node2labels=node2labels)

    # 5) Color map (override special labels)
    all_labels = sorted({lab for labs in node2labels.values() for lab in labs})
    label2color = color_map_for_labels(all_labels)

    # White for “minor”/missing buckets
    specials_present = {lab for labs in node2labels.values() for lab in labs}
    if "__MINOR_GROUP__" in specials_present: label2color["__MINOR_GROUP__"] = (1.0, 1.0, 1.0, 1.0)
    if "__NO_BLOCK__"   in specials_present: label2color["__NO_BLOCK__"]   = (1.0, 1.0, 1.0, 1.0)
    if "__MINOR_META__" in specials_present: label2color["__MINOR_META__"] = (1.0, 1.0, 1.0, 1.0)

    # 6) Plot for seeds
    for i in range(args.num_layouts):
        seed = args.seed_base + i
        suffix = f"seed{seed}"
        title = (args.title + " " if args.title else "") + f"(seed {seed}, color={'+'.join(color_sources)})"
        out_png = args.output_dir / f"network_{'+'.join(color_sources)}_{suffix}.png"
        plot_network(G=G, node2labels=node2labels, label2color=label2color,
                     node_size=args.node_size, label_font=args.label_size, show_labels=args.show_labels,
                     layout_seed=seed, output_path=out_png, title=title,
                     legend_max=args.legend_max, legend_ncol=args.legend_ncol)

if __name__ == "__main__":
    main()
