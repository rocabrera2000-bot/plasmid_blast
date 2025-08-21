#!/usr/bin/env python3
"""
Identify plasmid groups (maximal cliques) and export per-block CORE/ACCESSORY
segments.

CORE: anchor sweep + projection (unchanged from your accepted logic).
ACCESSORY: anchor-seeded clustering so each accessory ID names exactly the
seed fragment and fragments in *other plasmids* that directly align to it by
≥ --accessory-link-min-overlap bp. No transitive chaining across third samples.

Also:
- Split accessory only at nucmer row borders whose aligned length ≥ link threshold.
- Optional --accessory-merge-gap to collapse micro-gaps after splitting.

Recommended when grouping at 40kb but core shrinks a bit:
  --edge-min-length 40000 --report-min-length 30000
"""

from __future__ import annotations
import argparse, csv, logging, shutil, bisect
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# -------------------- logging --------------------

def setup_logging(verbose: bool, quiet: bool):
    level = logging.INFO
    if verbose: level = logging.DEBUG
    if quiet: level = logging.WARNING
    logging.basicConfig(
        level=level,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

# -------------------- file helpers --------------------

def find_fasta_files(directory: Path) -> List[Path]:
    exts = {".fa", ".fasta", ".fna", ".ffn"}
    return sorted(p for p in directory.iterdir() if p.suffix.lower() in exts and p.is_file())

def parse_genes(filename: str) -> List[str]:
    stem = Path(filename).stem
    if "bla" not in stem: return []
    gene_part = stem.split("bla", 1)[1]
    gene_part = gene_part.strip("_-. ")
    return [tok for tok in gene_part.split("and") if tok]

def canon_pair(a: str, b: str) -> Tuple[str,str]:
    return (a, b) if a <= b else (b, a)

# -------------------- intervals (1-based closed) --------------------

def merge_intervals(iv: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    if not iv: return []
    iv = sorted((min(a,b), max(a,b)) for a,b in iv)
    out = []
    s, e = iv[0]
    for a, b in iv[1:]:
        if a <= e + 1:
            e = max(e, b)
        else:
            out.append((s, e))
            s, e = a, b
    out.append((s, e))
    return out

def subtract_intervals(a: List[Tuple[int,int]], b: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    """Return a - b (parts of a not covered by b)."""
    a = merge_intervals(a); b = merge_intervals(b)
    out = []
    j = 0
    for s, e in a:
        cur = s
        while j < len(b) and b[j][1] < cur:
            j += 1
        k = j
        while k < len(b) and b[k][0] <= e:
            if b[k][0] > cur:
                out.append((cur, b[k][0]-1))
            cur = max(cur, b[k][1] + 1)
            k += 1
        if cur <= e:
            out.append((cur, e))
    return out

# -------------------- read nucmer coords --------------------
# Keep:
#   per_pair_rows[(u,v)] = normalized rows: {u,v,u_s,u_e,v_s,v_e}
#   per_seq_partner[S][J] = [(s,e)] intervals on S that align to J
#   per_side_with_id[(u,v)][u or v] = [(s,e,id), ...]  (for unbiased edge metrics)
#   row_spans_by_seq[S] = ALL nucmer row spans touching S, with aligned length
#                         (x, y, row_len_on_alignment) for splitting

def load_coords(coords_path: Path, min_id: float, log_interval: int):
    per_pair_rows: Dict[Tuple[str,str], List[Dict[str,int]]] = defaultdict(list)
    per_seq_partner: Dict[str, Dict[str, List[Tuple[int,int]]]] = defaultdict(lambda: defaultdict(list))
    per_side_with_id: Dict[Tuple[str,str], Dict[str, List[Tuple[int,int,float]]]] = defaultdict(lambda: defaultdict(list))
    row_spans_by_seq: Dict[str, List[Tuple[int,int,int]]] = defaultdict(list)
    all_seqs: Set[str] = set()

    total = kept = 0
    with coords_path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            total += 1
            if total % log_interval == 0:
                logging.info(f"[read] rows read {total:,}… kept {kept:,} (≥{min_id}%% id)")
            ident = float(row["%_IDY"])
            if ident < min_id:  # filter early
                continue
            kept += 1
            ref = row["ref_file"]; qry = row["query_file"]
            s1, e1 = int(row["S1"]), int(row["E1"])
            s2, e2 = int(row["S2"]), int(row["E2"])

            u, v = canon_pair(ref, qry)

            # for unbiased edges
            per_side_with_id[(u, v)][ref].append((s1, e1, ident))
            per_side_with_id[(u, v)][qry].append((s2, e2, ident))

            # per-seq partner coverages (for core/accessory)
            per_seq_partner[ref][qry].append((s1, e1))
            per_seq_partner[qry][ref].append((s2, e2))

            # normalized row for projection & accessory linking
            if ref == u:
                per_pair_rows[(u, v)].append({"u": u, "v": v, "u_s": s1, "u_e": e1, "v_s": s2, "v_e": e2})
            else:
                per_pair_rows[(u, v)].append({"u": u, "v": v, "u_s": s2, "u_e": e2, "v_s": s1, "v_e": e1})

            # aligned row length (shorter side)
            row_len = min(abs(e1 - s1) + 1, abs(e2 - s2) + 1)
            row_spans_by_seq[ref].append((min(s1,e1), max(s1,e1), row_len))
            row_spans_by_seq[qry].append((min(s2,e2), max(s2,e2), row_len))

            all_seqs.add(ref); all_seqs.add(qry)

    logging.info(f"[read] total rows: {total:,}; pass id: {kept:,}; pairs: {len(per_side_with_id):,}; unique sequences: {len(all_seqs):,}")
    return per_pair_rows, per_seq_partner, per_side_with_id, row_spans_by_seq, all_seqs

# -------------------- pair graph (bias-fixed) --------------------

def merge_union_weighted(intervals_with_id: List[Tuple[int,int,float]]) -> Tuple[int, float, int]:
    """
    De-overlap with identity: greedy fill uncovered bp with higher-identity first.
    Returns (union_len, id_len_sum, merged_piece_count).
    """
    if not intervals_with_id: return 0, 0.0, 0
    ivals = [(min(s,e), max(s,e), float(iden)) for s,e,iden in intervals_with_id]
    ivals.sort(key=lambda x: (-x[2], x[0], x[1]))  # high id first
    covered: List[Tuple[int,int,float]] = []
    for s, e, ident in ivals:
        a, b = s, e
        new_parts = []
        for cs, ce, _ in covered:
            if ce < a or cs > b:
                continue
            if a < cs:
                new_parts.append((a, cs - 1))
            a = max(a, ce + 1)
            if a > b: break
        if a <= b:
            new_parts.append((a, b))
        for x, y in new_parts:
            covered.append((x, y, ident))
    covered.sort(key=lambda x: (x[0], x[1]))
    union_len = 0; id_len_sum = 0.0; pieces = 0
    for s, e, ident in covered:
        seglen = e - s + 1
        union_len += seglen
        id_len_sum += ident * seglen
        pieces += 1
    return union_len, id_len_sum, pieces

def build_pair_graph(per_side_with_id, edge_min_len: int) -> nx.Graph:
    G = nx.Graph()
    kept = 0
    for (u, v), sides in per_side_with_id.items():
        u_len, u_idsum, u_pieces = merge_union_weighted(sides.get(u, []))
        v_len, v_idsum, v_pieces = merge_union_weighted(sides.get(v, []))
        shared_bp = min(u_len, v_len)
        if shared_bp < edge_min_len: continue
        if u_len == 0 and v_len == 0:
            mean_id, blocks = 0.0, 0
        elif u_len <= v_len:
            mean_id, blocks = (u_idsum / u_len), u_pieces
        else:
            mean_id, blocks = (v_idsum / v_len), v_pieces
        G.add_edge(u, v, shared_bp=shared_bp, mean_id=mean_id, blocks=blocks)
        kept += 1
    logging.info(f"[graph] nodes: {G.number_of_nodes():,}; edges kept (≥{edge_min_len} bp): {kept:,}")
    return G

# -------------------- groups (same semantics) --------------------

def find_groups(G: nx.Graph, max_missing: int, overlap: float):
    logging.info("[groups] finding maximal cliques…")
    cliques = [set(c) for c in nx.find_cliques(G) if len(c) > 1]
    cliques.sort(key=len, reverse=True)
    logging.info(f"[groups] cliques>1: {len(cliques):,}")
    groups = []; alts = 0
    for clique in cliques:
        placed = False
        for main in groups:
            members: Set[str] = main["members"]
            missing = len(members - clique)
            cov = len(members & clique) / len(members)
            if missing <= max_missing and cov >= overlap:
                main.setdefault("alts", []).append(clique)
                alts += 1; placed = True; break
        if not placed:
            groups.append({"members": clique, "alts": []})
    logging.info(f"[groups] main: {len(groups):,}; alts: {alts:,}")
    return groups

def make_group_names(groups, edge_min_len: int) -> Dict[frozenset, str]:
    kb = max(1, int(round(edge_min_len / 1000)))
    return {frozenset(g["members"]): f"group_{kb}kb_{i+1}" for i, g in enumerate(groups)}

# -------------------- copying & reports --------------------

def copy_files(members: Set[str], fasta_dir: Path, dest: Path) -> None:
    if dest.exists() and any(dest.iterdir()):
        logging.info(f"[copy] skip {dest} (exists & not empty)")
        return
    dest.mkdir(parents=True, exist_ok=True)
    for name in members:
        src = fasta_dir / name
        if src.exists():
            shutil.copy(src, dest / name)

def write_reports(outdir: Path, group_names, groups, genes, edge_min_len: int):
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
            for m in members: plasmid_to_groups[m].append(main_name)
            counts = Counter(g for m in members for g in genes.get(m, []))
            for gene, count in counts.items():
                gg.write(f"{main_name}\t{gene}\t{count}\n")
            for idx, alt in enumerate(grp["alts"], 1):
                alt_name = f"{main_name}_alt{idx}"
                alt_members = sorted(alt)
                gr.write(f"{alt_name}\t{len(alt_members)}\t{','.join(alt_members)}\n")
                for m in alt_members: plasmid_to_groups[m].append(alt_name)
                counts = Counter(g for m in alt_members for g in genes.get(m, []))
                for gene, count in counts.items():
                    gg.write(f"{alt_name}\t{gene}\t{count}\n")

    with multi_report.open("w") as mr:
        mr.write("plasmid\tgroups\n")
        for plasmid, grps in plasmid_to_groups.items():
            if len(grps) > 1:
                mr.write(f"{plasmid}\t{','.join(sorted(grps))}\n")
    logging.info("[write] group reports done")

# -------------------- plotting (optional) --------------------

def plot_graph(G, genes, group_color, group_to_color, group_members, args):
    H = G.subgraph([n for n, deg in G.degree() if deg > 0]).copy()
    if H.number_of_nodes() == 0:
        logging.info("[plot] no edges to plot"); return
    plt.figure(figsize=(12, 10), dpi=300)
    pos = nx.kamada_kawai_layout(H)
    if args.color_by == "gene":
        palette = plt.get_cmap("tab20")
        gene_names = sorted({(genes.get(n, ["unknown"])[0]) for n in H.nodes})
        denom = max(1, len(gene_names) - 1)
        cmap = {g: palette(i/denom) for i, g in enumerate(gene_names)}
        node_colors = [cmap.get(genes.get(n, ["unknown"])[0], "lightgrey") for n in H.nodes]
        legend_labels = gene_names
        legend_colors = [cmap[g] for g in legend_labels]
    else:
        node_colors = [group_color.get(n, "lightgrey") for n in H.nodes]
        visible_groups = [g for g, mem in group_members.items() if any(m in H for m in mem)]
        legend_labels = [g for g in group_to_color if g in visible_groups]
        legend_colors = [group_to_color[g] for g in legend_labels]
    nx.draw_networkx_edges(H, pos, alpha=0.5)
    nx.draw_networkx_nodes(H, pos, node_color=node_colors, node_size=args.node_size,
                           edgecolors="black", linewidths=0.8)
    if args.show_labels: nx.draw_networkx_labels(H, pos, font_size=args.label_size)
    handles = [mpatches.Patch(color=c, label=l) for l, c in zip(legend_labels, legend_colors)]
    if handles:
        plt.legend(handles=handles, title=args.color_by.capitalize(),
                   loc="center left", bbox_to_anchor=(1.02,0.5), borderaxespad=0)
    plt.axis("off")
    if args.output_figure: plt.savefig(args.output_figure, bbox_inches="tight")
    else: plt.show()

# -------------------- CORE: anchor sweep + projection (unchanged) --------------------

def sweep_on_anchor(anchor: str, partners: List[str], per_seq_partner, report_min_len: int) -> List[Tuple[int,int,int]]:
    """
    Minimal anchor segments with constant count of distinct partners covering.
    Returns (start, end, count) filtered by report_min_len; count>=1 only.
    """
    if not partners: return []
    per = []
    for j in partners:
        ivs = merge_intervals(per_seq_partner.get(anchor, {}).get(j, []))
        if ivs:
            per.append(ivs)
    if not per: return []
    # sweep events
    events = defaultdict(int)
    for ivs in per:
        for s,e in ivs:
            events[s] += 1
            events[e+1] -= 1
    segments = []
    run = 0
    last = None
    for pos in sorted(events.keys()):
        if last is not None and run > 0:
            s = last; e = pos - 1
            if e >= s and (e - s + 1) >= report_min_len:
                segments.append((s, e, run))
        run += events[pos]
        last = pos
    return segments

def map_interval(anchor: str, member: str, per_pair_rows, interval: Tuple[int,int]) -> List[Tuple[int,int]]:
    """
    Project [a,b] on 'anchor' to coordinates on 'member' via raw nucmer rows.
    Handles strand and trims to alignment borders. Merges outputs.
    """
    a, b = interval
    u, v = canon_pair(anchor, member)
    rows = per_pair_rows.get((u, v), [])
    out = []
    for r in rows:
        if anchor == r["u"]:
            u_s, u_e, v_s, v_e = r["u_s"], r["u_e"], r["v_s"], r["v_e"]
        else:
            u_s, u_e, v_s, v_e = r["v_s"], r["v_e"], r["u_s"], r["u_e"]
        # overlap on anchor side
        a1, a2 = min(u_s, u_e), max(u_s, u_e)
        s = max(a1, a); e = min(a2, b)
        if s > e: continue
        sign_u = 1 if u_e >= u_s else -1
        sign_v = 1 if v_e >= v_s else -1
        dist_s = (s - u_s) * sign_u
        dist_e = (e - u_s) * sign_u
        y1 = v_s + dist_s * sign_v
        y2 = v_s + dist_e * sign_v
        ys, ye = (min(y1, y2), max(y1, y2))
        out.append((int(ys), int(ye)))
    return merge_intervals(out)

# -------------------- ACCESSORY helpers --------------------

def build_accessory_union_by_seq(groups, per_seq_partner, per_pair_rows, report_min_len: int):
    """
    Compute accessory union per sequence:
      acc_by_seq[S] = union(all coverage from ANY partner) - union(all projected CORE on S)
    CORE is computed per group (anchor sweep + projection) and accumulated.
    """
    core_by_seq: Dict[str, List[Tuple[int,int]]] = defaultdict(list)
    acc_by_seq: Dict[str, List[Tuple[int,int]]] = defaultdict(list)

    for grp in groups:
        members = sorted(grp["members"])
        if not members: continue
        anchor = min(members)
        partners = [x for x in members if x != anchor]

        # sweep anchor -> segments; core = count == len(partners)
        segments = sweep_on_anchor(anchor, partners, per_seq_partner, report_min_len)
        core_anchor = [(s,e) for (s,e,c) in segments if c == len(partners)]

        # map anchor core to all members and collect
        per_member_core = {m: [] for m in members}
        per_member_core[anchor].extend(core_anchor)
        for m in partners:
            mapped = []
            for (a,b) in core_anchor:
                mapped.extend(map_interval(anchor, m, per_pair_rows, (a,b)))
            per_member_core[m].extend(merge_intervals(mapped))

        # accumulate cores globally (in case of multi-group memberships)
        for m in members:
            core_by_seq[m].extend(per_member_core[m])

    # merge global cores per seq
    for s in core_by_seq:
        core_by_seq[s] = merge_intervals(core_by_seq[s])

    # partner union per seq across ALL partners in coords
    for s, partners in per_seq_partner.items():
        all_cov = []
        for ivs in partners.values():
            all_cov.extend(ivs)
        partner_union = merge_intervals(all_cov)
        acc = subtract_intervals(partner_union, core_by_seq.get(s, []))
        # keep only ≥ report_min_len
        acc = [iv for iv in acc if (iv[1]-iv[0]+1) >= report_min_len]
        if acc:
            acc_by_seq[s] = merge_intervals(acc)

    return acc_by_seq, core_by_seq

def split_accessory_by_row_boundaries(acc_by_seq, row_spans_by_seq, min_row_len: int, merge_gap: int, report_min_len: int):
    """
    Split accessory intervals only at boundaries of rows whose aligned length >= min_row_len.
    Then optionally merge adjacent fragments separated by <= merge_gap.
    Keep only fragments >= report_min_len.
    Returns frags_by_seq[S] = [(s,e), ...].
    """
    frags_by_seq: Dict[str, List[Tuple[int,int]]] = {}

    for s, acc_list in acc_by_seq.items():
        if not acc_list:
            continue

        cuts = set()
        for a,b in acc_list:
            cuts.add(a); cuts.add(b+1)
        for (x,y,row_len) in row_spans_by_seq.get(s, []):
            if row_len >= min_row_len:
                cuts.add(x); cuts.add(y+1)

        cutpoints = sorted(cuts)
        frags = []
        for (a,b) in acc_list:
            idx_start = bisect.bisect_left(cutpoints, a+1)
            idx_end   = bisect.bisect_right(cutpoints, b)
            prev = a
            for i in range(idx_start, idx_end):
                cp = cutpoints[i]
                if prev <= cp-1:
                    frags.append((prev, cp-1))
                prev = cp
            if prev <= b:
                frags.append((prev, b))

        # Optional: merge tiny gaps
        if merge_gap > 0 and frags:
            merged = []
            cur_s, cur_e = frags[0]
            for s2, e2 in frags[1:]:
                if s2 - cur_e - 1 <= merge_gap:  # gap ≤ merge_gap
                    cur_e = max(cur_e, e2)
                else:
                    merged.append((cur_s, cur_e))
                    cur_s, cur_e = s2, e2
            merged.append((cur_s, cur_e))
            frags = merged

        # Enforce report_min_len
        frags = [iv for iv in frags if (iv[1] - iv[0] + 1) >= report_min_len]
        if frags:
            frags_by_seq[s] = frags

    return frags_by_seq

# Anchor-seeded accessory naming (NO transitive chaining)
def assign_accessory_ids(
    per_pair_rows: Dict[Tuple[str,str], List[Dict[str,int]]],
    frags_by_seq: Dict[str, List[Tuple[int,int]]],
    link_min_overlap: int
) -> List[Tuple[str, str, int, int]]:
    """
    Process accessory fragments in (seq, start) order. For each unassigned fragment F on seed
    sequence S, start a new ACC ID and include ONLY fragments on other sequences that have a
    DIRECT nucmer row to S overlapping F by ≥ link_min_overlap (on both sides).
    Returns a list of lines: (ACC_ID, seq, start, end)
    """

    # Build quick indices for fragments per sequence
    seq2frags = {s: sorted(frags) for s, frags in frags_by_seq.items()}
    seq2starts = {s: [x for x,_ in seq2frags[s]] for s in seq2frags}
    seq2ends   = {s: [y for _,y in seq2frags[s]] for s in seq2frags}
    assigned: Dict[str, List[bool]] = {s: [False]*len(seq2frags[s]) for s in seq2frags}

    # Map sequence -> list of pair keys it participates in (for speed)
    pairs_by_seq: Dict[str, List[Tuple[str,str]]] = defaultdict(list)
    for key in per_pair_rows.keys():
        u, v = key
        pairs_by_seq[u].append(key)
        pairs_by_seq[v].append(key)

    def frags_overlapping(seq: str, a: int, b: int, min_len: int) -> List[int]:
        """Return indices of fragments on seq overlapping [a,b] by at least min_len."""
        if seq not in seq2frags:
            return []
        frags = seq2frags[seq]; starts = seq2starts[seq]; ends = seq2ends[seq]
        i = bisect.bisect_left(ends, a)
        out_idx = []
        while i < len(frags) and frags[i][0] <= b:
            fa, fb = frags[i]
            ov_a = max(fa, a); ov_b = min(fb, b)
            if ov_b >= ov_a and (ov_b - ov_a + 1) >= min_len:
                out_idx.append(i)
            i += 1
        return out_idx

    def map_overlap_from_seed(seed_seq: str, row: Dict[str,int], seed_a: int, seed_b: int, other_seq: str) -> List[Tuple[int,int]]:
        """
        Given a row touching seed_seq, intersect with [seed_a, seed_b] on the seed side,
        and map that intersection to the other side. Return merged intervals on other_seq.
        """
        if seed_seq == row["u"]:
            s_s, s_e, o_s, o_e = row["u_s"], row["u_e"], row["v_s"], row["v_e"]
        else:
            s_s, s_e, o_s, o_e = row["v_s"], row["v_e"], row["u_s"], row["u_e"]
        sa, sb = min(s_s, s_e), max(s_s, s_e)
        ia, ib = max(sa, seed_a), min(sb, seed_b)
        if ib < ia:
            return []
        sign_s = 1 if s_e >= s_s else -1
        sign_o = 1 if o_e >= o_s else -1
        # distances from seed row start
        d1 = (ia - s_s) * sign_s
        d2 = (ib - s_s) * sign_s
        y1 = o_s + d1 * sign_o
        y2 = o_s + d2 * sign_o
        ys, ye = (min(y1, y2), max(y1, y2))
        return [(int(ys), int(ye))]

    acc_lines: List[Tuple[str, str, int, int]] = []
    acc_id_counter = 0

    for s in sorted(seq2frags.keys()):
        for idx, (a, b) in enumerate(seq2frags[s]):
            if assigned[s][idx]:
                continue
            # Start a new component, anchor = (s, a, b)
            acc_id_counter += 1
            acc_id = f"ACC_{acc_id_counter:06d}"
            # assign anchor
            assigned[s][idx] = True
            acc_lines.append((acc_id, s, a, b))

            # explore only DIRECT rows from s to others
            for key in pairs_by_seq.get(s, []):
                u, v = key
                other = v if s == u else u
                if other not in seq2frags:
                    continue
                for r in per_pair_rows[key]:
                    # choose row orientation so seed is on the "seed side"
                    seed_side = r["u"] if s == r["u"] else r["v"]
                    if seed_side != s:
                        continue  # (shouldn't happen, but guard)

                    # overlap of row with seed frag
                    seed_row_a = min(r["u_s"], r["u_e"]) if s == r["u"] else min(r["v_s"], r["v_e"])
                    seed_row_b = max(r["u_s"], r["u_e"]) if s == r["u"] else max(r["v_s"], r["v_e"])
                    ov_a = max(seed_row_a, a); ov_b = min(seed_row_b, b)
                    if ov_b < ov_a:
                        continue
                    if (ov_b - ov_a + 1) < link_min_overlap:
                        continue  # row doesn't give enough overlap to the seed

                    # map the overlapped slice to 'other'
                    mapped = map_overlap_from_seed(s, r, ov_a, ov_b, other)
                    for (ys, ye) in mapped:
                        # find 'other' fragments that overlap mapped region by ≥ link_min_overlap
                        for j in frags_overlapping(other, ys, ye, link_min_overlap):
                            if not assigned[other][j]:
                                assigned[other][j] = True
                                oa, ob = seq2frags[other][j]
                                acc_lines.append((acc_id, other, oa, ob))

    return acc_lines

# -------------------- CLI --------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("coords", type=Path, help="Path to all_coords.tsv")
    ap.add_argument("fasta_dir", type=Path, help="Directory with plasmid FASTA files")
    # Split thresholds
    ap.add_argument("--edge-min-length", type=int, required=True,
                    help="Minimum shared length (bp) to CREATE a graph edge (grouping threshold)")
    ap.add_argument("--report-min-length", type=int,
                    help="Minimum segment length (bp) to REPORT core/accessory (defaults to --edge-min-length)")
    ap.add_argument("--min-identity", type=float, default=90.0,
                    help="Minimum %% identity to accept nucmer rows (filter before grouping)")  # '%%' for argparse
    # Accessory linkage + splitting controls
    ap.add_argument("--accessory-link-min-overlap", type=int, default=1000,
                    help="Minimum aligned-overlap (bp) required to link two accessory fragments into the same block")
    ap.add_argument("--accessory-merge-gap", type=int, default=0,
                    help="Merge adjacent accessory fragments separated by ≤ this many bp (default: 0 = off)")
    ap.add_argument("--output-dir", type=Path, default=Path("."), help="Output directory")
    ap.add_argument("--max-missing", type=int, default=20,
                    help="Maximum members missing to call an alternative group")
    ap.add_argument("--overlap", type=float, default=0.90,
                    help="Fractional overlap (0–1) to attach a clique as alt of a larger group")
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--color-by", choices=["gene","group"], default="gene")
    ap.add_argument("--node-size", type=int, default=300)
    ap.add_argument("--label-size", type=int, default=8)
    ap.add_argument("--show-labels", action="store_true")
    ap.add_argument("--output-figure", type=Path)
    ap.add_argument("--log-interval", type=int, default=100_000)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    setup_logging(args.verbose, args.quiet)
    outdir = args.output_dir; outdir.mkdir(parents=True, exist_ok=True)
    report_min_len = args.report_min_length if args.report_min_length is not None else args.edge_min_length

    # Read coords
    (per_pair_rows, per_seq_partner, per_side_with_id,
     row_spans_by_seq, all_seqs) = load_coords(args.coords, args.min_identity, args.log_interval)

    # Pair graph + groups
    G = build_pair_graph(per_side_with_id, args.edge_min_length)
    for s in all_seqs: G.add_node(s)
    groups = find_groups(G, args.max_missing, args.overlap)
    group_names = make_group_names(groups, args.edge_min_length)

    # Copy FASTAs and prep colors
    fasta_files = find_fasta_files(args.fasta_dir)
    genes = {f.name: parse_genes(f.name) for f in fasta_files}
    length_kb = max(1, int(round(args.edge_min_length/1000)))
    group_members: Dict[str, Set[str]] = {}
    group_to_color: Dict[str, str] = {}
    group_color: Dict[str, str] = {}
    node_best_size = defaultdict(int)
    palette = plt.get_cmap("tab20")
    for idx, grp in enumerate(groups, 1):
        name = f"group_{length_kb}kb_{idx}"
        members = set(grp["members"])
        group_members[name] = members
        copy_files(members, args.fasta_dir, outdir / name)
        color = palette((idx - 1) % 20); group_to_color[name] = color
        size = len(members)
        for m in members:
            if size > node_best_size[m]:
                group_color[m] = color; node_best_size[m] = size
        for alt_idx, alt in enumerate(grp["alts"], 1):
            alt_name = f"{name}_alt{alt_idx}"
            alt_members = set(alt)
            copy_files(alt_members, args.fasta_dir, outdir / alt_name)
            alt_size = len(alt_members)
            for m in alt_members:
                if alt_size > node_best_size[m]:
                    group_color[m] = color; node_best_size[m] = alt_size

    write_reports(outdir, group_names, groups, genes, args.edge_min_length)

    # Export full graph table incl. singletons
    with (outdir / "network_graph.tsv").open("w") as out:
        out.write("node1\tnode2\tshared_bp\tblocks\tmean_id\tgroup1\tgroup2\n")
        for u in G.nodes:
            nbrs = list(G.neighbors(u))
            if nbrs:
                for v in nbrs:
                    if u < v:
                        d = G[u][v]
                        g1 = next((g for g, mem in group_names.items() if u in mem), "")
                        g2 = next((g for g, mem in group_names.items() if v in mem), "")
                        out.write(f"{u}\t{v}\t{d.get('shared_bp','')}\t{d.get('blocks','')}\t{d.get('mean_id','')}\t{g1}\t{g2}\n")
            else:
                g1 = next((g for g, mem in group_names.items() if u in mem), "")
                out.write(f"{u}\t\t\t\t\t{g1}\n")
    logging.info("[write] network_graph.tsv done")

    # ===================== CORE (unchanged behavior) =====================
    out_tsv = outdir / "blocks_core_accessory.tsv"
    with out_tsv.open("w") as out:
        out.write("Block_ID\tgroup_code\tsequence_name\tstart_coords\tend_coords\tcore_or_accessory\n")

        core_line_count = 0
        block_ctr = 0
        for grp in groups:
            members = sorted(grp["members"])
            if not members: continue
            gcode = group_names[frozenset(grp["members"])]
            anchor = min(members)
            partners = [x for x in members if x != anchor]

            segments = sweep_on_anchor(anchor, partners, per_seq_partner, report_min_len)
            core_anchor = [(s,e) for (s,e,c) in segments if c == len(partners)]

            for (a,b) in core_anchor:
                block_ctr += 1
                bid = f"{gcode}__B{block_ctr:04d}"
                out.write(f"{bid}\t{gcode}\t{anchor}\t{a}\t{b}\tcore\n"); core_line_count += 1
                for m in partners:
                    mapped = map_interval(anchor, m, per_pair_rows, (a,b))
                    for (x,y) in mapped:
                        out.write(f"{bid}\t{gcode}\t{m}\t{x}\t{y}\tcore\n"); core_line_count += 1
        logging.info(f"[core] lines written: {core_line_count:,}")

    # ===================== ACCESSORY (anchor-seeded, non-transitive) =====================

    # 1) Accessory union per seq (partner union minus core)
    acc_by_seq, _core_by_seq = build_accessory_union_by_seq(
        groups, per_seq_partner, per_pair_rows, report_min_len
    )
    logging.info(f"[accessory] sequences with accessory: {len(acc_by_seq):,}")

    # 2) Split accessory using only eligible rows; optional merge of tiny gaps
    frags_by_seq = split_accessory_by_row_boundaries(
        acc_by_seq,
        row_spans_by_seq,
        args.accessory_link_min_overlap,   # only rows that could be linked anyway
        args.accessory_merge_gap,
        report_min_len
    )
    total_frags = sum(len(v) for v in frags_by_seq.values())
    logging.info(f"[accessory] total accessory fragments after splitting: {total_frags:,}")

    # 3) Assign accessory IDs by seeding from each unassigned fragment; link only DIRECT overlaps to seed
    acc_lines = assign_accessory_ids(
        per_pair_rows,
        frags_by_seq,
        args.accessory_link_min_overlap
    )

    # 4) Write accessory lines (no group code)
    with out_tsv.open("a") as out:
        for acc_id, seq, a, b in acc_lines:
            out.write(f"{acc_id}\t\t{seq}\t{a}\t{b}\taccessory\n")

    unique_member_count = len(set().union(*[set(g["members"]) for g in groups])) if groups else 0
    summed_member_count = sum(len(g["members"]) for g in groups)
    logging.info(f"[write] {out_tsv}  (core lines: {core_line_count:,}; accessory frags: {len(acc_lines):,}; "
                 f"unique sequences in coords: {len(all_seqs):,}; summed group members: {summed_member_count:,}; "
                 f"link_min_overlap: {args.accessory_link_min_overlap} bp; merge_gap: {args.accessory_merge_gap} bp)")

    # Optional plot
    if args.plot:
        plot_graph(G, genes, group_color, group_to_color, group_members, args)
        
    if args.output_figure:
        plt.savefig(args.output_figure, bbox_inches="tight")
        logging.info(f"[plot] saved figure to {args.output_figure}")
    else:
        plt.show()


if __name__ == "__main__":
    main()

