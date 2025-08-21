#!/usr/bin/env python3
"""
Export per-sequence blocks as CORE (group backbone) or ACCESSORY.

Reads all_pairs nucmer coords TSV (same input as identify_plasmid_groups.py),
reconstructs block "components" (connected via pairwise alignments), merges
overlaps within each sequence, then labels blocks:

- CORE (by group): A block is CORE for group G if *every* main member of G
  has ≥ --min-length bp covered by that block.
- ACCESSORY: Present in ≥1 sequence but not CORE for that sequence's group
  (or the sequence has no group).

Emits:
  blocks_core_accessory.tsv with columns:
    Block_ID   group_code_or_empty   sequence_name   start   end   core_or_accessory
"""

from __future__ import annotations
import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

import networkx as nx

# ---------- Helpers copied/adapted from your grouping script ----------
# (Mirrors CLI + grouping behavior so results match your pipeline.)
# Ref: identify_plasmid_groups.py CLI & find_groups().

def _union_length(intervals: List[Tuple[int, int]]) -> int:
    """Length of union of 1-based closed intervals."""
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

def _merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Return merged non-overlapping 1-based closed intervals, sorted."""
    if not intervals:
        return []
    ivals = sorted((min(s, e), max(s, e)) for s, e in intervals)
    merged = []
    cs, ce = ivals[0]
    for s, e in ivals[1:]:
        if s <= ce + 1:
            ce = max(ce, e)
        else:
            merged.append((cs, ce))
            cs, ce = s, e
    merged.append((cs, ce))
    return merged

def build_graph(coords_file: Path, min_len: int, min_id: float) -> nx.Graph:
    """
    (From your script) Build graph of plasmids; edges connect pairs that share
    at least min_len bp of union coverage at >= min_id identity.
    """
    ref_intervals = defaultdict(list)
    qry_intervals = defaultdict(list)
    id_len_sum = defaultdict(float)
    len_sum = defaultdict(int)
    block_count = defaultdict(int)

    with coords_file.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ident = float(row["%_IDY"])
            if ident < min_id:
                continue
            ref = row["ref_file"]; qry = row["query_file"]
            S1, E1 = int(row["S1"]), int(row["E1"])
            S2, E2 = int(row["S2"]), int(row["E2"])
            seg_len = min(int(row["LEN_1"]), int(row["LEN_2"]))
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
            G.add_edge(ref, qry, shared_bp=shared_bp,
                       blocks=block_count[key], mean_id=mean_id)
    return G

def find_groups(G: nx.Graph, max_missing: int, overlap: float):
    """
    (From your script) Maximal cliques as main groups; alt groups are small
    near-duplicates of a larger group (<= max_missing missing AND >= overlap).
    Returns: [{"members": set(...), "alts": [set(...), ...]}, ...]
    """
    cliques = [set(c) for c in nx.find_cliques(G) if len(c) > 1]
    cliques.sort(key=len, reverse=True)
    groups = []
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

# ---------- New: build "block components" across all sequences ----------

def load_components(coords_file: Path, min_id: float) -> Tuple[List[List[Tuple[str,int,int]]], Dict[int, Dict[str, List[Tuple[int,int]]]]]:
    """
    Build connected components of intervals across all sequences.

    Returns:
      components: list of components; each component is a list of (seq, s, e) nodes
      comp_seq_ivals: component_id -> { seq_name -> [ (s,e), ... ] }  (raw, not merged)
    """
    # Build a graph whose nodes are per-row, per-side intervals; edges connect ref-side to qry-side for the same row.
    # This way, connected components capture transitive sharing across many pairs (also across different groups).
    G = nx.Graph()

    nodes = []  # (seq, s, e)
    # For mapping row to node indices
    with coords_file.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for idx, row in enumerate(reader):
            ident = float(row["%_IDY"])
            if ident < min_id:
                continue
            ref = row["ref_file"]; qry = row["query_file"]
            s1, e1 = int(row["S1"]), int(row["E1"])
            s2, e2 = int(row["S2"]), int(row["E2"])

            n_ref = ("ref", idx)  # unique per row-side
            n_qry = ("qry", idx)

            # attach interval metadata to node attributes
            G.add_node(n_ref, seq=ref, s=s1, e=e1)
            G.add_node(n_qry, seq=qry, s=s2, e=e2)
            G.add_edge(n_ref, n_qry)  # link the two sides of the same nucmer row

    components = []
    comp_seq_ivals: Dict[int, Dict[str, List[Tuple[int,int]]]] = {}

    for cid, comp in enumerate(nx.connected_components(G)):
        comp_nodes = list(comp)
        components.append([(G.nodes[n]["seq"], G.nodes[n]["s"], G.nodes[n]["e"]) for n in comp_nodes])

        seq2ivals: Dict[str, List[Tuple[int,int]]] = defaultdict(list)
        for n in comp_nodes:
            seq = G.nodes[n]["seq"]
            s = int(G.nodes[n]["s"]); e = int(G.nodes[n]["e"])
            seq2ivals[seq].append((s, e))
        comp_seq_ivals[cid] = seq2ivals

    return components, comp_seq_ivals

# ---------- Main ----------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coords", type=Path, help="Path to all_coords.tsv")
    parser.add_argument("fasta_dir", type=Path, help="Directory containing plasmid FASTA files")
    parser.add_argument("--min-length", type=int, required=True, help="Minimum block length to REPORT (bp)")
    parser.add_argument("--min-identity", type=float, default=90.0, help="Minimum percent identity to accept nucmer rows")
    parser.add_argument("--output-dir", type=Path, default=Path("."), help="Directory for results")
    parser.add_argument("--max-missing", type=int, default=20, help="Max members missing to call alt group (same behavior)")
    parser.add_argument("--overlap", type=float, default=0.90, help="Fractional overlap to call alt group (same behavior)")
    args = parser.parse_args()

    # 1) Rebuild your groupings so "core" matches your pipeline’s semantics
    Gpairs = build_graph(args.coords, args.min_length, args.min_identity)  # same as your grouping graph
    groups = find_groups(Gpairs, args.max_missing, args.overlap)           # main groups + alts

    # Map main groups to names like group_XXkb_N so we can print group codes
    length_kb = args.min_length // 1000
    group_names: Dict[frozenset, str] = {}
    main_groups: List[Set[str]] = []
    for idx, grp in enumerate(groups, 1):
        name = f"group_{length_kb}kb_{idx}"
        members = set(grp["members"])
        group_names[frozenset(members)] = name
        main_groups.append(members)

    # 2) Build block components across ALL sequences (cross-group allowed)
    _, comp_seq_ivals = load_components(args.coords, args.min_identity)

    # 3) For each component, merge per-sequence intervals and filter by --min-length
    outdir = args.output_dir
    outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / "blocks_core_accessory.tsv"

    with outfile.open("w") as out:
        out.write("Block_ID\tgroup_code\tsequence_name\tstart_coords\tend_coords\tcore_or_accessory\n")

        for cid, seq2ivals in comp_seq_ivals.items():
            # Merge per sequence
            merged_per_seq: Dict[str, List[Tuple[int,int]]] = {
                seq: _merge_intervals(ivals) for seq, ivals in seq2ivals.items()
            }
            # Only intervals >= min-length are reportable (per-sequence contribution)
            reportable: Dict[str, List[Tuple[int,int]]] = {}
            for seq, ivals in merged_per_seq.items():
                kept = [(s,e) for (s,e) in ivals if (e - s + 1) >= args.min_length]
                if kept:
                    reportable[seq] = kept

            if not reportable:
                continue

            block_id = f"B{cid:06d}"

            # Determine which groups (if any) consider this component a CORE block
            # A block is CORE for group G if ALL main members of G have a reportable interval.
            core_groups_for_component: List[str] = []
            for members in main_groups:
                if all(m in reportable for m in members):
                    gname = group_names[frozenset(members)]
                    core_groups_for_component.append(gname)

            # 4) Emit one line per (block, sequence, merged interval)
            for seq, ivals in reportable.items():
                # Determine label: 'core' if seq belongs to ANY group that is core for this component
                # Prefer to print the group code if core; otherwise group_code is empty.
                seq_is_core = False
                seq_group_code = ""
                for members in main_groups:
                    if seq in members:
                        gname = group_names[frozenset(members)]
                        if gname in core_groups_for_component:
                            seq_is_core = True
                            seq_group_code = gname
                            break  # one is enough

                label = "core" if seq_is_core else "accessory"
                group_code_to_print = seq_group_code if seq_is_core else ""

                for (s,e) in ivals:
                    out.write(f"{block_id}\t{group_code_to_print}\t{seq}\t{s}\t{e}\t{label}\n")

    print(f"[OK] Wrote {outfile}")

if __name__ == "__main__":
    main()

