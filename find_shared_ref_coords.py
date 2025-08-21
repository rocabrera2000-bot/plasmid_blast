#!/usr/bin/env python3
"""
Find reference coordinates shared by all members of a group.

Inputs
------
- group_report.tsv: columns: group, size, members (members are comma-separated filenames)
- all_coords.tsv  : nucmer coords table with header like:
    ref_file  query_file  S1  E1  S2  E2  LEN_1  LEN_2  %_IDY  LEN_R  LEN_Q  %COV_R  %COV_Q  REF_TAG  QRY_TAG
- group name      : e.g. group_40kb_47
- reference file  : a member filename to serve as the coordinate reference

What it does
------------
1) For each non-reference member, pulls all alignments to the reference from all_coords.tsv,
   normalizes them to reference coordinates (using S1/E1 if ref_file==REF, else S2/E2),
   applies min-identity and min-length filters, and merges overlaps per member.
2) Intersects these per-member merged intervals across all members to get regions on the
   REFERENCE that are present in *every* member.
3) Prints the resulting interval list (1-based inclusive), plus optional per-member coverage stats.

Notes
-----
- Orientation doesn’t matter for shared-ref coverage; we normalize to ref coordinates and merge.
- If any member has no qualifying intervals, the shared set will be empty (script will warn).
- If `all_coords.tsv` contains both REF→QRY and QRY→REF rows, both directions are supported.

Example
-------
python find_shared_ref_coords.py \
  --group-report group_report.tsv \
  --coords ../polished_plasmids/nucmer/all_coords.tsv \
  --group group_40kb_47 \
  --ref E1192_hybrid2025_ctg_2_len_240710_blaOXA181.fasta \
  --min-idy 95 --min-len 500 \
  --print-coverage
"""

import sys, csv, argparse
from pathlib import Path
from typing import List, Tuple, Dict, Optional

# ---------- basic interval helpers ----------

def _merge_intervals(ivals: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    """Merge 1-based inclusive intervals; returns sorted, non-overlapping."""
    if not ivals:
        return []
    ivals = sorted((min(a,b), max(a,b)) for a,b in ivals)
    out = [ivals[0]]
    for a,b in ivals[1:]:
        la, lb = out[-1]
        if a <= lb + 1:  # overlap or touch
            out[-1] = (la, max(lb, b))
        else:
            out.append((a,b))
    return out

def _intersect_two(a: List[Tuple[int,int]], b: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    """Intersect two sorted, non-overlapping interval lists (1-based inclusive)."""
    i=j=0; out=[]
    while i < len(a) and j < len(b):
        a1,a2 = a[i]; b1,b2 = b[j]
        s = max(a1,b1); e = min(a2,b2)
        if s <= e:
            out.append((s,e))
        if a2 < b2:
            i += 1
        else:
            j += 1
    return out

def _intersect_many(lists: List[List[Tuple[int,int]]]) -> List[Tuple[int,int]]:
    if not lists:
        return []
    cur = lists[0]
    for nxt in lists[1:]:
        cur = _intersect_two(cur, nxt)
        if not cur:
            break
    return cur

def _sum_len(ivals: List[Tuple[int,int]]) -> int:
    return sum(b-a+1 for a,b in ivals)

# ---------- IO helpers ----------

def read_group_members(path: Path) -> Dict[str, List[str]]:
    groups: Dict[str, List[str]] = {}
    with path.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        need = {"group","members"}
        if not need.issubset(set(r.fieldnames or [])):
            raise SystemExit(f"[ERROR] {path} missing required columns: {sorted(need)}")
        for row in r:
            name = row["group"].strip()
            members = [m.strip() for m in row["members"].split(",") if m.strip()]
            groups[name] = members
    return groups

def _to_int(row: dict, key: str, default: Optional[int]=None) -> Optional[int]:
    try:
        return int(float(row[key]))
    except Exception:
        return default

def _to_float(row: dict, key: str, default: Optional[float]=None) -> Optional[float]:
    try:
        return float(row[key])
    except Exception:
        return default

def load_coords(path: Path) -> List[dict]:
    with path.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        rows = list(r)
    # sanity columns
    need = {"ref_file","query_file","S1","E1","S2","E2"}
    if not need.issubset(set(rows[0].keys() if rows else [])):
        raise SystemExit(f"[ERROR] {path} missing columns: {sorted(need)}")
    return rows

# ---------- main logic ----------

def collect_member_ref_intervals(rows: List[dict],
                                 ref: str,
                                 member: str,
                                 min_idy: float,
                                 min_len: int,
                                 idy_col_candidates = ("%_IDY","%_idy","pct_idy","%IDY")) -> List[Tuple[int,int]]:
    """
    From all_coords rows, return merged intervals on REFERENCE covered by 'member'
    after applying min identity and min length thresholds (based on the coordinates
    on the REFERENCE).
    """
    ivals: List[Tuple[int,int]] = []

    # helper to pull identity with flexible header name
    def get_idy(row: dict) -> float:
        for k in idy_col_candidates:
            if k in row:
                try:
                    return float(row[k])
                except Exception:
                    pass
        return float("nan")

    for row in rows:
        rf = row["ref_file"].strip()
        qf = row["query_file"].strip()
        if rf == ref and qf == member:
            s = _to_int(row, "S1"); e = _to_int(row, "E1")
            if s is None or e is None: continue
            seg_len = abs(e - s) + 1
            idy = get_idy(row)
            if (min_len is None or seg_len >= min_len) and (min_idy is None or (idy == idy and idy >= min_idy)):
                ivals.append((min(s,e), max(s,e)))
        elif rf == member and qf == ref:
            s = _to_int(row, "S2"); e = _to_int(row, "E2")
            if s is None or e is None: continue
            seg_len = abs(e - s) + 1
            idy = get_idy(row)
            if (min_len is None or seg_len >= min_len) and (min_idy is None or (idy == idy and idy >= min_idy)):
                ivals.append((min(s,e), max(s,e)))
        else:
            continue

    return _merge_intervals(ivals)

def main():
    ap = argparse.ArgumentParser(description="Find reference coordinates shared by all members of a group.")
    ap.add_argument("--group-report", required=True, type=Path)
    ap.add_argument("--coords", required=True, type=Path)
    ap.add_argument("--group", required=True, help="Group name, e.g. group_40kb_47")
    ap.add_argument("--ref", required=True, help="Reference filename (must be a member of the group)")
    ap.add_argument("--min-idy", type=float, default=95.0, help="Min %% identity for segments (default 95.0). Use 0 to disable.")
    ap.add_argument("--min-len", type=int, default=200, help="Min segment length on reference (default 200). Use 0 to disable.")
    ap.add_argument("--print-coverage", action="store_true", help="Print per-member coverage stats on the reference.")
    ap.add_argument("--tsv-out", type=Path, default=None, help="Optional TSV to write shared intervals.")
    args = ap.parse_args()

    groups = read_group_members(args.group_report)
    if args.group not in groups:
        raise SystemExit(f"[ERROR] Group '{args.group}' not found in {args.group_report}")

    members = groups[args.group]
    if args.ref not in members:
        raise SystemExit(f"[ERROR] Reference '{args.ref}' is not in group '{args.group}'")

    rows = load_coords(args.coords)

    # thresholds
    min_idy = None if args.min_idy is not None and args.min_idy <= 0 else args.min_idy
    min_len = None if args.min_len is not None and args.min_len <= 0 else args.min_len

    # compute per-member coverage on reference
    per_member: Dict[str, List[Tuple[int,int]]] = {}
    for m in members:
        if m == args.ref:
            continue
        ivals = collect_member_ref_intervals(rows, ref=args.ref, member=m,
                                             min_idy=min_idy, min_len=min_len)
        per_member[m] = ivals

    # report per-member coverage
    if args.print_coverage:
        print(f"# Reference: {args.ref}")
        for m, ivals in per_member.items():
            covered = _sum_len(ivals)
            # Try to infer reference length from coords if available
            ref_len = None
            for row in rows:
                if row["ref_file"].strip() == args.ref:
                    ref_len = _to_int(row, "LEN_R", None)
                    if ref_len: break
                if row["query_file"].strip() == args.ref:
                    ref_len = _to_int(row, "LEN_Q", None)
                    if ref_len: break
            cov_pct = (100.0 * covered / ref_len) if ref_len else None
            cov_txt = f"{cov_pct:.2f}%" if cov_pct is not None else "NA"
            print(f"# {m}\tcovered_bp={covered}\tcovered_pct={cov_txt}\tblocks={len(ivals)}")

    # intersect across all members
    lists: List[List[Tuple[int,int]]] = list(per_member.values())
    shared = _intersect_many(lists) if lists else []
    total_bp = _sum_len(shared)

    # print to stdout
    print("ref_file\tgroup\tblock_idx\tstart\tend\tlength")
    for i,(a,b) in enumerate(shared, start=1):
        print(f"{args.ref}\t{args.group}\t{i}\t{a}\t{b}\t{b-a+1}")

    # optional TSV
    if args.tsv_out:
        args.tsv_out.parent.mkdir(parents=True, exist_ok=True)
        with args.tsv_out.open("w") as out:
            out.write("ref_file\tgroup\tblock_idx\tstart\tend\tlength\n")
            for i,(a,b) in enumerate(shared, start=1):
                out.write(f"{args.ref}\t{args.group}\t{i}\t{a}\t{b}\t{b-a+1}\n")

    # summary to stderr
    print(f"[INFO] Members (excluding ref): {len(per_member)}", file=sys.stderr)
    print(f"[INFO] Shared blocks: {len(shared)} ; total_bp={total_bp}", file=sys.stderr)
    if not shared:
        print("[WARN] No shared intervals found. Consider lowering --min-idy/--min-len.", file=sys.stderr)

if __name__ == "__main__":
    main()

