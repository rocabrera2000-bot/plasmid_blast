#!/usr/bin/env python3
"""
Build a reference-anchored MSA from NUCmer blocks and append accessory clusters.

Backbone: same as v1 (place all blocks that involve the chosen MSA reference).
Accessory: cluster non-reference-only blocks across plasmids, orient, optionally
MAFFT-align per cluster, and append after the backbone with separators.

Outputs:
  <out>.aln.fasta                 Reference+Accessory pan-MSA
  <out>_conflicts.tsv             Backbone overwrite conflicts
  <out>_absent_in_ref.tsv         Non-ref hits not touching ref-anchored regions
  <out>_sanity_missing.tsv        Expected A<->B links (given shared ref region) that are absent
  <out>_accessory_clusters.tsv    Cluster summary (size, members, appended columns)
  <out>_panmap.tsv                MSA segments (REF or ACC_C###) -> column ranges
  <out>_acc/                      Per-cluster FASTAs and (optionally) aligned FASTAs
"""

from __future__ import annotations
import argparse, csv, os, sys, subprocess, tempfile, shutil
from pathlib import Path
from typing import Dict, List, Tuple, Iterable, Optional, Set, Any
from collections import defaultdict

def die(msg: str, code: int = 2):
    print(f"[ERROR] {msg}", file=sys.stderr); sys.exit(code)

# --------------------------- FASTA utils ---------------------------

def read_fasta_one(path: Path) -> Tuple[str, str]:
    if not path.exists(): die(f"FASTA not found: {path}")
    header = None; chunks: List[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if header is None: header = line[1:].strip()
                else: break
            else: chunks.append(line.strip())
    if header is None: die(f"No FASTA header in {path}")
    seq = "".join(chunks).replace(" ", "").upper()
    if not seq: die(f"Empty sequence in {path}")
    return header, seq

_rc_map = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn",
                        "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")
def revcomp(s: str) -> str: return s.translate(_rc_map)[::-1]
def subseq_1based(seq: str, s: int, e: int) -> str:
    if s <= 0 or e <= 0: die(f"Invalid 1-based coords: {s}-{e}")
    return seq[s-1:e] if s <= e else seq[e-1:s]
def as_int(d: Dict[str, str], k: str) -> int:
    try: return int(d[k])
    except: die(f"Bad integer in '{k}': {d.get(k)}")
def as_float(d: Dict[str, str], k: str) -> float:
    try: return float(d[k])
    except: die(f"Bad float in '{k}': {d.get(k)}")

# --------------------------- intervals ---------------------------

def interval_overlap(a: Tuple[int,int], b: Tuple[int,int]) -> int:
    a1,a2 = min(a),max(a); b1,b2 = min(b),max(b)
    left, right = max(a1,b1), min(a2,b2)
    return max(0, right - left + 1)

def merge_intervals(ivs: List[Tuple[int,int]], slop: int = 0) -> List[Tuple[int,int]]:
    if not ivs: return []
    ivs = [tuple(sorted(iv)) for iv in ivs]
    ivs.sort()
    out = [ivs[0]]
    for s,e in ivs[1:]:
        ps,pe = out[-1]
        if s <= pe + 1 + slop:
            out[-1] = (ps, max(pe, e))
        else:
            out.append((s,e))
    return out

def interval_any_overlap(iv: Tuple[int,int], ivs: Iterable[Tuple[int,int]]) -> bool:
    for x in ivs:
        if interval_overlap(iv, x) > 0: return True
    return False

# --------------------------- IO ---------------------------

def load_group_members(group_report: Path, group_name: str) -> List[str]:
    members: List[str] = []
    with group_report.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "group" not in reader.fieldnames or "members" not in reader.fieldnames:
            die("group_report.tsv must have columns: group, size, members")
        for row in reader:
            if row["group"] == group_name:
                raw = row["members"].strip()
                if not raw: die(f"Group '{group_name}' has empty members list.")
                members = [m.strip() for m in raw.split(",") if m.strip()]
                break
    if not members: die(f"Group '{group_name}' not found in {group_report}")
    return members

def read_coords(coords_path: Path) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    with coords_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        needed = {"ref_file","query_file","S1","E1","S2","E2","%_IDY","LEN_1","LEN_2"}
        if not needed.issubset(reader.fieldnames or []):
            die(f"{coords_path} missing required columns: {sorted(needed)}")
        for r in reader:
            r["_S1"]=as_int(r,"S1"); r["_E1"]=as_int(r,"E1")
            r["_S2"]=as_int(r,"S2"); r["_E2"]=as_int(r,"E2")
            r["_LEN_1"]=as_int(r,"LEN_1"); r["_LEN_2"]=as_int(r,"LEN_2")
            r["_PIDY"]=as_float(r,"%_IDY")
            # convenience
            r["_A_iv"]=(min(r["_S1"],r["_E1"]),max(r["_S1"],r["_E1"]))
            r["_B_iv"]=(min(r["_S2"],r["_E2"]),max(r["_S2"],r["_E2"]))
            r["_ALEN"]=min(abs(r["_E1"]-r["_S1"])+1, abs(r["_E2"]-r["_S2"])+1)
            # orientation parity between sides (+1 same, -1 reverse)
            strand1 = 1 if r["_E1"]>=r["_S1"] else -1
            strand2 = 1 if r["_E2"]>=r["_S2"] else -1
            r["_PARITY"] = 1 if strand1==strand2 else -1
            rows.append(r)
    return rows

def filter_rows_for_group(rows: List[Dict[str,Any]],
                          members: Set[str],
                          min_len: int,
                          min_id: float) -> List[Dict[str,Any]]:
    out = []
    for r in rows:
        if r["ref_file"] not in members or r["query_file"] not in members: continue
        if r["_ALEN"] < min_len: continue
        if r["_PIDY"] < min_id: continue
        out.append(r)
    return out

# --------------------------- core MSA backbone ---------------------------

def build_backbone(rows: List[Dict[str,Any]],
                   members: List[str],
                   fasta_dir: Path,
                   msa_reference: str):

    seqs: Dict[str,str] = {}; headers: Dict[str,str] = {}
    for name in members:
        h,s = read_fasta_one(fasta_dir / name)
        seqs[name]=s; headers[name]=h

    if msa_reference not in seqs: die(f"--msa-reference '{msa_reference}' not in members")
    ref_seq = seqs[msa_reference]; ref_len = len(ref_seq)
    msa: Dict[str, List[str]] = {n: ["-"]*ref_len for n in members}
    msa[msa_reference] = list(ref_seq)

    involve_ref, non_ref_rows = [], []
    for r in rows:
        if r["ref_file"]==msa_reference or r["query_file"]==msa_reference:
            involve_ref.append(r)
        else:
            non_ref_rows.append(r)

    def ref_interval(r):
        if r["ref_file"]==msa_reference: a,b=r["_S1"],r["_E1"]
        else: a,b=r["_S2"],r["_E2"]
        return (min(a,b), max(a,b))

    involve_ref.sort(key=lambda r: ref_interval(r)[0])

    conflicts: List[Tuple[str,int,str,str,str,str,int,int,int,int,float]] = []
    placed_by_plasmid: Dict[str, List[Tuple[Tuple[int,int], Tuple[int,int]]]] = defaultdict(list)

    # place blocks
    for r in involve_ref:
        if r["ref_file"]==msa_reference:
            r1,r2 = r["_S1"], r["_E1"]
            other = r["query_file"]; p1,p2 = r["_S2"], r["_E2"]
        else:
            r1,r2 = r["_S2"], r["_E2"]
            other = r["ref_file"];  p1,p2 = r["_S1"], r["_E1"]

        R=(min(r1,r2), max(r1,r2)); P=(min(p1,p2), max(p1,p2))
        seg = subseq_1based(seqs[other], P[0], P[1])
        if p2 < p1: seg = revcomp(seg)
        if (R[1]-R[0]+1)!=len(seg):
            die(f"Length mismatch placing {other}: ref {R} vs seg len {len(seg)}")

        row_chars = msa[other]; s0 = R[0]-1
        for i,base in enumerate(seg):
            pos = s0 + i
            if row_chars[pos] in ("-", base):
                row_chars[pos]=base
            else:
                conflicts.append((other, pos+1, row_chars[pos], base,
                                  r["ref_file"], r["query_file"],
                                  r["_S1"], r["_E1"], r["_S2"], r["_E2"], r["_PIDY"]))
        placed_by_plasmid[other].append((P, R))

    placed_by_plasmid[msa_reference].append(((1, ref_len), (1, ref_len)))

    return msa, seqs, placed_by_plasmid, non_ref_rows, conflicts, ref_len

# --------------------------- accessory clustering ---------------------------

def path_which(cmd: str) -> Optional[str]:
    return shutil.which(cmd)

def choose_aligner(name: str) -> Optional[str]:
    if name == "none": return None
    if name == "mafft":
        exe = path_which("mafft")
        if exe: return exe
        print("[WARN] mafft not found; falling back to 'none'")
        return None
    die(f"Unknown aligner: {name}")

def organize_accessory(msa: Dict[str,List[str]],
                       seqs: Dict[str,str],
                       members: List[str],
                       non_ref_rows: List[Dict[str,Any]],
                       placed_by_plasmid: Dict[str, List[Tuple[Tuple[int,int], Tuple[int,int]]]],
                       out_prefix: Path,
                       acc_min_len: int,
                       acc_min_id: float,
                       acc_merge_slop: int,
                       acc_separator: int,
                       acc_aligner: Optional[str]):

    outdir = out_prefix.parent; outdir.mkdir(parents=True, exist_ok=True)
    acc_dir = outdir / (out_prefix.name + "_acc"); acc_dir.mkdir(exist_ok=True)

    # 1) "Covered" IVs per plasmid (those that touched the reference)
    covered: Dict[str,List[Tuple[int,int]]] = {}
    for p, lst in placed_by_plasmid.items():
        covered[p] = merge_intervals([piv for (piv,_miv) in lst], slop=0)

    # 2) Candidate accessory IVs per plasmid (from non_ref_rows, not overlapping covered)
    cand_per_plasmid: Dict[str,List[Tuple[int,int]]] = defaultdict(list)
    # collect raw
    for r in non_ref_rows:
        A,B = r["ref_file"], r["query_file"]
        A_iv, B_iv = r["_A_iv"], r["_B_iv"]
        if r["_ALEN"]>=acc_min_len and r["_PIDY"]>=acc_min_id:
            if not interval_any_overlap(A_iv, covered.get(A, [])):
                cand_per_plasmid[A].append(A_iv)
            if not interval_any_overlap(B_iv, covered.get(B, [])):
                cand_per_plasmid[B].append(B_iv)
    # merge per plasmid
    for p in list(cand_per_plasmid.keys()):
        cand_per_plasmid[p] = merge_intervals(cand_per_plasmid[p], slop=acc_merge_slop)

    # quick finder: which merged candidate IV contains/overlaps a given IV?
    def find_hit(plasmid: str, iv: Tuple[int,int]) -> Optional[int]:
        for i,(s,e) in enumerate(cand_per_plasmid.get(plasmid, [])):
            if interval_overlap((s,e), iv) > 0:
                return i
        return None

    # 3) Build graph of candidate nodes with parity edges from non_ref_rows
    # node key = (plasmid, idx)
    edges: Dict[Tuple[str,int], List[Tuple[Tuple[str,int], int]]] = defaultdict(list)
    support: Dict[Tuple[Tuple[str,int],Tuple[str,int]], int] = defaultdict(int)

    for r in non_ref_rows:
        if r["_ALEN"]<acc_min_len or r["_PIDY"]<acc_min_id: continue
        A,B = r["ref_file"], r["query_file"]
        ia = find_hit(A, r["_A_iv"])
        ib = find_hit(B, r["_B_iv"])
        if ia is None or ib is None: continue
        na, nb = (A,ia), (B,ib)
        parity = r["_PARITY"]  # +1 same, -1 reverse
        edges[na].append((nb, parity))
        edges[nb].append((na, parity))
        k = tuple(sorted((na,nb)))  # undirected pair
        support[k] += 1

    # 4) Connected components + orientation assignment by BFS
    visited: Set[Tuple[str,int]] = set()
    components: List[List[Tuple[str,int]]] = []
    orientation: Dict[Tuple[str,int], int] = {}  # +1 keep, -1 revcomp
    par_conflicts = 0

    for node in edges.keys():
        if node in visited: continue
        comp = []
        stack = [(node, 1)]
        orientation[node]=1; visited.add(node)
        while stack:
            u, uori = stack.pop()
            comp.append(u)
            for v, parity in edges[u]:
                if v not in orientation:
                    orientation[v] = uori * parity
                    visited.add(v)
                    stack.append((v, orientation[v]))
                else:
                    if orientation[v] != uori * parity:
                        par_conflicts += 1
        if comp: components.append(comp)

    # 5) Prepare to append clusters
    # Current MSA length is length of any row (backbone length)
    pan_len = len(next(iter(msa.values())))
    sep_cols = ["-"] * acc_separator

    # writers
    clusters_tsv = outdir / (out_prefix.name + "_accessory_clusters.tsv")
    panmap_tsv = outdir / (out_prefix.name + "_panmap.tsv")

    with clusters_tsv.open("w") as ct, panmap_tsv.open("w") as pm:
        ct.write("cluster_id\tsize\tmembers\taln_cols\n")
        pm.write("segment_id\tmsa_start\tmsa_end\tnote\n")
        pm.write(f"REF\t1\t{pan_len}\tbackbone\n")

        # iterate components by decreasing size
        comp_order = sorted(range(len(components)), key=lambda i: len(components[i]), reverse=True)
        for ci in comp_order:
            comp = components[ci]
            members_in_comp = sorted(set(p for (p,_i) in comp))
            cluster_id = f"ACC_C{ci+1:03d}"

            # extract oriented sequences per node
            seqs_in = []
            for (p, idx) in comp:
                s,e = cand_per_plasmid[p][idx]
                seg = subseq_1based(seqs[p], s, e)
                if orientation.get((p,idx),1) == -1:
                    seg = revcomp(seg)
                seqs_in.append((p, f"{p}:{s}-{e}", seg))

            # write raw cluster fasta
            raw_fa = acc_dir / f"{cluster_id}.raw.fasta"
            with raw_fa.open("w") as f:
                for p,lab,seg in seqs_in:
                    f.write(f">{p}|{lab}\n")
                    for i in range(0,len(seg),80): f.write(seg[i:i+80]+"\n")

            # align (optional)
            aligned: List[Tuple[str,str]] = []
            if acc_aligner is not None and len(seqs_in) >= 2:
                with tempfile.TemporaryDirectory() as td:
                    in_fa = Path(td) / "in.fa"
                    with in_fa.open("w") as f:
                        for p,lab,seg in seqs_in:
                            f.write(f">{p}\n{seg}\n")
                    try:
                        if "mafft" in acc_aligner:
                            # quiet, auto
                            result = subprocess.run([acc_aligner, "--auto", "--quiet", str(in_fa)],
                                                    check=True, capture_output=True, text=True)
                            out = result.stdout.strip().splitlines()
                            name=None; buf=[]
                            for line in out:
                                if line.startswith(">"):
                                    if name is not None:
                                        aligned.append((name,"".join(buf)))
                                    name = line[1:].strip(); buf=[]
                                else:
                                    buf.append(line.strip())
                            if name is not None: aligned.append((name,"".join(buf)))
                        else:
                            aligned = []
                    except Exception as e:
                        print(f"[WARN] MAFFT failed for {cluster_id}: {e}; using left-justify fallback.")
                        aligned = []
            if not aligned:
                # left-justify + pad to max
                maxlen = max(len(seg) for _p,_lab,seg in seqs_in)
                aligned = [(p, seg + "-"*(maxlen-len(seg))) for (p,_lab,seg) in seqs_in]

            # write aligned cluster fasta
            aln_fa = acc_dir / f"{cluster_id}.aln.fasta"
            with aln_fa.open("w") as f:
                for p,aseq in aligned:
                    f.write(f">{p}\n")
                    for i in range(0,len(aseq),80): f.write(aseq[i:i+80]+"\n")

            # append separator
            if acc_separator>0:
                for n in msa: msa[n].extend(sep_cols)
                pan_len += acc_separator

            # determine aligned width
            width = len(aligned[0][1])
            start = pan_len + 1
            end = pan_len + width

            # append cluster columns (fill gaps for non-members)
            present = {p for p,_a in aligned}
            for n in msa:
                if n in present:
                    aseq = next(a for (p,a) in aligned if p==n)
                    msa[n].extend(list(aseq))
                else:
                    msa[n].extend(["-"]*width)
            pan_len = end

            # summaries
            ct.write(f"{cluster_id}\t{len(members_in_comp)}\t{','.join(members_in_comp)}\t{width}\n")
            pm.write(f"{cluster_id}\t{start}\t{end}\taccessory_cluster\n")

    if par_conflicts:
        print(f"[WARN] Orientation parity conflicts detected in {par_conflicts} edges (cycles with mixed parity). "
              f"Cluster orientations chosen by BFS; inspect per-cluster raw/aln FASTAs if something looks off.")

    return msa  # extended to pan-MSA

# --------------------------- sanity check & writers ---------------------------

def write_outputs(msa: Dict[str,List[str]],
                  members: List[str],
                  out_prefix: Path,
                  conflicts: List[Tuple],
                  absent_records: List[Tuple],
                  sanity_missing: List[Tuple]):

    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    # FASTA
    msa_fa = out_prefix.with_suffix(".aln.fasta")
    with msa_fa.open("w") as out:
        for name in members:
            out.write(f">{name}\n")
            row = "".join(msa[name])
            for i in range(0, len(row), 80):
                out.write(row[i:i+80] + "\n")
    print(f"[INFO] Wrote pan-MSA: {msa_fa} (cols={len(next(iter(msa.values())))})")

    # Conflicts
    conflicts_tsv = out_prefix.parent / (out_prefix.name + "_conflicts.tsv")
    with conflicts_tsv.open("w") as out:
        out.write("plasmid\tmsa_pos\texisting\tnew\tref_file\tquery_file\tS1\tE1\tS2\tE2\tpct_idy\n")
        for rec in conflicts: out.write("\t".join(map(str, rec)) + "\n")
    print(f"[INFO] Conflicts: {len(conflicts)} -> {conflicts_tsv}")

    # Absent-in-reference
    absent_tsv = out_prefix.parent / (out_prefix.name + "_absent_in_ref.tsv")
    with absent_tsv.open("w") as out:
        out.write("A\tB\tA_start\tA_end\tB_start\tB_end\tpct_idy\talen\n")
        for rec in absent_records: out.write("\t".join(map(str, rec)) + "\n")
    print(f"[INFO] Absent-in-ref blocks: {len(absent_records)} -> {absent_tsv}")

    # Sanity missing
    sanity_tsv = out_prefix.parent / (out_prefix.name + "_sanity_missing.tsv")
    with sanity_tsv.open("w") as out:
        out.write("A\tB\tmsa_start\tmsa_end\toverlap_len\n")
        for rec in sanity_missing: out.write("\t".join(map(str, rec)) + "\n")
    print(f"[INFO] Sanity missing links: {len(sanity_missing)} -> {sanity_tsv}")

# --------------------------- main ---------------------------

def main():
    ap = argparse.ArgumentParser(description="Reference-anchored MSA + accessory clusters from NUCmer blocks")
    ap.add_argument("--coords", required=True, type=Path)
    ap.add_argument("--group-report", required=True, type=Path)
    ap.add_argument("--group", required=True)
    ap.add_argument("--fasta-dir", required=True, type=Path)
    ap.add_argument("--msa-reference", required=True)
    ap.add_argument("--out-prefix", required=True, type=Path)
    ap.add_argument("--min-len", default=1, type=int)
    ap.add_argument("--min-id", default=0.0, type=float)
    ap.add_argument("--sanity-min-overlap", default=200, type=int)

    # accessory opts
    ap.add_argument("--accessory", action="store_true", default=True,
                    help="Append accessory clusters (default: on)")
    ap.add_argument("--no-accessory", dest="accessory", action="store_false",
                    help="Disable accessory clustering/appending")
    ap.add_argument("--acc-aligner", choices=["mafft","none"], default="mafft")
    ap.add_argument("--acc-merge-slop", type=int, default=50,
                    help="bp to merge nearby candidate blocks per plasmid (default 50)")
    ap.add_argument("--acc-separator", type=int, default=50,
                    help="gap columns between segments (default 50)")
    ap.add_argument("--acc-min-len", type=int, default=200,
                    help="min block length to consider for accessory (default 200)")
    ap.add_argument("--acc-min-id", type=float, default=90.0,
                    help="min %%identity to consider for accessory (default 90)")

    args = ap.parse_args()

    members = load_group_members(args.group_report, args.group)
    rows = filter_rows_for_group(read_coords(args.coords), set(members), args.min_len, args.min_id)

    print(f"[INFO] Using {len(rows)} blocks within '{args.group}' (min_len={args.min_len}, min_id={args.min_id})")

    # backbone
    msa, seqs, placed_by_plasmid, non_ref_rows, conflicts, ref_len = build_backbone(
        rows, members, args.fasta_dir, args.msa_reference
    )

    # absent-in-ref list (for reporting)
    absent_records: List[Tuple[str,str,int,int,int,int,float,int]] = []
    for r in non_ref_rows:
        A,B = r["ref_file"], r["query_file"]
        A_iv, B_iv = r["_A_iv"], r["_B_iv"]
        A_over = interval_any_overlap(A_iv, (piv for (piv,_miv) in placed_by_plasmid.get(A, [])))
        B_over = interval_any_overlap(B_iv, (piv for (piv,_miv) in placed_by_plasmid.get(B, [])))
        if not A_over and not B_over:
            absent_records.append((A,B,A_iv[0],A_iv[1],B_iv[0],B_iv[1], r["_PIDY"], r["_ALEN"]))

    # sanity check
    pair_index: Dict[frozenset, List[Dict[str,Any]]] = defaultdict(list)
    for r in non_ref_rows: pair_index[frozenset((r["ref_file"], r["query_file"]))].append(r)

    sanity_missing: List[Tuple[str,str,int,int,int]] = []
    members_no_ref = [m for m in members if m != args.msa_reference]
    for i in range(len(members_no_ref)):
        A = members_no_ref[i]
        A_mivs = [miv for (_piv, miv) in placed_by_plasmid.get(A, [])]
        if not A_mivs: continue
        for j in range(i+1, len(members_no_ref)):
            B = members_no_ref[j]
            B_mivs = [miv for (_piv, miv) in placed_by_plasmid.get(B, [])]
            if not B_mivs: continue
            for a_m in A_mivs:
                for b_m in B_mivs:
                    ovl = interval_overlap(a_m, b_m)
                    if ovl >= args.sanity_min_overlap:
                        rows_ab = pair_index.get(frozenset((A,B)), [])
                        found=False
                        for r in rows_ab:
                            A_iv = (min(r["_S1"],r["_E1"]), max(r["_S1"],r["_E1"])) if r["ref_file"]==A else (min(r["_S2"],r["_E2"]), max(r["_S2"],r["_E2"]))
                            B_iv = (min(r["_S1"],r["_E1"]), max(r["_S1"],r["_E1"])) if r["ref_file"]==B else (min(r["_S2"],r["_E2"]), max(r["_S2"],r["_E2"]))
                            A_has = interval_any_overlap(A_iv, (piv for (piv,_miv) in placed_by_plasmid.get(A, [])))
                            B_has = interval_any_overlap(B_iv, (piv for (piv,_miv) in placed_by_plasmid.get(B, [])))
                            if A_has and B_has:
                                found=True; break
                        if not found:
                            sanity_missing.append((A,B, max(a_m[0],b_m[0]), min(a_m[1],b_m[1]), ovl))

    # accessory clustering / appending
    if args.accessory:
        aligner_exe = choose_aligner(args.acc_aligner)
        msa = organize_accessory(
            msa, seqs, members, non_ref_rows, placed_by_plasmid, args.out_prefix,
            args.acc_min_len, args.acc_min_id, args.acc_merge_slop, args.acc_separator,
            aligner_exe
        )

    # write everything
    write_outputs(msa, members, args.out_prefix, conflicts, absent_records, sanity_missing)

if __name__ == "__main__":
    main()
