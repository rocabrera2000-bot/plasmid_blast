#!/usr/bin/env python3
"""
Reference-anchored MSA from NUCmer blocks, with gapped block support.

- Filters all_coords.tsv to a chosen group.
- Builds an MSA canvas of length = MSA reference, fills reference with its sequence.
- Places all blocks involving the reference; reverse-complements query segments when needed.
- If a block has an indel (span lengths differ), performs a local column-wise paste using edlib
  (global alignment) so insertions/deletions are respected relative to the reference.
- Overwrites only gaps; logs conflicts when a non-gap differs from an incoming base.
- Reports:
  * *_conflicts.tsv
  * *_absent_in_ref.tsv (blocks seen only among non-ref pairs, not overlapping any region that mapped to ref)
  * *_sanity_missing.tsv (A and B overlap on the reference but no A↔B block supports that overlap)

Usage example:
python build_block_msa_v3.py \
  --coords ../polished_plasmids/nucmer/all_coords.tsv \
  --group-report results_40000/group_report.tsv \
  --group group_40kb_4 \
  --fasta-dir results_40000/group_40kb_4 \
  --msa-reference E1112_hybrid2025_ctg_3_len_89633_blaNDM1.fasta \
  --out-prefix out/group_40kb_4/E1112_ref \
  --min-len 200 --min-id 90

"""


from __future__ import annotations
import argparse, csv, sys
from pathlib import Path
from typing import Dict, List, Tuple, Iterable, Optional, Set, Any
from collections import defaultdict

# --------------------------- Optional edlib ---------------------------
_HAVE_EDLIB = False
try:
    import edlib  # pip install edlib
    _HAVE_EDLIB = True
except Exception:
    _HAVE_EDLIB = False

# --------------------------- utils ---------------------------

def die(msg: str, code: int = 2):
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(code)

def log(msg: str):
    print(f"[INFO] {msg}", file=sys.stderr)

def warn(msg: str):
    print(f"[WARN] {msg}", file=sys.stderr)

def load_group_members(group_report: Path, group_name: str) -> List[str]:
    members: List[str] = []
    with group_report.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        need = {"group", "members"}
        if not need.issubset(reader.fieldnames or []):
            die("group_report.tsv must have columns: group, size, members")
        for row in reader:
            if row["group"] == group_name:
                raw = (row.get("members") or "").strip()
                if not raw:
                    die(f"Group '{group_name}' found but members is empty.")
                members = [m.strip() for m in raw.split(",") if m.strip()]
                break
    if not members:
        die(f"Group '{group_name}' not found in {group_report}")
    return members

def read_fasta_one(path: Path) -> Tuple[str, str]:
    if not path.exists():
        die(f"FASTA not found: {path}")
    header = None
    chunks: List[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if header is None:
                    header = line[1:].strip()
                else:
                    break  # only first sequence
            else:
                chunks.append(line.strip())
    if header is None:
        die(f"No FASTA header found in {path}")
    seq = "".join(chunks).replace(" ", "").upper()
    if not seq:
        die(f"Empty sequence in {path}")
    return (header, seq)

_rc_map = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn",
                        "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")

def revcomp(s: str) -> str:
    return s.translate(_rc_map)[::-1]

def subseq_1based(seq: str, start: int, end: int) -> str:
    if start <= 0 or end <= 0:
        die(f"Invalid 1-based coords: {start}-{end}")
    if start <= end:
        return seq[start-1:end]
    else:
        return seq[end-1:start]  # caller handles orientation

def as_int(d: Dict[str, str], k: str) -> int:
    try:
        return int(d[k])
    except Exception:
        die(f"Bad integer in column '{k}': {d.get(k)}")

def as_float(d: Dict[str, str], k: str) -> float:
    try:
        return float(d[k])
    except Exception:
        die(f"Bad float in column '{k}': {d.get(k)}")

def interval_overlap(a: Tuple[int,int], b: Tuple[int,int]) -> int:
    a1, a2 = min(a), max(a)
    b1, b2 = min(b), max(b)
    left = max(a1, b1); right = min(a2, b2)
    return max(0, right - left + 1)

def interval_any_overlap(iv: Tuple[int,int], ivs: Iterable[Tuple[int,int]]) -> bool:
    for x in ivs:
        if interval_overlap(iv, x) > 0:
            return True
    return False


# --- conflict resolution helpers (ADD) ---
EPS = 1e-9

# per-plasmid per-MSA-position provenance of the base currently in the MSA
# winner_meta_by_pos[plasmid][pos] = {
#   "pct_idy": float, "block_len": int, "ref_start": int, "ref_end": int,
#   "plasmid_start": int, "plasmid_end": int, "ref_file": str, "query_file": str,
#   "S1": int, "E1": int, "S2": int, "E2": int, "role": str
# }
winner_meta_by_pos: Dict[str, Dict[int, Dict[str, Any]]] = defaultdict(dict)

def _should_overwrite(existing_meta: Optional[Dict[str,Any]],
                      new_meta: Dict[str,Any]) -> Tuple[bool, bool]:
    """
    Returns (overwrite?, special_warn?).
    Priority: higher identity wins; if equal (within EPS), longer block wins.
    special_warn=True when the winner has higher identity but *shorter* block.
    """
    if existing_meta is None:
        return True, False
    a = new_meta; b = existing_meta
    if a["pct_idy"] > b["pct_idy"] + EPS:
        return True, a["block_len"] < b["block_len"]
    if b["pct_idy"] > a["pct_idy"] + EPS:
        return False, False
    # tie on identity -> prefer longer
    if a["block_len"] > b["block_len"]:
        return True, False
    if b["block_len"] > a["block_len"]:
        return False, False
    return False, False  # exact tie -> keep existing

def _place_base_with_resolution(msa_row: List[str], plasmid: str, pos0: int, base: str,
                                conflict_sink: List[Tuple], new_meta: Dict[str,Any],
                                warn_fn=warn):
    """
    Resolve conflicts using identity, then length. Log conflicts either way.
    """
    existing = msa_row[pos0]
    if existing == '-' or existing == base:
        msa_row[pos0] = base
        winner_meta_by_pos[plasmid][pos0] = new_meta
        return

    ex_meta = winner_meta_by_pos[plasmid].get(pos0)
    overwrite, special = _should_overwrite(ex_meta, new_meta)

    if overwrite:
        # warn on smaller-but-higher-identity winner
        if special and ex_meta is not None:
            warn_fn(
                f"Choosing smaller but higher-identity block for {plasmid}: "
                f"new {new_meta['pct_idy']:.2f}% len={new_meta['block_len']} "
                f"(ref {new_meta['ref_start']}-{new_meta['ref_end']}, "
                f"qry {new_meta['plasmid_start']}-{new_meta['plasmid_end']}, {new_meta['role']}) "
                f"over old {ex_meta['pct_idy']:.2f}% len={ex_meta['block_len']}"
            )
        msa_row[pos0] = base
        winner_meta_by_pos[plasmid][pos0] = new_meta

    # In both cases, record the conflict attempt (kept base is whatever sits in MSA after this)
    conflict_sink.append((
        plasmid, pos0+1, existing, base,
        new_meta['ref_file'], new_meta['query_file'],
        new_meta['S1'], new_meta['E1'], new_meta['S2'], new_meta['E2'],
        new_meta['pct_idy'], new_meta['role']
    ))

# --------------------------- I/O ---------------------------

def read_coords(coords_path: Path) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    with coords_path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        need = {"ref_file","query_file","S1","E1","S2","E2","%_IDY","LEN_1","LEN_2"}
        if not need.issubset(reader.fieldnames or []):
            die(f"{coords_path} missing required columns: {sorted(need)}")
        for r in reader:
            r["_S1"] = as_int(r, "S1")
            r["_E1"] = as_int(r, "E1")
            r["_S2"] = as_int(r, "S2")
            r["_E2"] = as_int(r, "E2")
            r["_LEN_1"] = as_int(r, "LEN_1")
            r["_LEN_2"] = as_int(r, "LEN_2")
            r["_PIDY"] = as_float(r, "%_IDY")
            rows.append(r)
    return rows

def filter_rows_for_group(rows: List[Dict[str,Any]],
                          members: Set[str],
                          min_len: int,
                          min_id: float) -> List[Dict[str,Any]]:
    out = []
    for r in rows:
        rf = r["ref_file"]; qf = r["query_file"]
        if rf not in members or qf not in members:
            continue
        alen = min(abs(r["_E1"] - r["_S1"]) + 1, abs(r["_E2"] - r["_S2"]) + 1)
        if alen < min_len:
            continue
        if r["_PIDY"] < min_id:
            continue
        r["_ALEN"] = alen
        out.append(r)
    return out

# --------------------------- gapped pasting ---------------------------

def paste_ungapped(msa_row: List[str], ref_start_1based: int, segment: str,
                   conflict_sink: List[Tuple],
                   plasmid: str, ref_file: str, query_file: str,
                   S1: int, E1: int, S2: int, E2: int, pct_idy: float,
                   role_str: str, block_len: int,
                   # NEW params for equal-span edlib:
                   ref_seq_slice: Optional[str] = None,
                   plasmid_start: Optional[int] = None,
                   plasmid_end: Optional[int] = None) -> str:
    """
    Returns mode used: 'ungapped' or 'edlib'.
    If ref_seq_slice is provided and edlib finds internal gaps (even when spans equal),
    fall back to edlib path, warn, and report as edlib.
    """
    # Try edlib if we have the reference slice and edlib available
    if ref_seq_slice is not None and _HAVE_EDLIB:
        res = edlib.align(segment, ref_seq_slice, mode="NW", task="path")
        nice = edlib.getNiceAlignment(res, segment, ref_seq_slice)
        if ('-' in nice["query_aligned"]) or ('-' in nice["target_aligned"]):
            warn(f"Equal-span gaps detected -> using edlib for {plasmid} "
                 f"block ref {ref_start_1based}-{ref_start_1based+len(ref_seq_slice)-1} "
                 f"(qry {plasmid_start}-{plasmid_end}, {role_str})")
            # make sure role says '|edlib' (not '|ungapped')
            role_local = role_str.replace("ungapped", "edlib") if "|ungapped" in role_str else f"{role_str}|edlib"
            paste_gapped_edlib(msa_row, ref_seq_slice, segment, ref_start_1based,
                               plasmid, ref_file, query_file,
                               S1, E1, S2, E2, pct_idy,
                               role_local, block_len,
                               plasmid_start or -1, plasmid_end or -1,
                               conflict_sink=conflict_sink) 
            return "edlib"


    # Plain 1:1 pasting with conflict resolution
    s0 = ref_start_1based - 1
    new_meta_common = {
        "pct_idy": pct_idy, "block_len": block_len,
        "ref_start": ref_start_1based, "ref_end": ref_start_1based + len(segment) - 1,
        "plasmid_start": plasmid_start if plasmid_start is not None else -1,
        "plasmid_end": plasmid_end if plasmid_end is not None else -1,
        "ref_file": ref_file, "query_file": query_file,
        "S1": S1, "E1": E1, "S2": S2, "E2": E2,
        "role": role_str,
    }
    for i, base in enumerate(segment):
        pos = s0 + i
        _place_base_with_resolution(msa_row, plasmid, pos, base,
                                    conflict_sink=conflict_sink,
                                    new_meta=new_meta_common)
    return "ungapped"



def paste_gapped_edlib(msa_row: List[str],
                       ref_seq_slice: str,
                       qry_seq_slice: str,
                       ref_start_1based: int,
                       plasmid: str, ref_file: str, query_file: str,
                       S1: int, E1: int, S2: int, E2: int, pct_idy: float,
                       role_str: str, block_len: int,
                       plasmid_start: int, plasmid_end: int,
                       conflict_sink: List[Tuple]):   # <-- ADDED
    """
    Align qry vs ref globally (edlib mode "NW"), then iterate aligned columns.
    - Advance MSA position when REF column is a base.
    - Write QUERY base using conflict resolution.
    """
    if not _HAVE_EDLIB:
        raise RuntimeError("edlib not available")

    res = edlib.align(qry_seq_slice, ref_seq_slice, mode="NW", task="path")
    nice = edlib.getNiceAlignment(res, qry_seq_slice, ref_seq_slice)
    q_aln = nice["query_aligned"]; r_aln = nice["target_aligned"]

    pos = ref_start_1based - 1  # 0-based MSA cursor
    assert len(q_aln) == len(r_aln)
    new_meta_common = {
        "pct_idy": pct_idy, "block_len": block_len,
        "ref_start": ref_start_1based, "ref_end": ref_start_1based + len(ref_seq_slice) - 1,
        "plasmid_start": plasmid_start, "plasmid_end": plasmid_end,
        "ref_file": ref_file, "query_file": query_file,
        "S1": S1, "E1": E1, "S2": S2, "E2": E2,
        "role": role_str,
    }

    for qi, ri in zip(q_aln, r_aln):
        if ri != '-':
            if qi != '-':
                _place_base_with_resolution(msa_row, plasmid, pos, qi,
                                            conflict_sink=conflict_sink,
                                            new_meta=new_meta_common)
            pos += 1
        # if ri == '-' (insertion vs ref), we don't advance MSA and we don't add columns
# --------------------------- core ---------------------------

def build_msa(rows: List[Dict[str,Any]],
              members: List[str],
              fasta_dir: Path,
              msa_reference: str,
              out_prefix: Path,
              sanity_min_overlap: int,
              gap_mode: str) -> None:
    # local conflict sink for this build
    conflicts: List[Tuple] = []
    # reset per-run provenance store
    winner_meta_by_pos.clear()

    members_set = set(members)
    if msa_reference not in members_set:
        die(f"--msa-reference '{msa_reference}' not in group members ({len(members)} members).")

    # Load sequences
    seqs: Dict[str,str] = {}
    headers: Dict[str,str] = {}
    for name in members:
        header, seq = read_fasta_one(fasta_dir / name)
        seqs[name] = seq
        headers[name] = header

    ref_seq = seqs[msa_reference]
    ref_len = len(ref_seq)
    log(f"MSA reference: {msa_reference} (len={ref_len}); group_size={len(members)}; gap_mode={gap_mode}"
        + (" (edlib missing → skip)" if gap_mode == "edlib" and not _HAVE_EDLIB else ""))

    # MSA canvas
    msa: Dict[str, List[str]] = {name: ["-"] * ref_len for name in members}
    msa[msa_reference] = list(ref_seq)

    # Track where each plasmid mapped to the ref (for absent-in-ref & sanity)
    placed_by_plasmid: Dict[str, List[Tuple[Tuple[int,int], Tuple[int,int]]]] = defaultdict(list)

    # Placement index (all blocks we actually placed)
    placed_blocks: List[Tuple[Any, ...]] = []

    # Split rows
    involve_ref: List[Dict[str,Any]] = []
    non_ref_rows: List[Dict[str,Any]] = []
    for r in rows:
        if r["ref_file"] == msa_reference or r["query_file"] == msa_reference:
            involve_ref.append(r)
        else:
            non_ref_rows.append(r)

    # Sort placements by ref interval start
    def ref_interval(r: Dict[str,Any]) -> Tuple[int,int]:
        if r["ref_file"] == msa_reference:
            a, b = r["_S1"], r["_E1"]
        else:
            a, b = r["_S2"], r["_E2"]
        return (min(a,b), max(a,b))
    involve_ref.sort(key=lambda r: ref_interval(r)[0])

    # Counters
    skipped_gapped = 0
    gapped_used = 0              # edlib used for unequal spans
    eqspan_edlib_used = 0        # edlib used for equal spans (slippage)

    # Place blocks involving the reference
    for r in involve_ref:
        if r["ref_file"] == msa_reference:
            r1, r2 = r["_S1"], r["_E1"]          # ref coords (raw)
            other = r["query_file"]
            o1, o2 = r["_S2"], r["_E2"]          # other coords (raw)
            pair_role = "ref_as_ref"
        else:
            r1, r2 = r["_S2"], r["_E2"]
            other = r["ref_file"]
            o1, o2 = r["_S1"], r["_E1"]
            pair_role = "ref_as_query"

        # --- robust strand decision BEFORE sorting
        ref_dir = 1 if (r2 - r1) >= 0 else -1
        qry_dir = 1 if (o2 - o1) >= 0 else -1
        rev_needed = (ref_dir * qry_dir) < 0

        # Sort for slicing only
        R1, R2 = (r1, r2) if r1 <= r2 else (r2, r1)
        P1, P2 = (o1, o2) if o1 <= o2 else (o2, o1)

        # Slice
        ref_raw = subseq_1based(ref_seq, R1, R2)
        qry_raw = subseq_1based(seqs[other], P1, P2)
        if rev_needed:
            qry_raw = revcomp(qry_raw)

        len_ref = R2 - R1 + 1
        len_qry = P2 - P1 + 1
        strand = '+' if o2 >= o1 else '-'

        # Default role label (may be promoted to |edlib inside paste_ungapped)
        role_str = f"{pair_role}|ungapped"

        if len_ref == len_qry:
            # Try equal-span edlib detection inside paste_ungapped (warn & switch if needed)
            mode_used = paste_ungapped(
                msa_row=msa[other],
                ref_start_1based=R1,
                segment=qry_raw,
                conflict_sink=conflicts,
                plasmid=other, ref_file=r["ref_file"], query_file=r["query_file"],
                S1=r["_S1"], E1=r["_E1"], S2=r["_S2"], E2=r["_E2"], pct_idy=r["_PIDY"],
                role_str=role_str, block_len=len_ref,
                ref_seq_slice=ref_raw, plasmid_start=P1, plasmid_end=P2
            )
            if mode_used == "edlib":
                eqspan_edlib_used += 1
            placed_by_plasmid[other].append(((P1, P2), (R1, R2)))
        else:
            # unequal spans
            if gap_mode == "die":
                die(f"Length mismatch placing {other}: ref ({R1}, {R2}) len={len_ref} vs seg len {len_qry}")
            elif gap_mode == "skip" or (gap_mode == "edlib" and not _HAVE_EDLIB):
                skipped_gapped += 1
                continue
            else:
                # edlib path
                role_str = f"{pair_role}|edlib"
                paste_gapped_edlib(
                    msa_row=msa[other],
                    ref_seq_slice=ref_raw,
                    qry_seq_slice=qry_raw,
                    ref_start_1based=R1,
                    plasmid=other, ref_file=r["ref_file"], query_file=r["query_file"],
                    S1=r["_S1"], E1=r["_E1"], S2=r["_S2"], E2=r["_E2"], pct_idy=r["_PIDY"],
                    role_str=role_str, block_len=len_ref,
                    plasmid_start=P1, plasmid_end=P2,
                    conflict_sink=conflicts
                )
                placed_by_plasmid[other].append(((P1, P2), (R1, R2)))
                gapped_used += 1
            mode_used = "edlib"

        # Log placed block (index)
        placed_blocks.append((
            other,                 # plasmid
            pair_role,             # role (without mode)
            mode_used,             # mode: 'ungapped' or 'edlib'
            strand,
            R1, R2, len_ref,
            P1, P2, len_qry,
            r["_PIDY"],
            r["ref_file"], r["query_file"],
            r["_S1"], r["_E1"], r["_S2"], r["_E2"]
        ))

    # Reference "placed" fully (useful for sanity)
    placed_by_plasmid[msa_reference].append(((1, ref_len), (1, ref_len)))

    # Absent-in-reference: non-ref pair blocks where neither side overlaps any ref-anchored interval
    absent_records: List[Tuple[str,str,int,int,int,int,float,int]] = []
    for r in non_ref_rows:
        A, B = r["ref_file"], r["query_file"]
        A_iv = (min(r["_S1"], r["_E1"]), max(r["_S1"], r["_E1"]))
        B_iv = (min(r["_S2"], r["_E2"]), max(r["_S2"], r["_E2"]))
        A_over = interval_any_overlap(A_iv, (piv for (piv, _miv) in placed_by_plasmid.get(A, [])))
        B_over = interval_any_overlap(B_iv, (piv for (piv, _miv) in placed_by_plasmid.get(B, [])))
        if not A_over and not B_over:
            absent_records.append((A, B, A_iv[0], A_iv[1], B_iv[0], B_iv[1], r["_PIDY"], r["_ALEN"]))

    # Sanity check across non-reference pairs
    pair_index: Dict[frozenset, List[Dict[str,Any]]] = defaultdict(list)
    for r in non_ref_rows:
        pair_index[frozenset((r["ref_file"], r["query_file"]))].append(r)

    sanity_missing: List[Tuple[str,str,int,int,int]] = []
    members_no_ref = [m for m in members if m != msa_reference]
    for i in range(len(members_no_ref)):
        A = members_no_ref[i]
        A_mivs = [miv for (_piv, miv) in placed_by_plasmid.get(A, [])]
        if not A_mivs:
            continue
        for j in range(i+1, len(members_no_ref)):
            B = members_no_ref[j]
            B_mivs = [miv for (_piv, miv) in placed_by_plasmid.get(B, [])]
            if not B_mivs:
                continue
            for a_m in A_mivs:
                for b_m in B_mivs:
                    ovl = interval_overlap(a_m, b_m)
                    if ovl >= sanity_min_overlap:
                        rows_ab = pair_index.get(frozenset((A, B)), [])
                        found = False
                        for r in rows_ab:
                            A_iv = (min(r["_S1"], r["_E1"]), max(r["_S1"], r["_E1"])) if r["ref_file"] == A else (min(r["_S2"], r["_E2"]), max(r["_S2"], r["_E2"]))
                            B_iv = (min(r["_S1"], r["_E1"]), max(r["_S1"], r["_E1"])) if r["ref_file"] == B else (min(r["_S2"], r["_E2"]), max(r["_S2"], r["_E2"]))
                            A_has = interval_any_overlap(A_iv, (piv for (piv, _miv) in placed_by_plasmid.get(A, [])))
                            B_has = interval_any_overlap(B_iv, (piv for (piv, _miv) in placed_by_plasmid.get(B, [])))
                            if A_has and B_has:
                                found = True
                                break
                        if not found:
                            sanity_missing.append((A, B, max(a_m[0], b_m[0]), min(a_m[1], b_m[1]), ovl))

    # --------------------------- outputs ---------------------------

    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    # MSA FASTA
    msa_fa = out_prefix.with_suffix(".aln.fasta")
    with msa_fa.open("w") as out:
        for name in members:
            out.write(f">{name}\n")
            row = "".join(msa[name])
            for i in range(0, len(row), 80):
                out.write(row[i:i+80] + "\n")
    log(f"Wrote MSA: {msa_fa}")

    # Conflicts
    conflicts_tsv = out_prefix.parent / (out_prefix.name + "_conflicts.tsv")
    with conflicts_tsv.open("w") as out:
        out.write("plasmid\tmsa_pos\texisting\tnew\tref_file\tquery_file\tS1\tE1\tS2\tE2\tpct_idy\trole\n")
        for rec in conflicts:
            out.write("\t".join(map(str, rec)) + "\n")
    log(f"Conflicts: {len(conflicts)} rows -> {conflicts_tsv}")

    # Absent-in-reference
    absent_tsv = out_prefix.parent / (out_prefix.name + "_absent_in_ref.tsv")
    with absent_tsv.open("w") as out:
        out.write("A\tB\tA_start\tA_end\tB_start\tB_end\tpct_idy\talen\n")
        for rec in absent_records:
            out.write("\t".join(map(str, rec)) + "\n")
    log(f"Absent-in-ref blocks: {len(absent_records)} -> {absent_tsv}")

    # Sanity missing
    sanity_tsv = out_prefix.parent / (out_prefix.name + "_sanity_missing.tsv")
    with sanity_tsv.open("w") as out:
        out.write("A\tB\tmsa_start\tmsa_end\toverlap_len\n")
        for rec in sanity_missing:
            out.write("\t".join(map(str, rec)) + "\n")
    log(f"Sanity missing links: {len(sanity_missing)} -> {sanity_tsv}")

    # Placed blocks index
    placed_tsv = out_prefix.parent / (out_prefix.name + "_placed_blocks.tsv")
    with placed_tsv.open("w") as out:
        out.write(
            "plasmid\trole\tmode\tstrand\t"
            "ref_start\tref_end\tref_span_len\t"
            "plasmid_start\tplasmid_end\tplasmid_span_len\t"
            "pct_idy\tref_file\tquery_file\tS1\tE1\tS2\tE2\n"
        )
        for rec in placed_blocks:
            out.write("\t".join(map(str, rec)) + "\n")
    log(f"Placed blocks: {len(placed_blocks)} -> {placed_tsv}")

    # Counters summary
    if gap_mode == "edlib":
        if not _HAVE_EDLIB:
            warn("gap_mode=edlib requested but edlib not installed; gapped blocks were skipped.")
        else:
            log(f"Edlib (unequal spans): {gapped_used}; equal-span edlib (slippage): {eqspan_edlib_used}; skipped gapped: {skipped_gapped}")
    elif gap_mode == "skip":
        log(f"Gapped blocks skipped: {skipped_gapped}")


def main():
    ap = argparse.ArgumentParser(description="Reference-anchored MSA from NUCmer blocks (gapped aware).")
    ap.add_argument("--coords", required=True, type=Path, help="Path to all_coords.tsv")
    ap.add_argument("--group-report", required=True, type=Path, help="Path to group_report.tsv")
    ap.add_argument("--group", required=True, help="Group name (e.g., group_40kb_4)")
    ap.add_argument("--fasta-dir", required=True, type=Path, help="Dir with member FASTAs")
    ap.add_argument("--msa-reference", required=True, help="Filename (must be in group)")
    ap.add_argument("--out-prefix", required=True, type=Path, help="Output prefix (path)")
    ap.add_argument("--min-len", default=1, type=int, help="Min block length (default: 1)")
    ap.add_argument("--min-id", default=0.0, type=float, help="Min %% identity (default: 0.0)")
    ap.add_argument("--sanity-min-overlap", default=200, type=int,
                    help="Min MSA-overlap (bp) to require A<->B link (default: 200)")
    ap.add_argument("--gap-mode", choices=["edlib","skip","die"], default="edlib",
                    help="When ref/query spans differ: use edlib, skip block, or die (default: edlib)")
    args = ap.parse_args()

    members = load_group_members(args.group_report, args.group)
    rows = read_coords(args.coords)
    rows = filter_rows_for_group(rows, set(members), args.min_len, args.min_id)

    log(f"Using {len(rows)} blocks within '{args.group}' (min_len={args.min_len}, min_id={args.min_id})")

    build_msa(rows, members, args.fasta_dir, args.msa_reference,
              args.out_prefix, args.sanity_min_overlap, args.gap_mode)

if __name__ == "__main__":
    main()

