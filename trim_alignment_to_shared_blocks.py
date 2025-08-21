#!/usr/bin/env python3
"""
Trim an alignment to shared reference blocks, with optional N-spacers.

Inputs
------
- --alignment     : multi-FASTA MSA (headers must match the filenames used in TSV)
- --blocks-tsv    : TSV with columns: ref_file, group, block_idx, start, end, length
- --group         : group name to select in TSV (e.g., group_40kb_47)
- --ref           : reference header (must match 'ref_file' in TSV and a header in the MSA)
- --out           : output FASTA path
- --spacer-len    : number of 'N' characters to place between blocks (default 0)
- --wrap          : FASTA line wrap width (default 80, use 0 for no wrapping)

What it does
------------
1) Reads the MSA and finds the reference row.
2) Builds a map from reference *ungapped* positions (1-based) to alignment columns.
   (The reference counter increments on any non-'-' character in the ref row.)
3) Loads shared blocks from the TSV for (--group, --ref), sorts by start, and
   for each block collects all alignment columns whose ref-position falls in [start, end].
   (Columns where the ref row is '-' are excluded — i.e., insertions relative to the reference
    are *not* included.)
4) Extracts those columns from every sequence and concatenates them, inserting a spacer
   of 'N' * --spacer-len between blocks.
5) Writes the trimmed alignment.

Notes
-----
- If a block contains reference gaps inside the alignment (insertions in other samples),
  those insertion columns are skipped (since they do not map to a reference position).
- If no blocks are found for the given filters, an empty sequence will be produced for all rows.
"""

import sys
import csv
import argparse
from pathlib import Path
from collections import OrderedDict
from typing import List, Tuple, Dict, Optional

def read_fasta_dict(path: Path) -> "OrderedDict[str, str]":
    seqs = OrderedDict()
    name = None
    chunks = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs

def load_blocks_tsv(path: Path, group: str, ref_file: str) -> List[Tuple[int, int, int]]:
    """
    Returns list of (block_idx, start, end) for rows matching group & ref_file.
    """
    out: List[Tuple[int, int, int]] = []
    with path.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        need = {"ref_file", "group", "block_idx", "start", "end"}
        if not need.issubset(set(r.fieldnames or [])):
            raise SystemExit(f"[ERROR] {path} missing required columns: {sorted(need)}")
        for row in r:
            if row["group"].strip() != group:
                continue
            if row["ref_file"].strip() != ref_file:
                continue
            try:
                bidx = int(row["block_idx"])
                start = int(row["start"])
                end = int(row["end"])
            except Exception:
                continue
            if end < start:
                start, end = end, start
            out.append((bidx, start, end))
    # Sort by start (then block_idx)
    out.sort(key=lambda x: (x[1], x[0]))
    return out

def build_refpos_per_column(ref_row: str) -> List[Optional[int]]:
    """
    For each alignment column i (0-based), returns the 1-based reference position
    mapped to that column, or None if the reference has '-' at that column.
    Reference position counter increments for every non-'-' character in ref_row.
    """
    refpos_per_col: List[Optional[int]] = [None] * len(ref_row)
    ref_pos = 0
    for i, ch in enumerate(ref_row):
        if ch != "-":
            ref_pos += 1
            refpos_per_col[i] = ref_pos
        else:
            refpos_per_col[i] = None
    return refpos_per_col

def indices_for_block(refpos_per_col: List[Optional[int]], start: int, end: int) -> List[int]:
    """Return alignment column indices whose reference position is within [start, end]."""
    lo, hi = min(start, end), max(start, end)
    keep = []
    for i, rp in enumerate(refpos_per_col):
        if rp is None:
            continue
        if lo <= rp <= hi:
            keep.append(i)
    return keep

def wrap_fasta(seq: str, width: int) -> str:
    if width is None or width <= 0:
        return seq + ("\n" if not seq.endswith("\n") else "")
    out = []
    for i in range(0, len(seq), width):
        out.append(seq[i:i+width])
    return "\n".join(out) + ("\n" if out else "\n")

def main():
    ap = argparse.ArgumentParser(description="Trim an alignment to shared reference blocks with optional N-spacers.")
    ap.add_argument("--alignment", required=True, type=Path, help="Input MSA FASTA")
    ap.add_argument("--blocks-tsv", required=True, type=Path, help="Shared blocks TSV")
    ap.add_argument("--group", required=True, help="Group name in TSV (e.g. group_40kb_47)")
    ap.add_argument("--ref", required=True, help="Reference header (must match both MSA and TSV ref_file)")
    ap.add_argument("--out", required=True, type=Path, help="Output FASTA")
    ap.add_argument("--spacer-len", type=int, default=0, help="Spacer length (N's) between blocks")
    ap.add_argument("--wrap", type=int, default=80, help="FASTA wrap width (0 for no wrap)")
    args = ap.parse_args()

    # Read alignment
    msa = read_fasta_dict(args.alignment)
    if args.ref not in msa:
        raise SystemExit(f"[ERROR] Reference '{args.ref}' not found in alignment headers.")
    ref_row = msa[args.ref]
    aln_len = len(ref_row)
    if aln_len == 0:
        raise SystemExit("[ERROR] Alignment reference row is empty.")

    # Build per-column reference positions
    refpos_per_col = build_refpos_per_column(ref_row)

    # Load blocks
    blocks = load_blocks_tsv(args.blocks_tsv, args.group, args.ref)
    if not blocks:
        print(f"[WARN] No blocks found for group='{args.group}' ref='{args.ref}' in {args.blocks_tsv}", file=sys.stderr)

    # For each block, collect alignment columns to keep (skip columns where ref='-')
    block_col_indices: List[List[int]] = []
    for bidx, start, end in blocks:
        idxs = indices_for_block(refpos_per_col, start, end)
        if not idxs:
            print(f"[WARN] Block {bidx} [{start}-{end}] produced 0 columns (possible heavy ref gaps?)", file=sys.stderr)
        block_col_indices.append(idxs)

    # Prepare spacer
    spacer = "N" * max(0, args.spacer_len)

    # Build trimmed sequences
    trimmed: Dict[str, str] = {}
    for hdr, seq in msa.items():
        parts: List[str] = []
        for k, idxs in enumerate(block_col_indices):
            if not idxs:
                continue
            # Collect only the columns that map to reference positions
            part = "".join(seq[i] for i in idxs)
            parts.append(part)
            if args.spacer_len > 0 and k != len(block_col_indices) - 1:
                parts.append(spacer)
        trimmed_seq = "".join(parts)
        trimmed[hdr] = trimmed_seq

    # Write output FASTA
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as out:
        for hdr in trimmed:
            out.write(f">{hdr}\n")
            out.write(wrap_fasta(trimmed[hdr], args.wrap))

    # Report
    total_cols = sum(len(x) for x in block_col_indices)
    print(f"[INFO] Alignment length: {aln_len}", file=sys.stderr)
    print(f"[INFO] Blocks used: {len(blocks)} ; total kept ref-mapped columns: {total_cols}", file=sys.stderr)
    if args.spacer_len > 0:
        print(f"[INFO] Spacer of {args.spacer_len} N's inserted between blocks.", file=sys.stderr)

if __name__ == "__main__":
    main()

