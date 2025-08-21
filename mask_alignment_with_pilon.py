#!/usr/bin/env python3
"""
Mask an MSA using Pilon evidence (Amb / LowCov) and enrich conflicts with VCF support.

Changes in this version:
- LowCov rule now masks when GT is reference and DP < 3 (not only 0). Reason logged as LowCov_ref_DPX.
- QUAL/AF/DP-based polish/keep/mask retained (tunable via flags).
- If an MSA column maps to multiple blocks for a sample (overlap), evaluate all candidate mappings;
  if decisions differ, choose the decision from the record with higher DP (ties → higher QUAL → first).

Inputs:
  - --alignment            MSA fasta produced by build_block_msa (headers = plasmid filenames)
  - --placed-blocks        *_placed_blocks.tsv (from build_block_msa)
  - --conflicts            *_conflicts.tsv (from build_block_msa)
  - --pilon-dir            Directory with per-sample Pilon outputs:
                             <SAMPLE>_pilon.vcf.gz (+.tbi) and <SAMPLE>_pilon.changes
  - --fasta-dir            Directory with the per-plasmid FASTAs used for the MSA (to get lengths & slices)
  - --msa-reference        The MSA reference filename (same string used earlier)
  - --out-prefix           Output prefix

Notes:
  - We rebuild per-block edlib maps (only for blocks with mode='edlib') so ref columns inside those blocks map
    correctly to sample coordinates even with internal indels/slippage.
  - Mapping to VCF uses Pilon .changes to convert polished (post) → pre coordinates; VCF is on *pre* coords.
"""

import sys, os, re, csv, argparse, bisect
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from collections import defaultdict, OrderedDict

# ---------------- small helpers ----------------

BASES = set("ACGT")

def is_real_base(ch: str) -> bool:
    return ch in BASES

def _ctx(seq: str, pos1: int, w: int = 10) -> str:
    a = max(1, pos1 - w); b = min(len(seq), pos1 + w)
    left = "…" if a > 1 else ""
    right = "…" if b < len(seq) else ""
    return f"{left}{seq[a-1:b]}{right}"

_rc_map = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn",
                        "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")
def revcomp(s: str) -> str:
    return s.translate(_rc_map)[::-1]

def comp_base(b: str) -> str:
    return revcomp(b) if len(b) == 1 else b

# --- debug helper: pretty print a VariantRecord as a VCF line
def vcf_record_as_line(rec) -> str:
    alts = ",".join(rec.alts) if rec.alts else "."
    filt = ";".join(rec.filter.keys()) if rec.filter else "PASS"
    info_items = []
    for k, v in rec.info.items():
        if v is True:
            info_items.append(str(k))
        else:
            if isinstance(v, tuple):
                v = ",".join(map(str, v))
            info_items.append(f"{k}={v}")
    info = ";".join(info_items) if info_items else "."
    fmt_keys = list(rec.format.keys()) if rec.format else []
    fmt = ":".join(fmt_keys) if fmt_keys else "."
    smp = rec.samples[0] if rec.samples else None
    smp_vals = ":".join(str(smp.get(k, ".")) for k in fmt_keys) if smp else "."
    return (
        f"{rec.contig}\t{rec.pos}\t{rec.id or '.'}\t"
        f"{rec.ref}\t{alts}\t{rec.qual or '.'}\t{filt}\t{info}\t{fmt}\t{smp_vals}"
    )

# Extract contig num/len from names like "E1110_hybrid2025_ctg_7_len_71859_blaKPC2.fasta"
_ctg_re = re.compile(r"_ctg_(\d+)_len_(\d+)")
def parse_ctg_num_len_from_fname(fname: str) -> Tuple[Optional[int], Optional[int]]:
    m = _ctg_re.search(fname)
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))

# ---------------- FASTA I/O ----------------

def read_fasta_dict(path: Path) -> "OrderedDict[str,str]":
    seqs = OrderedDict()
    name = None
    buf: List[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs

def read_fasta_one(path: Path) -> Tuple[str, str]:
    header = None
    chunks: List[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if header is None:
                    header = line[1:].strip()
                else:
                    break
            else:
                chunks.append(line.strip())
    if header is None:
        raise ValueError(f"No header in {path}")
    return header, "".join(chunks).upper()

# ---------------- blocks (placed_blocks.tsv) ----------------

_HAVE_EDLIB = False
try:
    import edlib
    _HAVE_EDLIB = True
except Exception:
    _HAVE_EDLIB = False

class Block:
    __slots__ = ("plasmid","role","mode","strand","R1","R2","P1","P2","len_ref","len_qry","pidy",
                 "_ref_seq","_qry_seq","_ref2qry_off")
    def __init__(self, row: Dict[str,str]):
        self.plasmid = row["plasmid"]
        self.role = row["role"]
        self.mode = row["mode"]  # 'ungapped' or 'edlib'
        self.strand = row["strand"]  # '+' or '-'
        self.R1 = int(row["ref_start"]); self.R2 = int(row["ref_end"])
        self.P1 = int(row["plasmid_start"]); self.P2 = int(row["plasmid_end"])
        self.len_ref = int(row["ref_span_len"]); self.len_qry = int(row["plasmid_span_len"])
        self.pidy = float(row["pct_idy"])
        self._ref_seq = None
        self._qry_seq = None
        self._ref2qry_off: Optional[List[Optional[int]]] = None  # size len_ref

    def contains_ref(self, rpos: int) -> bool:
        return self.R1 <= rpos <= self.R2

    def load_seqs(self, ref_seq: str, sample_seq: str):
        if self._ref_seq is not None:
            return
        self._ref_seq = ref_seq[self.R1-1:self.R2]
        q = sample_seq[self.P1-1:self.P2]
        if self.strand == '-':
            q = revcomp(q)
        self._qry_seq = q

    def build_ref2qry_map(self):
        if self._ref2qry_off is not None:
            return
        self._ref2qry_off = [None] * self.len_ref
        if self.mode == "ungapped":
            for i in range(self.len_ref):
                self._ref2qry_off[i] = i if i < self.len_qry else None
            return
        if not _HAVE_EDLIB:
            raise RuntimeError("edlib not installed but block needs edlib mapping")
        res = edlib.align(self._qry_seq, self._ref_seq, mode="NW", task="path")
        nice = edlib.getNiceAlignment(res, self._qry_seq, self._ref_seq)
        q_aln = nice["query_aligned"]; r_aln = nice["target_aligned"]
        qi = -1; ri = -1
        for qa, ra in zip(q_aln, r_aln):
            if qa != '-': qi += 1
            if ra != '-': ri += 1
            if qa != '-' and ra != '-':
                if 0 <= ri < self.len_ref:
                    self._ref2qry_off[ri] = qi

    def refpos_to_sample_post(self, rpos: int) -> Optional[int]:
        """Map ref 1-based position to sample post-polish 1-based coordinate (contig-local).
           Returns None if the position is a deletion (query gap) inside edlib mapping."""
        if not self.contains_ref(rpos):
            return None
        off = rpos - self.R1
        self.build_ref2qry_map()
        qo = self._ref2qry_off[off]
        if qo is None:
            return None
        if self.strand == '+':
            return self.P1 + qo
        else:
            return self.P2 - qo

def load_blocks(placed_path: Path) -> Dict[str, List[Block]]:
    out: Dict[str, List[Block]] = defaultdict(list)
    with placed_path.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        need = {"plasmid","mode","strand","ref_start","ref_end","ref_span_len",
                "plasmid_start","plasmid_end","plasmid_span_len","pct_idy"}
        if not need.issubset(set(r.fieldnames or [])):
            raise SystemExit(f"[ERROR] {placed_path} missing columns {sorted(need)}")
        for row in r:
            b = Block(row)
            out[b.plasmid].append(b)
    for k in out:
        out[k].sort(key=lambda b: b.R1)
    return out

# ---------------- Pilon .changes post->pre ----------------

_change_re = re.compile(
    r"^(?P<pre_ctg>[^:]+):(?P<pre_span>\d+(?:-\d+)?)\s+"
    r"(?P<post_ctg>[^:]+):(?P<post_span>\d+(?:-\d+)?)\s+"
    r"(?P<pre_seq>\S+)\s+(?P<post_seq>\S+)$"
)

def _parse_span(s: str) -> tuple[int, int]:
    if "-" in s:
        a, b = s.split("-", 1)
        return int(a), int(b)
    x = int(s)
    return x, x

class PostToPreMap:
    """Piecewise-constant shift: pre = post + cumulative_shift(post)."""
    def __init__(self):
        self.cutoffs: List[int] = []   # post cutoffs where shift changes (start of 'after' region)
        self.shifts: List[int] = []    # cumulative shift to apply at >= cutoff

    def add_event(self, post_start: int, pre_len: int, post_len: int):
        cutoff = post_start + post_len
        delta = pre_len - post_len
        if self.cutoffs and cutoff < self.cutoffs[-1]:
            cutoff = self.cutoffs[-1]
        cum = (self.shifts[-1] if self.shifts else 0) + delta
        if self.cutoffs and cutoff == self.cutoffs[-1]:
            self.shifts[-1] = cum
        else:
            self.cutoffs.append(cutoff)
            self.shifts.append(cum)

    def finalize(self):
        pass

    def post_to_pre(self, pos_post: int) -> int:
        i = bisect.bisect_right(self.cutoffs, pos_post) - 1
        shift = self.shifts[i] if i >= 0 else 0
        return pos_post + shift

def build_post2pre_from_changes(changes_path: Path, target_post_ctg: str) -> PostToPreMap:
    """
    Build a piecewise-constant map: pre = post + shift(post).
    For each change line scoped to the *post* contig (RHS), we compute:
        delta = pre_len - post_len
    and apply this delta to all positions at/after the end of the POST edit
    (cutoff = post_start + post_len).
    """
    p2p = PostToPreMap()
    with changes_path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = _change_re.match(line)
            if not m:
                continue

            post_ctg = m.group("post_ctg")
            if post_ctg != target_post_ctg:
                continue

            pre_span = m.group("pre_span")
            post_span = m.group("post_span")
            pre_seq  = m.group("pre_seq")
            post_seq = m.group("post_seq")

            pre_s, pre_e   = _parse_span(pre_span)
            post_s, post_e = _parse_span(post_span)

            # Lengths: prefer spans; handle '.' (deletion/insertion)
            pre_len  = 0 if pre_seq  == "." else (pre_e  - pre_s  + 1)
            post_len = 0 if post_seq == "." else (post_e - post_s + 1)

            if pre_len != post_len:
                p2p.add_event(post_start=post_s, pre_len=pre_len, post_len=post_len)

    p2p.finalize()
    return p2p


# ---------------- VCF helpers ----------------

try:
    import pysam
except Exception:
    print("[ERROR] This script requires pysam (pip install pysam).", file=sys.stderr)
    raise

def vcf_open_and_pick_contig(vcf_path: Path,
                             contig_num: Optional[int],
                             contig_len_hint: Optional[int]) -> Tuple[pysam.VariantFile, str]:
    vf = pysam.VariantFile(vcf_path)
    ids = list(vf.header.contigs.keys())

    chosen = None
    if contig_num is not None:
        pref = [cid for cid in ids if cid.startswith(f"{contig_num}_")]
        if pref:
            if contig_len_hint is not None:
                def _hlen(cid):
                    return getattr(vf.header.contigs[cid], "length", None)
                pref_with_len = [(cid, _hlen(cid)) for cid in pref if _hlen(cid) is not None]
                if pref_with_len:
                    chosen = min(pref_with_len, key=lambda x: abs(x[1] - contig_len_hint))[0]
                else:
                    chosen = pref[0]
            else:
                chosen = pref[0]

    if chosen is None and contig_len_hint is not None:
        by_len_tag = [cid for cid in ids if f"length={contig_len_hint}" in cid]
        if by_len_tag:
            chosen = by_len_tag[0]

    if chosen is None:
        def _hlen(cid):
            return getattr(vf.header.contigs[cid], "length", None) or 0
        if contig_len_hint is not None:
            chosen = min(ids, key=lambda cid: abs(_hlen(cid) - contig_len_hint))
            print(f"[WARN] Could not find contig by number/length in {vcf_path.name}; using closest by length: {chosen}",
                  file=sys.stderr)
        else:
            chosen = ids[0]
            print(f"[WARN] Could not infer contig for {vcf_path.name}; using: {chosen}", file=sys.stderr)

    hdr_len = getattr(vf.header.contigs[chosen], "length", None)
    print(f"[INFO] Using VCF contig '{chosen}' (header length={hdr_len}) from {vcf_path.name}", file=sys.stderr)
    return vf, chosen

def vcf_fetch_site(vf: pysam.VariantFile, contig: str, pos1: int) -> Optional[Any]:
    for rec in vf.fetch(contig, pos1-1, pos1):
        if rec.pos == pos1:
            return rec
    return None

def vcf_basic_metrics(rec: Any) -> Tuple[str, Optional[int], Any, str, List[str], Optional[float], Optional[float]]:
    # FILTER
    filt = "PASS"
    if rec.filter and len(rec.filter.keys()) > 0:
        filt = ";".join(list(rec.filter.keys()))
    # DP
    dp = rec.info.get("DP", None)
    if not isinstance(dp, int) and dp is not None:
        try: dp = int(dp)
        except Exception: dp = None
    # GT
    smp = rec.samples[0] if rec.samples else {}
    gt = smp.get("GT", None)
    # REF / ALTS
    ref = rec.ref
    alts = list(rec.alts or [])
    # QUAL
    qual = float(rec.qual) if rec.qual is not None else None
    # AF: Pilon reports AF in INFO; can be tuple. For ALT='.' (no variant), treat AF=0.0
    af_raw = rec.info.get("AF", None)
    if af_raw is None:
        af = 0.0 if not alts else None
    else:
        af = float(af_raw[0] if isinstance(af_raw, (list, tuple)) else af_raw)
    return filt, dp, gt, ref, alts, qual, af

def gt_is_ref(gt: Any) -> bool:
    return gt in {(0,0), "0/0", "0|0", (0,), "0"}

def gt_is_hom_alt(gt: Any) -> bool:
    return gt in {(1,1), "1/1", "1|1", (1,), "1"}

def gt_is_het(gt: Any) -> bool:
    return gt in {(0,1), (1,0), "0/1", "1/0", "0|1", "1|0"}

def orient_allele(allele: Optional[str], strand: str) -> Optional[str]:
    if allele is None:
        return None
    if len(allele) != 1:
        return None  # only act on SNPs
    return allele if strand == '+' else revcomp(allele)

# ---------------- conflicts index ----------------

def load_conflicts_index(conflicts_path: Path) -> Dict[Tuple[str,int], set]:
    idx: Dict[Tuple[str,int], set] = defaultdict(set)
    if not conflicts_path.exists():
        return idx
    with conflicts_path.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        need = {"plasmid","msa_pos","existing","new"}
        if not need.issubset(set(r.fieldnames or [])):
            return idx
        for row in r:
            s = row["plasmid"]
            try:
                pos = int(row["msa_pos"])
            except Exception:
                continue
            ex = row["existing"].strip().upper()
            nw = row["new"].strip().upper()
            if ex in BASES: idx[(s, pos)].add(ex)
            if nw in BASES: idx[(s, pos)].add(nw)
    return idx

# ---------------- sample file lookup ----------------

def guess_sample_code(plasmid_name: str) -> str:
    # Use leading token before first underscore, e.g. E1049_...
    return plasmid_name.split("_", 1)[0]

def find_pilon_files(pilon_dir: Path, sample_code: str) -> Tuple[Optional[Path], Optional[Path]]:
    # Anchor with underscore to avoid E39 matching E390, etc.
    vcf = None; chg = None
    vcfs = sorted(pilon_dir.glob(f"{sample_code}_*pilon.vcf.gz"))
    if vcfs:
        vcf = vcfs[0]
    chgs = sorted(pilon_dir.glob(f"{sample_code}_*pilon.changes"))
    if chgs:
        chg = chgs[0]
    return vcf, chg

# ---------------- main ----------------

def main():
    ap = argparse.ArgumentParser(description="Mask MSA using Pilon evidence and enrich conflicts.")
    ap.add_argument("--alignment", required=True, type=Path)
    ap.add_argument("--placed-blocks", required=True, type=Path)
    ap.add_argument("--conflicts", required=True, type=Path)
    ap.add_argument("--pilon-dir", required=True, type=Path)
    ap.add_argument("--fasta-dir", required=True, type=Path)
    ap.add_argument("--msa-reference", required=True)
    ap.add_argument("--out-prefix", required=True, type=Path)
    ap.add_argument("--gap-dp-thresh", type=int, default=5, help="If both flanks < this DP, replace '-' with 'N'")
    ap.add_argument("--polish-qual", type=float, default=30.0, help="Minimum QUAL to consider polishing/keeping")
    ap.add_argument("--polish-af-high", type=float, default=0.90, help="AF threshold to polish to ALT")
    ap.add_argument("--polish-af-low", type=float, default=0.10, help="AF threshold to keep REF")
    ap.add_argument("--polish-dp-min", type=int, default=2, help="Minimum DP to consider polishing/keeping")
    ap.add_argument("--uncertain", choices=["mask","keep"], default="mask",
                    help="What to do when neither REF nor ALT rule passes")
    ap.add_argument("--debug", action="store_true")
    args = ap.parse_args()

    # Load MSA
    msa_dict = read_fasta_dict(args.alignment)
    samples = list(msa_dict.keys())
    if args.msa_reference not in samples:
        print(f"[ERROR] --msa-reference header not found in alignment: {args.msa_reference}", file=sys.stderr)
        sys.exit(2)
    ref = args.msa_reference
    ref_len = len(msa_dict[ref])
    print(f"[INFO] Alignment loaded: {len(samples)} sequences; length={ref_len}", file=sys.stderr)

    # Load placed blocks and build maps
    blocks_by_plasmid = load_blocks(args.placed_blocks)

    # Load per-sample raw polished sequences (the ones used in the MSA)
    sample_raw_seq: Dict[str, str] = {}
    plasmid_len_map: Dict[str, int] = {}
    for s in samples:
        fa = args.fasta_dir / s
        header, seq = read_fasta_one(fa)
        sample_raw_seq[s] = seq
        plasmid_len_map[s] = len(seq)

    # Pre-build edlib maps per block
    ref_seq_raw = sample_raw_seq[ref]
    for s, blist in blocks_by_plasmid.items():
        raw_seq = sample_raw_seq.get(s)
        if raw_seq is None:
            continue
        for b in blist:
            b.load_seqs(ref_seq_raw, raw_seq)
            b.build_ref2qry_map()

    # Prepare per-sample VCF + post→pre mapper (INCLUDE REF)
    vcf_by_sample: Dict[str, Tuple[pysam.VariantFile, str]] = {}
    post2pre_by_sample: Dict[str, PostToPreMap] = {}
    for s in samples:
        code = guess_sample_code(s)
        vcf_path, chg_path = find_pilon_files(args.pilon_dir, code)
        if not vcf_path or not vcf_path.exists():
            print(f"[WARN] No VCF for {s} (code {code}) in {args.pilon_dir}", file=sys.stderr)
            continue
        ctg_num, ctg_len_hint = parse_ctg_num_len_from_fname(s)
        if ctg_len_hint is None:
            ctg_len_hint = plasmid_len_map.get(s)
        vf, contig = vcf_open_and_pick_contig(vcf_path, ctg_num, ctg_len_hint)
        vcf_by_sample[s] = (vf, contig)

        if not chg_path or not chg_path.exists():
            print(f"[WARN] No .changes for {s} (code {code}); assuming no indel shifts.", file=sys.stderr)
            post2pre_by_sample[s] = PostToPreMap()
        else:
            # RHS contig may be "<contig>_pilon" or just "<contig>"
            candidates = [f"{contig}_pilon", contig]
            best_p2p = None
            best_n = -1
            best_label = None
            for cand in candidates:
                p2p = build_post2pre_from_changes(chg_path, cand)
                n = len(p2p.cutoffs)
                if n > best_n:
                    best_n = n
                    best_p2p = p2p
                    best_label = cand
            if best_p2p is None:
                best_p2p = PostToPreMap()
            print(f"[INFO] Using .changes contig '{best_label}' with {best_n} shift events for {s}", file=sys.stderr)
            post2pre_by_sample[s] = best_p2p

    # Conflicts index (for early conflict masking)
    conflicts_idx = load_conflicts_index(args.conflicts)

    # Identify variant columns (ignore '-' and 'N')
    var_cols: List[int] = []
    col_bases_tmp = set()
    for i in range(ref_len):
        col_bases_tmp.clear()
        for s in samples:
            ch = msa_dict[s][i]
            if ch in "-N":
                continue
            if is_real_base(ch):
                col_bases_tmp.add(ch)
        if len(col_bases_tmp) >= 2:
            var_cols.append(i+1)  # 1-based
    print(f"[INFO] Variant columns: {len(var_cols)}", file=sys.stderr)

    # Mutable rows for masking
    masked_rows: Dict[str, List[str]] = {s: list(msa_dict[s]) for s in samples}

    # Helper: all blocks covering rpos (for overlap cases)
    def find_blocks(sample: str, rpos: int) -> List[Block]:
        return [b for b in blocks_by_plasmid.get(sample, []) if b.contains_ref(rpos)]

    # Map ref column to *all* candidate (block, post) mappings
    def map_ref_to_sample_post_all(sample: str, rpos: int) -> List[Tuple[Optional[Block], int, str]]:
        # (block_or_None, post_pos, strand)
        if sample == ref:
            return [(None, rpos, '+')]
        blks = find_blocks(sample, rpos)
        if not blks:
            return []
        out = []
        for b in blks:
            post_pos = b.refpos_to_sample_post(rpos)
            if post_pos is not None:
                out.append((b, post_pos, b.strand))
        return out

    # Log masked positions
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    mask_log_path = args.out_prefix.parent / (args.out_prefix.name + "_mask_log.tsv")
    mask_log = open(mask_log_path, "w")
    mask_log.write("plasmid\tmsa_pos\taction\treason\tvcf_contig\tpre_pos\tpost_pos\tFILTER\tDP\tGT\n")

    # Decision helper for a single mapping (one block candidate)
    def decide_from_rec(ch: str, strand: str, vf_contig: str, pre_pos: int, post_pos: int,
                        rec: Optional[Any],
                        polish_qual: float, af_hi: float, af_lo: float, dp_min: int,
                        uncertain_policy: str) -> Tuple[str, str, Optional[int], Optional[float], Optional[float], Any, str]:
        """
        Returns: (action, reason, dp, qual, af, gt, vline)
                 action ∈ {'mask','polish','keep','none'}
        """
        if rec is None:
            return ("none", "no_vcf", None, None, None, None, "<no VCF record at site>")

        filt, dp, gt, vref, alts, qual, af = vcf_basic_metrics(rec)
        vline = vcf_record_as_line(rec)

        # --- Mandatory masks
        if "Amb" in filt:
            return ("mask", "Amb", dp, qual, af, gt, vline)

        # LowCov (reference) now masks when DP < 3
        if ("LowCov" in filt) and gt_is_ref(gt) and (dp is None or dp < 3):
            return ("mask", "LowCov_ref_DPX", dp, qual, af, gt, vline)

        # --- NEW: No-variant rows (ALT missing) → keep if decent evidence
        # pysam gives rec.alts == None when ALT is '.', so "not alts" detects this.
        if not alts:
            if (dp is not None and dp >= dp_min) and (qual is not None and qual >= polish_qual):
                return ("keep", "no_variant", dp, qual, af, gt, vline)
            # Insufficient evidence → uncertain policy
            return (("mask", "uncertain_variant") if uncertain_policy == "mask"
                    else ("keep", "uncertain_keep")) + (dp, qual, af, gt, vline)

        # --- SNP rules (only when REF and first ALT are single bases)
        # Orient alleles for '-' strand
        vref_o = orient_allele(vref, strand)
        alt0   = (alts[0] if alts else None)
        alt0_o = orient_allele(alt0, strand)

        if vref_o and alt0_o and len(vref) == 1 and len(alt0 or "") == 1:
            # Polish to ALT
            if gt_is_hom_alt(gt) and (af is not None and af >= af_hi) and \
               (qual is not None and qual >= polish_qual) and (dp is not None and dp >= dp_min):
                if ch != alt0_o:
                    return ("polish", "ALT", dp, qual, af, gt, vline)
                else:
                    return ("keep", "ALT_already", dp, qual, af, gt, vline)

            # Keep REF
            if gt_is_ref(gt) and (af is not None and af <= af_lo) and \
               (qual is not None and qual >= polish_qual) and (dp is not None and dp >= dp_min):
                return ("keep", "REF", dp, qual, af, gt, vline)

        # --- Everything else → uncertain policy
        return (("mask", "uncertain_variant") if uncertain_policy == "mask"
                else ("keep", "uncertain_keep")) + (dp, qual, af, gt, vline)


    # Pick best decision among multiple mapping candidates
    def pick_best_decision(decisions: List[Tuple[str, str, Optional[int], Optional[float], Optional[float], Any, str]]
                           ) -> Tuple[int, Tuple[str,str,Optional[int],Optional[float],Optional[float],Any,str]]:
        """
        If decisions disagree, choose the one with higher DP (ties → higher QUAL → index 0).
        Returns (chosen_index, chosen_decision_tuple)
        """
        if not decisions:
            return (-1, ("none","no_mapping",None,None,None,None,""))
        # If all actions+reasons identical, keep the first
        first = (decisions[0][0], decisions[0][1])
        same = all((d[0], d[1]) == first for d in decisions)
        if same:
            return (0, decisions[0])
        # Else choose by DP, then QUAL
        best_i = 0
        best_dp = decisions[0][2] if decisions[0][2] is not None else -1
        best_qual = decisions[0][3] if decisions[0][3] is not None else -1.0
        for i in range(1, len(decisions)):
            dp = decisions[i][2] if decisions[i][2] is not None else -1
            qual = decisions[i][3] if decisions[i][3] is not None else -1.0
            if dp > best_dp or (dp == best_dp and qual > best_qual):
                best_i = i
                best_dp = dp
                best_qual = qual
        return (best_i, decisions[best_i])

    # ---------- mask variant bases (conflict-first) ----------
    amb_sites = []; lowcov_sites = []; uncertain_sites = []; conflict_mask_sites = []; polished_sites = []

    for rpos in var_cols:
        for s in samples:
            if s not in vcf_by_sample:
                continue
            ch = masked_rows[s][rpos-1]

            # Conflict-first handling
            ckey = (s, rpos)
            if ckey in conflicts_idx:
                allowed = sorted(conflicts_idx[ckey])
                if ch in allowed:
                    masked_rows[s][rpos-1] = 'N'
                    mask_log.write(f"{s}\t{rpos}\tmask\tconflict_allowed\t.\t.\t.\t.\t.\t.\n")
                    if args.debug:
                        print(f"[WARN] Masked {s} at MSA {rpos} due to conflict; base {ch} is in allowed={allowed}",
                              file=sys.stderr)
                    conflict_mask_sites.append((s, rpos, ch, "conflict_allowed"))
                    continue
                else:
                    masked_rows[s][rpos-1] = 'N'
                    mask_log.write(f"{s}\t{rpos}\tmask\tconflict_inconsistent\t.\t.\t.\t.\t.\t.\n")
                    if args.debug:
                        print(f"[WARN] Masked {s} at MSA {rpos} due to conflict inconsistency; base {ch} not in allowed={allowed}",
                              file=sys.stderr)
                    conflict_mask_sites.append((s, rpos, ch, "conflict_inconsistent"))
                    continue

            # All candidate mappings (overlapping blocks possible)
            candidates = map_ref_to_sample_post_all(s, rpos)
            if not candidates:
                continue

            vf, contig = vcf_by_sample[s]
            p2p = post2pre_by_sample[s]

            # Evaluate each mapping
            decisions = []
            meta = []  # to keep (post_pos, pre_pos, strand)
            for b, post_pos, strand in candidates:
                pre_pos = p2p.post_to_pre(post_pos)
                rec = vcf_fetch_site(vf, contig, pre_pos)
                d = decide_from_rec(
                    ch=ch, strand=strand, vf_contig=contig, pre_pos=pre_pos, post_pos=post_pos, rec=rec,
                    polish_qual=args.polish_qual, af_hi=args.polish_af_high, af_lo=args.polish_af_low,
                    dp_min=args.polish_dp_min, uncertain_policy=args.uncertain
                )
                decisions.append(d)
                meta.append((post_pos, pre_pos, strand, rec))

            chosen_idx, chosen = pick_best_decision(decisions)
            if chosen_idx < 0:
                continue
            action, reason, dp, qual, af, gt, vline = chosen
            post_pos, pre_pos, strand, rec = meta[chosen_idx]

            # Apply action
            if action == "mask" and ch != 'N':
                masked_rows[s][rpos-1] = 'N'
                mask_log.write(f"{s}\t{rpos}\tmask\t{reason}\t{contig}\t{pre_pos}\t{post_pos}\t"
                               f"{( ';'.join(rec.filter.keys()) if rec and rec.filter else 'PASS' )}\t{dp}\t{gt}\n")
                if args.debug:
                    print(f"[DEBUG] Masked {s} at MSA {rpos} due to {reason}; "
                          f"VCF {contig}:{pre_pos} -> {vline}", file=sys.stderr)
                if reason == "Amb":
                    amb_sites.append((s, rpos))
                elif reason == "LowCov_ref_DPX":
                    lowcov_sites.append((s, rpos))
                else:
                    uncertain_sites.append((s, rpos))

            elif action == "polish":
                # Polish to oriented ALT
                if rec is not None:
                    _, _, _, vref, alts, _, _ = vcf_basic_metrics(rec)
                    alt0 = (alts[0] if alts else None)
                    alt0_o = orient_allele(alt0, strand)
                    if alt0_o and ch != alt0_o:
                        masked_rows[s][rpos-1] = alt0_o
                        mask_log.write(f"{s}\t{rpos}\tpolish\tALT\t{contig}\t{pre_pos}\t{post_pos}\t"
                                       f"{( ';'.join(rec.filter.keys()) if rec.filter else 'PASS' )}\t{dp}\t{gt}\n")
                        if args.debug:
                            print(f"[INFO] Polished {s} msa={rpos} to ALT({alt0}->{alt0_o}) "
                                  f"[AF={af},QUAL={qual},DP={dp}] | {vline}", file=sys.stderr)
                        polished_sites.append((s, rpos, "ALT"))

            # else keep / none → do nothing

    # ---------- gap → N (low-coverage flanks) ----------
    for s in samples:
        if s not in vcf_by_sample:
            continue
        vf, contig = vcf_by_sample[s]
        p2p = post2pre_by_sample[s]
        row = masked_rows[s]
        i = 0
        while i < ref_len:
            if row[i] != '-':
                i += 1; continue
            j = i
            while j < ref_len and row[j] == '-':
                j += 1
            left_r = i
            while left_r > 0 and row[left_r] == '-':
                left_r -= 1
            right_r = j
            while right_r < ref_len and row[right_r] == '-':
                right_r += 1

            left_ok = False; right_ok = False
            left_dp = None; right_dp = None

            # Left flank
            if 0 <= left_r < ref_len and row[left_r] != '-':
                # Use primary mapping only for gaps (keep simple)
                posts = map_ref_to_sample_post_all(s, left_r+1)
                if posts:
                    _, post_pos, _ = posts[0]
                    pre_pos = p2p.post_to_pre(post_pos)
                    rec = vcf_fetch_site(vf, contig, pre_pos)
                    if rec is not None:
                        _, dp, _, _, _, _, _ = vcf_basic_metrics(rec)
                        left_dp = dp
                        left_ok = (dp is not None and dp >= args.gap_dp_thresh)

            # Right flank
            if 0 <= right_r < ref_len and row[right_r] != '-':
                posts = map_ref_to_sample_post_all(s, right_r+1)
                if posts:
                    _, post_pos, _ = posts[0]
                    pre_pos = p2p.post_to_pre(post_pos)
                    rec = vcf_fetch_site(vf, contig, pre_pos)
                    if rec is not None:
                        _, dp, _, _, _, _, _ = vcf_basic_metrics(rec)
                        right_dp = dp
                        right_ok = (dp is not None and dp >= args.gap_dp_thresh)

            if not left_ok and not right_ok:
                for k in range(i, j):
                    if row[k] == '-':
                        row[k] = 'N'
                        mask_log.write(f"{s}\t{k+1}\tgap->N\tlow_flank_DP<{args.gap_dp_thresh}\t{contig}\t.\t.\t.\t{left_dp}|{right_dp}\t.\n")
                        if args.debug:
                            print(f"[DEBUG] gap->N {s} MSA {k+1} flanks DP: L={left_dp} R={right_dp}", file=sys.stderr)

            i = j

    # ---------- enrich conflicts with VCF ----------
    conflicts_with_vcf = args.out_prefix.parent / (args.out_prefix.name + "_conflicts_with_vcf.tsv")
    with open(args.conflicts) as fh, open(conflicts_with_vcf, "w") as out:
        r = csv.DictReader(fh, delimiter="\t")
        cols = r.fieldnames or []
        out.write("\t".join(cols + ["vcf_contig","vcf_pre_pos","FILTER","DP","GT","vcf_ref","vcf_alt","supports_existing","supports_new"]) + "\n")
        for row in r:
            s = row["plasmid"]
            rpos = int(row["msa_pos"])
            existing = row["existing"]; new = row["new"]
            vcf_contig = "."; v_pre = "."; filt="."; dp="."; gt="."; vref="."; valt="."
            sup_e = "."; sup_n = "."
            if s in vcf_by_sample:
                vf, contig = vcf_by_sample[s]
                posts = map_ref_to_sample_post_all(s, rpos)
                # Just use the first mapping for conflict enrichment (traceable and simple)
                if posts:
                    _, post_pos, strand = posts[0]
                    pre_pos = post2pre_by_sample[s].post_to_pre(post_pos)
                    rec = vcf_fetch_site(vf, contig, pre_pos)
                    if rec:
                        vcf_contig = contig; v_pre = str(pre_pos)
                        filt0, dpv, gtv, vref0, alts0, _, _ = vcf_basic_metrics(rec)
                        filt = filt0
                        dp = str(dpv) if dpv is not None else "."
                        gt = str(gtv)
                        valt = ",".join(alts0) if alts0 else "."
                        # Decide strand to compare alleles in MSA orientation
                        vref_cmp = vref0 if strand == '+' else comp_base(vref0) if vref0 and len(vref0)==1 else vref0
                        alt0_cmp = (alts0[0] if alts0 else None)
                        if alt0_cmp and len(alt0_cmp)==1 and strand == '-':
                            alt0_cmp = comp_base(alt0_cmp)
                        # support:
                        if alt0_cmp and alt0_cmp in BASES and not gt_is_ref(gtv):
                            sup_e = "1" if existing == alt0_cmp else "0"
                            sup_n = "1" if new == alt0_cmp else "0"
                        else:
                            sup_e = "1" if existing == vref_cmp else "0"
                            sup_n = "1" if new == vref_cmp else "0"
                        vref = vref0
            out.write("\t".join([row[c] for c in cols] + [vcf_contig, v_pre, filt, dp, gt, vref, valt, sup_e, sup_n]) + "\n")

    # ---------- write masked alignment ----------
    masked_path = args.out_prefix.with_suffix(".masked.fasta")
    with masked_path.open("w") as out:
        for s in samples:
            out.write(f">{s}\n")
            row = "".join(masked_rows[s])
            for i in range(0, len(row), 80):
                out.write(row[i:i+80] + "\n")

    mask_log.close()

    # Debug summaries
    if args.debug:
        print(f"[DEBUG] Masked Amb sites: {len(amb_sites)}", file=sys.stderr)
        print(f"[DEBUG] Masked LowCov_ref_DPX sites: {len(lowcov_sites)}", file=sys.stderr)
        print(f"[DEBUG] Masked uncertain/other sites: {len(uncertain_sites)}", file=sys.stderr)
        print(f"[DEBUG] Masked due to conflict: {len(conflict_mask_sites)}", file=sys.stderr)
        print(f"[DEBUG] Polished sites: {len(polished_sites)}", file=sys.stderr)

    print(f"[INFO] Wrote masked alignment: {masked_path}", file=sys.stderr)
    print(f"[INFO] Mask log: {mask_log_path}", file=sys.stderr)
    print(f"[INFO] Conflicts+VCF: {conflicts_with_vcf}", file=sys.stderr)

if __name__ == "__main__":
    main()

