#!/usr/bin/env python3
import argparse
import csv
import os
import re
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Iterable, Set, Optional

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

# =========================
# db_xref parsing (generic)
# =========================

_UNIREF_RX = re.compile(r'^UniRef:UniRef(100|90|50)_(.+)$')
_GO_DOUBLE_RX = re.compile(r'^GO:GO:(.+)$')  # handles GO:GO:0008150

def split_dbxref_to_typed_pairs(dbx: str):
    """
    Parse a raw db_xref string into (column_key, value).
    Examples:
      RefSeq:WP_013213988.1        -> ('db_xref:RefSeq', 'WP_013213988.1')
      UniParc:UPI0001B9A894        -> ('db_xref:UniParc', 'UPI0001B9A894')
      UniRef:UniRef100_A0A142EBB1  -> ('db_xref:UniRef100', 'A0A142EBB1')
      SO:0001217                   -> ('db_xref:SO', '0001217')
      GO:0008150 / GO:GO:0008150   -> ('db_xref:GO', '0008150')
      KEGG:K01814, COG:COG1629, EC:1.1.1.1, PFAM:PF00005, VFDB:..., IS:IS26, ...
    """
    dbx = dbx.strip()
    m = _UNIREF_RX.match(dbx)
    if m:
        depth, val = m.group(1), m.group(2)
        return (f'db_xref:UniRef{depth}', val)
    m = _GO_DOUBLE_RX.match(dbx)
    if m:
        return ('db_xref:GO', m.group(1))
    if ':' in dbx:
        ns, rest = dbx.split(':', 1)
        ns = ns.strip()
        rest = rest.strip()
        if ns and rest:
            return (f'db_xref:{ns}', rest)
    return None

# =========================
# GBFF scanning / helpers
# =========================

def feature_bounds_1based(feat: SeqFeature) -> Tuple[int, int]:
    """
    Return (start,end) as 1-based inclusive for any Biopython location,
    handling SimpleLocation and CompoundLocation.
    """
    loc = feat.location
    if isinstance(loc, CompoundLocation):
        starts = [int(part.start) + 1 for part in loc.parts]
        ends   = [int(part.end) for part in loc.parts]
        start, end = min(starts), max(ends)
    else:
        start = int(loc.start) + 1
        end   = int(loc.end)
    if start > end:
        start, end = end, start
    return start, end

def overlaps_1based(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    """1-based, inclusive intervals overlap?"""
    return not (a_end < b_start or b_end < a_start)

def fully_contained_1based(inner_start: int, inner_end: int, outer_start: int, outer_end: int) -> bool:
    """Is [inner_start, inner_end] fully inside [outer_start, outer_end]?"""
    return inner_start >= outer_start and inner_end <= outer_end

def discover_available_columns(gbff_dir: str, sample_limit: Optional[int]=None) -> Dict[str, object]:
    """
    Scan gbff files to discover qualifier keys, db_xref 'types', and feature types.
    Also tallies db_xref namespace counts.
    """
    qkeys: Set[str] = set()
    dbxref_types: Set[str] = set()
    ftypes: Set[str] = set()
    dbxref_counts: Counter = Counter()

    files = [f for f in os.listdir(gbff_dir) if f.lower().endswith(('.gbff', '.gbk', '.gb'))]
    if sample_limit:
        files = files[:sample_limit]

    for fname in files:
        path = os.path.join(gbff_dir, fname)
        try:
            for rec in SeqIO.parse(path, 'genbank'):
                for feat in rec.features:
                    ftypes.add(feat.type)
                    if feat.qualifiers:
                        qkeys.update(feat.qualifiers.keys())
                        for dbx in feat.qualifiers.get('db_xref', []):
                            typed = split_dbxref_to_typed_pairs(dbx)
                            if typed:
                                col, _ = typed
                                dbxref_types.add(col)
                                dbxref_counts[col] += 1
        except Exception as e:
            print(f"[WARN] Could not scan {fname}: {e}")

    return {
        'qualifier_keys': qkeys,
        'dbxref_types': dbxref_types,
        'feature_types': ftypes,
        'dbxref_counts': dbxref_counts,
    }

def load_features_index(gbff_path: str) -> List[SeqFeature]:
    """Load features from a .gbff/.gbk file (aggregate across records if any)."""
    feats: List[SeqFeature] = []
    for rec in SeqIO.parse(gbff_path, 'genbank'):
        feats.extend(rec.features)
    feats = [f for f in feats if isinstance(f.location, FeatureLocation) or isinstance(f.location, CompoundLocation)]
    return feats

# =========================
# Extraction per feature
# =========================

def extract_fields_from_feature(
    feat: SeqFeature,
    selected_fields: List[str],
    multi_sep: str
) -> Dict[str, str]:
    """
    From a feature, pull selected fields (qualifiers and typed db_xrefs).
    Returns a dict {field -> string value} where multiple values inside
    the same feature are joined using multi_sep (default '{}').
    """
    out: Dict[str, List[str]] = defaultdict(list)
    quals = feat.qualifiers or {}

    # Qualifier keys
    for field in selected_fields:
        if field.startswith('db_xref:'):
            continue
        if field in quals:
            for v in quals[field]:
                if isinstance(v, str):
                    out[field].append(v.strip())

    # db_xref typed fields
    if 'db_xref' in quals:
        for dbx in quals['db_xref']:
            typed = split_dbxref_to_typed_pairs(dbx)
            if typed:
                col, val = typed
                if col in selected_fields:
                    out[col].append(val)

    # Join multi-values for the same feature
    joined: Dict[str, str] = {}
    for k in selected_fields:
        if out.get(k):
            joined[k] = multi_sep.join(out[k])
        else:
            joined[k] = ''
    return joined

# =========================
# Main annotation
# =========================

def annotate_blocks(
    blocks_tsv: str,
    gbff_dir: str,
    out_tsv: str,
    summary_tsv: str,
    selected_fields: List[str],
    feature_types: List[str],
    entry_sep: str,
    multi_sep: str,
    include_truncated: bool,
    gbff_suffixes: Tuple[str, ...] = ('.gbff', '.gbk', '.gb')
):
    """
    Process the blocks table, annotate each row by overlapping GBFF features.
    - entry_sep: separator between different features (default '<>')
    - multi_sep: separator for multiple values inside the same feature (default '{}')
    - include_truncated: include features overlapping but not fully contained in the block
    """
    # Cache GBFF features per file
    gbff_cache: Dict[str, List[SeqFeature]] = {}

    # Read blocks rows
    with open(blocks_tsv, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)
    base_cols = reader.fieldnames if reader.fieldnames else []

    # Prepare annotation headers
    add_cols = list(selected_fields)
    out_cols = base_cols + add_cols + ['truncated', '_feature_types_found', '_locus_tags']

    # Aggregation for summary (one line per Block_ID)
    # For counts allowing >1 per row, we just sum occurrences directly.
    per_block_counts: Dict[str, Dict[str, Counter]] = defaultdict(lambda: defaultdict(Counter))
    per_block_counts_trunc: Dict[str, Dict[str, Counter]] = defaultdict(lambda: defaultdict(Counter))
    per_block_tot_rows: Counter = Counter()
    per_block_group_code: Dict[str, str] = {}
    # (We will emit group_code from first occurrence; other base columns left blank in summary.)

    def find_gbff_for_seqname(seqname: str) -> Optional[str]:
        base = os.path.splitext(seqname)[0]
        for suf in gbff_suffixes:
            candidate = os.path.join(gbff_dir, base + suf)
            if os.path.exists(candidate):
                return candidate
        cand2 = os.path.join(gbff_dir, seqname)
        return cand2 if os.path.exists(cand2) else None

    # ------- Write annotated rows -------
    with open(out_tsv, 'w', newline='') as outf:
        writer = csv.DictWriter(outf, delimiter='\t', fieldnames=out_cols, extrasaction='ignore')
        writer.writeheader()

        for row in rows:
            bid = row['Block_ID']
            per_block_tot_rows[bid] += 1
            if 'group_code' in row:
                per_block_group_code.setdefault(bid, row['group_code'])

            seqname = row['sequence_name']
            try:
                bstart = int(row['start_coords'])
                bend   = int(row['end_coords'])
            except Exception:
                raise RuntimeError("Could not parse 'start_coords'/'end_coords' as integers.")

            gbff_path = find_gbff_for_seqname(seqname)
            feat_types_found: Set[str] = set()
            locus_tags: List[str] = []
            truncated_flags: List[str] = []  # 'complete' or 'truncated' per feature in output order

            # For each field, collect per-feature values (keep alignment by index)
            per_field_feature_values: Dict[str, List[str]] = {k: [] for k in selected_fields}

            if gbff_path is not None:
                if gbff_path not in gbff_cache:
                    gbff_cache[gbff_path] = load_features_index(gbff_path)
                feats = gbff_cache[gbff_path]

                # Sort features by their leftmost coordinate for stable output
                def feat_key(ft: SeqFeature):
                    s, e = feature_bounds_1based(ft)
                    return (s, e, ft.type)
                feats = sorted(feats, key=feat_key)

                for feat in feats:
                    if feature_types and feat.type not in feature_types:
                        continue
                    fstart, fend = feature_bounds_1based(feat)
                    if not overlaps_1based(bstart, bend, fstart, fend):
                        continue

                    complete = fully_contained_1based(fstart, fend, bstart, bend)
                    if not complete and not include_truncated:
                        continue  # skip truncated features when not requested

                    feat_types_found.add(feat.type)
                    # locus_tag is handy for debugging/traceability
                    for lt in (feat.qualifiers or {}).get('locus_tag', []):
                        locus_tags.append(lt)

                    # Extract selected fields for this feature
                    joined_vals = extract_fields_from_feature(feat, selected_fields, multi_sep=multi_sep)
                    for field in selected_fields:
                        per_field_feature_values[field].append(joined_vals.get(field, ''))

                    truncated_flags.append('complete' if complete else 'truncated')

            # Build output row with joined values per field (per feature order)
            out_row = dict(row)
            for field in selected_fields:
                vals = per_field_feature_values[field]
                # Do not deduplicate; keep order and multiplicity
                out_row[field] = entry_sep.join(vals) if vals else ''

                # --- aggregate into per-block counts (including duplicates) ---
                # Count each non-empty value occurrence; if a feature had multi values inside (multi_sep),
                # split them so each inner item counts individually.
                for i, v in enumerate(vals):
                    if not v:
                        continue
                    # split inner multi-values
                    inner_vals = [x for x in v.split(multi_sep) if x]
                    is_trunc = (i < len(truncated_flags) and truncated_flags[i] == 'truncated')
                    for inner in inner_vals:
                        if is_trunc:
                            per_block_counts_trunc[bid][field][inner] += 1
                        else:
                            per_block_counts[bid][field][inner] += 1

            out_row['truncated'] = entry_sep.join(truncated_flags) if truncated_flags else ''
            out_row['_feature_types_found'] = entry_sep.join(sorted(feat_types_found)) if feat_types_found else ''
            out_row['_locus_tags'] = entry_sep.join(locus_tags) if locus_tags else ''
            writer.writerow(out_row)

    # ------- Write summary (one line per Block_ID) -------
    # Columns: Block_ID, group_code, sequence_name, start_coords, end_coords, core_or_accessory, then selected fields.
    # We fill sequence_name/start/end/core_or_accessory as empty (multi-valued across rows).
    s_base = ['Block_ID', 'group_code', 'sequence_name', 'start_coords', 'end_coords', 'core_or_accessory']
    s_cols = s_base + selected_fields
    with open(summary_tsv, 'w', newline='') as sf:
        sw = csv.DictWriter(sf, delimiter='\t', fieldnames=s_cols)
        sw.writeheader()
        for bid in sorted(per_block_tot_rows.keys()):
            total_rows = per_block_tot_rows[bid]
            gcode = per_block_group_code.get(bid, '')
            out_line = {
                'Block_ID': bid,
                'group_code': gcode,
                'sequence_name': '',
                'start_coords': '',
                'end_coords': '',
                'core_or_accessory': '',
            }
            # for each selected field, compose VALUE{n=COUNT|pct=PCT} items, joining by entry_sep
            pieces_by_field: Dict[str, List[str]] = {}
            for field in selected_fields:
                pieces: List[str] = []
                # normal entries
                for val, cnt in per_block_counts[bid][field].items():
                    pct = (100.0 * cnt / total_rows) if total_rows else 0.0
                    pieces.append(f"{val}" + "{" + f"n={cnt}|pct={pct:.2f}" + "}")
                # truncated entries (suffix)
                for val, cnt in per_block_counts_trunc[bid][field].items():
                    pct = (100.0 * cnt / total_rows) if total_rows else 0.0
                    pieces.append(f"{val}_truncated" + "{" + f"n={cnt}|pct={pct:.2f}" + "}")
                out_line[field] = entry_sep.join(pieces)
            sw.writerow(out_line)

# =========================
# CLI
# =========================

def parse_args():
    p = argparse.ArgumentParser(
        description="Annotate block ranges from Bakta GenBank (.gbff) files and summarize by Block_ID."
    )
    sub = p.add_subparsers(dest='cmd', required=True)

    # list-columns
    l = sub.add_parser('list-columns', help="Scan GBFF directory and list available qualifier keys, db_xref types (with counts), and feature types.")
    l.add_argument('--gbff-dir', required=True, help="Directory containing .gbff/.gbk files.")
    l.add_argument('--sample-limit', type=int, default=None, help="Optional: only scan the first N files.")

    # annotate
    a = sub.add_parser('annotate', help="Annotate a blocks TSV and produce a 1-line-per-block summary.")
    a.add_argument('--blocks-tsv', required=True, help="Path to blocks_core_accessory.tsv (tab-separated).")
    a.add_argument('--gbff-dir', required=True, help="Directory containing .gbff/.gbk files.")
    a.add_argument('--out-tsv', required=True, help="Output annotated TSV (row-per-block-instance).")
    a.add_argument('--summary-tsv', required=True, help="Output summary TSV (one line per Block_ID).")
    a.add_argument('--fields', required=False,
                   default="product,gene,locus_tag,db_xref:RefSeq,db_xref:UniRef100,db_xref:UniRef90,db_xref:UniRef50,db_xref:UniParc,db_xref:GO,db_xref:SO,db_xref:COG,db_xref:KEGG,db_xref:EC,db_xref:PFAM,db_xref:VFDB",
                   help="Comma-separated fields to add. Qualifiers (e.g., product,gene,locus_tag) and typed dbxrefs like db_xref:RefSeq, db_xref:GO, etc.")
    a.add_argument('--feature-types', required=False, default="CDS,gene,tRNA,rRNA,ncRNA,regulatory",
                   help="Comma-separated feature types to consider (e.g., CDS,gene,tRNA,rRNA). Use 'All' for no filtering.")
    a.add_argument('--entry-sep', required=False, default="<>",
                   help="Separator between different features in output fields (default '<>').")
    a.add_argument('--multi-sep', required=False, default="{}",
                   help="Separator for multiple values **inside the same feature** (default '{}').")
    a.add_argument('--include-truncated', action='store_true',
                   help="Include features that overlap the block but are not fully contained; adds truncated accounting.")
    return p.parse_args()

def main():
    args = parse_args()

    if args.cmd == 'list-columns':
        found = discover_available_columns(args.gbff_dir, args.sample_limit)
        print("# Qualifier keys (use directly in --fields):")
        print("\n".join(sorted(found['qualifier_keys'])))
        print("\n# db_xref types (use as db_xref:TYPE in --fields):")
        print("\n".join(sorted(found['dbxref_types'])))
        if found.get('dbxref_counts'):
            print("\n# db_xref namespaces found (feature count occurrences):")
            for k, v in sorted(found['dbxref_counts'].items(), key=lambda kv: -kv[1]):
                print(f"{k.replace('db_xref:','')}\t{v}")
        print("\n# Feature types (can be used in --feature-types):")
        print("\n".join(sorted(found['feature_types'])))
        return

    if args.cmd == 'annotate':
        fields = [x.strip() for x in args.fields.split(',') if x.strip()]
        ftypes = [x.strip() for x in args.feature_types.split(',') if x.strip()]
        if len(ftypes) == 1 and ftypes[0].lower() == 'all':
            ftypes = []  # no filtering
        annotate_blocks(
            blocks_tsv=args.blocks_tsv,
            gbff_dir=args.gbff_dir,
            out_tsv=args.out_tsv,
            summary_tsv=args.summary_tsv,
            selected_fields=fields,
            feature_types=ftypes,
            entry_sep=args.entry_sep,
            multi_sep=args.multi_sep,
            include_truncated=args.include_truncated
        )

if __name__ == '__main__':
    main()
