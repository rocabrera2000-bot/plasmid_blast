#!/usr/bin/env python3
import argparse
import sys
import re
from pathlib import Path
import pandas as pd
import plotly.graph_objects as go

UNKNOWN_SENTINELS = {"0", 0, "", "NA", "N/A", "#N/A", "NaN", "nan", None}

def parse_group_index(name: str) -> int | None:
    """
    Extracts the trailing integer index from group names like:
      group_40kb_1, group_40kb_25
    Returns None if it can't find a clean integer at the end.
    """
    # Accept names that may end exactly with _<int> or _<int>_alt1
    m = re.match(r"^(.+?)_(\d+)(?:_alt1)?$", name)
    if not m:
        return None
    try:
        return int(m.group(2))
    except Exception:
        return None

def parse_group_report(path: Path, include_alt1: bool, max_groups: int | None) -> pd.DataFrame:
    """
    Returns a long dataframe with columns: ['group', 'Filename'].
    - Skips groups ending with '_alt1' unless include_alt1=True.
    - If max_groups is given, only retains groups with trailing index <= max_groups,
      e.g. group_40kb_1..group_40kb_N.
    """
    df = pd.read_csv(path, sep="\t", header=0, dtype=str)
    expected = {"group", "size", "members"}
    if not expected.issubset(df.columns):
        raise ValueError(f"{path} must have columns: group, size, members")

    if not include_alt1:
        df = df[~df["group"].str.endswith("_alt1", na=False)].copy()

    if max_groups is not None:
        keep_mask = df["group"].apply(lambda g: (parse_group_index(str(g)) is not None) and (parse_group_index(str(g)) <= max_groups))
        df = df[keep_mask].copy()

    rows = []
    for _, r in df.iterrows():
        g = str(r["group"])
        members = (r["members"] or "").split(",")
        for m in members:
            m = m.strip()
            if m:
                rows.append((g, m))
    if not rows:
        return pd.DataFrame(columns=["group", "Filename"])
    return pd.DataFrame(rows, columns=["group", "Filename"])

def normalize_unknowns(val, include_zeros: bool, unknown_label: str):
    if pd.isna(val):
        v = None
    else:
        v = str(val)
    if v in UNKNOWN_SENTINELS:
        return unknown_label if include_zeros else None
    return v

def build_sankey_links(df_chain: pd.DataFrame, stages: list[str]) -> tuple[list[str], list[int], list[int], list[float], list[int]]:
    """
    df_chain columns: ['group'] + stages
    Returns:
      node_labels, src, dst, values, layer_offsets_per_stage
      (layer_offsets_per_stage holds the starting node index for each layer)
    """
    layer_nodes = []
    layer_maps = []
    layer_offsets = []

    # Layer 0: groups
    groups = sorted(df_chain['group'].dropna().unique().tolist())
    layer_nodes.append(groups)
    offset = 0
    layer_offsets.append(offset)
    layer_maps.append({name: offset + i for i, name in enumerate(groups)})
    offset += len(groups)

    # Subsequent layers
    for col in stages:
        vals = df_chain[col].dropna()
        order = vals.value_counts().sort_values(ascending=False).index.tolist()
        layer_nodes.append(order)
        layer_offsets.append(offset)
        layer_maps.append({name: offset + i for i, name in enumerate(order)})
        offset += len(order)

    # Links
    src, dst, val = [], [], []

    # group -> first stage
    first = stages[0]
    agg = df_chain.groupby(['group', first], dropna=True).size().reset_index(name='count')
    for _, r in agg.iterrows():
        g, vlabel, c = r['group'], r[first], r['count']
        if pd.isna(g) or pd.isna(vlabel):
            continue
        s = layer_maps[0].get(g)
        d = layer_maps[1].get(vlabel)
        if s is None or d is None:
            continue
        src.append(s); dst.append(d); val.append(float(c))

    # stage i -> stage i+1
    for i in range(1, len(stages)):
        a = stages[i-1]
        b = stages[i]
        agg = df_chain.groupby([a, b], dropna=True).size().reset_index(name='count')
        for _, r in agg.iterrows():
            aval, bval, c = r[a], r[b], r['count']
            if pd.isna(aval) or pd.isna(bval):
                continue
            s = layer_maps[i].get(aval)
            d = layer_maps[i+1].get(bval)
            if s is None or d is None:
                continue
            src.append(s); dst.append(d); val.append(float(c))

    node_labels = [n for layer in layer_nodes for n in layer]
    return node_labels, src, dst, val, layer_offsets

def make_link_colors(src, dst, node_labels, mode: str, stages: list[str], layer_offsets: list[int]) -> list[str]:
    """
    Returns hex colors for each link.
      mode:
        - 'none'   : all default (Plotly colors)
        - 'source' : color by source node
        - 'target' : color by target node
        - 'stage:<name>' : color by the specified metadata stage (e.g., 'stage:Capsule_type')
    Deterministic palette cycling.
    """
    if mode == "none":
        return None

    # Simple repeating palette (12 distinct colors)
    palette = [
        "#1f77b4","#ff7f0e","#2ca02c","#d62728",
        "#9467bd","#8c564b","#e377c2","#7f7f7f",
        "#bcbd22","#17becf","#9edae5","#c7c7c7"
    ]
    def color_for_index(idx: int) -> str:
        return palette[idx % len(palette)]

    if mode == "source":
        return [color_for_index(s) for s in src]
    if mode == "target":
        return [color_for_index(t) for t in dst]

    if mode.startswith("stage:"):
        stage_name = mode.split(":",1)[1]
        # Find which stage index this is
        if stage_name not in stages:
            # fallback to source
            return [color_for_index(s) for s in src]
        stage_i = stages.index(stage_name) + 1  # +1 because layer 0 is 'group'
        start = layer_offsets[stage_i]
        # Assign colors based on the target node when link enters that stage
        colors = []
        for s, d in zip(src, dst):
            # If link's target is within that stage's node index range, color by target
            colors.append(color_for_index(d - start) if d >= start else color_for_index(s))
        return colors

    # Fallback
    return [color_for_index(s) for s in src]

def main():
    ap = argparse.ArgumentParser(description="Build a Sankey diagram linking plasmid group(s) to metadata columns.")
    ap.add_argument("--group-report", required=True, type=Path)
    ap.add_argument("--metadata", required=True, type=Path, help="TSV containing a 'Filename' column")
    ap.add_argument("--columns", default="Capsule_type",
                    help="Comma-separated metadata columns for chained stages, e.g., 'Capsule_type' or 'Capsule_type,PC'")
    ap.add_argument("--include-alt1", action="store_true",
                    help="Include groups ending with '_alt1' (default: excluded)")
    ap.add_argument("--max-groups", type=int, default=None,
                    help="Keep only groups named like group_*_1..group_*_N (default: keep all)")
    ap.add_argument("--include-zeros", action="store_true",
                    help="Include zero/unknown values and label them as Unknown (default: exclude them)")
    ap.add_argument("--unknown-label", default="Unknown")
    ap.add_argument("--link-color-mode", default="none",
                    choices=["none", "source", "target"] + [f"stage:{k}" for k in ["Capsule_type","PC","any"]],
                    help="Color links by 'source', 'target', or a specific metadata stage (use 'stage:<ColumnName>').")
    ap.add_argument("--output", required=True, type=Path, help="Output HTML file")
    ap.add_argument("--title", default=None)
    args = ap.parse_args()

    # Read & filter
    df_groups = parse_group_report(args.group_report, include_alt1=args.include_alt1, max_groups=args.max_groups)
    if df_groups.empty:
        print("No (group, Filename) pairs after filtering groups. Exiting.", file=sys.stderr)
        sys.exit(1)

    df_meta = pd.read_csv(args.metadata, sep="\t", header=0, dtype=str)
    if "Filename" not in df_meta.columns:
        raise ValueError(f"{args.metadata} must contain a 'Filename' column")

    stages = [c.strip() for c in args.columns.split(",") if c.strip()]
    missing = [c for c in stages if c not in df_meta.columns]
    if missing:
        raise ValueError(f"Column(s) not found in metadata: {missing}")

    # Join (many groups -> one metadata row)
    df = df_groups.merge(df_meta, on="Filename", how="left", validate="many_to_one")

    # Unknown handling
    for c in stages:
        df[c] = df[c].apply(lambda v: normalize_unknowns(v, args.include_zeros, args.unknown_label))
    if not args.include_zeros:
        mask = ~pd.isna(df[stages]).any(axis=1)
        df = df[mask].copy()

    if df.empty:
        print("No data after unknown/zero filtering. Exiting.", file=sys.stderr)
        sys.exit(1)

    # Build Sankey
    node_labels, src, dst, values, layer_offsets = build_sankey_links(df_chain=df[["group"] + stages], stages=stages)
    if not src:
        print("No links remain after filtering; nothing to plot. Exiting.", file=sys.stderr)
        sys.exit(1)

    # Link colors
    link_colors = make_link_colors(src, dst, node_labels, args.link_color_mode, stages, layer_offsets)

    sankey = go.Sankey(
        arrangement="snap",
        node=dict(label=node_labels, pad=12, thickness=16),
        link=dict(source=src, target=dst, value=values, **({"color": link_colors} if link_colors else {}))
    )
    fig = go.Figure(data=[sankey])
    fig.update_layout(title_text=(args.title or f"Groups → {' → '.join(stages)}"), font_size=12)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(args.output), include_plotlyjs="cdn")
    print(f"Wrote {args.output}")

if __name__ == "__main__":
    main()
