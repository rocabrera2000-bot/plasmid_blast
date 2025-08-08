#!/usr/bin/env python3
"""Run nucmer pairwise on all plasmid FASTA files in a directory.

This script executes MUMmer's ``nucmer`` for every pair of FASTA files in a
user-supplied directory.  Each comparison generates a ``.delta`` file and a
corresponding ``.coords`` file summarised with ``show-coords``.  All coordinate
records from every pair are collated into ``all_coords.tsv`` for downstream
analysis.

Example
-------
```
python all_vs_all_nucmer.py plasmids/ -o results/ --nucmer-args "--maxmatch"
```
"""

from __future__ import annotations

import argparse
import csv
import itertools
import os
from pathlib import Path
import shlex
import subprocess
from typing import Iterable, List, Tuple


def find_fasta_files(directory: Path) -> List[Path]:
    """Return a sorted list of FASTA files within *directory*."""
    exts = {".fa", ".fasta", ".fna", ".ffn"}
    return sorted(
        [p for p in directory.iterdir() if p.suffix.lower() in exts and p.is_file()]
    )


def pairwise_files(files: List[Path], include_self: bool, both_dirs: bool) -> Iterable[Tuple[Path, Path]]:
    """Generate file pairs based on the user's options."""
    if both_dirs:
        return itertools.product(files, files) if include_self else (
            (a, b) for a, b in itertools.product(files, files) if a != b
        )
    if include_self:
        return itertools.combinations_with_replacement(files, 2)
    return itertools.combinations(files, 2)


def run_nucmer(ref: Path, qry: Path, outdir: Path, extra_args: str) -> Path:
    """Execute nucmer and return the resulting ``.delta`` path."""
    prefix = outdir / f"{ref.stem}_vs_{qry.stem}"
    cmd = ["nucmer", "-p", str(prefix)] + shlex.split(extra_args) + [str(ref), str(qry)]
    subprocess.run(cmd, check=True)
    return prefix.with_suffix(".delta")


def run_show_coords(delta: Path, outdir: Path, extra_args: str) -> Path:
    """Execute show-coords and return the path to the ``.coords`` file."""
    coords_path = delta.with_suffix(".coords")
    cmd = [
        "show-coords",
        "-rclTH",
    ] + shlex.split(extra_args) + [str(delta)]
    with coords_path.open("w") as fh:
        subprocess.run(cmd, check=True, stdout=fh)
    return coords_path


def parse_coords_file(coords_path: Path, ref_file: Path, qry_file: Path) -> List[List[str]]:
    """Read a ``.coords`` file and return rows annotated with file names."""
    rows: List[List[str]] = []
    with coords_path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue
            rows.append([ref_file.name, qry_file.name] + row)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path, help="Directory containing FASTA plasmid files")
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("nucmer_results"),
        help="Directory in which to store nucmer outputs",
    )
    parser.add_argument("--include-self", action="store_true", help="Include self comparisons")
    parser.add_argument(
        "--both-directions",
        action="store_true",
        help="Run both reference/query orientations for each pair",
    )
    parser.add_argument(
        "--nucmer-args", default="", help="Additional arguments to pass to nucmer"
    )
    parser.add_argument(
        "--coords-args", default="", help="Additional arguments to pass to show-coords"
    )
    args = parser.parse_args()

    files = find_fasta_files(args.input_dir)
    if len(files) < 2 and not args.include_self:
        raise ValueError("Need at least two FASTA files for pairwise comparison")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Run pairwise comparisons
    all_rows: List[List[str]] = []
    header = [
        "ref_file",
        "query_file",
        "S1",
        "E1",
        "S2",
        "E2",
        "LEN_1",
        "LEN_2",
        "%_IDY",
        "LEN_R",
        "LEN_Q",
        "%COV_R",
        "%COV_Q",
        "REF_TAG",
        "QRY_TAG",
    ]

    for ref, qry in pairwise_files(files, args.include_self, args.both_directions):
        delta = run_nucmer(ref, qry, args.output_dir, args.nucmer_args)
        coords = run_show_coords(delta, args.output_dir, args.coords_args)
        all_rows.extend(parse_coords_file(coords, ref, qry))

    # Write aggregated table
    agg_path = args.output_dir / "all_coords.tsv"
    with agg_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(all_rows)


if __name__ == "__main__":
    main()
