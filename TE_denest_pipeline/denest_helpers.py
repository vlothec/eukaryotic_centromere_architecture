#!/usr/bin/env python3
"""Small, dependency-free helpers used by denest_pipeline.sh.

All coordinates handled here are BED coordinates: zero-based and half-open.
Keeping the sequence editing and coordinate projection in one program prevents
the common off-by-one errors that arise when FASTA and BED conventions mix.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable


def read_fasta(path: str) -> tuple[str, str]:
    """Read exactly one FASTA record and return its header and sequence."""
    header = None
    sequence: list[str] = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    raise ValueError(f"{path} contains more than one FASTA record")
                header = line[1:].split()[0]
            else:
                sequence.append(line)
    if header is None:
        raise ValueError(f"{path} contains no FASTA record")
    return header, "".join(sequence)


def write_fasta(path: str, name: str, sequence: str) -> None:
    with open(path, "w") as handle:
        handle.write(f">{name}\n")
        for start in range(0, len(sequence), 60):
            handle.write(sequence[start:start + 60] + "\n")


def read_bed(path: str) -> list[tuple[int, int]]:
    intervals: list[tuple[int, int]] = []
    with open(path) as handle:
        for number, line in enumerate(handle, start=1):
            fields = line.split()
            if len(fields) < 3 or line.startswith("#"):
                continue
            start, end = int(fields[1]), int(fields[2])
            if start < 0 or end <= start:
                raise ValueError(f"Invalid BED interval at {path}:{number}")
            intervals.append((start, end))
    return intervals


def merge(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping and touching intervals."""
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def write_bed(path: str, name: str, intervals: Iterable[tuple[int, int]]) -> None:
    with open(path, "w") as handle:
        for start, end in intervals:
            handle.write(f"{name}\t{start}\t{end}\n")


def cmd_excise(args: argparse.Namespace) -> None:
    name, sequence = read_fasta(args.fasta)
    intervals = merge(read_bed(args.intervals))
    if intervals and intervals[-1][1] > len(sequence):
        raise ValueError("An interval extends beyond the current sequence")

    kept: list[str] = []
    position = 0
    for start, end in intervals:
        kept.append(sequence[position:start])
        position = end
    kept.append(sequence[position:])
    write_fasta(args.output, args.name or name, "".join(kept))


def cmd_project(args: argparse.Namespace) -> None:
    """Project current-sequence removals to original ROI coordinates.

    The map is a BED file whose intervals are the still-retained blocks in the
    original ROI.  The removal BED uses coordinates of their concatenation.
    """
    mapping = read_bed(args.mapping)
    removals = merge(read_bed(args.removals))
    next_mapping: list[tuple[int, int]] = []
    removed_original: list[tuple[int, int]] = []
    current_start = 0

    for original_start, original_end in mapping:
        current_end = current_start + original_end - original_start
        for removal_start, removal_end in removals:
            overlap_start = max(current_start, removal_start)
            overlap_end = min(current_end, removal_end)
            if overlap_start >= overlap_end:
                continue
            original_overlap_start = original_start + overlap_start - current_start
            original_overlap_end = original_start + overlap_end - current_start
            removed_original.append((original_overlap_start, original_overlap_end))

        cursor = current_start
        for removal_start, removal_end in removals:
            start = max(current_start, removal_start)
            end = min(current_end, removal_end)
            if start >= end:
                continue
            if cursor < start:
                next_mapping.append((
                    original_start + cursor - current_start,
                    original_start + start - current_start,
                ))
            cursor = max(cursor, end)
        if cursor < current_end:
            next_mapping.append((
                original_start + cursor - current_start,
                original_end,
            ))
        current_start = current_end

    write_bed(args.next_map, args.name, next_mapping)
    write_bed(args.removed_original, args.name, removed_original)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    excise = subparsers.add_parser("excise", help="Remove local BED intervals from one FASTA record")
    excise.add_argument("--fasta", required=True)
    excise.add_argument("--intervals", required=True)
    excise.add_argument("--output", required=True)
    excise.add_argument("--name")
    excise.set_defaults(func=cmd_excise)

    project = subparsers.add_parser("project", help="Project removals to original ROI coordinates")
    project.add_argument("--mapping", required=True)
    project.add_argument("--removals", required=True)
    project.add_argument("--next-map", required=True)
    project.add_argument("--removed-original", required=True)
    project.add_argument("--name", required=True)
    project.set_defaults(func=cmd_project)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
