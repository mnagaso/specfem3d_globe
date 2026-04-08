#!/usr/bin/env python3

"""Fix io_band_*.txt files produced before shard-local multi-IO accounting.

The buggy multi-IO implementation logged the full global byte count on every
I/O rank. For a run with ``N`` I/O ranks, each per-rank entry therefore needs
to be rescaled from ``global_bytes`` to the shard-local byte count.

By default this script assumes an even split across I/O ranks and divides each
rank line by ``N``. If the shard sizes were uneven, provide ``--weights`` as a
comma-separated list in local I/O-rank order. The weights are normalized and
applied to the per-rank lines; the total line is always divided by ``N``.
"""

import argparse
import pathlib
import re
import sys
from typing import Iterable, List, Match, Optional, Sequence, Set, Tuple


RANK_LINE_RE = re.compile(
    r"^(?P<prefix>\s*mygroup:\s*\d+\s*,\s*myrank:\s*(?P<rank>\d+)\s*,\s*bytes_written:\s*)"
    r"(?P<bytes>\d+)"
    r"(?P<middle>\s*,\s*elapsed_time \(s\):\s*[0-9.Ee+\-]+\s*,\s*bandwidth:\s*)"
    r"(?P<bandwidth>[0-9.Ee+\-]+)"
    r"(?P<suffix>\s+MB/s.*)$"
)

TOTAL_LINE_RE = re.compile(
    r"^(?P<prefix>\s*mygroup:\s*\d+\s*,\s*total_bytes_written:\s*)"
    r"(?P<bytes>[0-9.Ee+\-]+)"
    r"(?P<middle>\s*,\s*max_elapsed_time \(s\):\s*[0-9.Ee+\-]+\s*,\s*total_bandwidth:\s*)"
    r"(?P<bandwidth>[0-9.Ee+\-]+)"
    r"(?P<suffix>\s+MB/s.*)$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fix io_band_*.txt files generated before shard-local multi-IO byte accounting."
    )
    parser.add_argument(
        "paths",
        nargs="*",
        default=["."],
        help="Files or directories containing io_band_*.txt files. Defaults to the current directory.",
    )
    parser.add_argument(
        "--io-nodes",
        type=int,
        required=True,
        help="Number of dedicated I/O ranks used in the buggy run.",
    )
    parser.add_argument(
        "--weights",
        help=(
            "Optional comma-separated shard weights in local I/O-rank order. "
            "Use this when shard sizes were uneven. The list is normalized automatically."
        ),
    )
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Deprecated. In-place writes are disabled to preserve the original logs.",
    )
    parser.add_argument(
        "--suffix",
        default="_mod",
        help=(
            "Suffix inserted before the original file extension. "
            "Default: '_mod' (for example io_band_0.txt -> io_band_0_mod.txt)."
        ),
    )
    parser.add_argument(
        "--skip-empty",
        action="store_true",
        help=(
            "Skip empty or whitespace-only input logs instead of failing. "
            "Skipped files are left untouched."
        ),
    )
    return parser.parse_args()


def normalize_weights(io_nodes: int, weights_arg: Optional[str]) -> List[float]:
    if io_nodes <= 1:
        raise ValueError("--io-nodes must be greater than 1 for this fixer")

    if weights_arg is None:
        return [1.0 / io_nodes] * io_nodes

    raw = [float(part.strip()) for part in weights_arg.split(",") if part.strip()]
    if len(raw) != io_nodes:
        raise ValueError(
            f"expected {io_nodes} weights for --weights, received {len(raw)}"
        )

    total = sum(raw)
    if total <= 0.0:
        raise ValueError("--weights must sum to a positive value")

    return [value / total for value in raw]


def iter_target_files(paths: Sequence[str]) -> Iterable[pathlib.Path]:
    seen = set()  # type: Set[pathlib.Path]

    for raw_path in paths:
        path = pathlib.Path(raw_path)
        if path.is_dir():
            for candidate in sorted(path.rglob("io_band_*.txt")):
                resolved = candidate.resolve()
                if resolved not in seen:
                    seen.add(resolved)
                    yield candidate
        elif path.is_file():
            resolved = path.resolve()
            if resolved not in seen:
                seen.add(resolved)
                yield path
        else:
            raise FileNotFoundError(f"path does not exist: {path}")


def format_rank_line(match: Match[str], weights: Sequence[float]) -> str:
    rank = int(match.group("rank"))
    if rank >= len(weights):
        raise ValueError(
            f"rank {rank} found in log line but only {len(weights)} weights were provided"
        )

    original_bytes = int(match.group("bytes"))
    original_bandwidth = float(match.group("bandwidth"))
    factor = weights[rank]

    corrected_bytes = int(round(original_bytes * factor))
    corrected_bandwidth = original_bandwidth * factor

    return (
        f"{match.group('prefix')}{corrected_bytes}"
        f"{match.group('middle')}{corrected_bandwidth:.6f}"
        f"{match.group('suffix')}"
    )


def format_total_line(match: Match[str], io_nodes: int) -> str:
    original_bytes = float(match.group("bytes"))
    original_bandwidth = float(match.group("bandwidth"))

    corrected_bytes = original_bytes / io_nodes
    corrected_bandwidth = original_bandwidth / io_nodes

    return (
        f"{match.group('prefix')}{corrected_bytes:.0f}"
        f"{match.group('middle')}{corrected_bandwidth:.6f}"
        f"{match.group('suffix')}"
    )


def fix_log_contents(text: str, io_nodes: int, weights: Sequence[float]) -> str:
    fixed_lines = []  # type: List[str]

    for line in text.splitlines():
        rank_match = RANK_LINE_RE.match(line)
        if rank_match:
                        fixed_lines.append(format_rank_line(rank_match, weights))
                        continue

        total_match = TOTAL_LINE_RE.match(line)
        if total_match:
                        fixed_lines.append(format_total_line(total_match, io_nodes))
                        continue

        fixed_lines.append(line)

    trailing_newline = "\n" if text.endswith("\n") else ""
    return "\n".join(fixed_lines) + trailing_newline


def output_path(input_path: pathlib.Path, suffix: str) -> pathlib.Path:
    if input_path.suffix:
        return input_path.with_name(f"{input_path.stem}{suffix}{input_path.suffix}")
    return input_path.with_name(input_path.name + suffix)


def build_fix_plan(
    files: Sequence[pathlib.Path],
    io_nodes: int,
    weights: Sequence[float],
    skip_empty: bool,
) -> Tuple[List[Tuple[pathlib.Path, str]], List[pathlib.Path]]:
    plan = []  # type: List[Tuple[pathlib.Path, str]]
    skipped_empty = []  # type: List[pathlib.Path]

    for file_path in files:
        original_text = file_path.read_text(encoding="utf-8")
        if not original_text.strip():
            if skip_empty:
                skipped_empty.append(file_path)
                continue
            raise ValueError(
                f"encountered empty io_band log: {file_path}; "
                "rerun with --skip-empty to leave it untouched"
            )

        fixed_text = fix_log_contents(original_text, io_nodes, weights)
        plan.append((file_path, fixed_text))

    return plan, skipped_empty


def main() -> int:
    args = parse_args()

    if args.in_place:
        print(
            "error: in-place writes are disabled; use the generated *_mod.txt sibling files",
            file=sys.stderr,
        )
        return 2

    try:
        weights = normalize_weights(args.io_nodes, args.weights)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    try:
        files = list(iter_target_files(args.paths))
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    if not files:
        print("error: no io_band_*.txt files found", file=sys.stderr)
        return 2

    try:
        plan, skipped_empty = build_fix_plan(
            files, args.io_nodes, weights, args.skip_empty
        )
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    for file_path, fixed_text in plan:
        destination = output_path(file_path, args.suffix)
        destination.write_text(fixed_text, encoding="utf-8")
        print(f"wrote {destination}")

    for file_path in skipped_empty:
        print(f"skipped empty {file_path}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())