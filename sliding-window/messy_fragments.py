#!/usr/bin/env python3
"""Generate imperfect sliding-window fragments and compare them to a clean reference.

Run from any directory:

    python /path/to/sliding-window/messy_fragments.py

By default, the script reads the CDH1 and VIM FASTA files, uses non-overlapping
20-nt reference windows, makes length/offset/representative internal-edit
variants, and writes one CSV under sliding-window/data/exports/.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import RNA
from Bio import SeqIO


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_FASTA_FILES = (
    SCRIPT_DIR / "data" / "CDH1__NM_004360.5.fasta",
    SCRIPT_DIR / "data" / "VIM__NM_003380.5.fasta",
)
RNA_ALPHABET = "ACGU"
REVERSE_COMPLEMENT_TABLE = str.maketrans("AUCG", "UAGC")


@dataclass(frozen=True)
class Candidate:
    """A generated fragment and its coordinate mapping to its clean reference."""

    candidate_id: str
    variant_class: str
    description: str
    sequence: str
    source_start: int
    source_end: int
    offset_nt: int
    anchor: str
    # For each candidate nucleotide, the corresponding 0-based reference index.
    # Inserted or out-of-window nucleotides map to None.
    position_map: tuple[int | None, ...]


def reverse_complement(sequence: str) -> str:
    return sequence.translate(REVERSE_COMPLEMENT_TABLE)[::-1]


def levenshtein_distance(first: str, second: str) -> int:
    """Unit-cost global Levenshtein distance using O(min(n, m)) memory."""
    if len(first) < len(second):
        first, second = second, first
    previous = list(range(len(second) + 1))
    for i, first_base in enumerate(first, start=1):
        current = [i]
        for j, second_base in enumerate(second, start=1):
            current.append(
                min(
                    current[-1] + 1,
                    previous[j] + 1,
                    previous[j - 1] + (first_base != second_base),
                )
            )
        previous = current
    return previous[-1]


def make_candidates(
    transcript: str,
    reference_start: int,
    window_size: int,
    lengths: Iterable[int],
    offsets: Iterable[int],
    internal_edit_mode: str,
) -> list[Candidate]:
    """Create one clean window and separately labeled noisy versions of it."""
    reference_end = reference_start + window_size
    reference = transcript[reference_start:reference_end]
    candidates = [
        Candidate(
            candidate_id="clean_reference",
            variant_class="reference",
            description="Unmodified reference fragment",
            sequence=reference,
            source_start=reference_start,
            source_end=reference_end,
            offset_nt=0,
            anchor="reference",
            position_map=tuple(range(window_size)),
        )
    ]

    # Length variants use the same 5' coordinate as the reference. A longer
    # candidate may include source bases beyond the reference; those bases do
    # not map onto the reference structure.
    for length in sorted(set(lengths)):
        if length <= 0 or length == window_size:
            continue
        end = reference_start + length
        if end > len(transcript):
            continue
        candidates.append(
            Candidate(
                candidate_id=f"length_{length}_anchor_5p",
                variant_class="length",
                description=f"{length}-nt source window, anchored at the reference 5' end",
                sequence=transcript[reference_start:end],
                source_start=reference_start,
                source_end=end,
                offset_nt=0,
                anchor="5p",
                position_map=tuple(i if i < window_size else None for i in range(length)),
            )
        )

    # Offset variants keep the reference length and use transcript-derived
    # bases. The mapping records which positions still overlap the reference.
    for offset in sorted(set(offsets)):
        if offset == 0:
            continue
        start = reference_start + offset
        end = start + window_size
        if start < 0 or end > len(transcript):
            continue
        candidates.append(
            Candidate(
                candidate_id=f"offset_{offset:+d}",
                variant_class="offset",
                description=f"{window_size}-nt source window shifted {offset:+d} nt",
                sequence=transcript[start:end],
                source_start=start,
                source_end=end,
                offset_nt=offset,
                anchor="window",
                position_map=tuple(
                    i + offset if 0 <= i + offset < window_size else None
                    for i in range(window_size)
                ),
            )
        )

    if internal_edit_mode in {"center", "all"}:
        if internal_edit_mode == "center":
            substitution_positions = [window_size // 2]
            insertion_positions = [window_size // 2]
            deletion_positions = [window_size // 2]
        else:
            substitution_positions = list(range(window_size))
            insertion_positions = list(range(window_size + 1))
            deletion_positions = list(range(window_size))

        for position in substitution_positions:
            original = reference[position]
            for alternative in RNA_ALPHABET:
                if alternative == original:
                    continue
                edited = reference[:position] + alternative + reference[position + 1 :]
                candidates.append(
                    Candidate(
                        candidate_id=f"substitution_pos_{position + 1}_{original}_to_{alternative}",
                        variant_class="substitution",
                        description=(
                            f"Substitute reference position {position + 1}: "
                            f"{original}>{alternative}"
                        ),
                        sequence=edited,
                        source_start=reference_start,
                        source_end=reference_end,
                        offset_nt=0,
                        anchor="synthetic",
                        position_map=tuple(range(window_size)),
                    )
                )

        for position in insertion_positions:
            for inserted_base in RNA_ALPHABET:
                edited = reference[:position] + inserted_base + reference[position:]
                mapping = tuple(
                    i if i < position else (None if i == position else i - 1)
                    for i in range(window_size + 1)
                )
                candidates.append(
                    Candidate(
                        candidate_id=f"insertion_before_{position + 1}_{inserted_base}",
                        variant_class="insertion",
                        description=(
                            f"Insert {inserted_base} before reference position {position + 1}"
                        ),
                        sequence=edited,
                        source_start=reference_start,
                        source_end=reference_end,
                        offset_nt=0,
                        anchor="synthetic",
                        position_map=mapping,
                    )
                )

        for position in deletion_positions:
            edited = reference[:position] + reference[position + 1 :]
            mapping = tuple(i if i < position else i + 1 for i in range(window_size - 1))
            candidates.append(
                Candidate(
                    candidate_id=f"deletion_pos_{position + 1}_{reference[position]}",
                    variant_class="deletion",
                    description=f"Delete reference position {position + 1}",
                    sequence=edited,
                    source_start=reference_start,
                    source_end=reference_end,
                    offset_nt=0,
                    anchor="synthetic",
                    position_map=mapping,
                )
            )

    return candidates


def fold_mfe(sequence: str) -> float:
    fold_compound = RNA.fold_compound(sequence)
    _, mfe = fold_compound.mfe()
    return float(mfe)


def interstrand_pairs(structure: str, length_a: int, length_b: int) -> set[tuple[int, int]]:
    """Return (strand-A index, strand-B index) pairs from a cofold dot-bracket."""
    symbols: list[str] = []
    cut_point: int | None = None
    for symbol in structure:
        if symbol == "&":
            if cut_point is not None:
                raise ValueError(f"Expected one cofold cut point, got {structure!r}")
            cut_point = len(symbols)
        else:
            symbols.append(symbol)

    if cut_point is None:
        cut_point = length_a
    if cut_point != length_a or len(symbols) != length_a + length_b:
        raise ValueError(
            "Cofold structure length/cut point does not match the two input strands: "
            f"structure={structure!r}, lengths=({length_a}, {length_b})"
        )

    opening_to_closing = {"(": ")", "[": "]", "{": "}", "<": ">"}
    closing_to_opening = {value: key for key, value in opening_to_closing.items()}
    stacks: dict[str, list[int]] = {opening: [] for opening in opening_to_closing}
    all_pairs: list[tuple[int, int]] = []
    for position, symbol in enumerate(symbols):
        if symbol in opening_to_closing:
            stacks[symbol].append(position)
        elif symbol in closing_to_opening:
            opening = closing_to_opening[symbol]
            if not stacks[opening]:
                raise ValueError(f"Unbalanced cofold structure: {structure!r}")
            all_pairs.append((stacks[opening].pop(), position))

    if any(stack for stack in stacks.values()):
        raise ValueError(f"Unbalanced cofold structure: {structure!r}")

    result: set[tuple[int, int]] = set()
    for first, second in all_pairs:
        if first < cut_point <= second:
            result.add((first, second - cut_point))
        elif second < cut_point <= first:
            result.add((second, first - cut_point))
    return result


def longest_antiparallel_run(pairs: set[tuple[int, int]]) -> int:
    """Longest consecutive intermolecular helix in antiparallel orientation."""
    if not pairs:
        return 0
    ordered = sorted(pairs)
    longest = current = 1
    previous_a, previous_b = ordered[0]
    for strand_a, strand_b in ordered[1:]:
        if strand_a == previous_a + 1 and strand_b == previous_b - 1:
            current += 1
        else:
            current = 1
        longest = max(longest, current)
        previous_a, previous_b = strand_a, strand_b
    return longest


def normalized_distance(distance: int, first_length: int, second_length: int) -> float:
    denominator = max(first_length, second_length)
    return float(distance / denominator) if denominator else 0.0


def bp_set_distance(
    reference_pairs: set[tuple[int, int]],
    candidate_pairs: set[tuple[int, int]],
) -> tuple[int, float]:
    symmetric_difference = len(reference_pairs.symmetric_difference(candidate_pairs))
    union_size = len(reference_pairs.union(candidate_pairs))
    return symmetric_difference, (symmetric_difference / union_size if union_size else 0.0)


def shape_string(structure: str, level: int) -> str:
    """Abstract a cofold structure after removing its strand cut-point marker."""
    if not hasattr(RNA, "abstract_shapes"):
        raise RuntimeError(
            "This ViennaRNA installation does not expose RNA.abstract_shapes(). "
            "Install ViennaRNA 2.4.18 or newer."
        )
    return str(RNA.abstract_shapes(structure.replace("&", ""), level))


def cofold_metrics(
    candidate_sequence: str,
    partner_sequence: str,
    reference_pairs: set[tuple[int, int]],
    position_map: tuple[int | None, ...],
    reference_shape: str,
    shape_level: int,
    partner_mfe: float | None = None,
) -> dict[str, object]:
    mfe_candidate = fold_mfe(candidate_sequence)
    if partner_mfe is None:
        partner_mfe = fold_mfe(partner_sequence)

    compound = RNA.fold_compound(candidate_sequence + "&" + partner_sequence)
    structure, mfe_complex = compound.mfe()
    binding_dg = float(mfe_complex) - (mfe_candidate + partner_mfe)

    candidate_pairs = interstrand_pairs(
        str(structure), len(candidate_sequence), len(partner_sequence)
    )
    mapped_pairs = {
        (position_map[candidate_position], partner_position)
        for candidate_position, partner_position in candidate_pairs
        if candidate_position < len(position_map)
        and position_map[candidate_position] is not None
    }
    mapped_pairs = {(int(first), second) for first, second in mapped_pairs}
    bp_distance, bp_distance_norm = bp_set_distance(reference_pairs, mapped_pairs)

    candidate_shape = shape_string(str(structure), shape_level)
    shape_distance = levenshtein_distance(reference_shape, candidate_shape)
    shape_distance_norm = normalized_distance(
        shape_distance, len(reference_shape), len(candidate_shape)
    )

    max_possible_pairs = min(len(candidate_sequence), len(partner_sequence))
    binding_dg_per_possible_pair = binding_dg / max_possible_pairs if max_possible_pairs else 0.0

    return {
        "mfe_candidate": mfe_candidate,
        "mfe_partner": partner_mfe,
        "mfe_complex": float(mfe_complex),
        "binding_dg_mfe": binding_dg,
        "binding_dg_per_possible_pair": binding_dg_per_possible_pair,
        "structure_cofold": str(structure),
        "interstrand_pair_count": len(candidate_pairs),
        "longest_interstrand_run": longest_antiparallel_run(candidate_pairs),
        "mapped_bp_distance": bp_distance,
        "mapped_bp_distance_norm": bp_distance_norm,
        "mapped_candidate_pairs": sorted(mapped_pairs),
        "shape_level": shape_level,
        "reference_shape": reference_shape,
        "candidate_shape": candidate_shape,
        "shape_distance": shape_distance,
        "shape_distance_norm": shape_distance_norm,
        "shape_exact_match": int(reference_shape == candidate_shape),
    }


def lower_is_better_percentiles(values: list[float]) -> list[float]:
    """Map values to goodness percentiles (best=1, worst=0; tied ranks average)."""
    count = len(values)
    if count == 1:
        return [1.0]

    order = sorted(range(count), key=lambda index: values[index])
    percentiles = [0.0] * count
    start = 0
    while start < count:
        end = start + 1
        while end < count and values[order[end]] == values[order[start]]:
            end += 1
        average_zero_based_rank = (start + end - 1) / 2
        goodness = (count - 1 - average_zero_based_rank) / (count - 1)
        for index in order[start:end]:
            percentiles[index] = goodness
        start = end
    return percentiles


def average_ranks(values: list[float], higher_is_better: bool = False) -> list[float]:
    """Return 1-based average ranks, with ties receiving their mean rank."""
    count = len(values)
    order = sorted(
        range(count),
        key=lambda index: values[index],
        reverse=higher_is_better,
    )
    ranks = [0.0] * count
    start = 0
    while start < count:
        end = start + 1
        while end < count and values[order[end]] == values[order[start]]:
            end += 1
        average_one_based_rank = (start + 1 + end) / 2
        for index in order[start:end]:
            ranks[index] = average_one_based_rank
        start = end
    return ranks


def score_reference_group(rows: list[dict[str, object]]) -> None:
    sequence_percentiles = lower_is_better_percentiles(
        [float(row["edit_distance_norm"]) for row in rows]
    )
    energy_percentiles = lower_is_better_percentiles(
        [float(row["binding_dg_per_possible_pair"]) for row in rows]
    )
    structure_percentiles = lower_is_better_percentiles(
        [float(row["mapped_bp_distance_norm"]) for row in rows]
    )
    shape_percentiles = lower_is_better_percentiles(
        [float(row["shape_distance_norm"]) for row in rows]
    )

    for index, row in enumerate(rows):
        row["sequence_goodness_percentile"] = sequence_percentiles[index]
        row["energy_goodness_percentile"] = energy_percentiles[index]
        row["structure_goodness_percentile"] = structure_percentiles[index]
        row["shape_goodness_percentile"] = shape_percentiles[index]
        row["composite_score"] = (
            sequence_percentiles[index]
            + energy_percentiles[index]
            + structure_percentiles[index]
        ) / 3
        # Alternate structure ranking: swap abstract-shape distance in for
        # base-pair distance rather than counting both structure measures.
        row["shape_composite_score"] = (
            sequence_percentiles[index]
            + energy_percentiles[index]
            + shape_percentiles[index]
        ) / 3

    composite_ranks = average_ranks(
        [float(row["composite_score"]) for row in rows], higher_is_better=True
    )
    shape_composite_ranks = average_ranks(
        [float(row["shape_composite_score"]) for row in rows], higher_is_better=True
    )
    sequence_ranks = average_ranks(
        [float(row["edit_distance_norm"]) for row in rows]
    )
    energy_ranks = average_ranks(
        [float(row["binding_dg_per_possible_pair"]) for row in rows]
    )
    structure_ranks = average_ranks(
        [float(row["mapped_bp_distance_norm"]) for row in rows]
    )
    shape_ranks = average_ranks(
        [float(row["shape_distance_norm"]) for row in rows]
    )
    for index, row in enumerate(rows):
        row["sequence_rank"] = sequence_ranks[index]
        row["energy_rank"] = energy_ranks[index]
        row["structure_rank"] = structure_ranks[index]
        row["shape_rank"] = shape_ranks[index]
        row["composite_rank"] = composite_ranks[index]
        row["shape_composite_rank"] = shape_composite_ranks[index]

    reference_rank = next(
        float(row["composite_rank"])
        for row in rows
        if row["variant_class"] == "reference"
    )
    reference_shape_rank = next(
        float(row["shape_composite_rank"])
        for row in rows
        if row["variant_class"] == "reference"
    )
    for row in rows:
        row["clean_reference_composite_rank"] = reference_rank
        row["clean_reference_shape_composite_rank"] = reference_shape_rank


def analyze_reference(
    record_id: str,
    source_file: Path,
    transcript: str,
    reference_start: int,
    window_size: int,
    lengths: list[int],
    offsets: list[int],
    internal_edit_mode: str,
    shape_level: int,
    vienna_version: str,
    temperature_c: str | float,
) -> list[dict[str, object]]:
    reference_sequence = transcript[reference_start : reference_start + window_size]
    partner_sequence = reverse_complement(reference_sequence)
    reference_cofold = RNA.fold_compound(reference_sequence + "&" + partner_sequence)
    reference_structure, _ = reference_cofold.mfe()
    reference_pairs = interstrand_pairs(
        str(reference_structure), len(reference_sequence), len(partner_sequence)
    )
    reference_shape = shape_string(str(reference_structure), shape_level)
    partner_mfe = fold_mfe(partner_sequence)

    candidates = make_candidates(
        transcript,
        reference_start,
        window_size,
        lengths,
        offsets,
        internal_edit_mode,
    )
    reference_id = f"{record_id}:{reference_start + 1}-{reference_start + window_size}"
    rows: list[dict[str, object]] = []
    for candidate in candidates:
        edit_distance = levenshtein_distance(reference_sequence, candidate.sequence)
        metrics = cofold_metrics(
            candidate.sequence,
            partner_sequence,
            reference_pairs,
            candidate.position_map,
            reference_shape,
            shape_level,
            partner_mfe=partner_mfe,
        )
        rows.append(
            {
                "record_id": record_id,
                "source_file": source_file.name,
                "reference_id": reference_id,
                "candidate_id": candidate.candidate_id,
                "variant_class": candidate.variant_class,
                "description": candidate.description,
                "reference_start_1based": reference_start + 1,
                "reference_end_1based": reference_start + window_size,
                "candidate_source_start_1based": candidate.source_start + 1,
                "candidate_source_end_1based": candidate.source_end,
                "offset_nt": candidate.offset_nt,
                "anchor": candidate.anchor,
                "reference_sequence": reference_sequence,
                "candidate_sequence": candidate.sequence,
                "partner_sequence": partner_sequence,
                "candidate_length": len(candidate.sequence),
                "cofold_cut_point_after_nt": len(candidate.sequence),
                "candidate_to_reference_map_json": json.dumps(candidate.position_map),
                "edit_distance": edit_distance,
                "edit_distance_norm": normalized_distance(
                    edit_distance, len(reference_sequence), len(candidate.sequence)
                ),
                "viennarna_version": vienna_version,
                "thermodynamic_temperature_c": temperature_c,
                "thermodynamic_model": "ViennaRNA default RNA parameters",
                **{
                    **metrics,
                    "mapped_candidate_pairs": json.dumps(metrics["mapped_candidate_pairs"]),
                },
            }
        )

    score_reference_group(rows)
    return rows


CSV_FIELDS = [
    "record_id",
    "source_file",
    "reference_id",
    "candidate_id",
    "variant_class",
    "description",
    "reference_start_1based",
    "reference_end_1based",
    "candidate_source_start_1based",
    "candidate_source_end_1based",
    "offset_nt",
    "anchor",
    "reference_sequence",
    "candidate_sequence",
    "partner_sequence",
    "candidate_length",
    "cofold_cut_point_after_nt",
    "candidate_to_reference_map_json",
    "edit_distance",
    "edit_distance_norm",
    "mfe_candidate",
    "mfe_partner",
    "mfe_complex",
    "binding_dg_mfe",
    "binding_dg_per_possible_pair",
    "structure_cofold",
    "interstrand_pair_count",
    "longest_interstrand_run",
    "mapped_candidate_pairs",
    "mapped_bp_distance",
    "mapped_bp_distance_norm",
    "shape_level",
    "reference_shape",
    "candidate_shape",
    "shape_distance",
    "shape_distance_norm",
    "shape_exact_match",
    "sequence_goodness_percentile",
    "energy_goodness_percentile",
    "structure_goodness_percentile",
    "shape_goodness_percentile",
    "sequence_rank",
    "energy_rank",
    "structure_rank",
    "shape_rank",
    "composite_score",
    "composite_rank",
    "clean_reference_composite_rank",
    "viennarna_version",
    "thermodynamic_temperature_c",
    "thermodynamic_model",
    "shape_composite_score",
    "shape_composite_rank",
    "clean_reference_shape_composite_rank",
]


def normalized_rna_sequence(sequence: str) -> str:
    return "".join(str(sequence).upper().split()).replace("T", "U")


def safe_version() -> str:
    return str(
        getattr(RNA, "__version__", None)
        or getattr(RNA, "VERSION", None)
        or "unknown"
    )


def default_temperature() -> str | float:
    try:
        return float(RNA.md().temperature)
    except (AttributeError, TypeError, ValueError):
        return "unknown"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate length-, offset-, and edit-perturbed RNA fragments; "
            "compare sequence distance, cofold energy, and structure distance."
        )
    )
    parser.add_argument(
        "--fasta",
        nargs="+",
        type=Path,
        default=None,
        help="FASTA file(s). Defaults to the repository CDH1 and VIM inputs.",
    )
    parser.add_argument("--window-size", type=int, default=20)
    parser.add_argument(
        "--reference-start",
        type=int,
        default=None,
        help=(
            "Analyze exactly one reference window per FASTA record, beginning at "
            "this 1-based position within that record. Useful for buffered fragments."
        ),
    )
    parser.add_argument(
        "--step",
        type=int,
        default=20,
        help="Stride between clean reference windows (default: window size).",
    )
    parser.add_argument(
        "--lengths",
        nargs="+",
        type=int,
        default=[18, 19, 21],
        help="Additional natural fragment lengths, anchored at the reference 5' end.",
    )
    parser.add_argument(
        "--offsets",
        nargs="+",
        type=int,
        default=[-2, -1, 1, 2],
        help="Start offsets in nt for same-length source windows.",
    )
    parser.add_argument(
        "--internal-edits",
        choices=("none", "center", "all"),
        default="center",
        help=(
            "Synthetic internal edits: none; deterministic edits at the center "
            "(default); or every possible single-base insertion/deletion/substitution."
        ),
    )
    parser.add_argument(
        "--shape-level",
        type=int,
        choices=range(0, 6),
        default=3,
        help="ViennaRNA abstract-shape level (0-5; default 3).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output CSV path (defaults to data/exports/messy_fragments_k{window-size}.csv).",
    )
    args = parser.parse_args()
    if args.window_size <= 0:
        parser.error("--window-size must be positive")
    if args.step <= 0:
        parser.error("--step must be positive")
    if args.reference_start is not None and args.reference_start <= 0:
        parser.error("--reference-start must be a positive 1-based position")
    return args


def run() -> int:
    args = parse_args()
    if not hasattr(RNA, "abstract_shapes"):
        raise RuntimeError(
            "RNA.abstract_shapes() is required for the shape comparison. "
            "Install ViennaRNA 2.4.18 or newer."
        )

    fasta_files = args.fasta or list(DEFAULT_FASTA_FILES)
    missing = [path for path in fasta_files if not path.is_file()]
    if missing:
        raise FileNotFoundError("FASTA file(s) not found: " + ", ".join(map(str, missing)))

    records: list[tuple[str, Path, str]] = []
    for fasta_path in fasta_files:
        for record in SeqIO.parse(str(fasta_path), "fasta"):
            transcript = normalized_rna_sequence(record.seq)
            if transcript:
                records.append((record.id, fasta_path, transcript))
    if args.reference_start is None and not any(
        len(sequence) >= args.window_size for _, _, sequence in records
    ):
        raise ValueError(
            f"No reference windows of length {args.window_size} were found in the input FASTA file(s)."
        )
    if args.reference_start is not None:
        start_zero_based = args.reference_start - 1
        too_short = [
            record_id
            for record_id, _, sequence in records
            if start_zero_based + args.window_size > len(sequence)
        ]
        if too_short:
            raise ValueError(
                "The requested --reference-start plus --window-size exceeds the length of "
                "these FASTA records: "
                + ", ".join(too_short)
            )

    vienna_version = safe_version()
    temperature_c = default_temperature()
    output_path = args.output or (
        SCRIPT_DIR
        / "data"
        / "exports"
        / f"messy_fragments_k{args.window_size}.csv"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    reference_count = 0
    total_candidates = 0
    with output_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=CSV_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for record_id, fasta_path, transcript in records:
            starts = (
                [args.reference_start - 1]
                if args.reference_start is not None
                else range(0, len(transcript) - args.window_size + 1, args.step)
            )
            for start in starts:
                reference_count += 1
                rows = analyze_reference(
                    record_id=record_id,
                    source_file=fasta_path,
                    transcript=transcript,
                    reference_start=start,
                    window_size=args.window_size,
                    lengths=args.lengths,
                    offsets=args.offsets,
                    internal_edit_mode=args.internal_edits,
                    shape_level=args.shape_level,
                    vienna_version=vienna_version,
                    temperature_c=temperature_c,
                )
                total_candidates += len(rows)
                writer.writerows(rows)

    print(f"ViennaRNA version: {vienna_version}")
    print(f"Reference windows: {reference_count}")
    print(f"Candidate comparisons: {total_candidates}")
    print(f"CSV written to: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
