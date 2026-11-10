#!/usr/bin/env python3
"""Plot and summarize messy-fragment comparison CSVs.

Bars show the mean for each variant class, whiskers show sample standard
deviation when there is more than one candidate, and dots show individual
candidates. These are descriptive summaries of designed variants, not
independent biological replicates, so the script does not calculate p-values.

Example:
    python plot_messy_fragments.py \
        --input data/exports/my_fragments_with_shifts.csv
"""

from __future__ import annotations

import argparse
import csv
import os
import statistics
import tempfile
from collections import defaultdict
from pathlib import Path

# The shared environment may not allow writes to the user's home cache.
_CACHE_DIR = Path(tempfile.gettempdir()) / "messy-fragment-plot-cache"
_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE_DIR / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE_DIR / "xdg"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
VARIANT_ORDER = ("reference", "length", "offset", "substitution", "insertion", "deletion")
VARIANT_LABELS = {
    "reference": "Clean",
    "length": "Length",
    "offset": "Shift",
    "substitution": "Substitution",
    "insertion": "Insertion",
    "deletion": "Deletion",
}

SUMMARY_METRICS = {
    "binding_dg_per_possible_pair": "MFE binding ΔG per possible pair (kcal/mol)",
    "interstrand_pair_fraction": "Interstrand pair fraction",
    "continuous_run_fraction": "Longest continuous pair run fraction",
    "edit_distance_norm": "Normalized sequence edit distance",
    "mapped_bp_distance_norm": "Normalized base-pair distance",
    "shape_distance_norm": "Normalized abstract-shape distance",
    "composite_score": "Base-pair composite score",
    "shape_composite_score": "Abstract-shape composite score",
}
DERIVED_METRICS = {"interstrand_pair_fraction", "continuous_run_fraction"}


def read_candidates(path: Path) -> list[dict[str, str | float | int]]:
    with path.open("r", newline="", encoding="utf-8") as input_file:
        reader = csv.DictReader(input_file)
        required = {
            "record_id",
            "reference_id",
            "variant_class",
            "candidate_id",
            "candidate_length",
            "partner_sequence",
            "interstrand_pair_count",
            "longest_interstrand_run",
            *(set(SUMMARY_METRICS).difference(DERIVED_METRICS)),
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "Input CSV is missing required columns: " + ", ".join(sorted(missing))
            )

        candidates: list[dict[str, str | float | int]] = []
        for row in reader:
            max_possible_pairs = min(
                int(row["candidate_length"]), len(row["partner_sequence"])
            )
            if max_possible_pairs <= 0:
                raise ValueError(f"Invalid strand lengths in candidate {row['candidate_id']!r}")
            candidate: dict[str, str | float | int] = dict(row)
            candidate["candidate_length"] = int(row["candidate_length"])
            candidate["interstrand_pair_count"] = int(row["interstrand_pair_count"])
            candidate["longest_interstrand_run"] = int(row["longest_interstrand_run"])
            candidate["interstrand_pair_fraction"] = (
                candidate["interstrand_pair_count"] / max_possible_pairs
            )
            candidate["continuous_run_fraction"] = (
                candidate["longest_interstrand_run"] / max_possible_pairs
            )
            for metric in SUMMARY_METRICS:
                if metric not in DERIVED_METRICS:
                    candidate[metric] = float(row[metric])
            candidates.append(candidate)

    if not candidates:
        raise ValueError(f"No candidate rows found in {path}")
    return candidates


def grouped_candidates(
    candidates: list[dict[str, str | float | int]],
) -> dict[str, list[dict[str, str | float | int]]]:
    groups: dict[str, list[dict[str, str | float | int]]] = defaultdict(list)
    for candidate in candidates:
        groups[str(candidate["record_id"])].append(candidate)
    return dict(groups)


def variant_classes(rows: list[dict[str, str | float | int]]) -> list[str]:
    present = {str(row["variant_class"]) for row in rows}
    return [variant for variant in VARIANT_ORDER if variant in present]


def descriptive_summary(
    groups: dict[str, list[dict[str, str | float | int]]],
) -> list[dict[str, object]]:
    summary: list[dict[str, object]] = []
    for record_id, rows in groups.items():
        reference_window_counts = len({str(row["reference_id"]) for row in rows})
        for variant_class in variant_classes(rows):
            variant_rows = [row for row in rows if row["variant_class"] == variant_class]
            for metric, description in SUMMARY_METRICS.items():
                values = [float(row[metric]) for row in variant_rows]
                summary.append(
                    {
                        "record_id": record_id,
                        "variant_class": variant_class,
                        "metric": metric,
                        "description": description,
                        "candidate_count": len(values),
                        "reference_window_count": reference_window_counts,
                        "mean": statistics.mean(values),
                        "sample_sd": statistics.stdev(values) if len(values) > 1 else "",
                        "median": statistics.median(values),
                        "minimum": min(values),
                        "maximum": max(values),
                    }
                )
    return summary


def write_summary(summary: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "record_id",
        "variant_class",
        "metric",
        "description",
        "candidate_count",
        "reference_window_count",
        "mean",
        "sample_sd",
        "median",
        "minimum",
        "maximum",
    ]
    with path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(output_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary)


def jitter_offsets(count: int, span: float = 0.16) -> list[float]:
    if count <= 1:
        return [0.0] * count
    return [span * (index / (count - 1) - 0.5) for index in range(count)]


def mean_sd(rows: list[dict[str, str | float | int]], metric: str) -> tuple[float, float]:
    values = [float(row[metric]) for row in rows]
    return statistics.mean(values), (statistics.stdev(values) if len(values) > 1 else 0.0)


def tick_labels(
    rows: list[dict[str, str | float | int]], classes: list[str]
) -> list[str]:
    labels = []
    for variant_class in classes:
        count = sum(row["variant_class"] == variant_class for row in rows)
        labels.append(f"{VARIANT_LABELS.get(variant_class, variant_class)}\n(n={count})")
    return labels


def draw_single_metric(
    axis: plt.Axes,
    rows: list[dict[str, str | float | int]],
    classes: list[str],
    metric: str,
    title: str,
    ylabel: str,
    color: str,
    unit_interval: bool = False,
) -> None:
    x_positions = list(range(len(classes)))
    for x_position, variant_class in zip(x_positions, classes):
        group = [row for row in rows if row["variant_class"] == variant_class]
        mean, sd = mean_sd(group, metric)
        axis.bar(
            x_position,
            mean,
            width=0.58,
            color=color,
            alpha=0.72,
            edgecolor="black",
            linewidth=0.5,
            yerr=sd if len(group) > 1 else None,
            capsize=3,
            zorder=2,
        )
        axis.scatter(
            [x_position + offset for offset in jitter_offsets(len(group))],
            [float(row[metric]) for row in group],
            color="black",
            s=22,
            alpha=0.8,
            zorder=3,
        )

    axis.set_xticks(x_positions, tick_labels(rows, classes), rotation=0)
    axis.set_title(title)
    axis.set_ylabel(ylabel)
    axis.grid(axis="y", alpha=0.25, zorder=0)
    if unit_interval:
        axis.set_ylim(-0.05, 1.05)
    else:
        values = [float(row[metric]) for row in rows]
        minimum, maximum = min(values), max(values)
        span = max(maximum - minimum, 0.1)
        axis.set_ylim(minimum - 0.22 * span, max(0.0, maximum + 0.12 * span))
        axis.axhline(0, color="black", linewidth=0.7, alpha=0.6)


def draw_paired_metrics(
    axis: plt.Axes,
    rows: list[dict[str, str | float | int]],
    classes: list[str],
    first_metric: str,
    second_metric: str,
    title: str,
    ylabel: str,
) -> None:
    x_positions = list(range(len(classes)))
    bar_width = 0.32
    metric_specs = (
        (first_metric, "#4C78A8", -bar_width / 2),
        (second_metric, "#F58518", bar_width / 2),
    )
    for metric, color, offset in metric_specs:
        for x_position, variant_class in zip(x_positions, classes):
            group = [row for row in rows if row["variant_class"] == variant_class]
            mean, sd = mean_sd(group, metric)
            bar_position = x_position + offset
            axis.bar(
                bar_position,
                mean,
                width=bar_width,
                color=color,
                alpha=0.76,
                edgecolor="black",
                linewidth=0.5,
                yerr=sd if len(group) > 1 else None,
                capsize=3,
                zorder=2,
            )
            axis.scatter(
                [bar_position + jitter for jitter in jitter_offsets(len(group), span=0.1)],
                [float(row[metric]) for row in group],
                color="black",
                s=18,
                alpha=0.75,
                zorder=3,
            )

    axis.set_xticks(x_positions, tick_labels(rows, classes), rotation=0)
    axis.set_title(title)
    axis.set_ylabel(ylabel)
    axis.set_ylim(-0.05, 1.05)
    axis.grid(axis="y", alpha=0.25, zorder=0)


def make_figure(
    groups: dict[str, list[dict[str, str | float | int]]],
    path: Path,
    dpi: int,
) -> None:
    record_ids = list(groups)
    figure, axes = plt.subplots(
        nrows=len(record_ids),
        ncols=3,
        figsize=(16, max(4.4 * len(record_ids), 6.0)),
        squeeze=False,
    )

    for row_index, record_id in enumerate(record_ids):
        rows = groups[record_id]
        classes = variant_classes(rows)
        fragment_label = record_id.split("_")[0]
        draw_single_metric(
            axes[row_index, 0],
            rows,
            classes,
            "binding_dg_per_possible_pair",
            f"{fragment_label}: MFE binding energy per possible pair",
            "kcal/mol (more negative is more favorable)",
            "#4C78A8",
        )
        draw_paired_metrics(
            axes[row_index, 1],
            rows,
            classes,
            "interstrand_pair_fraction",
            "continuous_run_fraction",
            "Pairing retained (blue: all pairs; orange: continuous run)",
            "Fraction of maximum possible base pairs",
        )
        draw_paired_metrics(
            axes[row_index, 2],
            rows,
            classes,
            "composite_score",
            "shape_composite_score",
            "Combined scores (blue: base-pair; orange: shape)",
            "Percentile composite score (higher is better)",
        )

    figure.suptitle(
        "Messy-fragment variant summary\n"
        "Bars = mean; whiskers = ±1 sample SD; dots = individual candidates. "
        "Descriptive only; variants are not independent biological replicates.",
        fontsize=12,
    )
    figure.tight_layout(rect=(0.04, 0.04, 1, 0.93))
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot candidate-level metrics and descriptive statistics from a messy-fragment CSV."
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Result CSV produced by messy_fragments.py.",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=None,
        help="Output path prefix; defaults to the input filename without its suffix.",
    )
    parser.add_argument("--dpi", type=int, default=200)
    args = parser.parse_args()
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    return args


def main() -> int:
    args = parse_args()
    if not args.input.is_file():
        raise FileNotFoundError(f"Result CSV not found: {args.input}")

    prefix = args.output_prefix or args.input.with_suffix("")
    figure_path = Path(f"{prefix}_visualization.png")
    summary_path = Path(f"{prefix}_summary.csv")

    candidates = read_candidates(args.input)
    groups = grouped_candidates(candidates)
    summary = descriptive_summary(groups)
    make_figure(groups, figure_path, args.dpi)
    write_summary(summary, summary_path)

    print(f"Candidate rows summarized: {len(candidates)}")
    print(f"Reference records summarized: {len(groups)}")
    print(f"Figure written to: {figure_path}")
    print(f"Descriptive statistics written to: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
