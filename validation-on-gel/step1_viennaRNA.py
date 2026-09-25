from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Callable

@dataclass(frozen=True)
class FoldResult:
	"""The fold result for one construct."""

	name: str
	sequences: tuple[str, ...]
	structure: str
	energy: float


def validate_sequence(sequence: str) -> str:
	"""Normalize and validate one RNA sequence."""
	normalized = sequence.strip().upper().replace("T", "U")
	if not normalized:
		raise ValueError("RNA sequences must not be empty")
	invalid = set(normalized) - set("ACGU")
	if invalid:
		raise ValueError(
			f"Invalid RNA bases {sorted(invalid)} in sequence {sequence!r}"
		)
	return normalized


def fold_construct(
	name: str,
	sequences: list[str],
	fold: Callable[[str], tuple[str, float]],
) -> FoldResult:
	"""Fold a construct, preserving strand boundaries for RNA cofolding."""
	validated = tuple(validate_sequence(sequence) for sequence in sequences)
	input_sequence = "&".join(validated)
	structure, energy = fold(input_sequence)
	return FoldResult(name, validated, structure, float(energy))


def fold_all_constructs(constructs) -> list[FoldResult]:
	"""Fold all configured constructs with ViennaRNA."""
	try:
		import RNA
	except ImportError as error:
		raise RuntimeError(
			"ViennaRNA is required. Install it with `pip install ViennaRNA`."
		) from error

	return [
		fold_construct(name, sequences, RNA.cofold)
		for name, sequences in constructs.items()
	]


def write_results(results: list[FoldResult], output_path: Path) -> None:
	"""Write a human-readable, dated report for downstream modelling."""
	output_path.parent.mkdir(parents=True, exist_ok=True)
	with output_path.open("w", encoding="utf-8") as output:
		output.write(f"Validation-on-gel RNA folds ({date.today().isoformat()})\n")
		output.write("=" * 68 + "\n\n")
		for result in results:
			output.write(f"[{result.name}]\n")
			output.write(f"sequences: {' & '.join(result.sequences)}\n")
			output.write(f"lengths: {', '.join(str(len(sequence)) for sequence in result.sequences)}\n")
			output.write(f"structure: {result.structure}\n")
			output.write(f"free_energy_kcal_mol: {result.energy:.2f}\n\n")
