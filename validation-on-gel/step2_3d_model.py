"""Prepare and optionally run 3-D RNA modelling jobs.

The default output is a job directory that can be uploaded to RNAComposer,
3dRNA, or consumed by a local modelling command.  A local command can use
``{fasta}``, ``{dotbracket}``, and ``{output}`` placeholders.
"""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import webbrowser
from dataclasses import asdict, dataclass
from pathlib import Path


RNA_COMPOSER_URL = "https://rnacomposer.cs.put.poznan.pl/"


@dataclass(frozen=True)
class Fold:

    name: str
    sequences: tuple[str, ...]
    structure: str
    energy: float


def read_fold_report(path: Path) -> list[Fold]:
    """Read the report written by ``step1_viennaRNA.py``."""
    folds: list[Fold] = []
    current: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("[") and line.endswith("]"):
            if current:
                folds.append(_make_fold(current))
            current = {"name": line[1:-1]}
        elif ": " in line and current:
            key, value = line.split(": ", 1)
            current[key] = value
    if current:
        folds.append(_make_fold(current))
    if not folds:
        raise ValueError(f"No fold records found in {path}")
    return folds


def _make_fold(values: dict[str, str]) -> Fold:
    try:
        return Fold(
            values["name"],
            tuple(values["sequences"].split(" & ")),
            values["structure"],
            float(values["free_energy_kcal_mol"]),
        )
    except KeyError as error:
        raise ValueError(f"Incomplete fold record: missing {error.args[0]}") from error


def prepare_jobs(report: Path, output_dir: Path) -> list[dict[str, object]]:
    """Write one FASTA and one dot-bracket file per construct."""
    output_dir.mkdir(parents=True, exist_ok=True)
    jobs: list[dict[str, object]] = []
    for fold in read_fold_report(report):
        prefix = output_dir / fold.name
        fasta = prefix.with_suffix(".fasta")
        dotbracket = prefix.with_suffix(".db")
        fasta.write_text(
            ">" + fold.name + "\n" + "&".join(fold.sequences) + "\n",
            encoding="utf-8",
        )
        dotbracket.write_text(fold.structure + "\n", encoding="utf-8")
        jobs.append(
            {
                "name": fold.name,
                "fasta": str(fasta),
                "dotbracket": str(dotbracket),
                "output": str(prefix.with_suffix(".pdb")),
                "sequences": fold.sequences,
                "energy": fold.energy,
            }
        )
    (output_dir / "manifest.json").write_text(json.dumps(jobs, indent=2) + "\n", encoding="utf-8")
    return jobs


def open_rna_composer() -> None:
    """Open the hosted RNAComposer interactive modeller."""
    webbrowser.open(RNA_COMPOSER_URL)


def run_model_command(jobs: list[dict[str, object]], command: str) -> None:
    """Run a local modeller once per job using placeholder substitution."""
    for job in jobs:
        formatted = command.format(**job)
        subprocess.run(shlex.split(formatted), check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", type=Path, default=Path("date_results.txt"))
    parser.add_argument("-o", "--output", type=Path, default=Path("models"))
    parser.add_argument(
        "--command",
        help="optional local modeller command with {fasta}, {dotbracket}, and {output}",
    )
    parser.add_argument(
        "--open-rna-composer",
        action="store_true",
        help="open RNAComposer after preparing the input files",
    )
    args = parser.parse_args()
    jobs = prepare_jobs(args.input, args.output)
    if args.open_rna_composer:
        open_rna_composer()
    if args.command:
        run_model_command(jobs, args.command)
    print(f"Prepared {len(jobs)} modelling job(s) in {args.output}")


if __name__ == "__main__":
    main()