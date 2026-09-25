"""Run Rosetta FARFAR2 for the RNA jobs prepared by step 2."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
from pathlib import Path


def run_farfar2_model(
	model_dir: Path,
	output_dir: Path,
	executable: str,
	nstruct: int = 10,
	extra_args: list[str] | None = None,
) -> list[Path]:
	"""Run FARFAR2 for every step-2 FASTA/dot-bracket pair.

	FARFAR2 models the sequence and secondary structure. The step-2 PDB is
	kept alongside the inputs but is not used as a starting structure.
	"""
	if nstruct < 1:
		raise ValueError("nstruct must be at least 1")

	output_dir.mkdir(parents=True, exist_ok=True)
	models: list[Path] = []
	for fasta in sorted(model_dir.glob("*.fasta")):
		secstruct = fasta.with_suffix(".db")
		if not secstruct.exists():
			raise FileNotFoundError(f"Missing secondary structure file: {secstruct}")

		run_dir = output_dir / fasta.stem
		run_dir.mkdir(parents=True, exist_ok=True)
		silent = run_dir / f"{fasta.stem}.out"
		command = [
			*shlex.split(executable),
			"-fasta",
			str(fasta.resolve()),
			"-secstruct_file",
			str(secstruct.resolve()),
			"-nstruct",
			str(nstruct),
			"-out:path:all",
			str(run_dir.resolve()),
			"-out:file:silent",
			str(silent.resolve()),
			"-minimize_rna",
			"-dump",
		]
		command.extend(extra_args or [])
		subprocess.run(command, check=True, cwd=run_dir)
		models.extend(sorted(run_dir.glob("*.pdb")))

	if not models:
		raise RuntimeError(
			f"FARFAR2 produced no PDB files in {output_dir}; "
			"check that the Rosetta command includes -dump"
		)
	manifest = {
		"model_dir": str(model_dir.resolve()),
		"models": [str(path.resolve()) for path in models],
	}
	(output_dir / "manifest.json").write_text(
		json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
	)
	return models


def main() -> None:
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument("model_dir", type=Path)
	parser.add_argument("-o", "--output", type=Path, default=Path("farfar2"))
	parser.add_argument(
		"--executable",
		default="rna_denovo.default.linuxgccrelease",
		help="Rosetta FARFAR2/rna_denovo executable",
	)
	parser.add_argument("--nstruct", type=int, default=10)
	parser.add_argument(
		"--extra-arg",
		action="append",
		default=[],
		help="additional Rosetta argument; may be repeated",
	)
	args = parser.parse_args()
	models = run_farfar2_model(
		args.model_dir,
		args.output,
		args.executable,
		args.nstruct,
		args.extra_arg,
	)
	print(f"FARFAR2 produced {len(models)} PDB model(s) in {args.output}")


if __name__ == "__main__":
	main()