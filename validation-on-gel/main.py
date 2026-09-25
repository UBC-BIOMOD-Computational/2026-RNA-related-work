


"""Fold the validation-on-gel RNA constructs and export their structures."""

from __future__ import annotations

import argparse
from pathlib import Path

from step1_viennaRNA import fold_all_constructs, write_results
from step2_3d_model import open_rna_composer, prepare_jobs, run_model_command
from step3_oxrna import relax_models

# Sequences are kept in RNA alphabet because ViennaRNA expects U rather than T.
vim = "CGAGGAGGAGAUGCGGGAGC"
cdh1 = "CCGCCCUGCCAAUCCCGAUG"
helper = "CAUCGGGAUUGGGAAAGGCAGGGCGGGCUCCCGCAUCUCCUCCUCG"

X = [vim, helper]
Y = [vim, cdh1, helper]
Z = [cdh1, helper]



def main() -> None:
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument(
		"-o",
		"--output",
		type=Path,
		default=Path(__file__).with_name("date_results.txt"),
		help="report path (default: validation-on-gel/date_results.txt)",
	)
	parser.add_argument(
		"--models-dir",
		type=Path,
		default=Path(__file__).with_name("models"),
		help="step-2 job and PDB directory",
	)
	parser.add_argument(
		"--model-command",
		help="step-2 command with {fasta}, {dotbracket}, and {output} placeholders",
	)
	parser.add_argument(
		"--open-rna-composer",
		action="store_true",
		help="open RNAComposer after preparing step-2 input files",
	)
	parser.add_argument(
		"--snapshots-dir",
		type=Path,
		default=Path(__file__).with_name("snapshots"),
		help="step-3 output directory",
	)
	parser.add_argument(
		"--relax-command",
		help="step-3 command with {input}, {output}, {steps}, and {seed} placeholders",
	)
	parser.add_argument("--steps", type=int, default=1000)
	parser.add_argument("--seed", type=int)
	args = parser.parse_args()
    
	constructs = {"X": X, "Y": Y, "Z": Z}
	results = fold_all_constructs(constructs)
	write_results(results, args.output)
	for result in results:
		print(f"{result.name}: {result.energy:.2f} kcal/mol -> {result.structure}")
	print(f"Wrote {args.output}")

	jobs = prepare_jobs(args.output, args.models_dir)
	print(f"Prepared {len(jobs)} modelling job(s) in {args.models_dir}")
	if args.open_rna_composer:
		open_rna_composer()
	if args.model_command:
		run_model_command(jobs, args.model_command)
	if args.relax_command:
		manifest = relax_models(
			args.models_dir,
			args.snapshots_dir,
			args.relax_command,
			args.steps,
			args.seed,
		)
		print(f"Collected {sum(len(item['snapshots']) for item in manifest)} snapshot(s)")

if __name__ == "__main__":
	main()