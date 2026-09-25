"""Run an oxRNA/oxDNA relaxation command and collect PDB snapshots."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
from pathlib import Path


def relax_model(
    model: Path,
    output_dir: Path,
    command: str,
    steps: int,
    seed: int | None = None,
) -> list[Path]:
    """Run the supplied simulator command and collect generated PDB files.

    The command may use ``{input}``, ``{output}``, ``{steps}``, and ``{seed}``.
    The simulator is responsible for producing PDB snapshots in its output
    directory; this wrapper makes them available in a stable per-model folder.
    """
    run_dir = output_dir / model.stem
    run_dir.mkdir(parents=True, exist_ok=True)
    values = {
        "input": str(model.resolve()),
        "output": str(run_dir.resolve()),
        "steps": steps,
        "seed": "" if seed is None else seed,
    }
    subprocess.run(shlex.split(command.format(**values)), check=True, cwd=run_dir)
    snapshots = sorted(run_dir.rglob("*.pdb"))
    if not snapshots:
        raise RuntimeError(f"Relaxation produced no PDB snapshots in {run_dir}")
    return snapshots


def relax_models(
    model_dir: Path,
    output_dir: Path,
    command: str,
    steps: int = 1000,
    seed: int | None = None,
) -> list[dict[str, object]]:
    """Relax every PDB model in ``model_dir`` and write a snapshot manifest."""
    manifest: list[dict[str, object]] = []
    for model in sorted(model_dir.glob("*.pdb")):
        snapshots = relax_model(model, output_dir, command, steps, seed)
        manifest.append({"model": str(model), "snapshots": [str(path) for path in snapshots]})
    if not manifest:
        raise ValueError(f"No PDB models found in {model_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model_dir", type=Path)
    parser.add_argument("-o", "--output", type=Path, default=Path("snapshots"))
    parser.add_argument("--command", required=True, help="simulator command template")
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--seed", type=int)
    args = parser.parse_args()
    try:
        manifest = relax_models(args.model_dir, args.output, args.command, args.steps, args.seed)
    except ValueError as error:
        raise SystemExit(str(error)) from error
    print(f"Collected {sum(len(item['snapshots']) for item in manifest)} snapshot(s)")


if __name__ == "__main__":
    main()