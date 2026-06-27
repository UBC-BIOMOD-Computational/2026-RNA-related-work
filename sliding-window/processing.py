import os
import RNA
import matplotlib.pyplot as plt
from model import DuplexResult
from strategy.cofold import cofold_energies
from strategy.accessibility import mean_unpaired_probability

# ── pretty-print helper ───────────────────────────────────────────────────────
EXPORT_DIR = './data/exports'


def _sliding_window(seq: str, k: int, step: int):
    """Yield (start, window) over a single RNA sequence."""
    n = len(seq)
    for s in range(0, n - k + 1, step):
        yield s, seq[s:s + k]


def _analyze_pair(start_a: int, seq_a: str, start_b: int, seq_b: str) -> DuplexResult:
    """
    Run the 4-step duplex evaluation for one window pair from the same sequence.

    Steps:
      1. RNAfold on each window individually   (via cofold_energies)
      2. RNAcofold on the joined duplex        (via cofold_energies)
      3. ΔG_binding = ΔG_cofold − (ΔG_a + ΔG_b)
      4. Partition-function accessibility for each window
    """
    energies = cofold_energies(seq_a, seq_b)
    return DuplexResult(
        start_a=start_a,
        start_b=start_b,
        seq_a=seq_a,
        seq_b=seq_b,
        mfe_a=energies["mfe1"],
        mfe_b=energies["mfe2"],
        mfe_cofold=energies["mfe_cofold"],
        binding_dg=energies["binding_dg"],
        acc_a=mean_unpaired_probability(seq_a),
        acc_b=mean_unpaired_probability(seq_b),
        structure_cofold=energies["structure_cofold"],
    )


def _export_histogram(results: list[DuplexResult], k: int, label: str):
    """
    Save a two-panel histogram (binding ΔG and composite score) to data/exports/.

    Left panel:  distribution of ΔG_binding across all pairs; vertical lines mark
                 the −10 and −5 kcal/mol strength thresholds.
    Right panel: distribution of composite scores.

    Args:
        results: All DuplexResults for this window size (not just top_n).
        k:       Window length — used in the filename and title.
        label:   Sequence label (e.g. "CDH1") for the filename.
    """
    os.makedirs(EXPORT_DIR, exist_ok=True)

    binding_dgs = [r.binding_dg for r in results]
    scores      = [r.score      for r in results]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    fig.suptitle(f"{label}  |  window k={k}  |  {len(results)} pairs", fontsize=12)

    ax1.hist(binding_dgs, bins=40, color="steelblue", edgecolor="white")
    ax1.axvline(-10, color="red",    linestyle="--", linewidth=1, label="−10 kcal/mol (strong)")
    ax1.axvline(-5,  color="orange", linestyle="--", linewidth=1, label="−5 kcal/mol (moderate)")
    ax1.set_xlabel("ΔG_binding (kcal/mol)")
    ax1.set_ylabel("pair count")
    ax1.set_title("Binding ΔG distribution")
    ax1.legend(fontsize=8)

    ax2.hist(scores, bins=40, color="seagreen", edgecolor="white")
    ax2.set_xlabel("composite score")
    ax2.set_ylabel("pair count")
    ax2.set_title("Score distribution")

    plt.tight_layout()
    path = os.path.join(EXPORT_DIR, f"{label}_k{k}_histogram.png")
    plt.savefig(path, dpi=150)
    plt.close(fig)
    print(f"  histogram saved → {path}")


_RC_TABLE = str.maketrans("AUCG", "UAGC")

def _reverse_complement(seq: str) -> str:
    """Return the reverse complement of an RNA sequence (A↔U, G↔C, reversed)."""
    return seq.translate(_RC_TABLE)[::-1]


def _best_rc_match(rc: str, windows: list[tuple[int, str]]) -> tuple[int, str]:
    """
    Find the window whose sequence most closely matches the target reverse complement.

    Similarity is measured by the number of matching positions (Hamming-style).
    All windows have the same length as rc, so comparison is position-wise.

    Returns the (start, seq) of the best-matching window.
    """
    best_start, best_seq, best_matches = 0, windows[0][1], -1
    for start, w in windows:
        matches = sum(a == b for a, b in zip(rc, w))
        if matches > best_matches:
            best_matches = matches
            best_start, best_seq = start, w
    return best_start, best_seq


def scan_pairwise(
    seq: str,
    k: int,
    step: int = 10,
    top_n: int = 10,
    label: str = "seq",
) -> list[DuplexResult]:
    """
    For each window of length k, find its best reverse-complement partner among
    all other windows in the same sequence, then evaluate whether that pair
    preferentially forms a duplex.

    This gives exactly N results for N windows — one candidate partner per window —
    rather than the C(N,2) combinatorial explosion of comparing every pair.

    The partner is the window whose sequence most closely matches the reverse
    complement of the query window (highest position-wise base match count).
    RNAcofold then decides whether the thermodynamics actually favour binding.

    Args:
        seq:   Full transcript, RNA (U-alphabet).
        k:     Window length in nucleotides.
        step:  Stride between window starts.
               step=k → non-overlapping tiles (fastest).
               step=1 → single-nucleotide resolution (slowest).
        top_n: Number of top-ranked DuplexResults to return.
        label: Short name used in the exported histogram filename.

    Returns:
        List of up to top_n DuplexResults sorted by composite score (descending).
        Side-effect: saves a PNG histogram to data/exports/.
    """
    windows = list(_sliding_window(seq, k, step))
    print(f"  k={k}: {len(windows)} windows → {len(windows)} candidate pairs")

    results: list[DuplexResult] = []
    for i, (sa, wa) in enumerate(windows):
        rc = _reverse_complement(wa)
        # exclude the window itself from partner search
        others = [w for w in windows if w[0] != sa]
        if not others:
            continue
        sb, wb = _best_rc_match(rc, others)
        results.append(_analyze_pair(sa, wa, sb, wb))

    _export_histogram(results, k, label)

    results.sort(key=lambda r: r.score, reverse=True)
    return results[:top_n]
