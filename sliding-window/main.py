import RNA
import numpy as np
from dataclasses import dataclass, field
from Bio import SeqIO

from strategy.mfe import get_mfe
from strategy.cofold import cofold_energies
from strategy.accessibility import mean_unpaired_probability, opening_cost

# ── sequence sources ────────────────────────────────────────────────────────
CDH1_PATH = './data/CDH1__NM_004360.5.fasta'
VIM_PATH  = './data/VIM__NM_003380.5.fasta'

# ── window sizes to evaluate ────────────────────────────────────────────────
WINDOW_SIZES = [42, 60, 100, 200, 300]


# ── per-strand sliding window ────────────────────────────────────────────────
def sliding_window(sequence: str, k: int, step: int = 1):
    """
    Tile a nucleotide sequence into overlapping windows of length k.

    Converts DNA to RNA (T → U) before windowing so downstream ViennaRNA
    calls always receive RNA input.

    Args:
        sequence: Raw nucleotide string (DNA or RNA, any case).
        k:        Window length in nucleotides.
        step:     Stride between consecutive window starts.
                  step=1 is exhaustive; step=k gives non-overlapping tiles.

    Yields:
        (start, window_seq) — 0-based start index and the RNA window string.
    """
    seq = sequence.upper().replace("T", "U")
    n = len(seq)
    for s in range(0, n - k + 1, step):
        yield s, seq[s:s + k]


# ── single-strand analysis ───────────────────────────────────────────────────
@dataclass
class StrandProfile:
    start: int
    seq: str
    mfe: float
    accessibility: float   # mean P(unpaired)
    opening_cost: float    # kcal/mol estimated cost to open structure


def analyze_strand(start: int, seq: str) -> StrandProfile:
    """
    Run the single-strand arm of the workflow (step 1 of 4).

    Computes:
      - MFE secondary structure via RNAfold
      - Mean P(unpaired) across all nucleotides via partition function
      - Estimated opening cost (kcal/mol) — energy needed to melt local structure

    Args:
        start: 0-based position in the source sequence (for traceability).
        seq:   RNA window string (U-alphabet).

    Returns:
        StrandProfile populated with all three thermodynamic descriptors.
    """
    return StrandProfile(
        start=start,
        seq=seq,
        mfe=get_mfe(seq),
        accessibility=mean_unpaired_probability(seq),
        opening_cost=opening_cost(seq),
    )


# ── pairwise duplex analysis ─────────────────────────────────────────────────
@dataclass
class DuplexResult:
    cdh1_start: int
    vim_start: int
    cdh1_seq: str
    vim_seq: str
    mfe_cdh1: float
    mfe_vim: float
    mfe_cofold: float
    binding_dg: float           # ΔG_cofold − (ΔG1 + ΔG2)
    acc_cdh1: float
    acc_vim: float
    structure_cofold: str
    score: float = field(init=False)  # composite ranking score

    def __post_init__(self):
        """
        Compute composite ranking score after field initialisation.

        score = −ΔG_binding + mean(acc_CDH1, acc_VIM)

        Both terms are positive-good:
          −ΔG_binding  is large when duplex formation releases a lot of energy.
          mean_acc     is large when both binding regions are already open.
        Higher score → pair is both energetically and structurally favourable.
        """
        self.score = -self.binding_dg + (self.acc_cdh1 + self.acc_vim) / 2


def analyze_pair(
    cdh1_start: int, cdh1_seq: str,
    vim_start: int,  vim_seq: str,
) -> DuplexResult:
    """
    Run the full 4-step duplex evaluation for one CDH1/VIM window pair.

    Steps:
      1. RNAfold on each strand individually  (via cofold_energies)
      2. RNAcofold on the joined duplex       (via cofold_energies)
      3. Compute ΔG_binding = ΔG_cofold − (ΔG_CDH1 + ΔG_VIM)
      4. Partition-function accessibility for each strand

    Args:
        cdh1_start: 0-based start of the CDH1 window in the full transcript.
        cdh1_seq:   RNA window from CDH1.
        vim_start:  0-based start of the VIM window in the full transcript.
        vim_seq:    RNA window from VIM.

    Returns:
        DuplexResult with all energy terms, accessibilities, and composite score.
    """
    energies = cofold_energies(cdh1_seq, vim_seq)
    acc1 = mean_unpaired_probability(cdh1_seq)
    acc2 = mean_unpaired_probability(vim_seq)

    return DuplexResult(
        cdh1_start=cdh1_start,
        vim_start=vim_start,
        cdh1_seq=cdh1_seq,
        vim_seq=vim_seq,
        mfe_cdh1=energies["mfe1"],
        mfe_vim=energies["mfe2"],
        mfe_cofold=energies["mfe_cofold"],
        binding_dg=energies["binding_dg"],
        acc_cdh1=acc1,
        acc_vim=acc2,
        structure_cofold=energies["structure_cofold"],
    )


# ── full pairwise scan for one window size ───────────────────────────────────
def scan_pairwise(
    cdh1_seq: str,
    vim_seq: str,
    k: int,
    step: int = 10,
    top_n: int = 10,
) -> list[DuplexResult]:
    """
    Exhaustively compare every CDH1 window against every VIM window of length k.

    Complexity is O(N_CDH1 × N_VIM) cofold calls, where N = (len − k) / step + 1.
    Increasing step reduces cost linearly at the expense of positional resolution.

    Args:
        cdh1_seq: Full CDH1 transcript (RNA, U-alphabet).
        vim_seq:  Full VIM transcript (RNA, U-alphabet).
        k:        Window length in nucleotides.
        step:     Stride for both sequence windows.
                  step=k → non-overlapping tiles (fastest).
                  step=1 → single-nucleotide resolution (slowest).
        top_n:    Number of top-ranked DuplexResults to return.

    Returns:
        List of up to top_n DuplexResults sorted by composite score (descending).
    """
    results: list[DuplexResult] = []

    cdh1_windows = list(sliding_window(cdh1_seq, k, step))
    vim_windows  = list(sliding_window(vim_seq,  k, step))

    total = len(cdh1_windows) * len(vim_windows)
    print(f"  k={k}: {len(cdh1_windows)} × {len(vim_windows)} = {total} pairs")

    for i, (cs, cw) in enumerate(cdh1_windows):
        for vs, vw in vim_windows:
            results.append(analyze_pair(cs, cw, vs, vw))
        if (i + 1) % 10 == 0:
            print(f"    CDH1 window {i+1}/{len(cdh1_windows)} done …")

    results.sort(key=lambda r: r.score, reverse=True)
    return results[:top_n]


# ── pretty-print helper ───────────────────────────────────────────────────────
def print_result(r: DuplexResult, rank: int):
    """
    Print a single DuplexResult in a human-readable summary block.

    Binding strength labels follow the rough thresholds from the literature:
      STRONG   ΔG_binding < −10 kcal/mol
      moderate −10 ≤ ΔG_binding < −5 kcal/mol
      weak     ΔG_binding ≥ −5 kcal/mol

    Args:
        r:    DuplexResult to display.
        rank: 1-based rank within the current window-size results.
    """
    strength = (
        "STRONG"   if r.binding_dg < -10 else
        "moderate" if r.binding_dg < -5  else
        "weak"
    )
    print(
        f"  #{rank:>2}  CDH1[{r.cdh1_start}]  VIM[{r.vim_start}]\n"
        f"        ΔG₁={r.mfe_cdh1:+.2f}  ΔG₂={r.mfe_vim:+.2f}  "
        f"ΔG_cofold={r.mfe_cofold:+.2f}  ΔG_binding={r.binding_dg:+.2f} kcal/mol  [{strength}]\n"
        f"        acc_CDH1={r.acc_cdh1:.3f}  acc_VIM={r.acc_vim:.3f}  "
        f"score={r.score:.3f}\n"
        f"        structure: {r.structure_cofold}\n"
    )


# ── entry point ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    cdh1_record = next(SeqIO.parse(CDH1_PATH, "fasta"))
    vim_record  = next(SeqIO.parse(VIM_PATH,  "fasta"))

    cdh1_seq = str(cdh1_record.seq).upper().replace("T", "U")
    vim_seq  = str(vim_record.seq).upper().replace("T", "U")

    print(f"CDH1: {len(cdh1_seq)} nt   VIM: {len(vim_seq)} nt\n")

    for k in WINDOW_SIZES:
        print(f"=== Window size k={k} ===")
        # step=k gives non-overlapping tiles; lower values are more thorough but slower
        top = scan_pairwise(cdh1_seq, vim_seq, k=k, step=k, top_n=5)
        print(f"  Top {len(top)} duplex candidates:")
        for rank, r in enumerate(top, 1):
            print_result(r, rank)
