import RNA
import numpy as np


def mean_unpaired_probability(seq: str) -> float:
    """
    Accessibility score via partition function (RNAplfold-style).

    Computes the mean probability that each nucleotide is unpaired.
    Higher score → region is more open / accessible for intermolecular binding.

    Uses ViennaRNA partition function; bpp()[i][j] is P(i pairs with j).
    P(unpaired at i) = 1 − Σ_j P(i,j).
    """
    fc = RNA.fold_compound(seq)
    fc.pf()
    bpp = fc.bpp()  # bpp[i][j], 1-indexed, upper-triangle

    n = len(seq)
    paired_prob = np.zeros(n + 1)
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            paired_prob[i] += bpp[i][j]
            paired_prob[j] += bpp[i][j]

    unpaired = 1.0 - paired_prob[1:]  # drop index-0 placeholder
    return float(np.mean(unpaired))


def opening_cost(seq: str) -> float:
    """
    Approximate cost (kcal/mol) to make a sequence fully accessible.

    ΔG_open ≈ MFE × mean_paired_probability
    (structured windows have larger opening costs)
    """
    fc = RNA.fold_compound(seq)
    _, mfe = fc.mfe()
    acc = mean_unpaired_probability(seq)
    mean_paired = 1.0 - acc
    return mfe * mean_paired
