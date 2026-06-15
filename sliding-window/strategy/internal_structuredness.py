# Internal Structuredness Z-score

"""
Biological Rationale

Some regions appear stable merely because they are GC-rich.

A structuredness z-score asks:

Is this region more structured than expected given its composition?

Input
window_seq
n_shuffles

Output
z_score
Algorithm

Already implemented:

(native_mfe - mean(shuffled_mfes)) / std(shuffled_mfes)

using mononucleotide or dinucleotide preserving shuffles.

Interpretation

Z-score << 0	= Unusually structured
Z-score ≈ 0	= Expected
Z-score > 0 = 	Less structured than expected


Desired Direction

For accessibility:
- higher z-score preferred

For structure discovery:
- lower z-score preferred

"""


def structuredness_zscore(window_seq: str, n_shuffles: int = 50, seed: int | None = None) -> float:
    """
    A strongly negative `z` means the window is *more* stable than its base composition would predict by chance — a real candidate for "internally structured / locked up" rather than just GC-rich. Note: a plain `random.shuffle` only preserves mononucleotide composition; if you care about CpG-type biases, swap in a dinucleotide-preserving shuffle (Altschul–Erikson algorithm) instead.

    This is also the expensive step — `n_shuffles=50` multiplies your MFE-call count by 50, so run it on a coarser step (e.g. the `L/100` non-overlapping windows from the spec) rather than every single nucleotide offset.

    """

    rng = random.Random(seed)
    seq = window_seq.upper().replace("T", "U")

    native_mfe = RNA.fold_compound(seq).mfe()[1]

    shuffled_mfes = []
    chars = list(seq)
    for _ in range(n_shuffles):
        rng.shuffle(chars)
        shuffled_mfes.append(RNA.fold_compound("".join(chars)).mfe()[1])

    mu, sigma = np.mean(shuffled_mfes), np.std(shuffled_mfes)
    return (native_mfe - mu) / sigma if sigma > 0 else 0.0
