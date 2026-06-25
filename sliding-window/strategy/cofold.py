import RNA


def cofold_energies(seq1: str, seq2: str) -> dict:
    """
    Run RNAfold on each strand and RNAcofold on the duplex.

    Returns individual MFEs, cofold MFE, and the net binding ΔG:
        ΔG_binding = ΔG_cofold − (ΔG_strand1 + ΔG_strand2)

    More negative binding_dg → duplex formation more favored.
    Rough thresholds: < −10 strong, −5 to −10 moderate, > −5 weak.
    """
    fc1 = RNA.fold_compound(seq1)
    _, mfe1 = fc1.mfe()

    fc2 = RNA.fold_compound(seq2)
    _, mfe2 = fc2.mfe()

    fc_co = RNA.fold_compound(seq1 + "&" + seq2)
    structure_co, mfe_co = fc_co.mfe()

    binding_dg = mfe_co - (mfe1 + mfe2)

    return {
        "mfe1": mfe1,
        "mfe2": mfe2,
        "mfe_cofold": mfe_co,
        "binding_dg": binding_dg,
        "structure_cofold": structure_co,
    }
