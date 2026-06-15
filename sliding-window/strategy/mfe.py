import RNA


def get_mfe(seq: str):
    """
    Compute MFE per window with ViennaRNA.


    NOTE: 
    Raw MFE is strongly correlated with GC content — a high-GC window will look "structured" even if it's not unusually so for its composition.
    The standard fix (used by tools like RNAz) is to compare a window's MFE against the MFE of many randomized/shuffled versions of the same window,
    which preserves nucleotide composition but destroys any real structure.
    """
    fc = RNA.fold_compound(seq)
    _, mfe = fc.mfe()
    return mfe