import RNA
import random
import numpy as np
from strategy.mfe import get_mfe

strategies = {
    "mfe": get_mfe
    # TODO - please verify and check 
}

strategy = strategies['mfe']

def sliding_window(sequence: str, k: int, step: int = 1) -> tuple[np.ndarray, np.ndarray]:
    global strategy
    seq = sequence.upper().replace("T", "U")
    n = len(seq)
    starts = np.arange(0, n - k + 1, step)
    results = np.empty(len(starts))

    for i, s in enumerate(starts):
        results[i] = strategy(seq[s:s + k])

    return starts, results
