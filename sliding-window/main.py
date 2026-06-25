import os
import RNA
from Bio import SeqIO

from model import DuplexResult
from processing import scan_pairwise

# ── sequence sources ────────────────────────────────────────────────────────
CDH1_PATH = './data/CDH1__NM_004360.5.fasta'
VIM_PATH  = './data/VIM__NM_003380.5.fasta'

# ── window sizes to evaluate ────────────────────────────────────────────────
WINDOW_SIZES = [42, 60, 100, 200, 300]



def print_result(r: DuplexResult, rank: int):
    """
    Print a single DuplexResult in a human-readable summary block.

    Binding strength labels follow the rough thresholds from the literature:
      STRONG   ΔG_binding < −10 kcal/mol
      moderate −10 ≤ ΔG_binding < −5 kcal/mol
      weak     ΔG_binding ≥ −5 kcal/mol
    """
    strength = (
        "STRONG"   if r.binding_dg < -10 else
        "moderate" if r.binding_dg < -5  else
        "weak"
    )
    print(
        f"  #{rank:>2}  [{r.start_a}]  vs  [{r.start_b}]\n"
        f"        ΔG_a={r.mfe_a:+.2f}  ΔG_b={r.mfe_b:+.2f}  "
        f"ΔG_cofold={r.mfe_cofold:+.2f}  ΔG_binding={r.binding_dg:+.2f} kcal/mol  [{strength}]\n"
        f"        acc_a={r.acc_a:.3f}  acc_b={r.acc_b:.3f}  "
        f"score={r.score:.3f}\n"
        f"        structure: {r.structure_cofold}\n"
    )


def determine_best_strand_for_duplex(rna_seq: str, label: str = "seq"):
    """
    Run scan_pairwise across all WINDOW_SIZES for a single transcript and
    print the top candidates at each resolution.

    Args:
        rna_seq: Full transcript, RNA (U-alphabet).
        label:   Short name for this sequence (used in histogram filenames).
    """
    for k in WINDOW_SIZES:
        print(f"=== {label}  window k={k} ===")
        top = scan_pairwise(rna_seq, k=k, step=k, top_n=5, label=label)
        print(f"  Top {len(top)} duplex candidates:")
        for rank, r in enumerate(top, 1):
            print_result(r, rank)

# ── entry point ──────────────────────────────────────────────────────────────
if __name__ == "__main__":
    cdh1_record = next(SeqIO.parse(CDH1_PATH, "fasta"))
    vim_record  = next(SeqIO.parse(VIM_PATH,  "fasta"))

    cdh1_seq = str(cdh1_record.seq).upper().replace("T", "U")
    vim_seq  = str(vim_record.seq).upper().replace("T", "U")

    print(f"CDH1: {len(cdh1_seq)} nt   VIM: {len(vim_seq)} nt\n")
    determine_best_strand_for_duplex(cdh1_seq, label="CDH1")
    determine_best_strand_for_duplex(vim_seq,  label="VIM")


    

    