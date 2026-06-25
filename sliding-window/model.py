from dataclasses import dataclass, field


# ── single-strand analysis ───────────────────────────────────────────────────
@dataclass
class StrandProfile:
    start: int
    seq: str
    mfe: float
    accessibility: float   # mean P(unpaired)
    opening_cost: float    # kcal/mol estimated cost to open structure



# ── pairwise duplex analysis ─────────────────────────────────────────────────
@dataclass
class DuplexResult:
    start_a: int
    start_b: int
    seq_a: str
    seq_b: str
    mfe_a: float
    mfe_b: float
    mfe_cofold: float
    binding_dg: float           # ΔG_cofold − (ΔG_a + ΔG_b)
    acc_a: float
    acc_b: float
    structure_cofold: str
    score: float = field(init=False)  # composite ranking score

    def __post_init__(self):
        """
        Compute composite ranking score after field initialisation.

        score = −ΔG_binding + mean(acc_a, acc_b)

        Both terms are positive-good:
          −ΔG_binding  is large when duplex formation releases a lot of energy.
          mean_acc     is large when both binding regions are already open.
        Higher score → pair is both energetically and structurally favourable.
        """
        self.score = -self.binding_dg + (self.acc_a + self.acc_b) / 2
