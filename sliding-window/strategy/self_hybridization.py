# Self-Self Hybridization Potential
"""
Biological Rationale

A region should strongly hybridize to a complementary partner.

Regions capable of forming stable intermolecular duplexes have favorable hybridization energies.

Input
window_seq: str

Output
deltaG_hybridization: float
Algorithm

Generate reverse complement:

partner = reverse_complement(window_seq)

Compute duplex energy:

RNA.duplexfold(window_seq, partner)

Return:
ΔG_duplex

Interpretation
ΔG = Very negative	= Strong binding
ΔG = Near zero =	Weak binding
"""