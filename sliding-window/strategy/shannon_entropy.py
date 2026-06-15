#Shannon Entropy Score

"""
Biological Rationale

Low complexity sequences often indicate repeats and poor targeting regions.

Entropy measures sequence complexity.

Input
window_seq

Output
entropy
Algorithm

Compute nucleotide frequencies:
pA,pC,pG,pU

Then calculate H = - ∑p_i *   log_2(p_i)

H = entropy


Interpretation
Entropy = ~2.0	Highly diverse
Entropy Low	Repetitive

"""