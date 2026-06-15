# Codon Bias / Translational Structure Score
"""
Biological Rationale

Within coding regions, synonymous codon choice affects:

translation rate
ribosome occupancy
cotranslational folding
mRNA structure

Highly biased codon regions often coincide with structured mRNA segments.

Input
coding_sequence
organism_codon_table
Output
codon_bias_score
Algorithm

Compute one of:

CAI

Codon Adaptation Index

CAI(window)
tAI

tRNA Adaptation Index

tAI(window)
Relative Synonymous Codon Usage
RSCU(window)
Interpretation
score = High	Strong codon preference
score = Low	Neutral codon usage
"""
