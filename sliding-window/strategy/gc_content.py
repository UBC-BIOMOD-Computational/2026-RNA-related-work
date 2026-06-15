#GC Content Score
"""
Biological Rationale

GC content strongly influences duplex stability and secondary structure.

Input
window_seq
Output
gc_fraction
Algorithm
(G + C) / length
Interpretation
gc % = Very low	Weak duplex
gc % = Very high	Overstructured
gc % = Moderate	Usually optimal

Typical desirable range:

40–65%
"""
