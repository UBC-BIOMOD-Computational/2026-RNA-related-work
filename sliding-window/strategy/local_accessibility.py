# Local Accessibility Score

"""
Biological Rationale

A region can only participate in intermolecular hybridization if it is physically accessible. Highly structured regions are often buried inside stable intramolecular base pairs and are less available for binding.

This approximates accessibility under in vitro conditions.

Input
sequence: str
window_size: int

Output
accessibility_score: float

Higher = more accessible.

Algorithm
For each sliding window:

Predict secondary structure using ViennaRNA.
Compute MFE structure.

Calculate:
accessibility = fraction of nucleotides that are unpaired

Alternative:
Use partition function and base-pair probabilities.

accessibility = mean(P(unpaired)) across all nucleotides.

Interpretation
score = High = Open region
score = Low = Buried in local structure
"""