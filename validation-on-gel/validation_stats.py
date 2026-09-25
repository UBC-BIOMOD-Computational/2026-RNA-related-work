from collections import defaultdict

# ============================================================
# RNA native-PAGE structural analysis
# ============================================================

# Your sequences
vim = "CGAGGAGGAGAUGCGGGAGC"
cdh1 = "CCGCCCUGCCAAUCCCGAUG"
helper = "CAUCGGGAUUGGGAAAGGCAGGGCGGGCUCCCGCAUCUCCUCCUCG"

# X = vim + helper
X_sequences = [vim, helper]
X_structure = "((((((((((((((((((((..(((...(((.......)))..)))))))))))))))))))))"

# Z = cdh1 + helper
Z_sequences = [cdh1, helper]
Z_structure = "(((((((((((((((((((())))))))))......))))))))))((....))............"


# ============================================================
# Parsing
# ============================================================

def clean_structure(structure):
    """Remove whitespace from a dot-bracket structure."""
    return "".join(structure.split())


def combined_sequence(sequences):
    return "&".join(sequences)


def sequence_length(sequences):
    return sum(len(s) for s in sequences)


def chain_boundaries(sequences):
    """
    Return the zero-based nucleotide ranges occupied by each chain
    in the concatenated sequence.
    """
    boundaries = []
    start = 0

    for i, seq in enumerate(sequences):
        end = start + len(seq) - 1
        boundaries.append((i + 1, start, end))
        start = end + 1

    return boundaries


def parse_structure(structure, sequences):
    """
    Parse dot-bracket notation.

    Supports:
        ( )
        .
        &
    
    & separates chains and is NOT counted as a nucleotide.
    """
    structure = clean_structure(structure)

    # Remove chain separators for pairing analysis
    structure_no_amp = structure.replace("&", "")

    expected_length = sequence_length(sequences)

    if len(structure_no_amp) != expected_length:
        raise ValueError(
            f"Structure/sequence length mismatch:\n"
            f"  sequence length = {expected_length}\n"
            f"  structure length = {len(structure_no_amp)}\n"
        )

    # Map positions in the structure-with-& to nucleotide positions
    nt_positions = []
    nt_index = 0

    for char in structure:
        if char == "&":
            nt_positions.append(None)
        else:
            nt_positions.append(nt_index)
            nt_index += 1

    stack = []
    pairs = []

    for char, nt_pos in zip(structure, nt_positions):
        if char == "(":
            stack.append(nt_pos)

        elif char == ")":
            if not stack:
                raise ValueError("Unmatched ')' in structure.")

            partner = stack.pop()
            pairs.append((partner, nt_pos))

    if stack:
        raise ValueError("Unmatched '(' in structure.")

    return sorted(pairs)


# ============================================================
# Chain assignment
# ============================================================

def nucleotide_chain_map(sequences):
    """
    Map each nucleotide index to its chain number.
    """
    mapping = []

    for chain_id, seq in enumerate(sequences, start=1):
        mapping.extend([chain_id] * len(seq))

    return mapping


# ============================================================
# Base-pair statistics
# ============================================================

def pair_statistics(pairs, sequences):
    chain_map = nucleotide_chain_map(sequences)

    intramolecular = []
    intermolecular = []

    for i, j in pairs:
        if chain_map[i] == chain_map[j]:
            intramolecular.append((i, j))
        else:
            intermolecular.append((i, j))

    return intramolecular, intermolecular


def paired_unpaired(pairs, n):
    paired = set()

    for i, j in pairs:
        paired.add(i)
        paired.add(j)

    unpaired = set(range(n)) - paired

    return len(paired), len(unpaired), paired, unpaired


# ============================================================
# Stems
# ============================================================

def find_stems(pairs):
    """
    Identify contiguous helices where base pairs are stacked:

        (i,j)
        (i+1,j-1)
        (i+2,j-2)

    Returns stems as lists of base pairs.
    """

    pair_set = set(pairs)
    visited = set()
    stems = []

    for i, j in pairs:

        if (i, j) in visited:
            continue

        # This pair is the beginning of a stem if the previous
        # stacked pair does not exist.
        if (i - 1, j + 1) in pair_set:
            continue

        stem = []
        a, b = i, j

        while (a, b) in pair_set:
            stem.append((a, b))
            visited.add((a, b))
            a += 1
            b -= 1

        stems.append(stem)

    return stems


# ============================================================
# Loop / unpaired-region analysis
# ============================================================

def unpaired_runs(pairs, n):
    """
    Find consecutive runs of unpaired nucleotides.
    """
    paired = set()

    for i, j in pairs:
        paired.add(i)
        paired.add(j)

    runs = []

    start = None

    for i in range(n):

        if i not in paired:
            if start is None:
                start = i

        else:
            if start is not None:
                runs.append((start, i - 1))
                start = None

    if start is not None:
        runs.append((start, n - 1))

    return runs


# ============================================================
# Structural topology
# ============================================================

def topology_description(pairs, sequences):
    """
    Provide a simple qualitative topology description.
    """

    n = sequence_length(sequences)
    stems = find_stems(pairs)

    paired, unpaired, _, _ = paired_unpaired(pairs, n)

    # Count paired neighbors around each nucleotide
    pair_degree = defaultdict(int)

    for i, j in pairs:
        pair_degree[i] += 1
        pair_degree[j] += 1

    # Approximate branching:
    # number of distinct stem starts/endpoints
    stem_count = len(stems)

    if stem_count == 1:
        topo = "Predominantly single-stem / hairpin-like"
    elif stem_count == 2:
        topo = "Two-stem architecture"
    elif stem_count >= 3:
        topo = "Multi-stem / branched architecture"
    else:
        topo = "Largely unstructured"

    if unpaired / n > 0.45:
        topo += "; relatively extended due to substantial unpaired sequence"
    elif unpaired / n < 0.20:
        topo += "; relatively highly base-paired"

    return topo


# ============================================================
# Compactness heuristic
# ============================================================

def compactness_score(pairs, sequences):
    """
    HEURISTIC ONLY.

    Higher score = predicted more compact.
    This is NOT a physical prediction of hydrodynamic radius.

    Factors:
        + base-pairing fraction
        + longer stems
        + fewer long unpaired runs
        - large unpaired regions
        - many separate structural segments
    """

    n = sequence_length(sequences)

    paired, unpaired, _, _ = paired_unpaired(pairs, n)

    pairing_fraction = paired / n

    stems = find_stems(pairs)

    if stems:
        mean_stem = sum(len(s) for s in stems) / len(stems)
    else:
        mean_stem = 0

    runs = unpaired_runs(pairs, n)

    if runs:
        longest_unpaired = max(end - start + 1 for start, end in runs)
    else:
        longest_unpaired = 0

    # Base score
    score = 50

    # More paired bases generally favors compactness
    score += 35 * (pairing_fraction - 0.5)

    # Longer continuous stems favor a more organized structure
    score += min(mean_stem * 2.0, 12)

    # Long unpaired stretches favor extension
    score -= min(longest_unpaired * 1.0, 20)

    # Multiple independent stems can indicate branching/complexity
    if len(stems) >= 3:
        score -= 5

    return score


# ============================================================
# Native PAGE interpretation
# ============================================================

def mobility_prediction(score_x, score_z):
    difference = score_x - score_z

    if abs(difference) < 5:
        return (
            "X and Z have similar predicted compactness. "
            "A strong mobility difference would therefore be "
            "difficult to explain from this simple structural metric alone."
        )

    if difference > 0:
        return (
            "X is predicted to be more compact than Z. "
            "On a native gel, X would heuristically be expected "
            "to migrate faster than Z."
        )

    return (
        "Z is predicted to be more compact than X. "
        "On a native gel, Z would heuristically be expected "
        "to migrate faster than X."
    )


# ============================================================
# Report
# ============================================================

def analyze(name, sequences, structure):

    print("=" * 70)
    print(f"{name} ANALYSIS")
    print("=" * 70)

    print("\nChains:")
    for i, seq in enumerate(sequences, start=1):
        print(f"  Chain {i}: {len(seq)} nt")

    n = sequence_length(sequences)

    print(f"\nTotal nucleotides: {n}")

    print("\nSequence:")
    print(combined_sequence(sequences))

    print("\nSecondary structure:")
    print(structure)

    pairs = parse_structure(structure, sequences)

    paired, unpaired, _, _ = paired_unpaired(pairs, n)

    print("\n1. PAIRED vs UNPAIRED")
    print(f"  Paired nucleotides:   {paired}")
    print(f"  Unpaired nucleotides: {unpaired}")
    print(f"  Paired fraction:      {paired/n:.2%}")
    print(f"  Unpaired fraction:    {unpaired/n:.2%}")

    stems = find_stems(pairs)

    print("\n2. STEMS")

    if not stems:
        print("  No stems detected.")
    else:
        for i, stem in enumerate(stems, start=1):
            print(
                f"  Stem {i}: {len(stem)} bp "
                f"(positions {stem[0][0]+1}-{stem[0][1]+1} "
                f"through {stem[-1][0]+1}-{stem[-1][1]+1})"
            )

    runs = unpaired_runs(pairs, n)

    print("\n3. UNPAIRED REGIONS / LOOP-LIKE REGIONS")

    if not runs:
        print("  None.")
    else:
        for start, end in runs:
            print(
                f"  nt {start+1}-{end+1}: "
                f"{end-start+1} unpaired nt"
            )

    intra, inter = pair_statistics(pairs, sequences)

    print("\n4. BRANCHING / TOPOLOGY")
    print(f"  Number of stems: {len(stems)}")

    print("\n5. INTERMOLECULAR BASE PAIRS")

    if not inter:
        print("  None detected.")
    else:
        chain_map = nucleotide_chain_map(sequences)

        for i, j in inter:
            print(
                f"  Chain {chain_map[i]} nt {i+1} "
                f"<-> Chain {chain_map[j]} nt {j+1}"
            )

    print("\n6. PREDICTED SECONDARY-STRUCTURE TOPOLOGY")
    print(f"  {topology_description(pairs, sequences)}")

    score = compactness_score(pairs, sequences)

    print("\n7. COMPACTNESS HEURISTIC")
    print(f"  Score: {score:.1f}")
    print("  Higher score = predicted more compact.")
    print("  IMPORTANT: this is a qualitative heuristic, not a physical measurement.")

    return score


# ============================================================
# Run X and Z
# ============================================================

score_X = analyze(
    "X",
    X_sequences,
    X_structure
)

score_Z = analyze(
    "Z",
    Z_sequences,
    Z_structure
)

print("\n" + "=" * 70)
print("8. NATIVE PAGE MOBILITY INTERPRETATION")
print("=" * 70)

print(mobility_prediction(score_X, score_Z))

print("\nInterpretation:")
print(
    "Native PAGE mobility depends on conformation, compactness, "
    "hydrodynamic shape, gel conditions, and RNA-RNA interactions. "
    "Therefore this script provides a hypothesis, not a guaranteed "
    "prediction of which band will run farther."
)