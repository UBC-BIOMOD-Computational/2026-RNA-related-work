# Messy-fragment evaluation plan

**Status:** initial implementation is available in [`messy_fragments.py`](./messy_fragments.py). The buffered CDH1 and VIM input is in [`data/messy_fragment_inputs_with_buffer.fasta`](./data/messy_fragment_inputs_with_buffer.fasta). The two-base flanks allow transcript-derived shifts up to ±2 nt and a natural 21-nt extension. The generated CSV is a run output and is not committed.

## Goal

Reuse the existing sliding-window workflow to ask how robustly three kinds of measurements distinguish a clean RNA fragment from plausible messy versions:

1. **Sequence similarity:** how many edits separate a candidate from the expected fragment?
2. **Binding thermodynamics:** how favorably does the candidate bind its intended partner?
3. **Structure similarity:** does the predicted paired structure keep the same base-pair pattern or abstract shape?

The experiment should make length and position errors explicit. For the current `k=20` setting, examples include 18-, 19-, 20-, and 21-nt fragments and fragments whose start is shifted by one or more nucleotides.

## What the existing code already does

- [`main.py`](./main.py) currently sets `WINDOW_SIZES = [20]` and calls `scan_pairwise(..., step=k)`. That samples non-overlapping 20-nt windows.
- [`processing.py`](./processing.py) has the reusable window generator. Setting `step=1` would enumerate every possible start position; explicit start offsets would make the perturbation cases easier to label.
- [`strategy/cofold.py`](./strategy/cofold.py) folds each strand alone, folds the two strands together with ViennaRNA, and calculates `mfe_cofold - (mfe_1 + mfe_2)`. It returns the MFE structure and energies. It does **not** currently calculate a cofold partition function or equilibrium concentrations.
- `scan_pairwise` currently pairs each window with **its own reverse complement**. That is appropriate for asking whether a sequence can form a self-complementary duplex, but it changes both strands when the candidate changes. To isolate the effect of a messy fragment, the recommended first experiment keeps the partner fixed: derive it once from the clean reference fragment, then cofold every candidate against that same partner.
- The current composite score also includes single-strand accessibility. Keep accessibility as an optional diagnostic for this experiment; first compare the three requested measurements on their own so the source of a ranking is visible.

## Proposed experiment setup

### 1. Define a clean reference

Choose a 20-nt parent window from the transcript and record its sequence, transcript label, and start/end coordinates. Treat this as the expected clean fragment. The corresponding reference partner is its reverse complement, fixed for all candidate comparisons.

### 2. Generate labeled candidate variants

Keep the error types separate so we can tell what caused a score to change.

| Variant class | Example for a 20-nt reference | What it tests |
|---|---|---|
| Clean control | Exact 20-nt reference | Baseline and self-consistency |
| Length | 18, 19, 21, and optionally 22 nt, with a documented anchor | Truncation or extension |
| Offset | 20-nt windows starting 1 or 2 nt to the left/right | Boundary or coordinate error |
| Internal edit | One inserted, deleted, or substituted nucleotide | Sequence corruption independent of window location |

For length variants, record whether the 5′ end or 3′ end is held fixed. A 21-nt natural window and a 20-nt window with an inserted base are different cases, so label them separately. Begin with one change at a time; only add combined errors after the single-error behavior is understood. Use fixed random seeds if internal edits are generated randomly.

### 3. Compare each candidate against the same reference context

Use the candidate as strand A and the fixed clean-reference partner as strand B for cofolding. This lets the measurements respond to a truncation, extension, offset, or internal edit in strand A. Record the source coordinates even when the sequence itself was edited synthetically.

## The three measurements

### A. Sequence edit distance

Use global Levenshtein distance as the first, interpretable sequence baseline: the minimum number of substitutions, insertions, and deletions needed to turn the reference into the candidate. Save both:

- the raw edit count; and
- a length-normalized distance, `edit_count / max(reference_length, candidate_length)`.

The alphabet being RNA does not, by itself, make ordinary sequence edit distance different: normalize T/U consistently and use the same edit costs for A, C, G, and U. Also save the known genomic/transcript offset as its own field. A sequence edit metric measures sequence changes; the offset field measures where the fragment came from.

“RNA edit distance” can also mean a **structure-aware** distance. Those methods treat paired bases as arcs and can assign costs to disrupting or changing a base pair; that is a different comparison from sequence Levenshtein distance. For this first experiment, keep sequence edit distance simple and test structure-aware distances under the structure measurement below. The RNA-structure alignment literature describes these as distinct edit operations and alignment problems ([Chiu & Chen, 2017](https://doi.org/10.1093/bib/bbw009)).

### B. Cofold free energy

For each candidate and the fixed partner, record:

- each strand’s MFE;
- the two-strand cofold MFE and cofold dot-bracket structure;
- the current MFE-based binding difference, `ΔG_bind,MFE = MFE_complex - (MFE_A + MFE_B)`;
- the number of interstrand base pairs and the longest continuous interstrand run, as explanatory diagnostics.

Because longer fragments can make more base pairs and gain more total energy, retain the raw binding energy but also report it normalized by the maximum possible contact length, `min(candidate_length, partner_length)`. Do not use raw total ΔG alone to compare different lengths.

ViennaRNA cofolding can report more than one number: it predicts a joint secondary structure and its MFE; partition-function calculations add ensemble free energy and base-pair probabilities; concentration-aware calculations can estimate monomer/dimer equilibrium when starting concentrations are supplied. The current `cofold_energies` helper uses the MFE route only. If ensemble binding is important, include it as a separately labeled follow-up feature, using consistent conditions and parameters across candidates. The ViennaRNA documentation describes cofolding and its partition-function/concentration calculations ([RNA-RNA interaction documentation](https://viennarna.readthedocs.io/en/latest/cofold.html), [RNAcofold manual](https://www.tbi.univie.ac.at/RNA/usage/man/RNAcofold.html)).

### C. Structure distance and shape abstraction

Dot-bracket notation is a **way to write a predicted structure**, not a ranking algorithm by itself. ViennaRNA can compare structures with base-pair distance or RNAdistance’s string/tree-edit modes. A base-pair distance counts pairs present in one structure but absent from the other; compare the actual base-pair sets rather than applying plain character edit distance to `(`, `)`, and `.` because the two bracket positions jointly encode a pair. ViennaRNA warns that base-pair distance is not recommended for structures with different lengths ([RNAdistance manual](https://www.tbi.univie.ac.at/RNA/ViennaRNA/refman/usage/man/RNAdistance.html)).

For the first structure score:

1. Predict the reference and candidate cofold structures against the fixed partner.
2. Map positions back to the reference coordinates (and retain inserted/deleted positions separately).
3. Compare mapped interstrand base-pair sets; report raw symmetric-difference distance and normalize it by the size of the union of the two pair sets (define the distance as 0 if both sets are empty). Do not compare unrelated sequence positions as if they were aligned.

RNAshapes is a real, literature-described package for grouping structures into abstract shapes and estimating shape probabilities. It offers five abstraction levels: lower levels preserve more detail; higher levels collapse more details. It is useful for testing whether small changes preserve overall fold topology, but a shape label alone is not a distance or a rank. Compare the candidate and reference shape strings, or use the probability mass assigned to the reference shape when ensemble shape probabilities are available. The original work describes shape classes and representatives ([Giegerich et al., 2004](https://doi.org/10.1093/nar/gkh779)); the package paper describes shape probabilities and five levels ([Steffen et al., 2006](https://doi.org/10.1093/bioinformatics/btk010)).

Use mapped base-pair distance as the primary structure feature in the initial composite. Treat abstract-shape distance/probability as an exploratory alternative or sensitivity check, rather than adding both to the composite at first; both describe structure and could double-count the same signal.

**Cofold compatibility needs checking before choosing RNAshapes as the implementation.** The RNAshapes manual documents folding single RNA sequences and an `abstract` mode that maps a Vienna dot-bracket structure to a shape string, but it does not document a two-strand cofold mode. A practical route is to test whether it accepts ViennaRNA’s cofold dot-bracket output with the `&` strand separator. If it does not, ViennaRNA also provides abstract-shape functions; use one of those on the cofold structure while retaining the strand cut-point as metadata. Do not silently drop the strand boundary.

Literature comparisons are task-dependent. A benchmark of evolutionary structure conservation found simple energy-based and base-pair-distance measures strong in that setting, with more elaborate tree-edit measures not consistently improving results; it did not benchmark noisy-fragment recovery, so it motivates trying simple baselines rather than declaring a winner in advance ([Gruber et al., 2008](https://doi.org/10.1186/1471-2105-9-122)).

## Candidate ranking formula

Keep all raw measurements in the output. For a first composite ranking, convert each of the three measurements to a within-reference-set **goodness percentile** on `[0, 1]`, with 1 meaning best:

- `P_seq`: smaller normalized edit distance is better;
- `P_energy`: more favorable length-normalized binding energy is better;
- `P_shape`: smaller mapped interstrand base-pair distance is better. If the cofold-compatible shape method works, it can replace this feature in a separate comparison.

Then use an equal-weight exploratory score:

```text
composite_score = (P_seq + P_energy + P_shape) / 3
```

The script reports this first with mapped base-pair distance as the structure feature. It also reports an alternative `shape_composite_score` that substitutes the abstract-shape distance for base-pair distance; it does not add both structure measures to one score.

This avoids adding unlike raw units (nucleotide edits, kcal/mol, and structure distance). It is a transparent starting point, not a validated biological model. Report the three individual ranks beside the composite rank, and check whether conclusions change when weights vary. If experimental binding or known successful/failed fragments become available, use those labels to tune weights and evaluate on held-out fragments. Without such labels, the synthetic perturbation categories can test whether the clean reference is recovered and how scores react to each error type, but cannot establish that one ranking is biologically best.

## Results to save

Write one row per reference/variant/partner combination to a CSV or TSV. Include:

- reference ID, candidate ID, source coordinates, length, offset, perturbation type, and cofold cut-point position;
- reference sequence, candidate sequence, fixed partner sequence;
- raw and normalized edit distance;
- individual and cofold MFE, MFE-based binding difference, normalized binding difference, cofold structure, interstrand pair count, and longest paired run;
- mapped base-pair distance, shape level/string/distance, and shape probability if calculated;
- each method’s goodness percentile/rank, the base-pair composite score, and the alternative abstract-shape composite score;
- ViennaRNA/RNAshapes versions and thermodynamic settings.

Summarize per perturbation class: the rank of the clean control, how rankings change with edit/offset size, and where the methods disagree. Keep energy and structure plots length-stratified as well as pooled.

## Suggested first-pass choices

- Start from the existing 20-nt window size.
- Test 18–21 nt, offsets of ±1 and ±2 nt, and single internal edits as separate sets.
- Hold the clean reference partner fixed for every candidate.
- Use global normalized Levenshtein, MFE-based cofold binding energy plus a length-normalized view, and mapped interstrand base-pair distance as the three initial features.
- Add RNAshapes/shape-probability analysis after confirming how cofold structures and the `&` strand cut-point can be represented.
- Use the equal-weight percentile formula only as an exploratory first ranking.
