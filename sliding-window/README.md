# RNA Sliding Window Duplex Scanner

Basically the purpose of this function is to see if there are subsequences of window k in the overall total sequence  that will preferentially form complementary duplexes (ie: 1-to-1 double stranded RNA) instead of staying as two single strands

## Setup

```bash
python -m venv .venv

# Windows
.venv\Scripts\activate
# macOS / Linux
source .venv/bin/activate

pip install -r requirements.txt
```

> If `pip install ViennaRNA` fails, use conda:
> ```bash
> conda install -c bioconda viennarna
> ```

## Run

```bash
python main.py
```

Scans all window sizes (42, 60, 100, 200, 300 nt) and prints the top 5 CDH1/VIM pairs per size.

## Messy-fragment comparison

`messy_fragments.py` generates noisy variants around each clean reference window and compares them to the same fixed reference partner. By default it uses non-overlapping 20-nt windows from the bundled CDH1 and VIM FASTA files. For each reference it evaluates:

- length variants of 18, 19, and 21 nt, anchored at the 5′ end;
- 20-nt windows shifted by ±1 or ±2 nt where transcript boundaries allow;
- deterministic single-base insertions, deletions, and substitutions at the center;
- the clean reference as a control.

Run it from any directory:

```bash
python /path/to/sliding-window/messy_fragments.py
```

The CSV is written to `sliding-window/data/exports/messy_fragments_k20.csv`. Each row preserves the raw measurements and includes sequence edit distance, MFE-based cofold binding energy, mapped interstrand base-pair distance, ViennaRNA abstract-shape distance, and method-specific ranks. It reports two equal-weight percentile composites: one using base-pair distance and a second using abstract-shape distance in its place.

Options for finer scanning or exhaustive internal edits:

```bash
python messy_fragments.py --step 1 --output data/exports/messy_fragments_k20_step1.csv
python messy_fragments.py --internal-edits all --output data/exports/messy_fragments_all_edits.csv
```

For a specific reference fragment embedded in a short FASTA record with flanking bases, use `--reference-start` (1-based within each record). The supplied CDH1/VIM fragments have two bases of context on both sides, so the reference starts at position 3:

```bash
python messy_fragments.py \
  --fasta data/messy_fragment_inputs_with_buffer.fasta \
  --reference-start 3 \
  --output data/exports/my_fragments_with_shifts.csv
```

This analyzes just the selected 20-nt reference from each record while using the surrounding transcript bases for the ±1/±2 shifts and the natural 21-nt extension.

`--internal-edits all` evaluates every single-base substitution, insertion, and deletion for every reference window and can generate a large CSV. Use `--internal-edits none` to omit synthetic edits. The shape metric uses ViennaRNA's `RNA.abstract_shapes()` interface, available starting with ViennaRNA 2.4.18.

## Tuning

In `main.py`, inside the `for k in WINDOW_SIZES` loop:

| Parameter | Default | Effect |
|-----------|---------|--------|
| `step` | `k` (non-overlapping tiles) | Lower → finer resolution, slower |
| `top_n` | `5` | More candidates returned per window size |

## Interpreting results

```
ΔG₁ / ΔG₂      MFE of each strand alone — how stable its self-structure is
ΔG_cofold       MFE of the joined duplex
ΔG_binding      ΔG_cofold − (ΔG₁ + ΔG₂)
                Negative = energy released by binding; more negative = stronger duplex
acc_CDH1/VIM    Mean P(unpaired) ∈ [0, 1]
                Higher = binding region is open and accessible
score           −ΔG_binding + mean(acc) — composite rank; higher is better
```

**Binding strength labels**

| Label | ΔG_binding |
|-------|-----------|
| STRONG | < −10 kcal/mol |
| moderate | −5 to −10 kcal/mol |
| weak | > −5 kcal/mol |

A pair is a good duplex candidate when:
- `ΔG_binding` is strongly negative
- `acc_*` values are high (binding regions not buried in hairpins)
- The cofold structure shows a long uninterrupted duplex stretch
