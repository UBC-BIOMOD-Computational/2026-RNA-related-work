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
