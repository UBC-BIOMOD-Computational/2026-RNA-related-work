TODO later:
- reason more difficult since it requires external api


### Conservation Score

#### Purpose

Measures evolutionary conservation of the nucleotide sequence within a transcript region. Highly conserved regions are more likely to be functionally important.

#### Biological Interpretation

A high score indicates that the sequence has been preserved across multiple species.

#### Inputs

```python
class ConservationScoreInput:
    chromosome: str
    start: int
    end: int
```

or for transcript coordinates:

```python
class ConservationScoreInput:
    transcript_id: str
    transcript_start: int
    transcript_end: int
```

#### External Dependencies

Required:

* Precomputed conservation tracks:

  * phastCons100way
  * phyloP100way
  * GERP (optional)

Recommended source:

* UCSC Genome Browser downloads

Storage format:

```text
.bigWig
.bed
```

#### Algorithm

1. Map transcript coordinates to genomic coordinates.
2. Extract conservation values for all positions in the window.
3. Compute:

```python
score = np.mean(conservation_values)
```

#### Output

```python
float
```

Range depends on source track.

Higher values indicate stronger conservation.

#### Runtime Requirements

No online API required after initial data download.

Recommended architecture:

```python
ConservationTrackProvider
    └── pyBigWig
```

#### Complexity

```text
O(window_length)
```

---

### Structural Conservation Score

#### Purpose

Measures whether RNA secondary structure is conserved across homologous transcripts, even when the underlying sequence differs.

#### Biological Interpretation

A high score indicates evolutionary pressure to preserve RNA structure rather than nucleotide identity.

Example:

```text
Human: GGGAAAUCC
Mouse: GGCAAAUCG

Structure:
(((...)))
```

Sequence changed, structure preserved.

#### Inputs

```python
class StructuralConservationInput:
    target_sequence: str
    homolog_sequences: list[str]
```

Minimum recommended:

```python
[
    human,
    mouse,
    rat,
    dog
]
```

#### External Dependencies

Required:

* Ortholog retrieval source:

  * NCBI Orthologs
  * Ensembl Orthologs
  * UCSC Multiz Alignments

* Structural conservation tool:

  * RNAz
  * LocARNA
  * Infernal (advanced)

#### Algorithm

1. Retrieve homologous transcript sequences.
2. Generate multiple sequence alignment.
3. Evaluate structural conservation.

Possible implementations:

##### RNAz

Input:

```text
multiple sequence alignment
```

Output:

```python
probability_structured_rna
```

##### LocARNA

Input:

```text
multiple sequence alignment
```

Output:

```python
structural_similarity_score
```

##### Infernal

Input:

```text
alignment + covariance model
```

Output:

```python
covariation_score
```

#### Output

```python
float
```

Normalized:

```text
0.0 = no evidence of structural conservation
1.0 = highly conserved structure
```

#### Runtime Requirements

Requires external homologous sequence data.

Options:

| Approach                     | API Needed      |
| ---------------------------- | --------------- |
| NCBI Ortholog lookup         | Yes (initially) |
| Ensembl REST                 | Yes (initially) |
| Downloaded Multiz alignments | No              |
| Local homolog database       | No              |

#### Complexity

Approximately:

```text
O(N × alignment_cost × folding_cost)
```

where N is the number of homologs.

Significantly more expensive than Conservation Score.

---

### Recommended Implementation Order

#### Phase 1 (Sequence-only)

No external databases required.

```python
MFE
Accessibility
SelfSelfHybridization
StructurednessZ
GCContent
Entropy
Uniqueness
```

#### Phase 2 (Add conservation)

Requires:

```text
phastCons100way.bigWig
```

Implement:

```python
ConservationScore
```

#### Phase 3 (Comparative genomics)

Requires:

```text
Ortholog retrieval
Multiple sequence alignment
RNAz / LocARNA
```

Implement:

```python
StructuralConservationScore
```

For a CDH1/VIM target-discovery pipeline, Conservation Score is typically a core feature. Structural Conservation Score is an advanced feature and usually provides the most value when identifying conserved RNA regulatory elements rather than simply finding accessible hybridization sites.
