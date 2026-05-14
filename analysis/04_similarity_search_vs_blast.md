# Similarity search vs BLAST

## Experimental setup

**Query proteins:** 50 proteins selected by stratified sampling — each protein family contributes proportionally to the query set. All queries have at least one GO biological process term, which is required for precision evaluation.

**Embedding search:** brute-force nearest neighbour search using cosine, Euclidean, and Manhattan distance. Queries run against a 376-protein dataset (500 proteins filtered to those with both family labels and GO biological process annotations). Embeddings computed with ESM2 `esm2_t6_8M_UR50D` (8 million parameters, 6 transformer layers, 320-dimensional output).

**BLAST:** `blastp` against full human SwissProt (~20,000 reviewed human proteins) via NCBI web API. Top-10 hits retrieved per query, excluding the query protein itself.

**Evaluation metric:** precision@10 — the fraction of top-10 retrieved proteins that share at least one GO biological process term with the query. GO biological process was chosen over molecular function or cellular component because it captures functional context at the level most relevant for assessing whether two proteins do the same thing in the cell.

## Distance metric comparison — embedding search

Precision@10 was computed across all 50 query proteins for each of the three distance metrics:

| Metric | Mean precision@10 |
|--------|------------------|
| Cosine | 0.114 |
| Euclidean | 0.114 |
| Manhattan | 0.116 |

All three metrics perform identically. This confirms the Week 1 finding — the tight L2 norm distribution (norms between 4.5 and 6.0) means cosine and Euclidean distances are effectively equivalent for this dataset. When embedding magnitudes are uniform, the angle between vectors (cosine) and the straight-line distance (Euclidean) measure the same thing.

## BLAST vs embedding search — 10 selected proteins

A subset of 10 proteins was selected for the BLAST comparison, covering a range of functional types: enzymes, receptors, structural proteins, signalling proteins, and immune proteins.

| UniProt ID | Protein | Type | BLAST p@10 | Embedding p@10 |
|-----------|---------|------|-----------|----------------|
| P08922 | Proto-oncogene tyrosine kinase ROS | enzyme | 0.600 | 0.200 |
| P08172 | Muscarinic acetylcholine receptor M2 | receptor | 1.000 | 0.700 |
| O75369 | Filamin-B | structural | 0.700 | 0.100 |
| P05787 | Keratin type II cytoskeletal 8 | structural | 0.000 | 0.000 |
| P01571 | Interferon alpha-17 | signalling | 0.900 | 0.000 |
| O00746 | Nucleoside diphosphate kinase D | enzyme | 0.750 | 0.000 |
| P06865 | Beta-hexosaminidase subunit alpha | enzyme | 0.200 | 0.100 |
| P01909 | HLA class II histocompatibility antigen | immune | 1.000 | 0.100 |
| P05997 | Collagen alpha-2(V) chain | structural | 0.800 | 0.300 |
| P09683 | Secretin | signalling | 0.500 | 0.200 |
| **Mean** | | | **0.645** | **0.170** |

## Precision@10 by protein type

Averaging results across protein types reveals where each method performs relatively better or worse:

| Type | BLAST | Embedding (cosine) |
|------|-------|-------------------|
| Enzyme | 0.517 | 0.100 |
| Immune | 1.000 | 0.100 |
| Receptor | 1.000 | 0.700 |
| Signalling | 0.700 | 0.100 |
| Structural | 0.500 | 0.133 |

The receptor category is the most notable — embedding search achieves 0.700, much closer to BLAST's 1.000 than in any other category. This is the GPCR result: GPCRs are the most abundant protein family in the dataset (7 proteins), so when querying a GPCR the relevant proteins are present in the 376-protein search space and the embeddings find them.

## Why the comparison is confounded

The BLAST vs embedding comparison differs in two variables simultaneously:

1. **Similarity metric** — BLAST uses sequence alignment (BLOSUM62 scoring matrix), embedding search uses vector distance in ESM2 representation space
2. **Search space size** — BLAST searches ~20,000 human proteins, embedding search searches 376 proteins

These cannot be disentangled from the current results. A protein not found by embedding search may be absent because the embeddings fail to capture the similarity, or simply because the protein is not in the 376-protein dataset.

Evidence that dataset size is a major factor: filtering BLAST hits to only proteins present in the 376-protein dataset showed that 7 out of 10 query proteins have zero BLAST hits in the local dataset. BLAST and embedding search are effectively searching different spaces.

The GPCR result (embedding p@10=0.700) provides the clearest evidence that embeddings can work well when the search space is adequate — GPCRs are well-represented in the dataset, and the embeddings correctly group them together.

## What a fair comparison requires

To isolate the effect of the similarity metric from the effect of dataset size, a fair comparison would require:

- ESM2 embeddings computed for all ~20,000 reviewed human proteins in SwissProt
- The same query proteins and the same GO biological process evaluation criterion
- Both methods searching the same protein space

Additionally, the embedding quality is limited by the model size. The ESM2 8M parameter model (6 layers) is the smallest available variant. The largest ESM2 model has 650 million parameters and 33 layers, and has been shown in published benchmarks to significantly outperform the small model on remote homology detection. A GPU would be required to run the larger model.

## Key conclusions

BLAST on full human SwissProt outperforms ESM2 small model embedding search on a 376-protein subset (mean precision@10: 0.645 vs 0.170). However, this result reflects both the metric difference and the dataset size difference and cannot be attributed to either alone.

The receptor/GPCR case demonstrates that embedding similarity can reach BLAST-comparable precision (0.700 vs 1.000) when the relevant proteins are present in the search space. This is consistent with the clustering results from Week 3, where GPCRs consistently formed the most coherent and well-separated cluster across all algorithms.

The most important limitation is dataset scale. Running this comparison on the full human SwissProt with a larger ESM2 model is the natural next step and would provide a definitive answer to the core biological question: do protein language model embeddings capture functional relationships that sequence alignment misses?