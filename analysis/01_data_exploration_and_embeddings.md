# Week 1 — Data pipeline, embeddings, and first search

## Dataset

- Source: UniProt REST API, reviewed human proteins (SwissProt, organism_id:9606)
- Size: 500 proteins with columns: uniprot_id, protein_name, sequence, go_terms, family
- Embedding model: ESM2 `esm2_t6_8M_UR50D` — Meta's protein language model, smallest variant, 6 transformer layers, 8 million parameters
- Embedding dimensions: 320 per protein
- Final dataset shape: annotations (500, 5), embeddings (500, 320)

## Data quality

| Field | Missing values |
|-------|---------------|
| go_terms | 1 (0.2%) |
| sequence | 0 |
| family | 117 (23.4%) |

340 unique protein families across 500 proteins. The top families by count are G-protein coupled receptors (7), heparin-binding growth factors (5), and intermediate filaments (4). Most families have only 1-2 representatives in the dataset — a natural consequence of sampling 500 proteins from a diverse proteome.

117 proteins have no family label. These were kept in the full dataset for clustering and similarity search. For evaluation tasks requiring ground truth labels (ARI computation), only the 383 labeled proteins were used.

## Sequence length distribution

| Statistic | Value |
|-----------|-------|
| Min | 16 aa |
| Median | 435 aa |
| Mean | 573 aa |
| Max | 4834 aa |

The high standard deviation (530 aa) reflects the diversity of human proteins — from small peptide hormones like insulin (110 aa) to large structural proteins. ESM embeddings were computed one protein at a time to avoid memory errors from batching long sequences.

## Embedding computation

Initial batch computation caused memory errors on CPU. Fixed by processing proteins one at a time in a loop. Each protein's embedding is computed as the mean of per-residue representations across all amino acid positions (mean pooling), giving one 320-dimensional vector per protein regardless of sequence length.

## L2 norm distribution

The L2 norm of each embedding vector was computed as a proxy for embedding magnitude. Results:

- Most norms cluster tightly between 4.5 and 6.0
- Distribution is roughly bell-shaped, slightly left-skewed
- No extreme outliers — no proteins with norms near 0 or above 10

This tight distribution has a direct methodological implication: when embedding magnitudes are similar across all proteins, cosine distance (which measures the angle between vectors, ignoring magnitude) and Euclidean distance (which measures straight-line distance, sensitive to magnitude) become nearly equivalent. This was confirmed empirically in the nearest neighbour experiment below.

## Nearest neighbour search — three distance metrics

All three distance metrics were implemented: cosine similarity, Euclidean distance, and Manhattan distance. The search was tested on insulin (P01308, Insulin family).

All three metrics returned identical top-5 neighbours in the same order:

| Rank | Protein | Family |
|------|---------|--------|
| 1 | Secretin | Glucagon family |
| 2 | Left-right determination factor 1 | TGF-beta family |
| 3 | Gastrin-releasing peptide | Bombesin/neuromedin-B/ranatensin family |
| 4 | Acetylcholine receptor subunit gamma | Ligand-gated ion channel family |
| 5 | Fibroblast growth factor 19 | Heparin-binding growth factors family |

The identical results across all three metrics confirm the L2 norm finding — when magnitudes are uniform, the three metrics measure essentially the same thing.

The top results are biologically plausible. Secretin and gastrin-releasing peptide are both peptide hormones involved in signalling, the same broad functional category as insulin. The embeddings are capturing functional similarity even in this small dataset.

## Key decisions

- Dataset size set to 500 proteins — large enough for meaningful clustering and visualisation, small enough to compute embeddings on CPU in reasonable time
- ESM2 small model (8M parameters) used instead of larger variants (up to 650M parameters) due to CPU constraint — this limits embedding quality compared to GPU-based runs
- All 500 proteins retained regardless of sequence length — memory errors that initially occurred with batched processing were resolved by single-protein processing
- Proteins without family labels retained for clustering but excluded from ARI evaluation
- Data saved to disk immediately after computation to avoid re-fetching and re-computing