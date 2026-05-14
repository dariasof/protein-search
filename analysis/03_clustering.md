# Clustering analysis

## Algorithm assumptions

Three clustering algorithms were compared, each making fundamentally different assumptions about the structure of the data.

K-means assumes clusters are spherical and of roughly equal size. It requires specifying k in advance and assigns every protein to a cluster — there is no concept of noise or outliers. It is sensitive to random initialisation, though sklearn's k-means++ initialisation mitigates this somewhat.

HDBSCAN is density-based and makes no assumptions about cluster shape. It determines k automatically by finding regions of high density separated by regions of low density. Points that do not belong to any dense region are labelled as noise (-1). This is biologically meaningful — some proteins genuinely do not fit into any family.

Hierarchical clustering (Ward linkage) builds a tree of nested clusters by iteratively merging the pair of clusters whose merger causes the smallest increase in within-cluster variance. It requires no k in advance — you choose the resolution by cutting the dendrogram at a given height. Ward linkage tends to produce compact, evenly-sized clusters and is well suited to biological data.

## Elbow experiment — how many clusters?

Silhouette score was computed on both the UMAP 2D projection and the full 320-dimensional embeddings across k = 5, 10, 15, 20, 25, 30, 40, 50.

| Space | k=5 | k=10 | k=15 | k=20 | k=25 | k=30 | k=40 | k=50 |
|-------|-----|------|------|------|------|------|------|------|
| UMAP 2D | 0.415 | 0.416 | 0.462 | 0.464 | 0.451 | 0.476 | 0.522 | 0.525 |
| Full 320D | 0.060 | 0.070 | 0.083 | 0.087 | 0.072 | 0.073 | 0.072 | 0.079 |

The UMAP scores are inflated — UMAP compresses the data into tight local clusters by design, so points naturally appear well-separated in 2D. This is an artifact of the dimensionality reduction, not a property of the underlying data.

The full 320D scores are more honest. They peak at k=20 (silhouette=0.087) and flatten afterwards. This suggests the natural number of clusters in the embedding space is around 20. The low absolute values (0.06-0.09) indicate that protein families do not form perfectly spherical, well-separated clusters — they have fuzzy boundaries and overlapping functions, which is biologically expected.

**Lesson:** internal validation metrics should always be computed on the original embedding space, not on the dimensionality-reduced projection.

## HDBSCAN parameter sweep

HDBSCAN was run with min_cluster_size = 2, 3, 5, 10 on both the full 320D embeddings and the UMAP 2D projection.

On full 320D embeddings, HDBSCAN failed completely — labelling 99.8-100% of proteins as noise regardless of min_cluster_size. This is the curse of dimensionality: in 320 dimensions, all pairwise distances become similar, so HDBSCAN cannot identify meaningful density differences.

On UMAP 2D embeddings, HDBSCAN found meaningful structure:

| min_cluster_size | clusters found | noise points |
|-----------------|---------------|-------------|
| 2 | 82 | 59 (11.8%) |
| 3 | 58 | 74 (14.8%) |
| 5 | 26 | 118 (23.6%) |
| 10 | 2 | 4 (0.8%) |

min_cluster_size=5 was selected as the best setting — it finds 26 clusters, consistent with the k-means elbow estimate of k=20, and labels 23.6% of proteins as noise. This noise fraction is biologically plausible: many proteins are genuinely unusual and do not belong to a well-defined family.

## ARI comparison — all three algorithms

Adjusted Rand Index was computed against known protein family labels for each algorithm at its best settings.

| Method | k | ARI |
|--------|---|-----|
| K-means (k known) | 29 | 0.494 |
| K-means (k from elbow) | 26 | 0.463 |
| HDBSCAN (automatic) | 26 | 0.445 |
| Hierarchical (Ward) | 30 | 0.433 |

All three methods perform similarly, with ARI values between 0.43 and 0.49. This consistency across fundamentally different algorithms confirms that the ESM embeddings contain genuine biological structure — the clustering results are not an artifact of any single method's assumptions.

K-means with k=29 scores highest, but this comparison is unfair — it was given the true number of families. When k-means uses the same k=26 that HDBSCAN found automatically, the gap narrows to 0.463 vs 0.445.

## Stability analysis

Each algorithm was run 10 times on 80% subsamples of the data. Pairwise ARI between runs measures reproducibility.

| Method | Mean stability ARI | Std |
|--------|-------------------|-----|
| K-means (k=26) | 0.697 | 0.061 |
| HDBSCAN (min_size=5) | 0.417 | 0.209 |

K-means is significantly more stable. This is expected — k-means always produces exactly k clusters regardless of which proteins are in the sample. HDBSCAN is sensitive to the density structure of the subsample: removing 20% of proteins can cause clusters to appear, disappear, or merge, leading to inconsistent label assignments across runs.

The high standard deviation of HDBSCAN (0.209) means some pairs of runs agree moderately well while others are nearly random. This instability is a practical limitation when reproducibility matters.

Stability and biological coherence are different things. HDBSCAN may find more meaningful clusters on the full dataset, but k-means is more reproducible across subsamples. For a production similarity search system, k-means is preferable. For exploratory analysis where you do not know k, HDBSCAN is the right choice.

## GO term enrichment validation

GO term enrichment was computed for each of the 26 k-means clusters using a hypergeometric test. All 26 clusters showed at least one significantly enriched GO term (p < 0.05), and most were highly significant (p < 0.001).

Selected biologically meaningful clusters:

| Cluster | Size | Top GO term | Biological meaning | p-value |
|---------|------|-------------|-------------------|---------|
| 3 | 9 | GO:0004930 | G-protein coupled receptor activity | 3.99e-12 |
| 5 | 4 | GO:0019814 | Immunoglobulin complex | 3.89e-10 |
| 7 | 5 | GO:0016020 | Membrane | 2.04e-04 |
| 9 | 7 | GO:0000786 | Nucleosome | 4.92e-10 |
| 16 | 30 | GO:0000981 | DNA-binding transcription factor activity | 6.27e-16 |
| 22 | 35 | GO:0005524 | ATP binding | 6.36e-12 |
| 23 | 19 | GO:0031012 | Extracellular matrix | 6.45e-10 |

These are biologically coherent functional categories that domain experts would recognise. Cluster 5 is particularly striking — all 4 proteins are immunoglobulins, forming a perfectly pure cluster despite the algorithm having no knowledge of protein families during clustering.

Cluster 3 (GPCRs) is the most consistently distinct group across all analyses — visible as an isolated island in every dimensionality reduction method, the first branch to separate in the dendrogram, and the cluster with one of the strongest enrichment signals.

## Key conclusions

ESM embeddings capture genuine biological structure. Three fundamentally different clustering algorithms — k-means, HDBSCAN, and hierarchical — all find similar groupings with ARI values between 0.43 and 0.49, and GO enrichment confirms these groupings correspond to real functional categories.

The natural cluster count in this dataset is approximately 20-26, consistent across the k-means elbow, HDBSCAN automatic detection, and hierarchical dendrogram cuts.

K-means is recommended for reproducibility and interpretability. HDBSCAN is recommended for exploratory work when k is unknown, despite its lower stability on subsamples.

The most important limitation of this analysis is that HDBSCAN requires a low-dimensional projection (UMAP) to function, while k-means and hierarchical clustering can operate directly on the full 320D embeddings. Evaluating all methods on the same embedding space is therefore not straightforward.