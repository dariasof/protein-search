# Dimensionality reduction analysis
 
## What each method optimises
 
PCA finds the directions of maximum variance in the data through linear projection. It is deterministic, fast, and interpretable — each principal component is a linear combination of the original 320 embedding dimensions. The limitation is the linearity assumption: if protein families form curved or non-linear manifolds in embedding space, PCA cannot capture that structure.
 
t-SNE minimises the divergence between pairwise similarity distributions in high and low dimensions. It prioritises preserving local neighbourhood structure — points that are close in 320D stay close in 2D. The tradeoff is that global structure is not preserved: distances between clusters are not meaningful, and the layout changes with the perplexity hyperparameter. t-SNE is non-deterministic without a fixed random seed.
 
UMAP constructs a fuzzy topological representation of the high-dimensional data and optimises a low-dimensional embedding to match it. It preserves both local neighbourhood structure (like t-SNE) and some global structure (unlike t-SNE). It is faster than t-SNE and more stable across runs. The two key hyperparameters are n_neighbors (controls local vs global structure) and min_dist (controls how tightly points pack together).
 
## ARI comparison
 
Adjusted Rand Index was computed for each method by running k-means (k=29) on the 2D projection and comparing cluster assignments to true protein family labels. Higher ARI means the 2D projection better separates protein families.
 
| Method | ARI |
|--------|-----|
| PCA (2D) | 0.166 |
| t-SNE (perplexity=10) | 0.670 |
| t-SNE (perplexity=30) | 0.675 |
| t-SNE (perplexity=50) | 0.617 |
| UMAP (n_neighbors=3) | 0.566 |
| UMAP (n_neighbors=4) | 0.633 |
| UMAP (n_neighbors=10) | 0.611 |
 
## Why PCA scores so much lower
 
PCA scored 0.166 — four times lower than the best non-linear methods. This is expected: PCA assumes that protein family differences align with the directions of maximum variance in a linear sense. In practice, protein families form non-linear manifolds in embedding space. A linear projection into 2D collapses much of this structure.
 
Visually, PCA clearly separated only the most structurally distinct family — G-protein coupled receptors (GPCRs). 7 out of 12 proteins in the right cluster (PC1 > 1) belong to the GPCR family. GPCRs all have 7 transmembrane helices, a structural motif so distinctive it creates a strong linear signal in the high-variance directions. Intermediate filaments grouped partially in the bottom-left, but other families overlapped heavily in the centre of the PCA plot.
 
This confirms that linear dimensionality reduction is insufficient for protein embedding space — the biologically meaningful structure is non-linear.
 
## Why t-SNE slightly beats UMAP on ARI
 
t-SNE (perplexity=30) achieved ARI=0.675, marginally higher than the best UMAP configuration (n_neighbors=4, ARI=0.633). However, this does not mean t-SNE is a better method for this project.
 
t-SNE creates very tight, well-separated local clusters — this is what k-means detects, and why ARI is higher. But t-SNE distances between clusters are not meaningful: the global layout is arbitrary and changes with different perplexity values. Running t-SNE twice with different random seeds produces different global arrangements.
 
UMAP preserves more global structure — the relative positions of clusters carry information about which protein families are broadly similar. This matters for downstream tasks like similarity search and clustering on the full embedding space.
 
## Hyperparameter effects
 
For t-SNE, perplexity controls how many neighbours each point considers. Perplexity=30 gave the best ARI. At perplexity=10, clusters were tighter and more isolated; at perplexity=50, structure was more compressed and less separated.
 
For UMAP, n_neighbors controls the balance between local and global structure. Very low n_neighbors (n=3 or 4) caused extreme separation — GPCRs formed a completely isolated island far from all other proteins. This looks dramatic but reflects overfitting to local structure. n_neighbors=10 gave more balanced results.
 
## Intrinsic dimensionality
 
PCA was also run with all 320 components to measure the intrinsic dimensionality of the embedding space. The cumulative explained variance analysis showed:
 
- 49 dimensions explain 90% of variance
- 76 dimensions explain 95% of variance
- The remaining 244 dimensions contribute only 5% of variance
 
This means the ESM model encodes protein information compactly — most biological variation lives in approximately 50 directions, not 320. The first principal component alone explains ~18% of variance.
 
This finding has a direct practical consequence for clustering. Running k-means on PCA-50D embeddings (which retain 90% of signal) outperformed both full 320D and PCA-2D:
 
| Space | ARI (k-means, k=29) |
|-------|---------------------|
| PCA 2D | 0.204 |
| Full 320D | 0.381 |
| PCA 50D | 0.405 |
| UMAP 2D | 0.494 |
 
PCA-50D outperforms full 320D because dimensions 50-320 add noise rather than signal — they distort distance calculations in k-means without contributing meaningful information. This is a well-known effect called noise reduction through dimensionality reduction.
 
## Which method to use going forward
 
UMAP with n_neighbors=4 was selected as the primary 2D projection for visualisation and clustering throughout the project. The reasons:
 
- Higher ARI than PCA (0.633 vs 0.166) — better family separation
- More stable than t-SNE — global layout is consistent across runs
- Global structure is preserved — cluster positions relative to each other are meaningful
- Faster than t-SNE on this dataset size
 
For clustering and similarity search, PCA-50D embeddings are used instead of full 320D, as they retain 90% of variance while removing noisy dimensions and improving k-means performance.
 
t-SNE is retained for visual comparison but not used as the primary projection, because its global layout is arbitrary and its distance metric is not meaningful for downstream analysis.