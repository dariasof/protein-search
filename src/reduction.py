import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap import UMAP


def pca(embeddings, n_components=2, random_state=42):
    """
    Project embeddings using PCA.

    Args:
        embeddings: 2D array of shape (n_samples, n_features).
        n_components: Number of output dimensions. Defaults to 2.
        random_state: Random seed for reproducibility. Defaults to 42.

    Returns:
        Tuple of (projected array of shape (n_samples, n_components), fitted PCA object).
    """
    reducer = PCA(n_components=n_components, random_state=random_state)
    return reducer.fit_transform(embeddings), reducer


def tsne(embeddings, n_components=2, perplexity=30, random_state=42):
    """
    Project embeddings using t-SNE.

    Args:
        embeddings: 2D array of shape (n_samples, n_features).
        n_components: Number of output dimensions. Defaults to 2.
        perplexity: t-SNE perplexity — controls local vs global structure trade-off. Defaults to 30.
        random_state: Random seed for reproducibility. Defaults to 42.

    Returns:
        Projected array of shape (n_samples, n_components).
    """
    reducer = TSNE(n_components=n_components, perplexity=perplexity, random_state=random_state)
    return reducer.fit_transform(embeddings)


def umap(embeddings, n_components=2, n_neighbors=4, min_dist=0.1, random_state=42):
    """
    Project embeddings using UMAP.

    Args:
        embeddings: 2D array of shape (n_samples, n_features).
        n_components: Number of output dimensions. Defaults to 2.
        n_neighbors: Controls local vs global structure. Lower values emphasise local structure. Defaults to 4.
        min_dist: Minimum distance between points in the low-dimensional space. Defaults to 0.1.
        random_state: Random seed for reproducibility. Defaults to 42.

    Returns:
        Projected array of shape (n_samples, n_components).
    """
    reducer = UMAP(n_components=n_components, n_neighbors=n_neighbors, min_dist=min_dist, random_state=random_state)
    return reducer.fit_transform(embeddings)
