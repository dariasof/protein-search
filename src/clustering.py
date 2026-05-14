import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
import hdbscan as hdbscan_lib


def kmeans_clustering(embeddings, k, max_iters=100, random_state=42):
    """
    Run k-means clustering on a set of embeddings.

    Args:
        embeddings: 2D array of shape (n_samples, n_features).
        k: Number of clusters.
        max_iters: Maximum number of iterations before stopping. Defaults to 100.
        random_state: Seed for reproducible centroid initialisation. Defaults to 42.

    Returns:
        closest_centres: 1D array of cluster assignments for each sample.
        centres: 2D array of shape (k, n_features) with final cluster centroids.
    """
    np.random.seed(random_state)
    n_samples = embeddings.shape[0]
    centre_indexes = np.random.choice(n_samples, k, replace=False)
    centres = embeddings[centre_indexes]
    for _ in range(max_iters):
        old_centres = centres.copy()
        differences = embeddings[:, np.newaxis, :] - centres
        distances = np.sqrt(np.sum(differences ** 2, axis=2))
        closest_centres = np.argmin(distances, axis=1)
        for i in range(k):
            if np.any(closest_centres == i):
                centres[i] = embeddings[closest_centres == i].mean(axis=0)
        if np.allclose(old_centres, centres):
            break
    return closest_centres, centres


def hdbscan_clustering(embeddings, min_cluster_size=5):
    """
    Run HDBSCAN density-based clustering on a set of embeddings.

    Args:
        embeddings: 2D array of shape (n_samples, n_features).
        min_cluster_size: Minimum number of points to form a cluster. Defaults to 5.

    Returns:
        labels: 1D array of cluster assignments. Points labelled -1 are noise.
    """
    clusterer = hdbscan_lib.HDBSCAN(min_cluster_size=min_cluster_size)
    return clusterer.fit_predict(embeddings)


def hierarchical_clustering(embeddings, n_clusters, method='ward'):
    """
    Run hierarchical (agglomerative) clustering on a set of embeddings.

    Args:
        embeddings: 2D array of shape (n_samples, n_features).
        n_clusters: Number of clusters to cut the dendrogram at.
        method: Linkage method (e.g. 'ward', 'complete', 'average'). Defaults to 'ward'.

    Returns:
        labels: 1D array of cluster assignments (1-indexed, as returned by fcluster).
        linkage_matrix: The linkage matrix from scipy, useful for plotting dendrograms.
    """
    linkage_matrix = linkage(embeddings, method=method)
    labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
    return labels, linkage_matrix
