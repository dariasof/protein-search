import numpy as np
from sklearn.metrics import adjusted_rand_score, silhouette_score, davies_bouldin_score
from src.search import nearest_neighbours


def ari(true_labels, predicted_labels):
    """
    Compute Adjusted Rand Index between true and predicted cluster assignments.

    Args:
        true_labels: 1D array of ground-truth labels.
        predicted_labels: 1D array of predicted cluster assignments.

    Returns:
        ARI score as a float in [-1, 1]. Higher is better; 1.0 is perfect agreement.
    """
    return adjusted_rand_score(true_labels, predicted_labels)


def silhouette(embeddings, labels):
    """
    Compute mean silhouette score for a clustering.

    Args:
        embeddings: 2D array of shape (n_samples, n_features).
        labels: 1D array of cluster assignments. Points labelled -1 (noise) are excluded.

    Returns:
        Silhouette score as a float in [-1, 1]. Higher is better.
    """
    mask = labels != -1
    if mask.sum() < 2:
        return float("nan")
    return silhouette_score(embeddings[mask], labels[mask])


def davies_bouldin(embeddings, labels):
    """
    Compute Davies-Bouldin index for a clustering.

    Args:
        embeddings: 2D array of shape (n_samples, n_features).
        labels: 1D array of cluster assignments. Points labelled -1 (noise) are excluded.

    Returns:
        Davies-Bouldin score as a float. Lower is better; 0 is perfect separation.
    """
    mask = labels != -1
    if mask.sum() < 2:
        return float("nan")
    return davies_bouldin_score(embeddings[mask], labels[mask])


def precision_at_k(query_idx, embeddings, df, metric="cosine", k=10):
    """
    Compute precision@k for embedding-based nearest neighbour retrieval.

    A neighbour is a hit if it shares at least one GO term with the query.

    Args:
        query_idx: Row index of the query protein in df and embeddings.
        embeddings: 2D array of all protein embeddings.
        df: DataFrame containing a 'go_list' column with per-protein GO terms.
        metric: Distance metric passed to nearest_neighbours. Defaults to 'cosine'.
        k: Number of neighbours to evaluate. Defaults to 10.

    Returns:
        Precision score as a float between 0 and 1.
    """
    query_embedding = embeddings[query_idx]
    query_go = set(df.iloc[query_idx]["go_list"])

    _, neighbour_indexes = nearest_neighbours(query_embedding, embeddings, df, n=k, metric=metric)

    hits = sum(
        1 for idx in neighbour_indexes
        if len(query_go & set(df.loc[idx]["go_list"])) > 0
    )
    return hits / k
