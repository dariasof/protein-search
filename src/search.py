import numpy as np

VALID_METRICS = ("cosine", "euclidean", "manhattan")

def nearest_neighbours(query_embedding, all_embeddings, df, n, metric="cosine"):
    """
    Find the n nearest neighbours to a query embedding.

    Args:
        query_embedding: 1D array of the query protein's embedding.
        all_embeddings: 2D array of all protein embeddings to search against.
        df: DataFrame with protein metadata (must contain 'protein_name' and 'family' columns).
        n: Number of neighbours to return.
        metric: Distance metric to use. One of 'cosine', 'euclidean', or 'manhattan'. Defaults to 'cosine'.

    Returns:
        neighbour_names: DataFrame with 'protein_name' and 'family' of the n nearest neighbours.
        neighbour_indexes: Array of row indices of the neighbours in df.

    Raises:
        ValueError: If metric is not one of the supported options.
    """
    if metric not in VALID_METRICS:
        raise ValueError(f"Invalid metric {metric!r}. Must be one of: {VALID_METRICS}")

    from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances, manhattan_distances

    if metric == "cosine":
        similarities = cosine_similarity([query_embedding], all_embeddings)[0]
        neighbour_indexes = np.argsort(similarities)[::-1][1:n+1] # index of the most similar protein (excluding itself)
    elif metric == "euclidean":
        distances = euclidean_distances([query_embedding], all_embeddings)[0]
        neighbour_indexes = np.argsort(distances)[1:n+1] # index of the closest protein (excluding itself)
    elif metric == "manhattan":
        distances = manhattan_distances([query_embedding], all_embeddings)[0]
        neighbour_indexes = np.argsort(distances)[1:n+1] # index of the closest protein (excluding itself)
    else:
        raise ValueError("Unsupported metric. Use 'cosine' or 'euclidean' or 'manhattan'.")
    
    neighbour_names = df.iloc[neighbour_indexes][["protein_name",'family']]
    return neighbour_names, neighbour_indexes