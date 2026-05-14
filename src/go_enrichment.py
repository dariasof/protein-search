from scipy.stats import hypergeom
from collections import Counter

def go_enrichment(cluster_indices, df, all_go_counts, top_n=1):
    """
    Perform GO term enrichment analysis on a cluster using the hypergeometric test.

    Args:
        cluster_indices: List of row indices belonging to the cluster.
        df: DataFrame containing a 'go_list' column with per-protein GO terms.
        all_go_counts: Counter or dict mapping each GO term to its count across the full dataset.
        top_n: Number of most enriched GO terms to return. Defaults to 1.

    Returns:
        List of (go_term, k, n, p_value) tuples sorted by p-value (ascending),
        where k is the count in the cluster and n is the count in the full dataset.
    """
    M = len(df)  # total proteins
    N = len(cluster_indices)  # cluster size
    
    cluster_go = Counter()
    for idx in cluster_indices:
        cluster_go.update(df["go_list"].iloc[idx])
    
    results = []
    for go_term, k in cluster_go.items():
        n = all_go_counts[go_term]  # proteins with this term in whole dataset
        p_value = hypergeom.sf(k-1, M, n, N)
        results.append((go_term, k, n, p_value))
    
    results.sort(key=lambda x: x[3])  # sort by p-value
    return results[:top_n]