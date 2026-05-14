from src.search import  nearest_neighbours

from Bio.Blast import NCBIWWW, NCBIXML
def precision_at_10(query_idx, embeddings, df, metric, n=10):
    """
    Compute precision@n for embedding-based nearest neighbour retrieval.

    A neighbour is considered a hit if it shares at least one GO term with the query.

    Args:
        query_idx: Row index of the query protein in df and embeddings.
        embeddings: 2D array of all protein embeddings.
        df: DataFrame containing a 'go_list' column with per-protein GO terms.
        metric: Distance metric passed to nearest_neighbours ('cosine', 'euclidean', or 'manhattan').
        n: Number of neighbours to evaluate. Defaults to 10.

    Returns:
        Precision score as a float between 0 and 1.
    """
    from sklearn.metrics.pairwise import cosine_similarity, euclidean_distances, manhattan_distances
    query_embedding = embeddings[query_idx]
    query_go = set(df.iloc[query_idx]["go_list"])
    
    _, neighbours_indexes = nearest_neighbours(query_embedding, embeddings, df, n=10, metric=metric)
    hits = 0
    for idx in neighbours_indexes:
        neighbour_go = set(df.loc[idx]["go_list"])
        if len(query_go & neighbour_go) > 0:
            hits += 1
    
    return hits / n
import re

def extract_uniprot_id(title):
    """
    Extract a UniProt accession ID from a BLAST alignment title.

    Args:
        title: BLAST alignment title string (e.g. 'sp|P01308.1|INS_HUMAN Insulin...').

    Returns:
        UniProt accession string, or None if not found.
    """
    match = re.search(r'sp\|(\w+)\.\d+\|', title)
    if match:
        return match.group(1)
    return None
def blast_precision(sequence, query_go, query_uniprot_id, df_clean):
    """
    Compute precision for BLAST-based retrieval against SwissProt.

    A hit is considered correct if it shares at least one GO term with the query.

    Args:
        sequence: Amino acid sequence string of the query protein.
        query_go: Set of GO terms for the query protein.
        query_uniprot_id: UniProt ID of the query protein (excluded from hits).
        df_clean: DataFrame with 'uniprot_id' and 'go_list' columns for GO term lookup.

    Returns:
        Tuple of (precision, n_hits) where precision is a float between 0 and 1
        and n_hits is the number of BLAST hits found. Returns (None, 0) if no hits.
    """
    
    result_handle = NCBIWWW.qblast("blastp", "swissprot", sequence,
                                    entrez_query="Homo sapiens[organism]")
    blast_records = NCBIXML.parse(result_handle)
    record = next(blast_records)
    
    hit_ids = []
    for alignment in record.alignments:
        uid = extract_uniprot_id(alignment.title)
        if uid and uid != query_uniprot_id:
            hit_ids.append(uid)
    
    if len(hit_ids) == 0:
        return None, 0
    
    # fetch GO terms for hits
    hits_go = 0
    for uid in hit_ids:
        hit_row = df_clean[df_clean["uniprot_id"] == uid]
        if len(hit_row) > 0:
            hit_go = set(hit_row["go_list"].iloc[0])
            if len(query_go & hit_go) > 0:
                hits_go += 1
    
    return hits_go / len(hit_ids), len(hit_ids)