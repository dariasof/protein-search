from Bio.Blast import NCBIWWW, NCBIXML

def blast_search(sequence, n=10):
    """
    Run a remote BLASTp search against SwissProt for a given protein sequence.

    Args:
        sequence: Amino acid sequence string to query.
        n: Number of top hits to return. Defaults to 10.

    Returns:
        List of alignment title strings for the top n hits.
    """
    result_handle = NCBIWWW.qblast("blastp", "swissprot", sequence, entrez_query="Homo sapiens[organism]")
    blast_records = NCBIXML.parse(result_handle)
    record = next(blast_records)
    
    hits = []
    for alignment in record.alignments[:n]:
        # extract UniProt ID from alignment title
        title = alignment.title
        hits.append(title)
    return hits