import os
import requests
from tqdm import tqdm
import pandas as pd
from io import StringIO
import torch
import esm
import numpy as np

def fetch_uniprot(n):
    """
    Fetch reviewed human proteins from UniProt.

    Args:
        n: Number of proteins to fetch.

    Returns:
        DataFrame with columns: uniprot_id, protein_name, sequence, go_terms, family.
        Returns None if the request fails.
    """
    link = f"https://rest.uniprot.org/uniprotkb/search?query=organism_id:9606+reviewed:true&fields=accession,protein_name,sequence,go,protein_families&format=tsv&size={n}"
    response = requests.get(link)
    if response.status_code != 200:
        print("Error:", response.status_code)
        return None
    df = pd.read_csv(StringIO(response.text), sep="\t")
    df.columns = ["uniprot_id", "protein_name", "sequence", "go_terms", "family"]
    return df

def fetch_go_biological_process(n):
    """
    Fetch GO biological process annotations for reviewed human proteins from UniProt.

    Args:
        n: Number of proteins to fetch.

    Returns:
        DataFrame with columns: uniprot_id, go_biological_process.
        Returns None if the request fails.
    """
    url = f"https://rest.uniprot.org/uniprotkb/search?query=organism_id:9606+reviewed:true&fields=accession,go_p&format=tsv&size={n}"
    response = requests.get(url)
    if response.status_code != 200:
        print("Error:", response.status_code)
        return None
    df = pd.read_csv(StringIO(response.text), sep="\t")
    df.columns = ["uniprot_id", "go_biological_process"]
    return df


def compute_embeddings(df, model, alphabet):
    """
    Compute ESM2 embeddings for each protein in the DataFrame.

    Runs each sequence through the model individually and mean-pools the
    per-token representations (excluding BOS/EOS tokens) to produce a
    single fixed-size vector per protein.

    Args:
        df: DataFrame containing 'uniprot_id' and 'sequence' columns.
        model: Loaded ESM2 model (from esm.pretrained).
        alphabet: Corresponding ESM2 alphabet.

    Returns:
        2D numpy array of shape (n_proteins, embedding_dim).
    """
    batch_converter = alphabet.get_batch_converter()
    
    sequence_representations = []
    
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Computing embeddings"):
        data = [(row["uniprot_id"], row["sequence"])]
        batch_labels, batch_strs, batch_tokens = batch_converter(data)
        batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)
        
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[6])
        token_representations = results["representations"][6]

        for i, tokens_len in enumerate(batch_lens):
            sequence_representations.append(
                token_representations[i, 1:tokens_len-1].mean(0)
            )
    
    return torch.stack(sequence_representations).numpy()

def save_data(df, embeddings):
    """
    Save annotations and embeddings to disk.

    Args:
        df: DataFrame of protein annotations to save as data/raw/annotations.csv.
        embeddings: Numpy array of embeddings to save as data/raw/embeddings.npy.

    Returns:
        Confirmation string.
    """
    base = os.path.join(os.path.dirname(__file__), "..")
    df.to_csv(os.path.join(base, "data/raw/annotations.csv"), index=False)
    np.save(os.path.join(base, "data/raw/embeddings.npy"), embeddings)
    return "Data saved successfully!"
def load_data():
    """
    Load saved annotations and embeddings from disk.

    Returns:
        df: DataFrame of protein annotations from data/raw/annotations.csv.
        embeddings: Numpy array of embeddings from data/raw/embeddings.npy.
    """
    base = os.path.join(os.path.dirname(__file__), "..")
    df = pd.read_csv(os.path.join(base, "data/raw/annotations.csv"))
    embeddings = np.load(os.path.join(base, "data/raw/embeddings.npy"))
    return df, embeddings

