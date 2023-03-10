"""
Created on 03/03/23

@author: dimiboeckaerts

PhageDEPOdetection DATA PROCESSING
"""

# 0 - INPUTS
# --------------------------------------------------------------------------------
"""
Raw data:
    - protein sequences (FASTA) from MillardLab phage genomes that had a hit in the interpro domains
    - protein IDs & linked interpro domains (dict)
    - protein IDs & annotations (dict)
    - protein ESM2 embeddings for affinity propagation clustering
"""
absolute_path = '/Users/dimi/Documents/GitHub/PhageDEPOdetection/'
protein_domains_path = absolute_path+'data/dbsuite_results.v2.json'
protein_df_path = absolute_path+'data/df_sequences.index.v2.csv'
annotations_path = absolute_path+'data/proteinID_annotation.v2.json'
protein_embeddings_path = absolute_path+'data/embeddings.proteins.v2.csv'

# 1 - LIBRARIES
# --------------------------------------------------------------------------------
import ast
import json
import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.cluster import AffinityPropagation

# 2 - FUNCTIONS
# --------------------------------------------------------------------------------


# 3 - MAIN
# --------------------------------------------------------------------------------
"""
First filter: at the domain level
    - remove all interpro domains with not a single congruent phage depo annotation (manual)
    - add interpro domains to remove to a list, add questionable domains to a second list
    - remove the sequences linked to the interpro domains to remove
"""
# Read the data
with open(protein_domains_path) as f1:
    domains_dict = json.load(f1)
with open(annotations_path) as f2:
    annotations_dict = json.load(f2)
embeddings = pd.read_csv(protein_embeddings_path, header=None).iloc[:,:-1]
seq_df = pd.read_csv(protein_df_path, sep='\t', header=None)

# define the interpro domains to remove
to_remove = ['IPR000490', 'IPR000852', 'IPR001088', 'IPR001137', 'IPR001547', 'IPR004185', 'IPR004888', 
             'IPR006048', 'IPR006425', 'IPR008902', 'IPR011496']
to_check = ['IPR000165', 'IPR000757', 'IPR000922', 'IPR001139', 'IPR001329', 'IPR001371', 'IPR001439',
             'IPR001554', 'IPR005199', 'IPR006065', 'IPR007724', 'IPR007781', 'IPR008291', 'IPR008929',
             'IPR010702', 'IPR010905', 'IPR011613']

# loop over the domains_dict and collect the protein IDs to remove or to check
proteins_to_remove = []
proteins_to_check = []
for proteinid in domains_dict.keys():
    domains_list = [domain.split('__')[0] for domain in domains_dict[proteinid]]
    if any([domain in domains_list for domain in to_remove]):
        proteins_to_remove.append(proteinid)
    elif any([domain in domains_list for domain in to_check]):
        proteins_to_check.append(proteinid)

"""
Second filter: at the individual sequence level
    - doublecheck the sequences linked to the questionable domains
    - do affinity propagation clustering on the embeddings of all protein sequences
    - check the annotations of the quproteins in each of those clusters
    - if any of the annotations of the proteins in the cluster are in the list of annotations to keep, keep the cluster
"""

# do affinity propagation with all of the sequences and their embeddings
X = embeddings.iloc[:, 1:]
af = AffinityPropagation(damping=0.90, preference=None, random_state=123, max_iter=1000,verbose=True).fit(X)
cluster_centers_indices = af.cluster_centers_indices_
cluster_ids = af.labels_

# loop over the proteins in the list of proteins to check
protein_ids_to_check = []
for protein in proteins_to_check:
    # retrieve the embeddings indices (column 0 in seq_df) of the proteins to check
    protein_id = int(seq_df.loc[seq_df.iloc[:,1]== protein, 0])
    protein_ids_to_check.append(int(protein_id))

    # retrieve the cluster_ids using these embeddings indices
    embedding_id = list(embeddings.iloc[:,0]).index(protein_id)
    cluster_id = cluster_ids[embedding_id]

    # retrieve all of the protein_ids that are in that same cluster
    cluster_protein_ids = embeddings.loc[cluster_ids == cluster_id,0]

    # get the protein names for the proteins in the cluster
    # NEEDS TO BE CHECKED! DOES NOT WORK PROPERLY
    protein_names = [str(seq_df.loc[seq_df.iloc[:,0] == this_id, 1]) for this_id in cluster_protein_ids]

    # get the annotations for the proteins in the cluster
    annotations = [annotations_dict[protein] for protein in protein_names]






"""
Ideas:
    - pre processing: get everything in one dataframe already so we only have
    one object to work with (i.e. annotation & linked domains)
    - better decision filter 2: keep the ambiguous individual sequences if they appear in a cluster in 
    which approved interpro domains are present?!

Questions:
    - recluster individual clusters to do what again?
"""