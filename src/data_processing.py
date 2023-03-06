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
    - interpro domains & linked protein IDs (dict)
    - protein ESM2 embeddings for affinity propagation clustering
"""
absolute_path = '/Users/dimi/Documents/GitHub/PhageDEPOdetection/'
interpro_domains_path = ...
protein_fasta_path = absolute_path+'data/millard_depos.fasta'
protein_embeddings_path = absolute_path+'data/protein_embeddings.csv'

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
with open(interpro_domains_path) as f:
    data = f.read()
    # data = json.load(f)
interpro_dict = ast.literal_eval(data)

fasta_dict = {}
for record in SeqIO.parse(protein_fasta_path,'fasta'):
    name, sequence = record.id, str(record.seq)
    fasta_dict[name] = sequence

embeddings = pd.read_csv(protein_embeddings_path, header=None).iloc[:,:-1]

# define the interpro domains to remove
to_remove = ['IPR000490', 'IPR000852', 'IPR001088', 'IPR001137', 'IPR001547', 'IPR004185', 'IPR004888', 
             'IPR006048', 'IPR006425', 'IPR008902', 'IPR011496']

# loop over the interpro domains to remove and remove the corresponding sequences in the fasta dict
# and the keys in the interpro dict
for domain in to_remove:
    for protein in interpro_dict[domain]:
        fasta_dict.pop(protein)
    interpro_dict.pop(domain)

"""
Second filter: at the individual sequence level
    - doublecheck the sequences linked to the questionable domains
    - do affinity propagation clustering on the embeddings of all protein sequences
    - check the annotations of the quproteins in each of those clusters
    - if any of the annotations of the proteins in the cluster are in the list of annotations to keep, keep the cluster
"""
# define the interpro domains to check
ambiguous = ['IPR000165', 'IPR000757', 'IPR000922', 'IPR001139', 'IPR001329', 'IPR001371', 'IPR001439',
             'IPR001554', 'IPR005199', 'IPR006065', 'IPR007724', 'IPR007781', 'IPR008291', 'IPR008929',
             'IPR010702', 'IPR010905', 'IPR011613']

# collect all of the sequences to check
fasta_dict_sub = {}
for domain in ambiguous:
    for protein in interpro_dict[domain]:
        fasta_dict_sub[protein] = fasta_dict[protein]

# do affinity propagation with all of the sequences and their embeddings
X = embeddings.iloc[:, 1:]
af = AffinityPropagation(damping=0.90, preference=None, random_state=123, max_iter=1000,verbose=True).fit(X)
cluster_centers_indices = af.cluster_centers_indices_
cluster_ids = af.labels_

# for each protein to check, retrieve the cluster_id it belongs to
cluster_dict = {}
for protein in fasta_dict_sub.keys():
    protein_index = list(fasta_dict.keys()).index(protein)
    embeddings_index = list(embeddings.iloc[:, 0]).index(protein_index)
    cluster_dict[protein] = cluster_ids[protein_index]

# for each cluster, check the annotations of the proteins in the cluster
to_keep = []
to_keep_annotations = [...] # TO FILL IN
for (key, cluster) in cluster_dict.items():
    # retrieve embeddings indices of the proteins in the cluster
    embeddings_indices = [embeddings.iloc[i,0] for i, x in enumerate(cluster_ids) if x == cluster]
    # retrieve the protein names of the proteins in the cluster
    annotations = [list(fasta_dict.keys())[i].split('__')[1] for i in embeddings_indices]
    # check if any of the annotations of the proteins in the cluster are in the list of annotations to keep
    if any([item in to_keep_annotations for item in annotations]):
        to_keep.append(key)

# better decision: keep the ambiguous individual sequences if they appear in a cluster in
# which approved interpro domains are present?!

# remove the sequences in fasta_dict that are not in to_keep
for protein in fasta_dict.keys():
    if protein not in to_keep:
        fasta_dict.pop(protein)
