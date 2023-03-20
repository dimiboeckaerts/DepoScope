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
# data paths
absolute_path = '/Users/dimi/Documents/GitHub/PhageDEPOdetection/'
protein_domains_path = absolute_path+'data/dbsuite_results.v2.json'
protein_df_path = absolute_path+'data/df_sequences.index.v2.csv'
annotations_path = absolute_path+'data/proteinID_annotation.v2.json'
protein_embeddings_path = absolute_path+'data/embeddings.proteins.v2.csv'

# domains and annotations to remove or keep
domains_to_remove = ['IPR000490', 'IPR000852', 'IPR001088', 'IPR001137', 
                     'IPR001547', 'IPR004185', 'IPR004888', 'IPR006048', 
                     'IPR006425', 'IPR008902', 'IPR011496']
domains_to_check = ['IPR000165', 'IPR000757', 'IPR000922', 'IPR001139', 
                    'IPR001329', 'IPR001371', 'IPR001439', 'IPR001554', 
                    'IPR005199', 'IPR006065', 'IPR007724', 'IPR007781', 
                    'IPR008291', 'IPR008929', 'IPR010702', 'IPR010905', 
                    'IPR011613']
annotations_to_remove = ['ribonucleoside', 'diphosphate', 'reductase', 'endolysin',
                         'RecA', 'UvsX-like', 'DNA helicase', 'Hoc-like head',
                         'epimerase', 'DNA repair', 'DNA annealing']

# 1 - LIBRARIES
# --------------------------------------------------------------------------------
import json
import numpy as np
import pandas as pd
from Bio import SeqIO
from sklearn.cluster import AffinityPropagation

# 2 - MAIN
# --------------------------------------------------------------------------------
"""
STEP 1: filter at the domain level
    - remove all interpro domains with not a single congruent phage depo annotation (manual)
    - add interpro domains to remove to a list, add questionable domains to a second list
    - remove the sequences linked to the interpro domains to remove
"""
print('Reading the data...')
# Read the data
with open(protein_domains_path) as f1:
    domains_dict = json.load(f1)
with open(annotations_path) as f2:
    annotations_dict = json.load(f2)
embeddings = pd.read_csv(protein_embeddings_path, header=None).iloc[:,:-1]
seq_df = pd.read_csv(protein_df_path, sep='\t', header=None)

# loop over the domains_dict and collect the protein IDs to remove or to check
proteins_to_remove = []
proteins_to_check = []
for proteinid in domains_dict.keys():
    domains_list = [domain.split('__')[0] for domain in domains_dict[proteinid]]
    if any([domain in domains_list for domain in domains_to_remove]):
        proteins_to_remove.append(proteinid)
    elif any([domain in domains_list for domain in domains_to_check]):
        proteins_to_check.append(proteinid)

"""
STEP 2: filter at the individual sequence level
    - doublecheck the sequences linked to the questionable domains
    - do affinity propagation clustering on the embeddings of all protein sequences
    - check the clusters in which questionable proteins are present 
    - get all of the annotations of the proteins in those clusters
    - if any of the annotations of the proteins in the cluster are in the 
        list of annotations to remove, than we remove the questionable protein
"""
print('Clustering with affinity propagation...')
# do affinity propagation with all of the sequences and their embeddings
X = embeddings.iloc[:, 1:]
af = AffinityPropagation(damping=0.90, preference=None, random_state=123, max_iter=1000,verbose=True).fit(X)
cluster_centers_indices = af.cluster_centers_indices_
cluster_ids = af.labels_

# loop over the proteins in the list of proteins to check
for protein in proteins_to_check:
    # retrieve the embeddings indices (column 0 in seq_df) of the proteins to check
    protein_id = int(seq_df.loc[seq_df.iloc[:,1]== protein, 0])

    # retrieve the cluster_ids using these embeddings indices
    embedding_index = list(embeddings.iloc[:,0]).index(protein_id)
    cluster_id = cluster_ids[embedding_index]

    # retrieve all of the protein_ids that are in that same cluster
    cluster_protein_ids = embeddings.loc[cluster_ids == cluster_id,0]

    # get the protein names for the proteins in the cluster
    protein_names = [list(seq_df.loc[seq_df.iloc[:,0]== this_id, 1]) for this_id in cluster_protein_ids]
    protein_names = [item for sublist in protein_names for item in sublist] # flatten the list

    # get the annotations for the proteins in the cluster
    annotations = [annotations_dict[protein] for protein in protein_names]

    # check if any of the annotations are in the list of annotations to remove
    # by checking the words of annotations_to_remove in each of the indiv annotations
    if any([annotation_rm in anno for anno in annotations for annotation_rm 
            in annotations_to_remove]):
        proteins_to_remove.append(protein)
    # if only 'hypothetical' or 'unknown' is present as annotation in the cluster, remove the protein
    elif all([(anno == 'hypothetical protein' or anno == 'unknown function') for anno in annotations]):
        proteins_to_remove.append(protein)

"""
STEP 3: constructing the final database

IDs - Interpro entry - sequences - domains - annotations
"""
final_ids = []
final_IPRs = []
final_sequences = []
final_annotations = []
for proteinid in domains_dict.keys():
    if proteinid not in proteins_to_remove:
        final_ids.append(proteinid)
        final_IPRs.append(domains_dict[proteinid])
        final_sequence = seq_df.loc[seq_df.iloc[:,1]== proteinid, 2].values[0]
        final_sequences.append(final_sequence)
        final_annotations.append(annotations_dict[proteinid])

final_database = {'protein_ID': final_ids, 'IPRs': final_IPRs, 
                  'sequence': final_sequences, 'annotation': final_annotations}

"""
STEP 4: Annotate domains and save dataframe for further use

Ideally, we do this with ResDom, but let's just keep it simple for now and follow the
simple rule that if a protein < 250 AAs, it is considered as a depo domain in its entirety,
and if it is larger than 250 AAs, the C-terminal part is considered as the depo domain.
"""
# identify domains with simple 200 AA cutoff
print('Annotating the domains...')
depo_domains = []
token_labels = []
count = 0
for sequence in final_database['sequence']:
    if len(sequence) < 250:
        depo_domains.append(sequence)
        token_labels.append([1]*len(sequence))
    else:
        depo_domains.append(sequence[200:])
        token_labels.append([0]*200 + [1]*(len(sequence)-200))

final_database['depo_domain'] = depo_domains
final_database['token_labels'] = token_labels
dbdf = pd.DataFrame(final_database)

# save the dataframe as json
print('Saving the final database...')
final_database_path = absolute_path+'data/final_database.json'
with open(final_database_path, 'w') as f:
    json.dump(final_database, f)

# save the database in format for finetuning
"""
https://huggingface.co/docs/datasets/loading#specify-features
{"version": "0.1.0",
 "train": [{"id": 1, "tokens":['M', 'L', 'P', ...], "labels": [0, 0, 1, 1, ...]},
          {"id": 322, "tokens": [...], "labels": [...]}]
 "test": ...
}
"""
max_length = 1024 # max length for ESM2
finetune_data = {'version': '0.1.0', 'train': [], 'test': []}
dbdf = dbdf.sample(frac=1).reset_index(drop=True) # shuffle
train_cutoff = int(dbdf.shape[0]*0.75)
for i, pid in enumerate(dbdf['protein_ID']):
    tokens = list(dbdf['sequence'][i])[:max_length]
    labels = dbdf['token_labels'][i][:max_length]
    this_data = {'id': pid, 'tokens': tokens, 'labels': labels}
    if i <= train_cutoff:
        finetune_data['train'].append(this_data)
    else:
        finetune_data['test'].append(this_data)

finetune_data_path = absolute_path+'data/finetune_data.json'
with open(finetune_data_path, 'w') as f:
    json.dump(finetune_data, f)

"""
Ideas & Remarks:
    - pre processing: get everything in one dataframe already so we only have
    one object to work with (i.e. annotation & linked domains)
    - better decision filter 2: keep the ambiguous individual sequences if they appear in a cluster in 
    which approved interpro domains are present?!
    - for now, the reclustering is skipped. Looking at all the annotation within one cluster,
    it seems that the clusters are already quite specific, thus might be overkill to recluster?
    - due to the max_length cutoff, we're actually losing a part of the data.
    we could also add the remaining part as an extra data point to avoid throuwing away data
"""
