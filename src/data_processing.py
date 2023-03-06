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
"""
absolute_path = '/Users/dimi/Documents/GitHub/PhageDEPOdetection/'
interpro_domains_path = ...
protein_fasta_path = absolute_path+'data/millard_depos.fasta'

# 1 - LIBRARIES
# --------------------------------------------------------------------------------
import numpy as np
import pandas as pd
from Bio import SeqIO
import json
import ast

# 2 - FUNCTIONS
# --------------------------------------------------------------------------------


# 3 - MAIN
# --------------------------------------------------------------------------------
"""
First filter:
    - remove all interpro domains with not a single congruent phage depo annotation (manual)
    - add interpro domains to remove to a list, add questionable domains to a second list
    - remove the sequences linked to the interpro domains to remove
    - doublecheck the sequences linked to the questionable domains
"""

# Read the data
with open(interpro_domains_path) as f:
    data = f.read()
    # data = json.load(f)
interpro_dict = ast.literal_eval(data)

# read the fasta file with BioPython and convert it to a dict {protein_id: sequence}
fasta_dict = {}
for record in SeqIO.parse(protein_fasta_path,'fasta'):
    name, sequence = record.id, str(record.seq)
    fasta_dict[name] = sequence

# define the interpro domains to remove
to_remove = ['IPR000490', 'IPR000852', 'IPR001088', 'IPR001137', 'IPR001547', 'IPR004185', 'IPR004888', 
             'IPR006048', 'IPR006425', 'IPR008902', 'IPR011496']

# loop over the interpro domains to remove and remove the corresponding sequences in the fasta dict
# and the keys in the interpro dict
for domain in to_remove:
    for protein in interpro_dict[domain]:
        fasta_dict.pop(protein)
    interpro_dict.pop(domain)

# define the interpro domains to check
ambiguous = ['IPR000165', 'IPR000757', 'IPR000922', 'IPR001139', 'IPR001329', 'IPR001371', 'IPR001439',
             'IPR001554', 'IPR005199', 'IPR006065', 'IPR007724', 'IPR007781', 'IPR008291', 'IPR008929',
             'IPR010702', 'IPR010905', 'IPR011613']

# collect all of the sequences to check
fasta_dict_sub = {}
for domain in ambiguous:
    for protein in interpro_dict[domain]:
        fasta_dict_sub[protein] = fasta_dict[protein]

# cluster the sequences to check with affinity propagation
