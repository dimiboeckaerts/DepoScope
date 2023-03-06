"""
Created on 03/03/23

@author: dimiboeckaerts

PhageDEPOdetection DATA PROCESSING
"""

# 0 - LIBRARIES
# --------------------------------------------------------------------------------


# 1 - FUNCTIONS
# --------------------------------------------------------------------------------


# 2 - MAIN
# --------------------------------------------------------------------------------
"""
Raw data:
    - interpro domains
    - protein sequences from MillardLab phage genomes that had a hit in the interpro domains
    - protein annotations linked to each of the interpro domains

First filter: 
    - remove all interpro domains with not a single congruent phage depo annotation (manual)
    - add interpro domains to remove to a list, add questionable domains to a second list
    - remove the sequences linked to the interpro domains to remove
    - doublecheck the sequences linked to the questionable domains
"""

to_remove = ['IPR000490', 'IPR000852', 'IPR001088', 'IPR001137', 'IPR001547', 'IPR004185', 'IPR004888', 
             'IPR006048', 'IPR006425', 'IPR008902', 'IPR011496']
ambiguous = ['IPR000165', 'IPR000757', 'IPR000922', 'IPR001139', 'IPR001329', 'IPR001371', 'IPR001439',
             'IPR001554', 'IPR005199', 'IPR006065', 'IPR007724', 'IPR007781', 'IPR008291', 'IPR008929',
             'IPR010702', 'IPR010905', 'IPR011613']