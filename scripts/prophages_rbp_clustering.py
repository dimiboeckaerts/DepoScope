"""
Created on 03/03/23

@author: dimiboeckaerts

PROPHAGE RBP DETECTION & CLUSTERING
"""

# 0 - LIBRARIES
# --------------------------------------------------
import os
import json
import subprocess
import numpy as np
import pandas as pd


# 1 - FUNCTIONS
# --------------------------------------------------
def cdhit_python(cdhit_path, input_file, output_file, c=0.50, n=3):
    """
    This function executes CD-HIT clustering commands from within Python. To install
    CD-HIT, do so via conda: conda install -c bioconda cd-hit. By default, CD-HIT
    works via a global alignment approach, which is good for our application as
    we cut the sequences to 'one unknown domain' beforehand.
    
    Input:
        - cdhit_path: path to CD-HIT software
        - input_file: FASTA file with protein sequences
        - output file: path to output (will be one FASTA file and one .clstr file)
        - c: threshold on identity for clustering
        - n: word length (3 for thresholds between 0.5 and 0.6)
    """
    
    cd_str = 'cd ' + cdhit_path # change directory
    raw_str = './cd-hit -i ' + input_file + ' -o ' + output_file + ' -c ' + str(c) + ' -n ' + str(n) + ' -d 0'
    command = cd_str+'; '+ raw_str
    #cd_process = subprocess.Popen(cd_str, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #cd_out, cd_err = cd_process.communicate()
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = process.communicate()
    
    return stdout, stderr


# 2 - SCRIPT
# --------------------------------------------------
rbps = pd.read_csv('/Users/dimi/GoogleDrive/PhD/3_PHAGEBASE/32_DATA/RBP_detection/annotated_RBPs_2022-01.csv')
prophage_rbps = pd.read_csv('/Users/dimi/GoogleDrive/PhD/0_INFO_OTHER/2023_02/valencia_compute/xgb_combined_detections.csv')

# drop entries in prophage_rbps that have the label 'nonRBP' in the column 'prediction'
prophage_rbps = prophage_rbps[prophage_rbps['prediction'] != 'nonRBP']
prophage_rbps = prophage_rbps.reset_index(drop=True)

# write one combined fasta file with all of the sequences in
with open('/Users/dimi/GoogleDrive/PhD/0_INFO_OTHER/2023_02/valencia_compute/prophage_annotated_rbps.fasta', 'w') as f:
    for i, sequence in enumerate(prophage_rbps['ProteinSeq']):
        name = 'PROPHAGE_'+prophage_rbps['name'][i]
        f.write('>'+name+'\n'+sequence+'\n')
    for i, sequence in enumerate(rbps['ProteinSeq']):
        name = 'RBP_'+rbps['protein_id'][i]
        f.write('>'+name+'\n'+sequence+'\n')

# run CD-HIT clustering
cdpath = '/Users/dimi/cd-hit-v4.8.1-2019-0228'
input_file = '/Users/dimi/GoogleDrive/PhD/0_INFO_OTHER/2023_02/valencia_compute/prophage_annotated_rbps.fasta'
output_file = '/Users/dimi/GoogleDrive/PhD/0_INFO_OTHER/2023_02/valencia_compute/CD_cluster'
stdout, stderr = cdhit_python(cdpath, input_file, output_file, c=0.50, n=3)

# explore the output
# How big are the clusters?
# How many clusters only contain one sequence?
# How is the division between prophage and annotated RBPs for each of the clusters?
# How diverse are sequences within each cluster and across the clusters? (e.g. pw align all the representatives)
# basically answering the question: how different or similar are prophage RBPs and annotated RBPs?


# 3 - FINDINGS
# --------------------------------------------------
"""
- In total, 9319 detected prophage RBPs and 6176 annotated RBPs were clustered (= 15495 total).
- These resulted in 2409 unique clusters at 50% identity threshold.

"""