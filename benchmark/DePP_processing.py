from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np
import pandas as pd
from math import pow
from itertools import product

# Read a multi-FASTA file and return a list of protein SeqRecords
def read_fasta(fasta_file):
    return list(SeqIO.parse(fasta_file, format = "fasta"))

# Calculate protein parameters using ProteinAnalysis and return a list of parameters
def calculate_protein_parameters(prot_record):
    def prot_params(seq):
        x = ProteinAnalysis(seq)
        aa_percents = x.get_amino_acids_percent()
        params = [
            x.molecular_weight(),
            x.aromaticity(),
            x.instability_index(),
            x.isoelectric_point(),
            *x.secondary_structure_fraction(),
            *x.molar_extinction_coefficient(),
            x.gravy(),
            sum(x.flexibility()) / len(x.flexibility()) if x.flexibility() else 0,
        ]
        return [prot_record.name, *params, *[aa_percents[aa] for aa in 'ACDEFGHIKLMNPQRSTVWY'], x.length]

    def k_mer(seq, k):
        groups = {'A': '1', 'V': '1', 'G': '1', 'I': '1', 'L': '1', 
                  'F': '2', 'W': '2', 'Y': '2','M': '3', 'C': '3', 
                  'H': '4', 'K': '4', 'R': '4', 'D': '5', 'E': '5',
                  'S': '6', 'T': '6', 'N': '6', 'Q': '6', 'P': '7'}

        # Generate all possible combinations of physicochemical groups
        iteratives = [''.join(i) for i in product("1234567", repeat=k)]

        # Initialize a dictionary to count the occurrences of physicochemical combinations
        counts = {combination: 0 for combination in iteratives}

        # Iterate through the sequence and count the occurrences of physicochemical combinations
        for j in range(0, len(seq) - k + 1):
            kmer = seq[j:j + k]
            c = ''.join(groups[aa] for aa in kmer)
            counts[c] += 1

        # Calculate the proportions of physicochemical combinations
        total_combinations = len(seq) - k + 1
        proportions = {combination: count / total_combinations for combination, count in counts.items()}

        proportions_array = np.array(list(proportions.values()))
        return proportions_array

    # Calculate protein parameters, dipeptides, and tripeptides for the given sequence  
    seq = str(prot_record.seq).rstrip('*')

    protein_params = prot_params(seq)
    di_mer_params = k_mer(seq, 2)
    tri_mer_params = k_mer(seq, 3)

    # Concatenate the results and return them as a single array
    return np.concatenate((protein_params, di_mer_params, tri_mer_params))
    # return protein_params


# Generate a pandas dataframe from a list of numpy arrays containing protein parameters
def generate_dataframe(protein_parameters):
    # Create a list of column names for the DataFrame
    column_names = ['name',
                    'molecular_weight', 'aromaticity', 'instability_index', 'isoelectric_point',
                    'sec_struc_helix', 'sec_struc_turn', 'sec_struc_sheet',
                    'molar_extinction_coefficient_reduced', 'molar_extinction_coefficient_cysteines',
                    'gravy', 'flexibility', 
                    *['aa_{}'.format(aa) for aa in 'ACDEFGHIKLMNPQRSTVWY'],
                    'length',
                    *['di_mer_{}'.format(i) for i in range(1, 50)],
                    *['tri_mer_{}'.format(i) for i in range(1, 344)]
                   ]

    # Convert the list of numpy arrays to a list of lists
    protein_data = [list(params) for params in protein_parameters]

    # Create a pandas DataFrame from the list of lists
    df = pd.DataFrame(protein_data, columns=column_names)

    return df



def save_dataframe(dataframe, output_file):
    dataframe.to_csv(output_file, index=False)