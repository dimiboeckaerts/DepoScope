import argparse
import subprocess
from io import StringIO

import pandas as pd
from tqdm import tqdm

import torch
import torch.nn.functional as F
from torch import nn

from transformers import AutoModelForTokenClassification, AutoTokenizer

from Bio import SeqIO
from Bio.Seq import Seq



def phanotate_processing(input_file_path: str) -> pd.DataFrame:
    """
    Function to process the output of PHANOTATE and return a dataframe with the genes and proteins.
    Input:
    - `input_file_path`, the path to the input file

    Output:
    - a dataframe with the following columns:
        - `gene_ID`: the gene ID
        - `start_position`: the start position of the gene
        - `stop_position`: the stop position of the gene
        - `gene_sequence`: the gene sequence
        - `protein_sequence`: the protein sequence

    Example:
    >>> phage_genes = phanotate_processing('phage.fasta')
    >>> print(phage_genes)
    gene_ID  start_position  stop_position  gene_sequence  protein_sequence
    0     gp-0            1021           1234          ATG...           M...
    1     gp-1            1244           2345          ATG...           M...
    2     gp-2            2365           3456          ATG...           M...
    """
    # get dataframe out of phanotate output
    phanotate_out = subprocess.run(["phanotate.py", input_file_path], capture_output=True)
    phanotate_out = phanotate_out.stdout[phanotate_out.stdout.find(b"#START")+1:]
    results_orfs = pd.read_csv(StringIO(str(phanotate_out, 'utf-8')), sep='\t')
    # fill up lists accordingly
    sequence = str(SeqIO.read(input_file_path, 'fasta').seq)
    gene_list = []
    protein_list = []
    gene_ids = []
    for j, strand in enumerate(results_orfs['FRAME']):
        start = results_orfs['START'][j]
        stop = results_orfs['STOP'][j]
        if strand == '+':
            gene = sequence[start-1:stop]
        else:
            sequence_part = sequence[stop-1:start]
            gene = str(Seq(sequence_part).reverse_complement())
        gene_list.append(gene)
        gene_ids.append(f'gp-{j}')
        protein_list.append(str(Seq(gene).translate())[:-1])
    # Export final genes database
    genebase = pd.DataFrame(
        list(zip(gene_ids, results_orfs['START'], results_orfs['STOP'], gene_list, protein_list)), 
        columns=['gene_ID', 'start_position', 'stop_position', 'gene_sequence', 'protein_sequence'])
    return genebase


class Dpo_classifier(nn.Module):
    def __init__(self, pretrained_model, tokenizer):
        super(Dpo_classifier, self).__init__()
        self.max_length = 1024
        self.pretrained_model = pretrained_model
        self.conv1 = nn.Conv1d(1, 64, kernel_size=5, stride=1)  # Convolutional layer
        self.conv2 = nn.Conv1d(64, 128, kernel_size=5, stride=1)  # Convolutional layer
        self.fc1 = nn.Linear(128 * (self.max_length - 2 * (5 - 1)), 32)  # calculate the output shape after 2 conv layers
        self.classifier = nn.Linear(32, 1)  # Binary classification
        self.tokenizer = tokenizer

    def make_prediction(self, fasta_txt):
        input_ids = self.tokenizer.encode(fasta_txt, truncation=True, return_tensors='pt')
        with torch.no_grad():
            outputs = self.pretrained_model(input_ids)
            probs = torch.nn.functional.softmax(outputs.logits, dim=-1)
            token_probs, token_ids = torch.max(probs, dim=-1)
            tokens = token_ids.view(1, -1) # ensure 2D shape
            return tokens

    def pad_or_truncate(self, tokens):
        if tokens.size(1) < self.max_length:
            tokens = F.pad(tokens, (0, self.max_length - tokens.size(1)))
        elif tokens.size(1) > self.max_length:
            tokens = tokens[:, :self.max_length]
        return tokens

    def forward(self, sequences):
        batch_size = len(sequences)
        tokens_batch = []
        for seq in sequences:
            tokens = self.make_prediction(seq)
            tokens = self.pad_or_truncate(tokens)
            tokens_batch.append(tokens)

        outputs = torch.cat(tokens_batch).view(batch_size, 1, self.max_length)  # ensure 3D shape
        outputs = outputs.float()  # Convert to float

        out = F.relu(self.conv1(outputs))
        out = F.relu(self.conv2(out))
        out = out.view(batch_size, -1)  # Flatten the tensor
        out = F.relu(self.fc1(out))
        out = self.classifier(out)
        return out, outputs


# define a helper sequence to make predictions with and plot the tokens if requested
def predict_sequence(model, sequence):
    model.eval()
    with torch.no_grad():
        sequence = [sequence]  # Wrap the sequence in a list to match the model's input format
        outputs, sequence_outputs = model(sequence)
        probas = torch.sigmoid(outputs)  # Apply sigmoid activation for binary classification
        predictions = (probas > 0.5).float()  # Convert probabilities to binary predictions
        sequence_outputs_list = sequence_outputs.cpu().numpy().tolist()[0][0]
        prob_predicted = probas[0].item()
        return (predictions.item(), prob_predicted), sequence_outputs_list


def parse_args():
    parser = argparse.ArgumentParser(description='DepoScope prediction script')
    parser.add_argument('-i', '--fasta_in', type=str, help='Input fasta file (single fasta)')
    parser.add_argument('-o', '--csv_out', type=str, help='Output CSV file')
    parser.add_argument('--esm2', type=str, help="""Path to the fine-tuned ESM-2 model checkpoint
                                                    (download with `wget https://zenodo.org/records/10957073/files/esm2_t12_finetuned_depolymerases.zip`,
                                                    this path must point to the folder containing the checkpoint files)""")
    parser.add_argument('--Dpo', type=str, help="""Path to the DepoDetection model file
                                                    (download with `https://zenodo.org/records/10957073/files/Deposcope.esm2_t12_35M_UR50D.2203.full.model`,
                                                    this path must point to the `.model` file)""")
    return parser.parse_args()

def main():
    # read input and output file names from command line
    args = parse_args()
    fasta_in = args.fasta_in
    csv_out = args.csv_out
    esm2_model_path = args.esm2
    DpoDetection_path = args.Dpo

    # check that it is a single fasta file
    if len(list(SeqIO.parse(fasta_in, 'fasta'))) != 1:
        raise ValueError('Input fasta file must contain a single sequence')

    # process the input file with PHANOTATE
    phage_genes = phanotate_processing(fasta_in)

    # load pre-trained models
    tokenizer = AutoTokenizer.from_pretrained(esm2_model_path, clean_up_tokenization_spaces=False) 
    esm2_finetuned = AutoModelForTokenClassification.from_pretrained(esm2_model_path)
    model_classifier = Dpo_classifier(esm2_finetuned, tokenizer) # Create an instance of Dpo_classifier
    model_classifier.load_state_dict(torch.load(DpoDetection_path, weights_only=True), strict = False) # Load the saved weights
    model_classifier.eval() # Set the model to evaluation mode for inference


    # run predictions
    scores_DepoScope = []
    token_scores_DepoScope = []
    for protein_seq in tqdm(phage_genes['protein_sequence']):
        prediction, sequence_outputs = predict_sequence(model_classifier, protein_seq)
        scores_DepoScope.append(prediction[1])
        token_scores_DepoScope.append(sequence_outputs)


    # save predictions
    phage_genes = pd.concat([phage_genes, pd.DataFrame({'scores_DepoScope': scores_DepoScope, 'token_labels' : token_scores_DepoScope})], axis=1)
    phage_genes.to_csv(csv_out, index=False)

if __name__ == "__main__":
    main()
