# DepoScope - Predict and annotate phage depolymerases

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1A2XJ_oUtlmIfU3XXmev5dzJUNxqR6VV9?usp=sharing)

## General information

This is the repository related to our manuscript *DepoScope: accurate phage depolymerase annotation and domain delineation using large language models*, currently in submission at PLOS Computational Biology.

## Quick start

The easiest way to get started with DepoScope to predict depolymerases from your phage genomes or genes is to run the Google Colab that we provide [here](https://colab.research.google.com/drive/1A2XJ_oUtlmIfU3XXmev5dzJUNxqR6VV9?usp=sharing). You only need a zip file of your phage genomes to get started! Alternatively, go to the the `scripts_clean` folder and run the `VII.DpoDetectionTool.ipynb` notebook.

To run the benchmarking against other depolymerase detection tools, go to the `benchmark` folder and run the `benchmark_notebook.ipynb`notebook.

## Running DepoScope as a script

To run DepoScope as a script, a few steps are required. 
The following instructions are for a Unix-based system, but they can be adapted to Windows with minimal changes. 
First, install [`uv`](https://docs.astral.sh/uv/) with:
```bash
curl -LsSf https://astral.sh/uv/install.sh | less
```
Then, make sure that `clang` is installed (it is required by `phanotate`) using `which clang`. If it is not installed, install it with `sudo apt install clang`.

Now we need to download the pre-trained model weights. First we download and extract the fine-tuned ESM-2 model weights:
```bash
wget https://zenodo.org/records/10957073/files/esm2_t12_finetuned_depolymerases.zip
unzip esm2_t12_finetuned_depolymerases.zip
```
Then we download the pre-trained DepoScope model:
```bash
wget https://zenodo.org/records/10957073/files/Deposcope.esm2_t12_35M_UR50D.2203.full.model
```

Finally, clone the repository and run `deposcope-predict.py` through `uv` with the following command:
```bash
uv run deposcope-predict.py -i <input_fasta_file> -o <output_file_name> --esm2 <path_to_esm2_checkpoint> --Dpo <path_to_Dpo_model>
```
where:
- `<input_fasta_file>` is the path to the input fasta file with *a single phage genome* to be annotated.
- `<output_file_name>` is the desired name for the output file.
- `<path_to_esm2_checkpoint>` is the path to the folder containing the fine-tuned ESM-2 model checkpoint files.
- `<path_to_Dpo_model>` is the path to the pre-trained DepoScope `.model` file.
