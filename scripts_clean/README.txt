Folder comprising the scripts used for the project.

The notebooks can be arranged according to the MM section of the manuscript.

I. Data Collection
    Collect_InterProEntries_PD_DB : 
### The cookbook for building the depolymerase database :
### A. Get the E.C associated with the depolymerase activity 
### B. Collect the IPR entries associated with the E.C of interest
### C. Scan the interproscan descriptions ; keep those that are relevant in our case
### D. Download the protein sequences, make a filtered hmm profile
### E. Make the DB

***
    
II. Generation of Polysaccharides degrading domains database 
    Make_PD_DB :
### The multi fasta files of the intersting entries have already been downloaded. Here are the steps to follow :
### I. Create a rep for each IPR entry
### II. Cluster each entry using mmseqs
### III. Generate the required files for the final DB
### IV. Make the DB 

***

III. Screening the proteins against the HMM profile Database 
    PT1_Build_MSA_candidates :
### Goal : Generate a MSA for each candidate protein
### A. Run mmseqs on the fasta database vs Uniref90
### B. Generate an MSA for each proteins
### C. Realign the MSAs with clustalo

    PT2_Scan_MSA_candidates :
### Goal : scan the (filtered) MSA with the depolymerase DB 
### 1. hmmcan the MSA
### 2. Scan the results

***

IV. Screening the 3D structure proteins against the depolymerase fold database
    PT1_Filter_PD_hits
### Goal : Save the Safada (or the ambiguous ones)
### I. Loading the dataframes
### II. Transform the features : PCA
### III. Clustering
### IV. Scan the results

    PT2_ESMFold_job
### GOAL : generate ESMfold prediction for the selected proteins after V3
### I. Generate the sequence file
### II. The prediction script

    PT3_Foldseek_job
### Goal : Scan the 3D predcitions with the Dpo fold database
### A.Run Foldseek
### B. Get the hits with good proba
### C. Report the folds

***

V. Determining the delimitation of the depolymerase domain
    PT1 : SWORD2_computation
### Make the SWORD2 prediction 
### I. The prediction 

    PT2 : Foldseek_Protein_Units
***

VI. Modelization
    PT1_esm2_FineTuning_Token_class
### Fine-tuning the model for token classification

    
    PT2_CNV_binary_class
### The Depo Detection tool
###    We finetuned the ESM2 model successfully (92% accuracy)
###    The goal now is to stack a RNN layer for a binary classification into Dpo or Not Dpo categories

### I. Load prebuilt model
### II. Stack RNN layer
### III. Train Eval
### IV. Metrics

***

VII. The DpoDetection Tool

***
VIII. Benchmarking