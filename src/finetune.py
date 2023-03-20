"""
Created on 14/03/23

@author: dimiboeckaerts

PhageDEPOdetection ESM FINETUNE
"""
# pip install datasets
# pip install transformers
absolute_path = '/Users/dimi/Documents/GitHub/PhageDEPOdetection/'

# 1 - LIBRARIES
# --------------------------------------------------------------------------------
import torch
import numpy as np
from datasets import load_dataset, load_metric, DatasetDict
from transformers import AutoTokenizer, AutoModelForTokenClassification
from transformers import TrainingArguments, Trainer, DataCollatorForTokenClassification
from seqeval.metrics import accuracy_score, f1_score, precision_score, recall_score

# 2 - FUNCTIONS
# --------------------------------------------------------------------------------
def tokenize_function(inputs):
    """
    ESM finetuning needs input_ids, attention_mask, labels, position_ids.
    If the dataset does not contain these columns, we need to add them w this function.
    """
    return tokenizer(inputs['tokens'], add_special_tokens=False, 
                     truncation=True, is_split_into_words=True)

def postprocess_tokenize(tokenized_data):
    """
    adjust label lengths if they dont match.
    """
    count = 0
    if len(tokenized_data['input_ids']) < len(tokenized_data['labels']):
        count += 1
        new_labels = tokenized_data['labels'][:len(tokenized_data['input_ids'])]
        tokenized_data["labels"] = new_labels
    print(count)

    return tokenized_data

def compute_metrics(p):
    predictions, labels = p
    labels = list(labels)
    predictions = list(np.argmax(predictions, axis=2))
    # Remove ignored index (special tokens)
    #true_predictions = [
    #    [label_list[p] for (p, l) in zip(prediction, label) if l != -100]
    #    for prediction, label in zip(predictions, labels)]
    #true_labels = [
    #    [label_list[l] for (p, l) in zip(prediction, label) if l != -100]
    #    for prediction, label in zip(predictions, labels)]
    return {"accuracy": accuracy_score(labels, predictions),
            "precision": precision_score(labels, predictions),
            "recall": recall_score(labels, predictions),
            "f1": f1_score(labels, predictions)}

# 3 - MAIN
# --------------------------------------------------------------------------------
# initialize the tokenizer and model
tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
#model = EsmForTokenClassification.from_pretrained("facebook/esm2_t6_8M_UR50D")
data_collator = DataCollatorForTokenClassification(tokenizer)
model = AutoModelForTokenClassification.from_pretrained("facebook/esm2_t6_8M_UR50D", num_labels=2)

# load the dataset
finetune_data_path = absolute_path+'data/finetune_data.json'
finetune_train = load_dataset("json", data_files=finetune_data_path, field="train")
finetune_test = load_dataset("json", data_files=finetune_data_path, field="test")
finetune_data = DatasetDict()
finetune_data['train'] = finetune_train['train']
finetune_data['test'] = finetune_test['train']

# tokenize the data
tokenized_data = finetune_data.map(tokenize_function, batched=True)
#tokenized_data = tokenized_data.map(postprocess_tokenize, batched=True)

# set the training arguments
training_args = TrainingArguments(
    output_dir=absolute_path+'data/finetune',
    evaluation_strategy="epoch",
    learning_rate=2e-5,
    per_device_train_batch_size=4, #16
    per_device_eval_batch_size=8, #16
    num_train_epochs=1,
    weight_decay=0.01,
)

# set the trainer and train
trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_data["train"],
    eval_dataset=tokenized_data["test"],
    data_collator=data_collator,
    tokenizer=tokenizer,
    compute_metrics=compute_metrics,
)
trainer.train()

# save the model & tokenizer
trainer.save_model('data/finetune')
tokenizer.save_pretrained('data/finetune')

# make predictions & evaluate
predictions, label_ids, metrics = trainer.predict("test_dataset")
predictions, labels, _ = trainer.predict(tokenized_data["validation"])
predictions = np.argmax(predictions, axis=2)
trainer.evaluate()

# load for later use
tokenizer = AutoTokenizer.from_pretrained("data/finetune")
finetuned = AutoModelForTokenClassification.from_pretrained("data/finetune")

"""
Ideas & remarks:
- use a bigger ESM model
- use a bigger batch size
- increase number of training epochs?
https://github.com/agemagician/ProtTrans/blob/master/Fine-Tuning/ProtBert-BFD-FineTune-SS3.ipynb
"""