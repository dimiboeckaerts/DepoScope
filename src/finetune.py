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
import json
import torch
import numpy as np
from datasets import load_dataset
from transformers import AutoTokenizer, EsmForTokenClassification, AutoModelForTokenClassification, TrainingArguments, Trainer, DataCollatorForTokenClassification
from seqeval.metrics import accuracy_score, f1_score, precision_score, recall_score


# 2 - FUNCTIONS
# --------------------------------------------------------------------------------
def tokenize_function(inputs):
    """
    ESM finetuning needs input_ids, attention_mask, labels, position_ids.
    If the dataset does not contain these columns, we need to add them w this function.
    """
    return tokenizer(inputs['tokens'], add_special_tokens=False, truncation=True, is_split_into_words=True)
    #return tokenizer(inputs['sequence'], padding="max_length", truncation=True)

def compute_metrics(p):
    predictions, labels = p
    predictions = np.argmax(predictions, axis=2)

    # Remove ignored index (special tokens)
    #true_predictions = [
    #    [label_list[p] for (p, l) in zip(prediction, label) if l != -100]
    #    for prediction, label in zip(predictions, labels)]
    #true_labels = [
    #    [label_list[l] for (p, l) in zip(prediction, label) if l != -100]
    #    for prediction, label in zip(predictions, labels)]

    return {
        "accuracy": accuracy_score(labels, predictions),
        "precision": precision_score(labels, predictions),
        "recall": recall_score(labels, predictions),
        "f1": f1_score(labels, predictions),
    }

# 3 - MAIN
# --------------------------------------------------------------------------------
# initialize the tokenizer and model
tokenizer = AutoTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
#model = EsmForTokenClassification.from_pretrained("facebook/esm2_t6_8M_UR50D")
data_collator = DataCollatorForTokenClassification(tokenizer)
model = AutoModelForTokenClassification.from_pretrained("facebook/esm2_t6_8M_UR50D", num_labels=2)

# load & tokenize the dataset
#ds = load_dataset('csv', data_files='path/to/local/my_dataset.csv')
path = absolute_path + 'data/processed_database.json'
inputs = load_dataset('json', data_files='path/to/local/my_dataset.json')
tokenized_inputs = inputs.map(tokenize_function, batched=True)

# set the training arguments
training_args = TrainingArguments(
    output_dir=absolute_path+'data/finetune',
    evaluation_strategy="epoch",
    learning_rate=2e-5,
    per_device_train_batch_size=1, #16
    per_device_eval_batch_size=8, #16
    num_train_epochs=3,
    weight_decay=0.01,
)

# set the trainer and train
trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_inputs["train"],
    eval_dataset=tokenized_inputs["test"],
    data_collator=data_collator,
    tokenizer=tokenizer,
    compute_metrics=compute_metrics,
)
trainer.train()

# save the model
trainer.save_model('data/finetune')


# make predictions & evaluate
predictions, label_ids, metrics = trainer.predict("test_dataset")
predictions, labels, _ = trainer.predict(tokenized_inputs["validation"])
predictions = np.argmax(predictions, axis=2)

trainer.evaluate()

# 4 - FLOW 3
# --------------------------------------------------------------------------------
"""
https://github.com/agemagician/ProtTrans/blob/master/Fine-Tuning/ProtBert-BFD-FineTune-SS3.ipynb
"""