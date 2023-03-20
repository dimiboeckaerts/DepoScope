"""
15/03/23: getting finetuning to work
"""
import random
# small test dataset

def tokenize_data(inputs):
    """
    ESM needs input_ids, attention_mask, labels, position_ids
    """
    tinputs = tokenizer(inputs['tokens'], add_special_tokens=False, truncation=True, is_split_into_words=True)
    tinputs['labels'] = inputs['ner_tags']
    #return tokenizer(inputs['sequence'], padding="max_length", truncation=True)
    return tinputs

datasets = load_dataset("conll2003")
tokenized_inputs = datasets.map(tokenize_data, batched=True)

tinput = tokenizer(datasets['train'][0]['tokens'], add_special_tokens=False, is_split_into_words=True, truncation=True)

# set the trainer and train
trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=tokenized_inputs["train"],
    eval_dataset=tokenized_inputs["test"],
    data_collator=data_collator,
    tokenizer=tokenizer,
)
trainer.train()


"""
20/03/23: splitting the dataset

SKLEARN does not work, or is confusing, found other solution
"""
from sklearn.model_selection import train_test_split
y = range(len(dataset['train']))
X = dataset['train']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)

count = 0
count_small = 0
for item in tokenized_data['train']:
    if len(item['input_ids']) != len(item['labels']):
        count += 1
    if len(item['input_ids']) < len(item['labels']):
        count_small += 1
print(count, count_small)

# adjust labels after mapping (optional)
new_labels_train = []
for item in tokenized_data['train']:
    if len(item['input_ids']) < len(item['labels']):
        this_new_labels = item['labels'][:len(item['input_ids'])]
        new_labels_train.append(this_new_labels)
tokenized_data['train']['labels'] = new_labels_train

new_labels_test = []
for item in tokenized_data['test']:
    if len(item['input_ids']) < len(item['labels']):
        this_new_labels = item['labels'][:len(item['input_ids'])]
        new_labels_test.append(this_new_labels)
tokenized_data['test']['labels'] = new_labels_test

new_labs = [labels[:1024] for labels in tokenized_data['train']['labels']]
tokenized_data['train']['labels'] = new_labs