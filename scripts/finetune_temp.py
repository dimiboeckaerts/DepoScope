test = 'MLPKKPPABBHEBHFF'
tokenizer(test, add_special_tokens=False, padding="max_length", truncation=True)

import random
# small test dataset

# make a dataset like inputs that contains 10 different entries (each with labels and tokens)
# pick a random number of labels and tokens for each entry
inputs = {'train': [], 'test': []}
for i in range(100):
    labels = []
    tokens = []
    for j in range(250):
        labels.append(random.randint(0, 1))
        tokens.append(random.choice(['M', 'L', 'P', 'K', 'K', 'P', 'P', 'A', 'B']))
    inputs['train'].append({'labels': labels, 'tokens': tokens})
for i in range(50):
    labels = []
    tokens = []
    for j in range(250):
        labels.append(random.randint(0, 1))
        tokens.append(random.choice(['M', 'L', 'P', 'K', 'K', 'P', 'P', 'A', 'B']))
    inputs['test'].append({'labels': labels, 'tokens': tokens})

def tokenize_function(inputs):
    return tokenizer(inputs['tokens'], add_special_tokens=False, truncation=True, is_split_into_words=True)
    #return tokenizer(inputs['sequence'], padding="max_length", truncation=True)


tokenized_inputs = inputs.map(tokenize_function, batched=True)


tinput = tokenizer(inputs['train'][0]['tokens'], add_special_tokens=False, is_split_into_words=True, truncation=True)
def tokenize_and_align_labels(examples):
    tokenized_inputs = tokenizer(examples["tokens"], truncation=True, is_split_into_words=True)

    labels = []
    for i, label in enumerate(examples['labels']):
        word_ids = tokenized_inputs.word_ids(batch_index=i)
        previous_word_idx = None
        label_ids = []
        for word_idx in word_ids:
            # Special tokens have a word id that is None. We set the label to -100 so they are automatically
            # ignored in the loss function.
            if word_idx is None:
                label_ids.append(-100)
            # We set the label for the first token of each word.
            elif word_idx != previous_word_idx:
                label_ids.append(label[word_idx])
            # For the other tokens in a word, we set the label to either the current label or -100, depending on
            # the label_all_tokens flag.
            else:
                label_ids.append(label[word_idx])
            previous_word_idx = word_idx

        labels.append(label_ids)

    tokenized_inputs["labels"] = labels
    return tokenized_inputs

tokenized_inputs = tokenize_and_align_labels(inputs['train'][0])