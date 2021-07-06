# preps data for task two, uses phil's preprocessing magic

# import stuff
import sys
import csv
import random
import numpy as np
import pandas as pd


#Defining the one hot encoding
alphabet = 'ACDEFGHIKLMNPQRSTVWY'
# define a mapping of chars to integers
char_to_int = dict((c, i) for i, c in enumerate(alphabet))
int_to_char = dict((i, c) for i, c in enumerate(alphabet))

def hotEncodingAAString(myString):
    # Defining the one hot encoding
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    # define a mapping of chars to integers
    char_to_int = dict((c, i) for i, c in enumerate(alphabet))
    int_to_char = dict((i, c) for i, c in enumerate(alphabet))
    onehot_encoded = list()
    integer_encoded = [char_to_int[char] for char in myString]
    for value in integer_encoded:
        letter = [0 for _ in range(len(alphabet))]
        letter[value] = 1
        # onehot_encoded.append(letter)
        onehot_encoded = onehot_encoded + letter
    #print(onehot_encoded)
    return[onehot_encoded]

def tokenize_labels(labels):
    '''
    tokenize the labels
    :param labels:
    :return:
    '''
    ulabels = sorted(set(labels))
    token_dict = {}
    print(len(ulabels))
    for i, item in enumerate(ulabels):
        token_dict[item] = i
    label_tokens = [token_dict[item] for item in labels]
    return label_tokens

def hotEncodingAAString(myString):
    # Defining the one hot encoding
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    # define a mapping of chars to integers
    char_to_int = dict((c, i) for i, c in enumerate(alphabet))
    int_to_char = dict((i, c) for i, c in enumerate(alphabet))
    onehot_encoded = list()
    integer_encoded = [char_to_int[char] for char in myString]
    for value in integer_encoded:
        letter = [0 for _ in range(len(alphabet))]
        letter[value] = 1
        # onehot_encoded.append(letter)
        onehot_encoded = onehot_encoded + letter
    #print(onehot_encoded)
    return[onehot_encoded]

def hotEncodingAAString_batch(myStrings, maxlen):
    # Defining the one hot encoding
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    # define a mapping of chars to integers
    char_to_int = dict((c, i) for i, c in enumerate(alphabet))
    int_to_char = dict((i, c) for i, c in enumerate(alphabet))
    X = []
    for myString in myStrings:
        onehot_encoded = list()
        integer_encoded = [char_to_int[char] for char in myString]
        for value in integer_encoded:
            letter = [0 for _ in range(len(alphabet))]
            letter[value] = 1
            # onehot_encoded.append(letter)
            onehot_encoded = onehot_encoded + letter
        padding_size = maxlen-len(myString)
        padding_vect= [0]*padding_size*len(alphabet)
        onehot_encoded = onehot_encoded + padding_vect
        X.append(onehot_encoded)
    return X

def prep_data_task4(infile, nseqtotal):
    '''
    preps data for taks4: multiclass paratope/epitope prediction
    :param infile:
    :return:
    '''
    df = pd.read_csv(infile, sep='\t')
    df2 = df[df.duplicates==1]
    df = df2
    max_input_len = max([len(item) for item in df.iloc[:,0]])
    print(max_input_len)
    label_tokens = tokenize_labels(df.iloc[:,1])
    df['labels'] = label_tokens
    nSeqTot = nseqtotal
    refMaxSeqFor20pctsTest = 1000000
    nAvailable = df.shape[0]
    nSeqTot = min(nSeqTot, nAvailable)
    train_size = int(0.8 * nSeqTot)
    test_size = int(0.2 * refMaxSeqFor20pctsTest)
    if test_size + train_size > nAvailable:
        test_size = min(test_size, int(0.2 * nAvailable))
    # shuffle the df prior to subsetting train, test datasets
    df = df.sample(frac=1).reset_index(drop=True)
    traindf = df.iloc[:train_size,:]
    testdf = df.iloc[nSeqTot:nSeqTot+test_size,:]
    print(traindf)
    train_x  =  np.array(traindf.iloc[:,0].tolist())
    test_x  = np.array(testdf.iloc[:,0].tolist())
    train_X = hotEncodingAAString_batch(traindf.iloc[:,0], max_input_len)
    train_labels = traindf.iloc[:,1]
    test_X = hotEncodingAAString_batch(testdf.iloc[:,0], max_input_len)
    test_labels = testdf.iloc[:,1]
    print(len(train_X), len(test_X), len(train_labels), len(test_labels))
    train_X = np.array(train_X)
    train_labels = np.array(train_labels)
    test_X = np.array(test_X)
    test_labels = np.array(test_labels)
    #print([len(item) for item in train_X])
    return train_X, train_labels, test_X, test_labels, train_x, test_x


# run stuff
# prep_data_task4('dataset/Task4_E_EpiSeq_ParaSeq_NoDeg_NoStar_n10000.txt', 100)
# prep_data_task4('dataset/Task4_E_EpiSeq_ParaSeq_NoDeg_NoStar.txt', 100)
