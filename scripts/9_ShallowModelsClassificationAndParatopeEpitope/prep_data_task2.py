# preps data for task two, uses phil's preprocessing magic

# import stuff
import sys
import csv
import random
import numpy as np


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

#input: [0,0,...,1,0,0]
def retrieveAA(onehot_encodedAA):
    # Defining the one hot encoding
    alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    # define a mapping of chars to integers
    char_to_int = dict((c, i) for i, c in enumerate(alphabet))
    int_to_char = dict((i, c) for i, c in enumerate(alphabet))
    #print("AA", onehot_encodedAA)
    #if(sum(onehot_encodedAA) != 1):
    #    return '?'
    foundAA = '?'
    for i in range(len(onehot_encodedAA)):
        if(onehot_encodedAA[i] == 1):
            foundAA = int_to_char[i]
    return foundAA

#input: [[[0,0,...,1,0,0] , ... , [0,1,0,0... 0]]]
def retrieveString(onehot_encodedString):
    #print("String", onehot_encodedString)
    foundString = ""
    for k in range(len(onehot_encodedString[0])):
        foundString = foundString + retrieveAA(onehot_encodedString[0][k])
    return foundString


def extractNDimLabels(sequence, listSelectedAntigens):
    binary = {"1": 1, "0": 0}
    justPositions = ''.join(sequence[x] for x in listSelectedAntigens)
    listLabels = [binary[x] for x in justPositions]
    return listLabels

def prep_data_task2(infile, nseqtotal):
    '''
    preps data for taks2
    :param infile:
    :return:
    '''

    balancedDataset = open(infile, newline = '')   #one line is a text with \t and \n
    data_reader = csv.reader(balancedDataset, delimiter='\t')  # transform lines into lists
    sequences = []
    labels = []
    for line in data_reader:
        if(not (line[0].startswith("#") or line[0].startswith("Slice"))):
            sequences.append(line[0])
            labels.append(line[2])
    ## phil params
    nSelectAntigens = 50
    nNeurons = 50
    # nSeqTot = 10000
    nSeqTot = nseqtotal
    refMaxSeqFor20pctsTest = 500000
    nRepeats = 10
    layers = 1
    condition = 1
    ## end phil params
    nAntigens = len(labels[0])
    possibleIDs = [*range(0, nAntigens)]
    random.shuffle(possibleIDs)
    selected = possibleIDs[0:min(nAntigens, nSelectAntigens)]
    multiDimensionalLabels = np.array([extractNDimLabels(code, selected) for code in labels])
    print(multiDimensionalLabels)
    hotEncodedKeys = np.array(list(map(hotEncodingAAString, sequences)))
    lenHotEncodedKeys = hotEncodedKeys.shape[0]
    print(lenHotEncodedKeys)
    # Now we will need to take only indices that do bind
    listSumTargets = list(map(sum, multiDimensionalLabels))
    print("Total Sequences=", len(listSumTargets), "Bining none of selected antigens:", listSumTargets.count(0))
    print("Binding 1:", listSumTargets.count(1), " 2:", listSumTargets.count(2), " 3:", listSumTargets.count(3), " 4:",
          listSumTargets.count(4), " 5:", listSumTargets.count(5))

    listBinders = [i for i, value in enumerate(listSumTargets) if value >= 1]
    print(len(listBinders), listBinders[0:10])
    listExactBinders = [i for i, value in enumerate(listSumTargets) if value == 1]
    print(len(listExactBinders), listExactBinders[0:10])
    listNonBinders = [i for i, value in enumerate(listSumTargets) if value == 0]
    print(len(listNonBinders), listNonBinders[0:10])
    random.shuffle(listNonBinders)
    # Here, we take only eact binders, this is multi-class (exclusive)
    possibleDataIDs = listExactBinders
    # + listNonBinders[0:min(len(listNonBinders), len(listExactBinders))]

    # Here, we could also include multibinders only, this is multi-label
    # possibleDataIDs = listBinders + listNonBinders[0:min(len(listNonBinders), len(listBinders))]

    # Will separate train and test datasets now, before stupid tensorflow datasets.
    # possibleDataIDs = [*range(0,len(hotEncodedKeys))]
    random.shuffle(possibleDataIDs)
    # If too many sequences are requested, will downscale, just to have this info in the output
    nSeqTot = min(nSeqTot, len(possibleDataIDs))
    # refMaxSeqFor20pctsTest = min(refMaxSeqFor20pctsTest, len(possibleDataIDs))

    # In[23]:

    nSeqTot = min(nSeqTot, len(possibleDataIDs))
    train_size = int(0.8 * nSeqTot)
    test_size = int(0.2 * refMaxSeqFor20pctsTest)
    if test_size + train_size > len(possibleDataIDs):
        test_size = min(test_size, int(0.2 * len(possibleDataIDs)))
    print('Requested: dataset=', nSeqTot, 'train=', train_size, 'test=', test_size, 'totUsed=', train_size + test_size)
    trainKeys = hotEncodedKeys[possibleDataIDs[0:train_size]].astype(float);
    trainLabels = multiDimensionalLabels[possibleDataIDs[0:train_size]];
    testKeys = hotEncodedKeys[possibleDataIDs[train_size: train_size + test_size]].astype(float);
    testLabels = multiDimensionalLabels[possibleDataIDs[train_size: train_size + test_size]];
    # print(len(trainKeys[1][0][0]),len(trainKeys))
    print(testLabels[0], len(testLabels))
    trainKeys = trainKeys[:,0,:]
    trainLabels = np.argmax(trainLabels,1)
    testKeys = testKeys[:,0,:]
    testLabels = np.argmax(testLabels,1)
    return trainKeys, trainLabels, testKeys, testLabels


def prep_data_task2_nantigen(infile, nseqtotal, nantigen):
    '''
    preps data for taks2
    :param infile:
    :return:
    '''

    balancedDataset = open(infile, newline = '')   #one line is a text with \t and \n
    data_reader = csv.reader(balancedDataset, delimiter='\t')  # transform lines into lists
    sequences = []
    labels = []
    for line in data_reader:
        if(not (line[0].startswith("#") or line[0].startswith("Slice"))):
            sequences.append(line[0])
            labels.append(line[2])
    ## phil params
    #nSelectAntigens = 50
    nSelectAntigens = nantigen
    nNeurons = 50
    # nSeqTot = 10000
    nSeqTot = nseqtotal
    refMaxSeqFor20pctsTest = 500000
    nRepeats = 10
    layers = 1
    condition = 1
    ## end phil params
    nAntigens = len(labels[0])
    possibleIDs = [*range(0, nAntigens)]
    random.shuffle(possibleIDs)
    selected = possibleIDs[0:min(nAntigens, nSelectAntigens)]
    multiDimensionalLabels = np.array([extractNDimLabels(code, selected) for code in labels])
    print(multiDimensionalLabels)
    hotEncodedKeys = np.array(list(map(hotEncodingAAString, sequences)))
    lenHotEncodedKeys = hotEncodedKeys.shape[0]
    print(lenHotEncodedKeys)
    # Now we will need to take only indices that do bind
    listSumTargets = list(map(sum, multiDimensionalLabels))
    print("Total Sequences=", len(listSumTargets), "Bining none of selected antigens:", listSumTargets.count(0))
    print("Binding 1:", listSumTargets.count(1), " 2:", listSumTargets.count(2), " 3:", listSumTargets.count(3), " 4:",
          listSumTargets.count(4), " 5:", listSumTargets.count(5))

    listBinders = [i for i, value in enumerate(listSumTargets) if value >= 1]
    print(len(listBinders), listBinders[0:10])
    listExactBinders = [i for i, value in enumerate(listSumTargets) if value == 1]
    print(len(listExactBinders), listExactBinders[0:10])
    listNonBinders = [i for i, value in enumerate(listSumTargets) if value == 0]
    print(len(listNonBinders), listNonBinders[0:10])
    random.shuffle(listNonBinders)
    # Here, we take only eact binders, this is multi-class (exclusive)
    possibleDataIDs = listExactBinders
    # + listNonBinders[0:min(len(listNonBinders), len(listExactBinders))]

    # Here, we could also include multibinders only, this is multi-label
    # possibleDataIDs = listBinders + listNonBinders[0:min(len(listNonBinders), len(listBinders))]

    # Will separate train and test datasets now, before stupid tensorflow datasets.
    # possibleDataIDs = [*range(0,len(hotEncodedKeys))]
    random.shuffle(possibleDataIDs)
    # If too many sequences are requested, will downscale, just to have this info in the output
    nSeqTot = min(nSeqTot, len(possibleDataIDs))
    # refMaxSeqFor20pctsTest = min(refMaxSeqFor20pctsTest, len(possibleDataIDs))

    # In[23]:

    nSeqTot = min(nSeqTot, len(possibleDataIDs))
    train_size = int(0.8 * nSeqTot)
    test_size = int(0.2 * refMaxSeqFor20pctsTest)
    if test_size + train_size > len(possibleDataIDs):
        test_size = min(test_size, int(0.2 * len(possibleDataIDs)))
    print('Requested: dataset=', nSeqTot, 'train=', train_size, 'test=', test_size, 'totUsed=', train_size + test_size)
    trainKeys = hotEncodedKeys[possibleDataIDs[0:train_size]].astype(float);
    trainLabels = multiDimensionalLabels[possibleDataIDs[0:train_size]];
    testKeys = hotEncodedKeys[possibleDataIDs[train_size: train_size + test_size]].astype(float);
    testLabels = multiDimensionalLabels[possibleDataIDs[train_size: train_size + test_size]];
    # print(len(trainKeys[1][0][0]),len(trainKeys))
    print(testLabels[0], len(testLabels))
    trainKeys = trainKeys[:,0,:]
    trainLabels = np.argmax(trainLabels,1)
    testKeys = testKeys[:,0,:]
    testLabels = np.argmax(testLabels,1)
    return trainKeys, trainLabels, testKeys, testLabels



def prep_data_task2_nantigen_aacomp(infile, nseqtotal, nantigen):
    '''
    preps data for taks2
    :param infile:
    :return:
    '''

    balancedDataset = open(infile, newline = '')   #one line is a text with \t and \n
    data_reader = csv.reader(balancedDataset, delimiter='\t')  # transform lines into lists
    sequences = []
    labels = []
    for line in data_reader:
        if(not (line[0].startswith("#") or line[0].startswith("Slice"))):
            sequences.append(line[0])
            labels.append(line[2])
    ## phil params
    #nSelectAntigens = 50
    nSelectAntigens = nantigen
    nNeurons = 50
    # nSeqTot = 10000
    nSeqTot = nseqtotal
    refMaxSeqFor20pctsTest = 500000
    nRepeats = 10
    layers = 1
    condition = 1
    ## end phil params
    nAntigens = len(labels[0])
    possibleIDs = [*range(0, nAntigens)]
    random.shuffle(possibleIDs)
    selected = possibleIDs[0:min(nAntigens, nSelectAntigens)]
    multiDimensionalLabels = np.array([extractNDimLabels(code, selected) for code in labels])
    print(multiDimensionalLabels)
    hotEncodedKeys = np.array(list(map(hotEncodingAAString, sequences)))
    lenHotEncodedKeys = hotEncodedKeys.shape[0]
    #hotEncodedKeys = np.sum(hotEncodedKeys, axis = 2)
    print(hotEncodedKeys.shape)
    print(hotEncodedKeys[1])
    sys.exit()
    # Now we will need to take only indices that do bind
    listSumTargets = list(map(sum, multiDimensionalLabels))
    print("Total Sequences=", len(listSumTargets), "Bining none of selected antigens:", listSumTargets.count(0))
    print("Binding 1:", listSumTargets.count(1), " 2:", listSumTargets.count(2), " 3:", listSumTargets.count(3), " 4:",
          listSumTargets.count(4), " 5:", listSumTargets.count(5))

    listBinders = [i for i, value in enumerate(listSumTargets) if value >= 1]
    print(len(listBinders), listBinders[0:10])
    listExactBinders = [i for i, value in enumerate(listSumTargets) if value == 1]
    print(len(listExactBinders), listExactBinders[0:10])
    listNonBinders = [i for i, value in enumerate(listSumTargets) if value == 0]
    print(len(listNonBinders), listNonBinders[0:10])
    random.shuffle(listNonBinders)
    # Here, we take only eact binders, this is multi-class (exclusive)
    possibleDataIDs = listExactBinders
    # + listNonBinders[0:min(len(listNonBinders), len(listExactBinders))]

    # Here, we could also include multibinders only, this is multi-label
    # possibleDataIDs = listBinders + listNonBinders[0:min(len(listNonBinders), len(listBinders))]

    # Will separate train and test datasets now, before stupid tensorflow datasets.
    # possibleDataIDs = [*range(0,len(hotEncodedKeys))]
    random.shuffle(possibleDataIDs)
    # If too many sequences are requested, will downscale, just to have this info in the output
    nSeqTot = min(nSeqTot, len(possibleDataIDs))
    # refMaxSeqFor20pctsTest = min(refMaxSeqFor20pctsTest, len(possibleDataIDs))

    # In[23]:

    nSeqTot = min(nSeqTot, len(possibleDataIDs))
    train_size = int(0.8 * nSeqTot)
    test_size = int(0.2 * refMaxSeqFor20pctsTest)
    if test_size + train_size > len(possibleDataIDs):
        test_size = min(test_size, int(0.2 * len(possibleDataIDs)))
    print('Requested: dataset=', nSeqTot, 'train=', train_size, 'test=', test_size, 'totUsed=', train_size + test_size)
    trainKeys = hotEncodedKeys[possibleDataIDs[0:train_size]].astype(float);
    trainLabels = multiDimensionalLabels[possibleDataIDs[0:train_size]];
    testKeys = hotEncodedKeys[possibleDataIDs[train_size: train_size + test_size]].astype(float);
    testLabels = multiDimensionalLabels[possibleDataIDs[train_size: train_size + test_size]];
    # print(len(trainKeys[1][0][0]),len(trainKeys))
    print(testLabels[0], len(testLabels))
    trainKeys = trainKeys[:,0,:]
    trainLabels = np.argmax(trainLabels,1)
    testKeys = testKeys[:,0,:]
    testLabels = np.argmax(testLabels,1)
    return trainKeys, trainLabels, testKeys, testLabels


# run stuff
# prep_data_task2('pyPhilippe/Treated142_sampled.txt')
# prep_data_task2('pyPhilippe/Treated142_sampled.txt')
#prep_data_task2_nantigen_aacomp('dataset/Treated142_sampled.txt', 250,5)
