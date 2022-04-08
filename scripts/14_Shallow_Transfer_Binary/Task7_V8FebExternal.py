#!/usr/bin/env python
# coding: utf-8
#!/usr/bin/env python3

#Script that performs ML Task 7 for an antigen. 
#Usage in command line:
#python ThisFile.py fileInput condition=1 removeDummies=false externalTestFile

from __future__ import division, print_function, absolute_import
#tf.keras.backend.set_floatx('float64') problems with float32 versus float64. Don't get.

runningInCommandLine = True   # write false for jupyter, then more plotting available

import tensorflow as tf
import collections
import os
import random
import urllib
import zipfile
import numpy as np
import csv
import sys
import math

from tensorflow.keras import regularizers

from collections import Counter

from sklearn import svm
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_score, recall_score, f1_score, accuracy_score, matthews_corrcoef

if(not(runningInCommandLine)):
    import matplotlib.pyplot as plt

#Hides the presence of the GPU, to avoid problems on immunhub
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

print(np.__version__) #should be 1.18.5



fileInput = "T7_1ADQ_A_Task1_SlicesBalancedData.txt_A.txt"
condition=1                #1: Normal, 3: Shuffled
removeDummies = "false"
externalTestFile = "TestDatasetFor1ADQ_A_Task1_SlicesBalancedData.txt"

if(runningInCommandLine):
    if(len(sys.argv) > 1):
        fileInput = sys.argv[1]

    if(len(sys.argv) > 2):
        condition = int(sys.argv[2])

    if(len(sys.argv) > 3):
        architecture = sys.argv[3]
        
    if(len(sys.argv) > 4):
        removeDummies = sys.argv[4]
        if(removeDummies == "True" or removeDummies == "TRUE"):
            removeDummies = "true"
        if(removeDummies == "False" or removeDummies == "FALSE"):
            removeDummies = "false"

    if(len(sys.argv) > 5):
        externalTestFile = sys.argv[4]
		
		
		
#Defining the one hot encoding functions and reverse
alphabet = 'ACDEFGHIKLMNPQRSTVWY'
# define a mapping of chars to integers
char_to_int = dict((c, i) for i, c in enumerate(alphabet))
int_to_char = dict((i, c) for i, c in enumerate(alphabet))

def hotEncodingAAString(myString):
    onehot_encoded = list()
    integer_encoded = [char_to_int[char] for char in myString]
    for value in integer_encoded:
        letter = [0 for _ in range(len(alphabet))]
        letter[value] = 1
        if(removeDummies == "true"):
            onehot_encoded.append(letter[0:len(alphabet)-2])
        else:
            onehot_encoded.append(letter)
    return onehot_encoded

#input: [0,0,...,1,0,0]
def retrieveAA(onehot_encodedAA):
    #print("AA", onehot_encodedAA)
    #if(sum(onehot_encodedAA) != 1):
    #    return '?'
    foundAA = '?'
    for i in range(len(onehot_encodedAA)):
        if(onehot_encodedAA[i] == 1):
            foundAA = int_to_char[i]
    if(removeDummies == "true" and foundAA == '?'):
        return alphabet[len(alphabet)-1]
    return foundAA

#input: [[[0,0,...,1,0,0] , ... , [0,1,0,0... 0]]]
def retrieveString(onehot_encodedString):
    #print("String", onehot_encodedString)
    foundString = ""
    for k in range(len(onehot_encodedString)):
        foundString = foundString + retrieveAA(onehot_encodedString[k])
    return foundString

def flattenedHotEncodingAAString(myString):
    return np.array(hotEncodingAAString(myString)).flatten()

#print(hotEncodingAAString("ACYD")); 
#print(retrieveString(hotEncodingAAString("ACYD")))
#print(flattenedHotEncodingAAString("ACYD"));
#print(np.array(hotEncodingAAString("ACYD")).shape)



def readTwoColumnsFile(fileName):
    balancedDataset = open(fileName, newline = '')   #one line is a text with \t and \n                                                                              
    data_reader = csv.reader(balancedDataset, delimiter='\t') #transform lines into lists 
    sequences = []
    labels = []
    for line in data_reader:
        if(not (line[0].startswith("CDR3") or line[0].startswith("Slide"))):
            sequences.append(line[0])
            labels.append(line[1])
    
    return (sequences, labels)
	
	
	
# Train-Val dataset
(sequences, labels) = readTwoColumnsFile(fileInput)

# For shallow and ANN: flat one-hot      
hotEncodedKeys = np.array(list(map(flattenedHotEncodingAAString, sequences)))
# For CNN and LSTM: 2D one-hot        
hotEncodedKeys2D = np.array(list(map(hotEncodingAAString, sequences)))

binary = {"Binder": 1,"NonBinder": 0}
binaryLabels = np.array([binary[item] for item in labels])

if(condition == 3):
    random.shuffle(binaryLabels)

nSeqTot = len(hotEncodedKeys)
train_size = int(0.8 * nSeqTot) 
test_size = int(0.2 * nSeqTot)

possibleDataIDs = [*range(0,len(hotEncodedKeys))]
random.shuffle(possibleDataIDs)



#External test dataset
(sequencesExt, labelsExt) = readTwoColumnsFile(fileInput)

hotEncodedKeysExt = np.array(list(map(flattenedHotEncodingAAString, sequencesExt)))      
hotEncodedKeys2DExt = np.array(list(map(hotEncodingAAString, sequencesExt)))
binaryLabelsExt = np.array([binary[item] for item in labelsExt])

#if(condition == 3):
#    random.shuffle(binaryLabelsTest)



# First, 1D encodings for shallow models
trainKeys = hotEncodedKeys[possibleDataIDs[0:train_size]].astype(float);
trainLabels = binaryLabels[possibleDataIDs[0:train_size]];
testKeys = hotEncodedKeys[possibleDataIDs[train_size: train_size + test_size]].astype(float);
testLabels = binaryLabels[possibleDataIDs[train_size: train_size + test_size]];

print(trainKeys.shape)
print('ML Task 7, fileInput ', fileInput, 'condition=', condition)
print("Data Ready: tot=", nSeqTot,  " train=", len(trainKeys), str(Counter(trainLabels)), ' test=', len(testKeys), str(Counter(testLabels)))



# Shallow architectures
for architecture in ["LR", "LSVM", "SVM", "RF", "AbsLR", "AbsRF"]:

    if(architecture == "LR"):
        #L2 penalty
        #tolerance 1e-4
        #solver lbfgs
        model = LogisticRegression(solver="lbfgs", C=1e-4, penalty='l2')    #   class_weight='balanced')

    if(architecture == "LSVM"):
        #class sklearn.svm.LinearSVC(penalty='l2', loss='squared_hinge', 
        #Loss function = Squared hinge
        #Tolerance = 1e-4
        #C = Regularization float = 1.0
        model = svm.LinearSVC(tol=1e-4, C = 1, loss='squared_hinge')

    if(architecture == "SVM"):
        #Kernel: Gaussian RBF
        #Degree = 3
        #Gamma = Scaled
        #1/n_features * X.var()
        #Tolerance = 1e-3
        model = svm.SVC(kernel='rbf', degree=3, gamma='scale', tol=1e-3)

    if(architecture == "RF"):
        #Number trees = 150
        #Criterion = Gini = mean decrease impurity
        model = RandomForestClassifier(n_estimators=150, criterion="gini")

    if(architecture == "AbsLR"):
        model = LogisticRegression(random_state=0)

    if(architecture == "AbsRF"):
        model = RandomForestClassifier(max_depth=2, random_state=0)


    model.fit(trainKeys, trainLabels)

    y_pred = model.predict(trainKeys)
    metricMatTrain = ["NoLoss", accuracy_score(trainLabels, y_pred), precision_score(trainLabels, y_pred), recall_score(trainLabels, y_pred), 
                 matthews_corrcoef(trainLabels, y_pred), f1_score(trainLabels, y_pred)]#, model.score(testKeys, trainLabels)]

    y_pred = model.predict(testKeys)
    metricMat = ["NoLoss", accuracy_score(testLabels, y_pred), precision_score(testLabels, y_pred), recall_score(testLabels, y_pred), 
                 matthews_corrcoef(testLabels, y_pred), f1_score(testLabels, y_pred)]#, model.score(testKeys, testLabels)]

    y_pred = model.predict(hotEncodedKeysExt)
    metricMatExt = ["NoLoss", accuracy_score(binaryLabelsExt, y_pred), precision_score(binaryLabelsExt, y_pred), recall_score(binaryLabelsExt, y_pred), 
                 matthews_corrcoef(binaryLabelsExt, y_pred), f1_score(binaryLabelsExt, y_pred)]#, model.score(testKeys, binaryLabelsExt)]

    
    print(architecture, metricMatTrain, metricMat, metricMatExt)
    
    file_object = open('HistoryTask7.txt', 'a')
    file_object.write(fileInput + "\t" + architecture + "\t" + str(condition) + "\t" + str(metricMatTrain) + "\t" + str(metricMat) + "\t" + str(metricMatExt) + "\n")
    file_object.close()
	
	
	
for architecture in ["ANN", "CNN", "LSTM"]:

    #Lstm and CNN require 2D inputs 
    if(architecture == "CNN" or architecture == "LSTM"):
        trainKeys = hotEncodedKeys2D[possibleDataIDs[0:train_size]].astype(float);
        trainLabels = binaryLabels[possibleDataIDs[0:train_size]];
        testKeys = hotEncodedKeys2D[possibleDataIDs[train_size: train_size + test_size]].astype(float);
        testLabels = binaryLabels[possibleDataIDs[train_size: train_size + test_size]];
        extKeys = hotEncodedKeys2DExt
        extLabels = binaryLabelsExt
    else: 
        extKeys = hotEncodedKeysExt
        extLabels = binaryLabelsExt
        
    #Transform into TF data structures
    TFkeysTrain = tf.data.Dataset.from_tensor_slices(trainKeys)
    TFvalTrain = tf.data.Dataset.from_tensor_slices(trainLabels)
    TFkeysTest = tf.data.Dataset.from_tensor_slices(testKeys)
    TFvalTest = tf.data.Dataset.from_tensor_slices(testLabels)
    TFkeysExt = tf.data.Dataset.from_tensor_slices(extKeys)
    TFvalExt = tf.data.Dataset.from_tensor_slices(extLabels)

    tfTrainElements = tf.data.Dataset.zip((TFkeysTrain, TFvalTrain))
    tfTestElements = tf.data.Dataset.zip((TFkeysTest, TFvalTest))
    tfExtElements = tf.data.Dataset.zip((TFkeysExt, TFvalExt))

    input_shape = trainKeys[0].shape
    print(architecture, "Input shape is", input_shape)
        
    if(architecture == "ANN"):
        #Hidden layers: 3
        #Nodes per layer = 70
        #Activation = ReLu
        #Dropout = 0.1
        #Output = Sigmoid
        #Optimizer = Adam
        #Loss = Binary cross entropy
        #Batch size = 16
        #Epochs = 20
        #in github, adam learning_rate=0.0001

        model = tf.keras.Sequential();
        model.add(tf.keras.layers.Dense(units=70,kernel_initializer='uniform',activation='relu',input_shape=input_shape))    
        model.add(tf.keras.layers.Dropout(rate=0.1))
        model.add(tf.keras.layers.Dense(units=70,kernel_initializer='uniform',activation='relu'))    
        model.add(tf.keras.layers.Dropout(rate=0.1))
        model.add(tf.keras.layers.Dense(units=70,kernel_initializer='uniform',activation='relu'))    
        model.add(tf.keras.layers.Dropout(rate=0.1))
        model.add(tf.keras.layers.Dense(1, kernel_initializer='uniform', activation='sigmoid'))
        optimizer = tf.keras.optimizers.Adam()
        loss = 'binary_crossentropy'
        nEpochs = 20
        batch_size = 16

    #For CNN, should not flatten (takes 10x20 input)
    if(architecture == "CNN"):
        #Kernel size = 3
        #Number filters = 400
        #Pool size = 2
        #Dense layer nodes = 50
        #Activation = ReLu
        #Dropout = 0.5
        #Output = Sigmoid
        #Optimizer = Adam
        #Loss = Binary cross entropy
        #Batch size = 16
        #Epochs = 20
        # In github, adam learning_rate=0.000075
        regularizer = None

        model = tf.keras.Sequential()
        model.add(tf.keras.layers.InputLayer(input_shape))    
        model.add(tf.keras.layers.Conv1D(filters=400, kernel_size=3, input_shape=input_shape))    
        model.add(tf.keras.layers.Conv1D(filters=400, kernel_size=3, strides=1, activation='relu', kernel_regularizer=regularizer, bias_regularizer=regularizer, padding='same'))    
        model.add(tf.keras.layers.Dropout(rate=0.5))        
        model.add(tf.keras.layers.MaxPool1D(pool_size=2, strides=1))
        model.add(tf.keras.layers.Flatten()) 
        model.add(tf.keras.layers.Dense(units=50, activation='relu', kernel_regularizer=regularizer, bias_regularizer=regularizer))   
        model.add(tf.keras.layers.Dense(1, activation='sigmoid'))

        optimizer = tf.keras.optimizers.Adam()
        loss = 'binary_crossentropy'
        nEpochs = 20
        batch_size = 16

    # The LSTM also should not be flattened
    if(architecture == "LSTM"):
        #l2_regularization=1e-4
        regularizer=tf.keras.regularizers.l2(1e-4)
        #regularizer=regularizers.L1L2(l1=0.0, l2=1e-4),

        model = tf.keras.Sequential()
        model.add(tf.keras.layers.LSTM(units=40, return_sequences=True, bias_regularizer=regularizer, input_shape=input_shape))
        model.add(tf.keras.layers.Dropout(rate=0.1))
        model.add(tf.keras.layers.LSTM(units=40, bias_regularizer=regularizer, return_sequences=True))
        model.add(tf.keras.layers.Dropout(rate=0.1))
        model.add(tf.keras.layers.LSTM(units=40, bias_regularizer=regularizer, return_sequences=False))
        model.add(tf.keras.layers.Dropout(rate=0.1))
        model.add(tf.keras.layers.Dense(units=1, activation='sigmoid'))

        optimizer = tf.keras.optimizers.Adam()
        loss = 'mean_squared_error'
        nEpochs = 20
        batch_size = 16

    #https://stackoverflow.com/questions/56865344/how-do-i-calculate-the-matthews-correlation-coefficient-in-tensorflow
    def mcc_metric(y_true, y_pred):
        predicted = tf.cast(tf.greater(y_pred, 0.5), tf.float32)
        true_pos = tf.math.count_nonzero(predicted * y_true)
        true_neg = tf.math.count_nonzero((predicted - 1) * (y_true - 1))
        false_pos = tf.math.count_nonzero(predicted * (y_true - 1))
        false_neg = tf.math.count_nonzero((predicted - 1) * y_true)
        x = tf.cast((true_pos + false_pos) * (true_pos + false_neg) 
          * (true_neg + false_pos) * (true_neg + false_neg), tf.float32)
        return tf.cast((true_pos * true_neg) - (false_pos * false_neg), tf.float32) / (tf.sqrt(x) + 0.00000001)

    metrics = ['accuracy','Precision', 'Recall', mcc_metric] #, mcc_metric]
    #'FalseNegatives', 'FalsePositives', 'TrueNegatives', 'TruePositives',
    # tfa.metrics.F1Score(num_classes=1, threshold=0.5, average='macro')
    model.build(input_shape)
    model.compile(loss=loss, optimizer=optimizer, metrics=metrics)
    model.summary()
    
    train_dataset = tfTrainElements.shuffle(train_size).batch(batch_size, drop_remainder=True)
    test_dataset = tfTestElements.shuffle(test_size).batch(batch_size, drop_remainder=True)
    ext_dataset = tfExtElements.shuffle(test_size).batch(batch_size, drop_remainder=True)
    
    history = model.fit(train_dataset, epochs=nEpochs, validation_data=test_dataset, validation_steps=1)
    print(architecture, str(model.evaluate(train_dataset)), "noF1", str(model.evaluate(test_dataset)), "noF1", str(model.evaluate(ext_dataset)), "noF1")
                  
    file_object = open('HistoryTask7.txt', 'a')
    file_object.write(fileInput + "\t" + architecture + "\t" + str(condition) + "\t" + str(model.evaluate(train_dataset))+ "\tnoF1\t" + str(model.evaluate(test_dataset)) + "\t" + "noF1" + "\t" + str(model.evaluate(ext_dataset)) + "noF1\n")
    file_object.close()
	
	
	
