#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python3

#Script that performs ML Task 1 for an antigen. 
#Usage in command line:
#python ThisFile.py antigenID nNeurons=5 nSeqTot=1e7 nRepeats=10 condition=1 layers= batch_size=min(nSeqTot, sequences in dfrom __future__ import division, print_function, absolute_importataset)


# In[2]:


#Should be the first command
from __future__ import division, print_function, absolute_import

from sklearn.metrics import roc_curve
from sklearn.metrics import auc
#tf.keras.backend.set_floatx('float64') problems with float32 versus float64. Don't get.


# In[3]:


# Local settings - Files necessary in that folder: antigenID + "_Task1_Slices[Imb/B]alancedData.txt"
# LocationDataFiles = "C:/Users/pprobert/Desktop/Etulos/ML1_ResultsTask1/PythonNotebook/"
LocationDataFiles = "/storage/pprobert/Task1/"
runningInCommandLine = True
nEpochs = 150
refMaxSeqFor20pctsTest = 100000  # test will be 20% of that at max.
nRepeatst = 1;      #For the jupyter version


# In[4]:


#Default parameters (changeable by command line)
antigenID = "1ADQ_A"
nNeurons = 10
nSeqTot = 50_000
nRepeats = 10
condition = 1
layers = 1


# In[5]:


import tensorflow as tf
import collections
import os
import random
import urllib
import zipfile

import numpy as np
import tensorflow as tf
import csv

import sys
import math
#pip install tensorflow_datasets
#import tensorflow_datasets as tfds

if(not(runningInCommandLine)):
    import matplotlib.pyplot as plt


# In[6]:


#Hides the presence of the GPU, to avoid problems on immunhub
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"


# In[7]:


#Parameters from command line. This cell does not work from Python notebook
#print(len(sys.argv), str(sys.argv))
# if(runningInCommandLine):
#     if(len(sys.argv) > 1):
#         antigenID = sys.argv[1]


#     if(len(sys.argv) > 2):
#         nNeurons = int(sys.argv[2])


#     if(len(sys.argv) > 3):
#         nSeqTot = int(sys.argv[3])


#     if(len(sys.argv) > 4):
#         nRepeats = int(sys.argv[4])

    # 1 = balanced dataset + sigmoid
    # 2 = imbalanced dataset + sigmoid (a bit useless)
    # 3 = shuffling labels + sigmoid 
    
    # First hard case is Mascotte (top 1%) versus loosers exclusive (top 5% to 1%)
    # 4 = HardV1 balanced + relu 
    # 5 = HardV1 shuffling + relu

    # Second hard case is Mascotte exclusive (top 1% to 0.1%) versus Heroes (top 0.1%)
    # 6 = HardV2 balanced + relu 
    # 7 = HardV2 shuffling + relu

    # if(len(sys.argv) > 5):
    #     condition = int(sys.argv[5])
        
    # if(len(sys.argv) > 6):
    #     layers = int(sys.argv[6])
        
#see below for argument 7


# 

# In[8]:


#Reading the data files and loading as two vectors, depending on antigenID and on the condition (balanced/imbalanced) Shuffling is made later
fileToOpen = LocationDataFiles + antigenID + "_Task1_SlicesBalancedData.txt"

if(condition == 2):
    fileToOpen = LocationDataFiles + antigenID + "_Task1_SlicesImbalancedData.txt"
    
if((condition == 4) or (condition == 5)):
   fileToOpen = LocationDataFiles + antigenID + "_Task1aMvsL_SlicesBalancedData.txt"

if((condition == 6) or (condition == 7)):
   fileToOpen = LocationDataFiles + antigenID + "_Task1bHvsM_SlicesBalancedData.txt"

fileToOpen = 'train1ADQ_D1_balanced.csv'
# In[9]:


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
        onehot_encoded.append(letter)
    #print(onehot_encoded)
    return[onehot_encoded]

#input: [0,0,...,1,0,0]
def retrieveAA(onehot_encodedAA):
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

# print(hotEncodingAAString("ACD")); print(retrieveString(hotEncodingAAString("ACD")))


# In[10]:


balancedDataset = open(fileToOpen, newline = '')   #one line is a text with \t and \n                                                                              
data_reader = csv.reader(balancedDataset, delimiter='\t') #transform lines into lists 
sequences = []
labels = []
for line in data_reader:
    if(not (line[0].startswith("#") or line[0].startswith("Slide"))):
        sequences.append(line[0])
        labels.append(line[1])


# In[11]:


#Conversion of the dataset into One Hot Encoding- takes 20 sec
hotEncodedKeys = np.array(list(map(hotEncodingAAString, sequences)))


# In[12]:


#Conversion of the Lavels into 0 or 1 (binary classification)
binary = {"Binder": 1,"NonBinder": 0}
binaryLabels = np.array([binary[item] for item in labels])


# In[13]:


#Shuffle the labels
if((condition == 3) or (condition == 5) or (condition == 7)):
    random.shuffle(binaryLabels)


# In[14]:


#If too many sequences are requested, will downscale, just to have this info in the output
nSeqTot = min(nSeqTot,len(hotEncodedKeys))
train_size = int(0.8 * nSeqTot) 
test_size = int(0.2 * refMaxSeqFor20pctsTest)
if test_size + train_size > len(hotEncodedKeys):
   test_size = min(test_size, int(0.2*len(hotEncodedKeys)))

print('Requested: dataset=', nSeqTot, 'train=', train_size, 'test=', test_size, 'totUsed=', train_size + test_size)


# In[15]:


#By default, we will take 50 batches, so the size of a batch is nr sequences / 50
batch_size = math.floor(nSeqTot / 50)
# if(runningInCommandLine):
#     if(len(sys.argv) > 7):
#         batch_size = int(sys.argv[7])


# In[16]:


possibleDataIDs = [*range(0,len(hotEncodedKeys))]
random.shuffle(possibleDataIDs)


# In[17]:


#Transform the lists into tensorflow datasets (not yet batched)

# First, we separate the test and train raw data manually because if we do it with datasets.take() or skip(), 
# there is a problem with shuffling (test and train might overlap), and then very difficult to get back
# what was the train and test datasets. 
# Also, for getting multilabel confusion matrix later, better to keep the raw data aside and recreate 
# a non-shuffled dataset to do model.predict in the good order.
trainKeys = hotEncodedKeys[possibleDataIDs[0:train_size]].astype(float);
trainLabels = binaryLabels[possibleDataIDs[0:train_size]];
testKeys = hotEncodedKeys[possibleDataIDs[train_size: train_size + test_size]].astype(float);
testLabels = binaryLabels[possibleDataIDs[train_size: train_size + test_size]];

#Transform into tensors
TFkeysTrain = tf.data.Dataset.from_tensor_slices(trainKeys)
TFvalTrain = tf.data.Dataset.from_tensor_slices(trainLabels)
TFkeysTest = tf.data.Dataset.from_tensor_slices(testKeys)
TFvalTest = tf.data.Dataset.from_tensor_slices(testLabels)

#Transforms them into input and output elements (together) 
tfTrainElements = tf.data.Dataset.zip((TFkeysTrain, TFvalTrain))
tfTestElements = tf.data.Dataset.zip((TFkeysTest, TFvalTest))



#Finally, this is batched so can be used inside the training. 
train_dataset = tfTrainElements.shuffle(train_size).batch(batch_size, drop_remainder=True)
test_dataset = tfTestElements.shuffle(test_size).batch(batch_size, drop_remainder=True)

print('Obtained: train=', len(list(train_dataset)), ' test=', len(list(test_dataset)))



# In[19]:


#This is creating an additional test dataset, comparing Mascotte versus Loosers...
#Should make a function at some point...
fileToOpen2 = LocationDataFiles + antigenID + "_Task1aMvsL_SlicesBalancedData.txt"

hardTestDataset = open(fileToOpen2, newline = '')   #one line is a text with \t and \n                                                                              
data_reader = csv.reader(hardTestDataset, delimiter='\t') #transform lines into lists 
sequencesExternalTest = []
labelsExternalTest = []
for line in data_reader:
    if(not (line[0].startswith("#") or line[0].startswith("Slide"))):
        sequencesExternalTest.append(line[0])
        labelsExternalTest.append(line[1])


# In[20]:


#print(sequencesExternalTest[0:5])
#print(labelsExternalTest[0:5])


# In[21]:


#This step is also a little long
hotEncodedKeysExternal = np.array(list(map(hotEncodingAAString, sequencesExternalTest)))
binaryLabelsExternal = np.array([binary[item] for item in labelsExternalTest])
TFkeysExt = tf.data.Dataset.from_tensor_slices(hotEncodedKeysExternal)
TFvalExt = tf.data.Dataset.from_tensor_slices(binaryLabelsExternal)
tfElementsExt = tf.data.Dataset.zip((TFkeysExt, TFvalExt))
external_dataset = tfElementsExt.shuffle(len(sequencesExternalTest)).batch(batch_size, drop_remainder=True)


# In[22]:


if(not(runningInCommandLine)):
    print(sequences[0:10])
    print(labels[0:10], binaryLabels[0:10])
    print(len(sequences), len(labels))
    print(hotEncodingAAString(sequences[0]))


# In[23]:


print('ML Task 1, Antigen ', antigenID, ' nNeurons=', nNeurons, 'nSeqTot=', nSeqTot, 'nRepeats=', nRepeats, 'condition=', condition, 'layers=', layers, 'batch_size=', batch_size, 'external_dataset', len(sequencesExternalTest), 'test_max_20pct_of', refMaxSeqFor20pctsTest)


# In[24]:


#Finally, this is batched so can be used inside the training. 
train_dataset = tfTrainElements.shuffle(train_size).batch(batch_size, drop_remainder=True)
test_dataset = tfTestElements.shuffle(test_size).batch(batch_size, drop_remainder=True)

# Rebuild the model from scratch [not sure how initialized]
model = tf.keras.Sequential();
model.add(tf.keras.layers.Flatten());

if(nNeurons > 1):
	model.add(tf.keras.layers.Dense(nNeurons, activation='relu'))

if(layers > 1):
	model.add(tf.keras.layers.Dense(nNeurons, activation='relu'))

if(layers > 2):
	model.add(tf.keras.layers.Dense(nNeurons, activation='relu'))

if(layers > 3):
	model.add(tf.keras.layers.Dense(nNeurons, activation='relu'))

model.add(tf.keras.layers.Dense(1, activation='sigmoid'))
print('Sigmoid activation function')


# In[25]:


#Now try with more binary classification stuff
#nEpochs = 40
optimizer = tf.keras.optimizers.Adam()
loss = 'binary_crossentropy'
metrics = ['accuracy', 'FalseNegatives', 'FalsePositives', 'Precision', 'Recall', 'TrueNegatives', 'TruePositives']
model.compile(loss=loss, optimizer=optimizer, metrics=metrics)


# In[26]:


input_shape = [1, 11, 20]
model.build(input_shape)
model.summary()


# In[27]:


history = model.fit(train_dataset, epochs=nEpochs, validation_data=test_dataset, validation_steps=1)


# In[28]:


history_dict = history.history
print(history_dict.keys())


# In[29]:


loss_values = history_dict['loss']
val_loss_values = history_dict['val_loss']
accuracy = history_dict['val_accuracy']
epochs = range(1, nEpochs+1)

if(not(runningInCommandLine)):
	plt.plot(epochs, loss_values, 'bo', label='Training loss')
	plt.plot(epochs, val_loss_values, 'b', label='Validation loss')
	plt.plot(epochs, accuracy, 'r', label='Accuracy')
	plt.title('Training and validation loss + validation accuracy')
	plt.xlabel('Epochs')
	plt.ylabel('Loss')
	plt.legend()
	plt.show()


# In[30]:


print('Evaluation:', model.evaluate(test_dataset))


# In[31]:


print('Evaluation External:', model.evaluate(external_dataset))


# In[32]:


#file_object = open('History.txt', 'a')
#file_object.write(antigenID + "\t" + str(nEpochs) + "\t" + str(nNeurons) + "\t" + str(nSeqTot) + "\t" + str(repeat) + "\t" + str(nRepeats)+ "\t" + str(condition) + "\t" + str(batch_size) + "\t" + str(layers) + '\t' + str(model.evaluate(test_dataset)) + "\n")
#file_object.close()


# In[33]:


def getDataset(dataset):
	testLabels = []
	testFeatures = []
	#i = 0
	for (features, labels) in dataset:
		for (seq , en1) in zip(features , labels):   # for each batch
			testFeatures.append(seq)
			testLabels.append(en1.numpy().astype(float)); # tf.dtypes.cast(en, tf.int32))
	return (testFeatures, testLabels)
			


# In[34]:


#if you want to retrieve the datasets
#(trainX, trainY) = getDataset(train_dataset)
#(testX, testY) = getDataset(test_dataset)
#decodedTrainX = list(map(retrieveString, trainX))
#decodedTestX = list(map(retrieveString, testX))


# In[35]:


#print(decodedTrainX[0:10])
#print(decodedTestX[0:10])


# In[36]:


#print(len(decodedTrainX), len(decodedTestX))
#a = np.ndarray.flatten(np.array(decodedTrainX))
#b = np.ndarray.flatten(np.array(decodedTestX))
#common_elements = list(set(a).intersection(set(b)))
#print(len(common_elements))
#common_elements[0:10]


# In[37]:


def myRound(logit):
	if(logit < 0.5):
		return 0
	return 1


# In[38]:


# Note: the lists of labels to build these datasets were:
#   trainLabels
#   testLabels
train_dataset = tfTrainElements.batch(batch_size, drop_remainder=True)
test_dataset = tfTestElements.batch(batch_size, drop_remainder=True)

pred_train = model.predict(train_dataset)
pred_test = model.predict(test_dataset)

classes_train = list(map(myRound, pred_train))
classes_test = list(map(myRound, pred_test))



# In[39]:


con_mat_train = tf.math.confusion_matrix(labels=trainLabels[0:len(classes_train)], predictions=classes_train).numpy()
print(con_mat_train)


# In[40]:


con_mat_test = tf.math.confusion_matrix(labels=testLabels[0:len(classes_test)], predictions=classes_test).numpy()
print(con_mat_test)




# In[42]:


fpr, tpr, _ = roc_curve(testLabels, pred_test)
roc_auc = auc(fpr, tpr)
print(roc_auc)


# In[43]:

# file_object = open('History.txt', 'a')
# file_object.write(antigenID + "\t" + str(nEpochs) + "\t" + str(nNeurons) + "\t" + str(nSeqTot) + "\t" + str(repeat) + "\t" + str(nRepeats)+ "\t" + str(condition) + "\t" + str(batch_size) + "\t" + str(layers) + '\t' + str(model.evaluate(test_dataset)) + '\t' + str(roc_auc) + '\t' + str(model.evaluate(train_dataset)) + "\t" + str(model.evaluate(external_dataset)) + "\n")
# file_object.close()

if(not(runningInCommandLine)):
	plt.figure()
	lw = 2
	plt.plot(fpr, tpr, color='darkorange',lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
	plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
	plt.xlim([0.0, 1.0])
	plt.ylim([0.0, 1.05])
	plt.xlabel('False Positive Rate')
	plt.ylabel('True Positive Rate')
	plt.title('Receiver operating characteristic example')
	plt.legend(loc="lower right")
	plt.show()



# convert to pyTorch
# model.get_config()
# model.summary()
# model.trainable_weights[i].numpy()
model.save('./model_d1_balanced/')
# loaded_model = tf.keras.models.load_model('/home/roberfra/test1/')
# In[ ]:




# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




