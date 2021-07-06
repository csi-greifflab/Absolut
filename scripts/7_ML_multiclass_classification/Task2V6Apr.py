#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python3

#Script that performs ML Task 1 for an antigen. 
#Usage in command line:
#python ThisFile.py nSelectAntigens=5 nNeurons=5 nSeqTot=1e7 nRepeats=10 condition=1 layers=1 batch_size=min(nSeqTot, sequences in dataset)

from __future__ import division, print_function, absolute_import


# In[2]:


#Local settings - Files necessary in that folder: "
LocationDataFiles = ""
runningInCommandLine = True
nEpochs = 150
repeat = 1   #for command line


# In[3]:


#Default parameters (changeable by command line)
nSelectAntigens = 50
nNeurons = 50
nSeqTot = 10000
refMaxSeqFor20pctsTest = 500000
nRepeats = 10
layers = 1
condition = 1
useAACompo = False

# In[4]:



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


# In[5]:


#Hides the presence of the GPU, to avoid problems on immunhub
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"


# In[6]:


#Parameters from command line. This cell does not work from Python notebook
#print(len(sys.argv), str(sys.argv))
if(runningInCommandLine):
    if(len(sys.argv) > 1):
        nSelectAntigens = int(sys.argv[1])

    if(len(sys.argv) > 2):
        nNeurons = int(sys.argv[2])


    if(len(sys.argv) > 3):
        nSeqTot = int(sys.argv[3])


    if(len(sys.argv) > 4):
        nRepeats = int(sys.argv[4])

    # 1 = only binders to exactly one of the selected antigen (maybe need to balance)
    # 3 = same + shuffling labels
    
    # 21 = multi-labels allowed
    # 23 = same + shuffling labels

    # 31 = only binders to exactly 1 + include 50% non-binders (but to other antigens) in the training
    # 33 = same + shuffling labels


    if(len(sys.argv) > 5):
        condition = int(sys.argv[5])

    if(condition > 10):
        useAACompo = True
        condition = condition - 10

    if(len(sys.argv) > 6):
        layers = int(sys.argv[6])
		
#see below for argument 7 and 8


# In[7]:


#Defining the one hot encoding
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


# In[8]:


if(not(runningInCommandLine)):
    test = hotEncodingAAString("ACD")
    print(test)
    got = retrieveString(test)
    print(got)


# In[9]:


#Reading the data files and loading as two vectors, depending on antigenID and on the condition (balanced/imbalanced) Shuffling is made later
fileToOpen = LocationDataFiles + "Treated142.txt" #"Task2_nBind1OrMore.txt"
#Will do the shuffling
#if((condition == 3) or (condition == 13)):
#   fileToOpen = LocationDataFiles + antigenID + "Task2_nBind1OrMore_Shuffled.txt"

balancedDataset = open(fileToOpen, newline = '')   #one line is a text with \t and \n                                                                              
data_reader = csv.reader(balancedDataset, delimiter='\t') #transform lines into lists 
sequences = []
labels = []
for line in data_reader:
    if(not (line[0].startswith("#") or line[0].startswith("Slice"))):
        sequences.append(line[0])
        labels.append(line[2])

nAntigens = len(labels[0])

# In[11]:


#Note: the IDs start at 0
possibleIDs = [*range(0,nAntigens)]
random.shuffle(possibleIDs)

print(nAntigens, nSelectAntigens)
selected = possibleIDs[0:min(nAntigens, nSelectAntigens)]
print(selected)


# In[12]:


#note: sequence should be of size nAntigens!! will not check
binary = {"1": 1,"0": 0}
def extractNDimLabels(sequence, listSelectedAntigens):
    justPositions = ''.join(sequence[x] for x in listSelectedAntigens)
    listLabels = [binary[x] for x in justPositions]
    return listLabels


# In[13]:


#This is just a test that the function extractNDimLabels works
if(not(runningInCommandLine)):
    print(extractNDimLabels("100001000010000100001", [0,5,10,1,6,11]))


# In[14]:


multiDimensionalLabels = np.array([extractNDimLabels(code, selected) for code in labels])


# In[15]:


multiDimensionalLabels.shape[0]


# In[16]:


#Conversion of the dataset into One Hot Encoding- takes 20 sec
hotEncodedKeys = np.array(list(map(hotEncodingAAString, sequences)))

if(useAACompo == True):
    hotEncodedKeys = np.sum(hotEncodedKeys, axis = 2)


# In[17]:


#len(hotEncodedKeys)
#hotEncodedKeys.shape()
lenHotEncodedKeys = hotEncodedKeys.shape[0]



# In[19]:


if(not(runningInCommandLine)):
    print(sequences[0:10])
    print(labels[0:10], multiDimensionalLabels[0:10])
    print(len(sequences), len(labels))
    print(hotEncodingAAString(sequences[0]))


# In[20]:


#Now we will need to take only indices that do bind
listSumTargets = list(map(sum, multiDimensionalLabels))
print("Total Sequences=", len(listSumTargets), "Bining none of selected antigens:", listSumTargets.count(0))
print("Binding 1:", listSumTargets.count(1), " 2:", listSumTargets.count(2), " 3:", listSumTargets.count(3), " 4:", listSumTargets.count(4), " 5:", listSumTargets.count(5))

listBinders = [i for i, value in enumerate(listSumTargets) if value >= 1] 
print(len(listBinders), listBinders[0:10])
listExactBinders = [i for i, value in enumerate(listSumTargets) if value == 1] 
print(len(listExactBinders), listExactBinders[0:10])
listNonBinders = [i for i, value in enumerate(listSumTargets) if value == 0] 
print(len(listNonBinders), listNonBinders[0:10])


# In[21]:


random.shuffle(listNonBinders)
#Here, we take only eact binders, this is multi-class (exclusive)
possibleDataIDs = listExactBinders
#+ listNonBinders[0:min(len(listNonBinders), len(listExactBinders))]

#Here, we could also include multibinders only, this is multi-label
#possibleDataIDs = listBinders + listNonBinders[0:min(len(listNonBinders), len(listBinders))]



#Will separate train and test datasets now, before stupid tensorflow datasets.
#possibleDataIDs = [*range(0,len(hotEncodedKeys))]
random.shuffle(possibleDataIDs)
print("total available sequences after balancing:", len(possibleDataIDs))


# In[ ]:





# In[22]:


#If too many sequences are requested, will downscale, just to have this info in the output
nSeqTot = min(nSeqTot,len(possibleDataIDs))
#refMaxSeqFor20pctsTest = min(refMaxSeqFor20pctsTest, len(possibleDataIDs))


# In[23]:


nSeqTot = min(nSeqTot,len(possibleDataIDs))
train_size = int(0.8 * nSeqTot) 
test_size = int(0.2 * refMaxSeqFor20pctsTest)
if test_size + train_size > len(possibleDataIDs):
   test_size = min(test_size, int(0.2*len(possibleDataIDs)))
print('Requested: dataset=', nSeqTot, 'train=', train_size, 'test=', test_size, 'totUsed=', train_size + test_size)


# In[24]:


#By default, we will take 50 batches, so the size of a batch is nr sequences / 50
batch_size = math.floor(nSeqTot / 50)
if(runningInCommandLine):
    if(len(sys.argv) > 7):
        batch_size = int(sys.argv[7])


# In[25]:


print('ML Task 2, Antigens ', selected, ' nNeurons=', nNeurons, 'nSeqTot=', nSeqTot, 'nRepeats=', nRepeats, 'condition=', condition, 'batch_size=', batch_size)


# In[26]:


# At this stage, all the data is prepared, and the requested size of batches is defined. 


# In[27]:


#In command line, here would be the repeat loop
#for repeat in range(1,nRepeats+1): #the last number is NOT reached

#Shuffle the labels

for repeat in range(nRepeats):
	if((condition == 3) or (condition == 13)):
		random.shuffle(multiDimensionalLabels)




	#Transform the lists into tensorflow datasets (not yet batched) - this step takes time so I do it only once (2 mins for 100 000 sequences)
	trainKeys = hotEncodedKeys[possibleDataIDs[0:train_size]].astype(float);
	trainLabels = multiDimensionalLabels[possibleDataIDs[0:train_size]];
	testKeys = hotEncodedKeys[possibleDataIDs[train_size: train_size + test_size]].astype(float);
	testLabels = multiDimensionalLabels[possibleDataIDs[train_size: train_size + test_size]];

	#Transform into tensors
	TFkeysTrain = tf.data.Dataset.from_tensor_slices(trainKeys)
	TFvalTrain = tf.data.Dataset.from_tensor_slices(trainLabels)
	TFkeysTest = tf.data.Dataset.from_tensor_slices(testKeys)
	TFvalTest = tf.data.Dataset.from_tensor_slices(testLabels)

	#Transforms them into input and output elements (together) 
	tfTrainElements = tf.data.Dataset.zip((TFkeysTrain, TFvalTrain))
	tfTestElements = tf.data.Dataset.zip((TFkeysTest, TFvalTest))

	train_dataset = tfTrainElements.shuffle(train_size).batch(batch_size, drop_remainder=True)
	test_dataset = tfTestElements.shuffle(test_size).batch(batch_size, drop_remainder=True)

	print('Obtained: train=', len(list(train_dataset)), ' test=', len(list(test_dataset)))


	# Rebuild the model from scratch [not sure how initialized]
	model = tf.keras.Sequential();
	model.add(tf.keras.layers.Flatten());
	if(nNeurons > 1):
		model.add(tf.keras.layers.Dense(nNeurons, activation='relu'))
		
	if( (nNeurons > 1) and (layers > 1)):
		model.add(tf.keras.layers.Dense(nNeurons, activation='relu'))
	if( (nNeurons > 1) and (layers > 2)):
		model.add(tf.keras.layers.Dense(nNeurons, activation='relu'))
	if( (nNeurons > 1) and (layers > 3)):
		model.add(tf.keras.layers.Dense(nNeurons, activation='relu'))
		
	model.add(tf.keras.layers.Dense(len(selected), activation='softmax'))

	#Epochs = 20
	optimizer = tf.keras.optimizers.Adam()
	#opt = SGD(lr=0.01, momentum=0.9)
	loss = 'categorical_crossentropy'
	#loss = 'binary_crossentropy'
	metrics = ['accuracy', 'FalseNegatives', 'FalsePositives', 'Precision', 'Recall', 'TrueNegatives', 'TruePositives']
	model.compile(loss=loss, optimizer=optimizer, metrics=metrics)


	input_shape = [1, 11, 20]
	if useAACompo == True:
		input_shape = [1, 20]
    
	model.build(input_shape)
	model.summary()


	history = model.fit(train_dataset, epochs=nEpochs, validation_data=test_dataset, validation_steps=10)


	history_dict = history.history
	#print(history_dict.keys())




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


	print('Evaluation:', model.evaluate(test_dataset))


	train_dataset = tfTrainElements.batch(batch_size, drop_remainder=True)
	test_dataset = tfTestElements.batch(batch_size, drop_remainder=True)

	pred_train = model.predict(train_dataset)
	pred_test = model.predict(test_dataset)


	def myRound(logit):
		if(logit <= 0.5):
			return 0
		return 1

	def getIDClass(floatingMultidimentionalLabel):
		multidimentionalLabel = list(map(myRound, floatingMultidimentionalLabel))
		dim = len(multidimentionalLabel)
		if(sum(multidimentionalLabel) > 1):
			return dim + 1
			print("Error, multilabel:" + str(multidimentionalLabel))
		return sum(multidimentionalLabel * np.array(range(1,dim+1)))    # starts at 0

	def getIDClassBiggestVal(floatingMultidimentionalLabel):
		return np.argmax(floatingMultidimentionalLabel, axis=0)


	# In[40]:


	# reminder how they were generated
	#trainKeys = hotEncodedKeys[possibleDataIDs[0:train_size]].astype(float);
	#trainLabels = multiDimensionalLabels[possibleDataIDs[0:train_size]];
	#testKeys = hotEncodedKeys[possibleDataIDs[train_size: train_size + test_size]].astype(float);
	#testLabels = multiDimensionalLabels[possibleDataIDs[train_size: train_size + test_size]];


	# In[41]:


	takeMoreHalf = False
	if(takeMoreHalf):
		testLabels1D = np.array(list(map(getIDClass, testLabels)))
		testPredicted1D = np.array(list(map(getIDClass, model.predict(test_dataset))))
		print("Test datasets: labels")
		print(testLabels[0:25])
		print(testLabels1D[0:50])
		print("Predicted from the model from the test dataset")
		print(testPredicted1D[0:50])
		con_mat = tf.math.confusion_matrix(labels=testLabels1D, predictions=testPredicted1D).numpy()
		print(con_mat)
		
		print('Evaluation:', model.evaluate(test_dataset))
		print(accuracy_score(testLabels1D[0:len(testPredicted1D)], testPredicted1D))

	from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score






# In[43]:


	if(not takeMoreHalf):
		testLabels1D = np.array(list(map(getIDClassBiggestVal, 0.5*testLabels)))
		testPredicted1D = np.array(list(map(getIDClassBiggestVal, model.predict(test_dataset))))
		
		print("Test datasets: labels")
		#print(testLabels[0:25])
		print(testLabels1D[0:50])
		print("Predicted from the model from the test dataset")
		print(testPredicted1D[0:50])
		con_mat = tf.math.confusion_matrix(labels=testLabels1D[0:len(testPredicted1D)], predictions=testPredicted1D).numpy()
		print(con_mat)
		
		print('Evaluation:', model.evaluate(test_dataset))
		print(accuracy_score(testLabels1D[0:len(testPredicted1D)], testPredicted1D))


	# In[44]:


	if(not(runningInCommandLine)):
		#https://deeplizard.com/learn/video/km7pxKy4UHU
		import itertools
		def plot_confusion_matrix(cm, classes,
								normalize=False,
								title='Confusion matrix',
								cmap=plt.cm.Blues):
			"""
			This function prints and plots the confusion matrix.
			Normalization can be applied by setting `normalize=True`.
			"""
			plt.figure(figsize = (25,25))
			plt.imshow(cm, interpolation='nearest', cmap=cmap)
			plt.title(title)
			plt.colorbar()
			tick_marks = np.arange(len(classes))
			plt.xticks(tick_marks, classes, rotation=45)
			plt.yticks(tick_marks, classes)

			if normalize:
				cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
				print("Normalized confusion matrix")
			else:
				print('Confusion matrix, without normalization')

			print(cm)

			thresh = cm.max() / 2.
			for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
				if(cm[i, j] >= 0.01):
					if normalize:
						plt.text(j, i, int(100*np.around(cm[i, j], decimals=2)),
							horizontalalignment="center",
							color="white" if cm[i, j] > thresh else "black")
					else:
						plt.text(j, i, cm[i, j],
							horizontalalignment="center",
							color="white" if cm[i, j] > thresh else "black")

			plt.tight_layout()
			plt.ylabel('True label')
			plt.xlabel('Predicted label')


	# In[45]:


	if(not(runningInCommandLine)):
		plot_confusion_matrix(con_mat, classes = [*range(0,50)])


	# In[46]:


	if(not(runningInCommandLine)):
		plot_confusion_matrix(con_mat, classes = [*range(0,50)], normalize=True)


	#https://gist.github.com/RyanAkilos/3808c17f79e77c4117de35aa68447045
	y_test = testLabels1D[0:len(testPredicted1D)]
	y_pred = testPredicted1D

	MiPrec = precision_score(y_test, y_pred, average='micro')
	MiRec = recall_score(y_test, y_pred, average='micro')
	MiF1 = f1_score(y_test, y_pred, average='micro')
	print('Micro Precision: {:.2f}'.format(MiPrec))
	print('Micro Recall: {:.2f}'.format(MiRec))
	print('Micro F1-score: {:.2f}\n'.format(MiF1))

	MaPrec = precision_score(y_test, y_pred, average='macro')
	MaRec = recall_score(y_test, y_pred, average='macro')
	MaF1 = f1_score(y_test, y_pred, average='macro')
	print('Macro Precision: {:.2f}'.format(MaPrec))
	print('Macro Recall: {:.2f}'.format(MaRec))
	print('Macro F1-score: {:.2f}\n'.format(MaF1))

	WeiPrec = precision_score(y_test, y_pred, average='weighted')
	WeiRec = recall_score(y_test, y_pred, average='weighted')
	WeiF1 = f1_score(y_test, y_pred, average='weighted')
	print('Weighted Precision: {:.2f}'.format(WeiPrec))
	print('Weighted Recall: {:.2f}'.format(WeiRec))
	print('Weighted F1-score: {:.2f}'.format(WeiF1))

	from sklearn.metrics import classification_report
	print('\nClassification Report\n')
	print(classification_report(y_test, y_pred))


	# In[48]:


	#nSelectAntigens = 50
	#nNeurons = 50
	#nSeqTot = 10000
	#refMaxSeqFor20pctsTest = 500000
	#nRepeats = 10
	#layers = 1
	#condition = 1
	file_object = open('HistoryTask2.txt', 'a')
	file_object.write(str(nSelectAntigens) + "\t" + str(nEpochs) + "\t" + str(nNeurons) + "\t" + str(nSeqTot) + "\t" + str(repeat) + "\t" + str(nRepeats)+ "\t" + str(condition) + "\t" + str(layers) + '\t' + str(batch_size) + "\t" + str(model.evaluate(test_dataset)) + "\t" + str(model.evaluate(train_dataset)) + '\t' + str(MiPrec)  + '\t' + str(MiRec) + '\t' + str(MiF1) + '\t' + str(MaPrec) + '\t' + str(MaRec) + '\t' + str(MaF1) + '\t' + str(WeiPrec) + '\t' + str(WeiRec) + '\t' + str(WeiF1) + "\n")
	file_object.close()
	
	file_object = open('Confusion.txt', 'a')
	file_object.write('#' + str(nEpochs) + "\t" + str(nNeurons) + "\t" + str(nSeqTot) + "\t" + str(repeat) + "\t" + str(nRepeats)+ "\t" + str(condition) + "\t" + str(layers) + '\t' + str(batch_size) + '\t' + str(con_mat.shape) + str(selected) + '\n' + str(con_mat) + '\n');
	file_object.close()
	# In[ ]:





	# In[49]:


	#Multiclass ROC is not supported
	#from sklearn.metrics import roc_curve
	#from sklearn.metrics import auc

	#fpr, tpr, _ = roc_curve(testLabels1D[0:len(testPredicted1D)], testPredicted1D)
	#roc_auc = auc(fpr, tpr)
	#print(roc_auc)

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



