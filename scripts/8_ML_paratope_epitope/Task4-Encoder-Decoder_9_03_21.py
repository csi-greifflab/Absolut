#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python3
# coding: utf-8

from __future__ import division, print_function, absolute_import


# In[ ]:


# Usage: 
# ThisScript.py embedding_dim=6  units=500  nSeqTot=1000  nRepeats=10  condition=1  wordSizeIn=1  wordSizeOut=1  
#             fileInput=seeBelow  colInput=0   colOutput=1  simSuffix=""   batch_size=nSeqTot/50  
# 
# Requires: 
# an input file with two columns: input and output. 
# the location is defined just below:
#
# Outputs:
#    creates a folder unique to this simulation
#    puts the weights of the fitted model inside


# In[ ]:


# local settings
runningInCommandLine = True


# In[ ]:


import tensorflow as tf
import collections
import os
import random
import urllib
import zipfile

import numpy as np
import tensorflow as tf
import csv

import time

import sys
import math


# In[ ]:


if(not(runningInCommandLine)):
    import matplotlib.pyplot as plt

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

#to get an ID of simulation
from datetime import datetime
from random import seed
import random
import jellyfish


# In[ ]:


#Default parameters (changeable by command line)
#in the same order of command line arguments

embedding_dim = 6
units = 500            #This is the hidden dimension, both encoder and decoder
nSeqTot = 10000
nRepeats = 10
condition = 1
wordSizeIn = 1
wordSizeOut = 1
FileInput = "Task4_E_EpiSeq_ParaSeq_NoDeg.txt"
colInput = 0
colOutput = 1
simSuffix = ""
batch_size = 200


# In[ ]:


#Non-command line options
refMaxSeqFor20pctsTest = 1000000
paddingInputs = 'post'
paddingOutputs = 'post'
EPOCHS = 250
removeStar = True   #Paratopes are distinguished by ending by a star. It means we will remove ending star from input and ouptuts sequences if they have it 


# In[ ]:


#get a simID to save the model weights
random.seed(datetime.now())
simID = random.randint(1,9999) #both numbers included
print("This run has ID",str(simID))
#Note: the repeats will be called ID.1, etc...


# In[ ]:


if(runningInCommandLine):
    if(len(sys.argv) > 1):
        embedding_dim = int(sys.argv[1])

    if(len(sys.argv) > 2):
        units = int(sys.argv[2])


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

    if(len(sys.argv) > 6):
        wordSizeIn = int(sys.argv[6])
        
    if(len(sys.argv) > 7):
        wordSizeOut = int(sys.argv[7])
            
    if(len(sys.argv) > 8):
        FileInput = sys.argv[8]
        print("Input file has been manually given to ", FileInput)

    if(len(sys.argv) > 9):
        colInput = int(sys.argv[9])
        print("Input column assign to be ", colInput)

    if(len(sys.argv) > 10):
        colOutput = int(sys.argv[10])
        print("Input column assign to be ", colOutput)

    if(len(sys.argv) > 11):
        simSuffix = sys.argv[11]
        print("Input column assign to be ", simSuffix)

    # argument 12 is the batch_size, see below


# In[7]:


def RemoveEndingStar(text):
    if(len(text) > 0):
        if(text[-1] == '*'):
            text = text[:-1]
    return(text)


# In[ ]:


def getFeaturesLabels(fileName, colI=0, colO=1):
    epiParaSeq = open(fileName, newline = '')   #one line is a text with \t and \n                                                                              
    data_reader = csv.reader(epiParaSeq, delimiter='\t') #transform lines into lists 
    sequences = []
    labels = []
    for line in data_reader:
        #print(line)
        if(not (line[0].startswith("#") or line[0].startswith("seqAGEpitope"))):
            if(removeStar):
                sequences.append(RemoveEndingStar(line[colI]))
                labels.append(RemoveEndingStar(line[colO]))
            else:
                sequences.append(line[colI])
                labels.append(line[colO])
    
    print("Got ", len(sequences), " sequences and ", len(labels), " from ", fileName)
    return (sequences, labels)


# In[ ]:


#Input: array of sequences, array of sequences (same size), size of words in sequences, size of words in labels
def createSentencePairs(sequences, labels, wordSizeSeq, wordSizeLabels):
    if(len(sequences) != len(labels)):
        print("ERR: createSentencePairs, different size labels and sequences")
        return()
    
    print("First 5 sentence pairs generated:")
    pairs = []
    for i in range(0,len(sequences)):
        inp = sequences[i]
        outp= labels[i]
        wordsIn = [inp[i:i+wordSizeIn] for i in range(0, len(inp), wordSizeIn)]
        if(wordSizeIn < 0):
            wordsIn = inp.split(' ')
        wordsOut = [outp[i:i+wordSizeOut] for i in range(0, len(outp), wordSizeOut)]
        if(wordSizeOut < 0):
            wordsOut = outp.split(' ')
        pair = ['<start> '+' '.join(list(wordsIn))+' <end>', '<start> '+' '.join(list(wordsOut))+' <end>']
        if(i < 5):
            print(pair)
            
        pairs.append(pair)
    
    return pairs


# In[ ]:


#Generates the list of words in sequences and labels
def getVocabulary(pairs):
    vocabIn = set()
    vocabOut = set()
    for couple in pairs:
        vocabIn.update(couple[0].split(' '))
        vocabOut.update(couple[1].split(' '))

    vocabIn = sorted(vocabIn)
    vocabOut = sorted(vocabOut)

    print("Feature vocabulary: ", len(vocabIn), "words: ", vocabIn)
    print("Labels  vocabulary: ", len(vocabOut), "words: ", vocabOut)
    
    return (vocabIn, vocabOut)


# In[ ]:


#Note: padding is included in the dictionnary manually.
def generateDico(vocab, addPadding = True):
    word2idx = {}
    idx2word = {}
    if(addPadding):
        word2idx['<pad>'] = 0

    for index, word in enumerate(vocab):
        word2idx[word] = index + 1

    for word, index in word2idx.items():
        idx2word[index] = word

    return (word2idx, idx2word) 


# In[ ]:


def processFile(fileName, wordSizeSeq = 1, wordSizeLabels = 1, colI = 0, colO = 1, paddingOptionInputs = 'post', paddingOptionOutputs = 'post'):
    (sequences, labels) = getFeaturesLabels(fileName, colI, colO)
    featurePairs = createSentencePairs(sequences, labels, wordSizeSeq, wordSizeLabels) 
    (vocabIn, vocabOut) = getVocabulary(featurePairs)
    
    (word2idxIN, idx2wordIN) = generateDico(vocabIn)
    (word2idxOUT, idx2wordOUT) = generateDico(vocabOut)
    
    input_tensor = [[word2idxIN[s] for s in first.split(' ')] for first, second in featurePairs]
    output_tensor = [[word2idxOUT[s] for s in second.split(' ')] for first, second in featurePairs]

    NI = max(len(t) for t in input_tensor)
    NO = max(len(t) for t in output_tensor)
    print("Max length input: ", NI, "Max length output:", NO)

    paddedInputTensor = tf.keras.preprocessing.sequence.pad_sequences(input_tensor,maxlen=NI,padding=paddingOptionInputs)
    paddedOutputTensor = tf.keras.preprocessing.sequence.pad_sequences(output_tensor,maxlen=NO,padding=paddingOptionOutputs)

    return(paddedInputTensor, paddedOutputTensor, word2idxIN, idx2wordIN, word2idxOUT, idx2wordOUT, NI, NO)


# In[ ]:


def padNewSentences(sentences, maxLength, word2IdDictionary, paddingOption):
    input_tensor = [[word2IdDictionary[s] for s in sentence.split(' ')] for sentence in sentences]
    paddedSequences = tf.keras.preprocessing.sequence.pad_sequences(input_tensor,maxlen=maxLength,padding=paddingOption)
    return paddedSequences


# In[ ]:


class Encoder(tf.keras.Model):
    def __init__(self, vocab_size, embedding_dim, enc_units, batch_sz):
        super(Encoder, self).__init__()
        self.batch_sz = batch_sz
        self.enc_units = enc_units
        self.embedding = tf.keras.layers.Embedding(vocab_size, embedding_dim)
        self.gru = tf.keras.layers.GRU(self.enc_units,
                                       return_sequences=True,
                                       return_state=True,
                                       recurrent_initializer='glorot_uniform')
        print("Initializing Encoder. Embedding into DE=", embedding_dim, ", GRU with hidden=", self.enc_units)

    def call(self, x, hidden):
        x = self.embedding(x)
        output, state = self.gru(x, initial_state = hidden)
        return output, state

    def initialize_hidden_state(self):
        return tf.zeros((self.batch_sz, self.enc_units))


# In[ ]:


class BahdanauAttention(tf.keras.layers.Layer):
    def __init__(self, units):
        super(BahdanauAttention, self).__init__()
        self.W1 = tf.keras.layers.Dense(units)
        self.W2 = tf.keras.layers.Dense(units)
        self.V = tf.keras.layers.Dense(1)
        print("Initializing BahdanauAttention for ", units, " hidden dimensions")

    def call(self, query, values):
        # query hidden state shape == (batch_size, hidden size)
        # query_with_time_axis shape == (batch_size, 1, hidden size)
        # values shape == (batch_size, max_len, hidden size)
        # we are doing this to broadcast addition along the time axis to calculate the score
        query_with_time_axis = tf.expand_dims(query, 1)

        # score shape == (batch_size, max_length, 1)
        # we get 1 at the last axis because we are applying score to self.V
        # the shape of the tensor before applying self.V is (batch_size, max_length, units)
        score = self.V(tf.nn.tanh(
            self.W1(query_with_time_axis) + self.W2(values)))

        # attention_weights shape == (batch_size, max_length, 1)
        attention_weights = tf.nn.softmax(score, axis=1)

        # context_vector shape after sum == (batch_size, hidden_size)
        context_vector = attention_weights * values
        context_vector = tf.reduce_sum(context_vector, axis=1)

        return context_vector, attention_weights


# In[ ]:


class Decoder(tf.keras.Model):
    def __init__(self, vocab_size, embedding_dim, dec_units, batch_sz):
        super(Decoder, self).__init__()
        self.batch_sz = batch_sz
        self.dec_units = dec_units
        self.embedding = tf.keras.layers.Embedding(vocab_size, embedding_dim)
        self.gru = tf.keras.layers.GRU(self.dec_units,
                                   return_sequences=True,
                                   return_state=True,
                                   recurrent_initializer='glorot_uniform')
        self.fc = tf.keras.layers.Dense(vocab_size)

        # used for attention
        self.attention = BahdanauAttention(self.dec_units)
        print("Initializing Decoder. Embedding from vocab size ", vocab_size, " into DE=", embedding_dim, ", GRU with hidden=", self.dec_units,
          " hidden dimensions; final layer outputs O=", vocab_size, " dimensions (alphabet)")

    def call(self, x, hidden, enc_output):
        # enc_output shape == (batch_size, max_length, hidden_size)
        context_vector, attention_weights = self.attention(hidden, enc_output)

        # x shape after passing through embedding == (batch_size, 1, embedding_dim)
        x = self.embedding(x)

        # x shape after concatenation == (batch_size, 1, embedding_dim + hidden_size)
        x = tf.concat([tf.expand_dims(context_vector, 1), x], axis=-1)
        #print(tf.reduce_max(x))

        # passing the concatenated vector to the GRU
        output, state = self.gru(x)

        # output shape == (batch_size * 1, hidden_size)
        output = tf.reshape(output, (-1, output.shape[2]))

        # output shape == (batch_size, vocab)
        x = self.fc(output)

        return x, state, attention_weights


# In[ ]:


def vocalise(tensor, dico, wordsize = 1):
    pulled = ''.join(dico[i] for i in tensor)
    if(wordsize < 0):
            pulled = ' '.join(dico[i] for i in tensor)
    return pulled.replace("<pad>", "").replace("<end>", "").replace("<start>", "")


# In[ ]:


def vocaliseBatch(batched_tensor, dico, wordsize = 1):
    result = []
    for tensor in batched_tensor:
        pulled = ''.join(dico[i] for i in tensor)
        if(wordsize < 0):
            pulled = ' '.join(dico[i] for i in tensor)
        text = pulled.replace("<pad>", "").replace("<start>", "")
        #Now, removes anything happening after the <end> 
        sep = '<end>'
        stripped = text.split(sep, 1)[0]
        result.append(stripped)
    
    return result


# In[ ]:


#Returns a list of 
def evaluateBatch(input_tensor, output_tensor, encoder, decoder, max_length_inp, max_length_targ, word2idIN, word2idOUT):
    input_batch_size = input_tensor.shape[0]
    
    #print("Evaluate ", input_tensor)
    attention_plot = np.zeros((max_length_targ, max_length_inp))

    # Calculate the encoder hidden states. Note: it doesn't take a python array but only a np array with a type (float)
    hidden = tf.zeros((input_batch_size, encoder.enc_units))
    
    #print("Input", input_tensor.shape)
    #print("Hidden", hidden.shape)
    
    enc_out, enc_hidden = encoder(input_tensor, hidden)

    dec_hidden = enc_hidden
    #dec_input = tf.expand_dims(tf.expand_dims([word2idOUT['<start>'], 1], 0))
    dec_input = tf.repeat(tf.expand_dims([word2idxOUT['<start>']], 0), input_batch_size, axis=0)
    #from the train function, dec_input = tf.expand_dims([word2idxOUT['<start>']] * batch_size, 1)
  
    
    
    #print("DecIn", dec_input.shape)
    
    result = tf.zeros((input_batch_size, 1), tf.int64) #np.empty([1, len(word2idOUT)]) 
    predictionEachPos = tf.zeros((input_batch_size, len(word2idOUT)))
    #print("PredEachPos", predictionEachPos.shape)
    
    loss = np.zeros((input_batch_size, 1))
    
    for t in range(max_length_targ):
        predictions, dec_hidden, attention_weights = decoder(dec_input, dec_hidden, enc_out)
        #print(dec_hidden)
        
        #print("predictions", predictions.shape)
        #print("dec_hidden", dec_hidden.shape)
        #print("attention_weights", attention_weights.shape)
    
        loss_t = loss_function(output_tensor[:,t], predictions)
        #print(output_tensor[:,t])
        #print("loss_t", loss_t.shape)
        
        loss += loss_t
        #print(loss)
        #print(loss_t.numpy())
        
        # storing the attention weights to plot later on
        ###attention_weights = tf.reshape(attention_weights, (-1,))
        ###attention_plot[t] = attention_weights.numpy()

        #print(predictions)
        predicted_id = tf.argmax(predictions, axis=1).numpy()
        #print("predicted_id", predicted_id.shape)
        #print(predicted_id)
        
        #print(result.shape)
        #print(tf.expand_dims(predicted_id, 1).shape)
        result = tf.concat([result, tf.expand_dims(predicted_id, 1)], axis = 1)
        #predictionEachPos = np.append(predictionEachPos, [predictions[0].numpy()])
        #tf.concat([predictionEachPos, tf.expand_dims(predictions[:,0].numpy(), 0)], axis = 0)

        #Here, we will predict after stop, and cut later. How do we do for loss, don't know.
        #if predicted_id == word2idOUT['<end>']:
        #    return result, loss.numpy() / max_length_targ

        # the predicted ID is fed back into the model
        #dec_input = tf.expand_dims([predicted_id], 0)
        dec_input = tf.expand_dims(predicted_id, 1)
        #print(dec_input.shape)
        #dec_input = tf.expand_dims(output_tensor[:,t], 1)
        #dec_input = tf.repeat(tf.expand_dims([word2idxOUT['<start>']], 0), input_batch_size, axis=0)
   
    return(result, attention_plot, np.array(predictionEachPos), loss.numpy() / max_length_targ) 


# In[ ]:


# All functions defined, now the script is starting


# In[ ]:


# For testing, different commands:

if(False):
    embedding_dim = 6
    units = 500            #This is the hidden dimension, both encoder and decoder
    nSeqTot = 10000
    nRepeats = 10
    condition = 1
    wordSizeIn = 1
    wordSizeOut = 1
    FileInput = "Task4_E_EpiSeq_ParaSeq_NoDeg.txt"
    colInput = 0
    colOutput = 1
    simSuffix = "TestSim"
    #batch_size = 200 => so it's calcualted as 1/50 of the usable data size
    
    refMaxSeqFor20pctsTest = 1000000
    paddingInputs = 'post'    #'pre' or 'post'
    paddingOutputs = 'post'
    EPOCHS = 250


# In[ ]:


# ============== Now opens the file and preprocess the data. ==================


#This step takes long time.
(paddedInputTensor, paddedOutputTensor, word2idxIN, idx2wordIN, word2idxOUT, idx2wordOUT, NI, NO) = processFile(FileInput, wordSizeIn, wordSizeOut, colInput, colOutput, paddingInputs, paddingOutputs)

#Note, we are actually working on np arrays. Tensorflow doesn't like the arrays from python.
type(paddedInputTensor)


# In[ ]:


# ============== Calculates the size of the test and training dataset ==================
nAvailable = len(paddedInputTensor)
nSeqTot = min(nSeqTot,nAvailable)
train_size = int(0.8 * nSeqTot) 
test_size = int(0.2 * refMaxSeqFor20pctsTest)
if test_size + train_size > nAvailable:
   test_size = min(test_size, int(0.2*nAvailable))
print('Requested sequences numbers: dataset=', nSeqTot, 'train=', train_size, 'test=', test_size, 'totUsed=', train_size + test_size)


# In[ ]:


#Calculates batch sizes.
#By default, we will take 50 batches, so the size of a batch is nr sequences / 50
batch_size = math.floor(nSeqTot / 50)
if(runningInCommandLine):
    if(len(sys.argv) > 12):
        batch_size = int(sys.argv[12])


# In[ ]:


#Defines the vocabulary sizes (I have already included pad as a word in the dictionary)
vocab_inp_size = len(word2idxIN)
vocab_tar_size = len(word2idxOUT)
print("Sizes including <pad> word: vocab_inp_size=", vocab_inp_size, "vocab_tar_size=", vocab_tar_size)
    
    
#Sum up of requested task
print('ML Task 4, suffix', simSuffix, 'embedding_dim=', embedding_dim, 'units=', units, 'nSeqTot=', nSeqTot, 'nRepeats=', nRepeats, 'condition=', condition, 'wordSizeIn=', wordSizeIn, 'wordSizeOut', wordSizeOut, 'FileInput', FileInput, 'colInput', colInput, 'colOutput', colOutput, 'batch_size=', batch_size, 'refMaxSeqFor20pctsTest', refMaxSeqFor20pctsTest)


# In[ ]:


#This part can be repeated in the python script file

for repeat in range(nRepeats):

# In[ ]:


	#Creates code and folder: ex: T4_E6_H500_S10.0k_C1_W1_1_CI0CO1_B200_ID1.0
	simCode = "T4_" + simSuffix + "_E" + str(embedding_dim) + "_H" + str(units) + "_S" + str(nSeqTot/1000) + "k_C" + str(condition) + "_W" + str(wordSizeIn) + "_" + str(wordSizeOut) + "_CI" + str(colInput) + "CO" + str(colOutput) + "_B" + str(batch_size) + "_ID" + str(simID) + "." + str(repeat)   
	print(simCode)

	if(not(os.path.exists(simCode))):
		try:
			os.mkdir(simCode)
		except OSError:
			print ("Creation of the directory %s failed" % simCode)
			sys.exit("Could not create folder => abort")
		else:
			print ("Successfully created the directory %s " % simCode)
	else:
		print("WRN: Simulation folder", simCode, " already existed...")


	# In[ ]:


	#Now randomly gets a subsampling 
	possibleDataIDs = [*range(0, nAvailable)]
	random.shuffle(possibleDataIDs)


	# In[ ]:


	#Transform the lists into tensorflow datasets (not yet batched) - this step takes time so I do it only once (2 mins for 100 000 sequences)
	trainKeys = paddedInputTensor[possibleDataIDs[0:train_size]].astype(float);
	testKeys = paddedInputTensor[possibleDataIDs[train_size: train_size + test_size]].astype(float);

	if((condition == 3) or (condition == 13)):
		random.shuffle(possibleDataIDs)

	trainLabels = paddedOutputTensor[possibleDataIDs[0:train_size]];
	testLabels = paddedOutputTensor[possibleDataIDs[train_size: train_size + test_size]];


	
	# In[ ]:


	#Transform into tensors
	TFkeysTrain = tf.data.Dataset.from_tensor_slices(trainKeys)
	TFvalTrain = tf.data.Dataset.from_tensor_slices(trainLabels)
	TFkeysTest = tf.data.Dataset.from_tensor_slices(testKeys)
	TFvalTest = tf.data.Dataset.from_tensor_slices(testLabels)


	# In[ ]:


	#Transforms them into input and output elements (together) 
	tfTrainElements = tf.data.Dataset.zip((TFkeysTrain, TFvalTrain))
	tfTestElements = tf.data.Dataset.zip((TFkeysTest, TFvalTest))

	train_dataset = tfTrainElements.shuffle(train_size).batch(batch_size, drop_remainder=True)
	test_dataset = tfTestElements.shuffle(test_size).batch(batch_size, drop_remainder=True)

	encoder = Encoder(vocab_inp_size, embedding_dim, units, batch_size)
	decoder = Decoder(vocab_tar_size, embedding_dim, units, batch_size)


	# In[ ]:


	optimizer = tf.keras.optimizers.Adam()
	loss_object = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True, reduction='none')


	# In[ ]:


	def loss_function(real, pred):
		mask = tf.math.logical_not(tf.math.equal(real, 0))
		loss_ = loss_object(real, pred)

		mask = tf.cast(mask, dtype=loss_.dtype)
		loss_ *= mask

		return tf.reduce_mean(loss_)


	# In[ ]:


	#This function implicitly needs:
	# - an encoder and decoder
	# - word2idxOUT
	# - batch_size
	# - an optimizer


	# In[ ]:


	@tf.function
	def train_step(inp, targ, enc_hidden, teacher_forcing = True):
		loss = 0

		with tf.GradientTape() as tape:
			enc_output, enc_hidden = encoder(inp, enc_hidden)

			dec_hidden = enc_hidden

			dec_input = tf.expand_dims([word2idxOUT['<start>']] * batch_size, 1)

			# Teacher forcing - feeding the target as the next input
			for t in range(1, targ.shape[1]):
				# passing enc_output to the decoder
				predictions, dec_hidden, _ = decoder(dec_input, dec_hidden, enc_output)

				loss_step = loss_function(targ[:, t], predictions)
				loss += loss_step
				#Note: impossible to get content of tensors if in mode "@tf.function" => remove if want to print
				#print(tf.keras.eval(loss).numpy())
				#print("foo = " + str(loss.numpy() / int(targ.shape[1])) + "Batch _size=" + str(targ.shape[1]))
				#note: this doesn't help either, tf.enable_eager_execution()  - In TensorFlow 2.0 Eager execution is enabled by default. No need to set it up.

				# using teacher forcing
				if(teacher_forcing):
					dec_input = tf.expand_dims(targ[:, t], 1)

				#This would be the non-teacher-forcing:
				else:
					predicted_id = tf.argmax(predictions, axis=1)
					dec_input = tf.expand_dims(predicted_id, 1)

		batch_loss = (loss / int(targ.shape[1]))

		variables = encoder.trainable_variables + decoder.trainable_variables

		gradients = tape.gradient(loss, variables)

		optimizer.apply_gradients(zip(gradients, variables))

		return batch_loss


	# In[ ]:


	#EPOCHS = 200
	steps_per_epoch = len(trainKeys)//batch_size
	print("The nr of steps per epoch will be", steps_per_epoch)


	# In[ ]:


	for epoch in range(EPOCHS):
		start = time.time()

		enc_hidden = encoder.initialize_hidden_state()
		total_loss = 0

		for (batch, (inp, targ)) in enumerate(train_dataset.take(steps_per_epoch)):
			batch_loss = train_step(inp, targ, enc_hidden)
			total_loss += batch_loss

			#if batch % 100 == 0:
			#  print('Epoch {} Batch {} Loss {:.4f}'.format(epoch + 1, batch, batch_loss.numpy()))

			# saving (checkpoint) the model every 2 epochs
			#if (epoch + 1) % 2 == 0:
			#  checkpoint.save(file_prefix = checkpoint_prefix)

		print('Forced Epoch {} Loss {:.4f}'.format(epoch + 1, total_loss / steps_per_epoch))
		print('Time taken for 1 epoch {} sec\n'.format(time.time() - start))


	# In[ ]:


	encoder.save_weights(simCode + "/" + 'encoder')
	decoder.save_weights(simCode + "/" + 'decoder')


	# In[ ]:


	for epoch in range(EPOCHS//10):
		start = time.time()

		enc_hidden = encoder.initialize_hidden_state()
		total_loss = 0

		for (batch, (inp, targ)) in enumerate(train_dataset.take(steps_per_epoch)):
			batch_loss = train_step(inp, targ, enc_hidden, False)
			total_loss += batch_loss

		print('Unforced Epoch {} Loss {:.4f}'.format(epoch + 1, total_loss / steps_per_epoch))
		print('Time taken for 1 epoch {} sec\n'.format(time.time() - start))


	# In[ ]:


	encoder.save_weights(simCode + "/" + 'encoderUnforced')
	decoder.save_weights(simCode + "/" + 'decoderUnforced')


	# In[ ]:


	#Now displaying the predictions for all sequences, train or test

	list_batched_inputs = np.array_split(trainKeys, 100, 0)
	list_batched_outputs =  np.array_split(trainLabels, 100, 0)
	
	for index in range(0,100):
		batched_inputs = list_batched_inputs[index]
		batched_outputs =  list_batched_outputs[index]
		(result, attention_plot, arrayPred, loss) = evaluateBatch(batched_inputs, batched_outputs, encoder, decoder, NI, NO,  word2idxIN, word2idxOUT)
		batched_predictions = result.numpy()

		textIn = vocaliseBatch(batched_inputs, idx2wordIN, wordSizeIn)
		textPred = vocaliseBatch(result.numpy(),  idx2wordOUT, wordSizeOut)
		textExpected = vocaliseBatch(batched_outputs,  idx2wordOUT, wordSizeOut)


		# In[ ]:


		AllResults = []
		f = open(simCode + "/" + "TrainResults.txt", "a")
		if(index == 0):
			f.write("lenout\tlenpred\tlenin\tld\tldnorm\tpred\tout\tin\n")

		# In[ ]:


		for i in range(0,min(len(batched_predictions),100000)):
			ld = jellyfish.levenshtein_distance(textPred[i],textExpected[i])
			ldnorm = ld/len(textPred[i])

			lenpara = len(textExpected[i])
			lenpredpara = len(textPred[i])
			lenepi = len(textIn[i])

			#AllResults.append([lenpredpara, lenpara, lenepi, ld, ldnorm, textPred[i], textExpected[i], textIn[i]])
			f.write(str(lenpredpara) + "\t" + str(lenpara) + "\t" + str(lenepi) + "\t" + str(ld) + "\t" + str(ldnorm) + "\t" + str(textPred[i]) + "\t" + str(textExpected[i]) + "\t" + str(textIn[i]) + "\n")

		f.close()


	list_batched_inputs = np.array_split(testKeys, 100, 0)
	list_batched_outputs =  np.array_split(testLabels, 100, 0)
	
	if(len(testKeys) > 100000):
		list_batched_inputs = np.array_split(testKeys[0:100000], 100, 0)
		list_batched_outputs =  np.array_split(testLabels[0:100000], 100, 0)
	
	for index in range(0,100):
		batched_inputs = list_batched_inputs[index]
		batched_outputs =  list_batched_outputs[index]

		(result, attention_plot, arrayPred, loss) = evaluateBatch(batched_inputs, batched_outputs, encoder, decoder, NI, NO,  word2idxIN, word2idxOUT)
		batched_predictions = result.numpy()


		# In[ ]:


		textIn = vocaliseBatch(batched_inputs, idx2wordIN, wordSizeIn)
		textPred = vocaliseBatch(result.numpy(),  idx2wordOUT, wordSizeOut)
		textExpected = vocaliseBatch(batched_outputs,  idx2wordOUT, wordSizeOut)


		# In[ ]:


		AllResults = []
		f = open(simCode + "/" + "TestResults.txt", "a")
		if(index == 0):
			f.write("lenout\tlenpred\tlenin\tld\tldnorm\tpred\tout\tin\n")


		# In[ ]:


		for i in range(0,min(len(batched_predictions),100000)):
			ld = jellyfish.levenshtein_distance(textPred[i],textExpected[i])
			ldnorm = ld/len(textPred[i])

			lenpara = len(textExpected[i])
			lenpredpara = len(textPred[i])
			lenepi = len(textIn[i])

			#AllResults.append([lenpredpara, lenpara, lenepi, ld, ldnorm, textPred[i], textExpected[i], textIn[i]])
			f.write(str(lenpredpara) + "\t" + str(lenpara) + "\t" + str(lenepi) + "\t" + str(ld) + "\t" + str(ldnorm) + "\t" + str(textPred[i]) + "\t" + str(textExpected[i]) + "\t" + str(textIn[i]) + "\n")

		f.close()


		# In[ ]:





	# In[ ]:




