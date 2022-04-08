#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python3
# coding: utf-8

from __future__ import division, print_function, absolute_import


# In[2]:


# Usage: 
# ThisScript.py d_model=16  dff=512  nSeqTot=1000  nRepeats=10  condition=1  wordSizeIn=1  wordSizeOut=1  
#             fileInput=seeBelow  colInput=0   colOutput=1  simSuffix=""   batch_size=nSeqTot/50  
# 
# Requires: 
# an input file with two columns: input and output. 
# the location is defined just below:
#
# Outputs:
#    creates a folder unique to this simulation
#    puts the weights of the fitted model inside
#
# Parameters of the transformer
# d_model = 16        # This is actually the embedding from vocab into this dimension.
# NEED d_model % self.num_heads == 0
# dff = 512           # Number of neurons in the dense layers


# In[3]:


# local settings
runningInCommandLine = True


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

import time

import sys
import math


# In[5]:


if(not(runningInCommandLine)):
    import matplotlib.pyplot as plt

os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

#to get an ID of simulation
from datetime import datetime
from random import seed
import random
import jellyfish


# In[6]:


#Default parameters (changeable by command line)
#in the same order of command line arguments

d_model = 12
dff = 512            
nSeqTot = 10000
nRepeats = 10
condition = 1
wordSizeIn = 1
wordSizeOut = 1
FileInput = "Task4_D_EpiChem_ParaChemD2.txt"
FileExternal = ""
colInput = 0
colOutput = 1
simSuffix = ""
batch_size = 200


# In[7]:


# Non-command line options
refMaxSeqFor20pctsTest = 4000 #1000000
paddingInputs = 'post'
paddingOutputs = 'post'
EPOCHS = 250

# Transformer hyperparameters
num_layers = 4      # Paper says 6, good range => we keep 4
num_heads = 8       # Number of parallel heads inside an encoder layer, I guess (multi-head attention).
dropout_rate = 0.1

removeStar = True   #Paratopes are distinguished by ending by a star. It means we will remove ending star from input and ouptuts sequences if they have it 
maxOutput = 100000  #For printint predictions


# In[8]:


#get a simID to save the model weights
random.seed(datetime.now())
simID = random.randint(1,9999) #both numbers included
print("This run has ID",str(simID))
#Note: the repeats will be called ID.1, etc...


# In[9]:


if(runningInCommandLine):
    if(len(sys.argv) > 1):
        d_model = int(sys.argv[1])

    if(len(sys.argv) > 2):
        dff = int(sys.argv[2])


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
        
    if(len(sys.argv) > 12):
        FileExternal = sys.argv[12]
        print("External test file assign to be ", FileExternal)


    # argument 13 is the batch_size, see below


# In[10]:


def RemoveEndingStar(text):
    if(len(text) > 0):
        if(text[-1] == '*'):
            text = text[:-1]
    return(text)


# In[11]:


def getFeaturesLabels(fileName, colI=0, colO=1):
    epiParaSeq = open(fileName, newline = '')   #one line is a text with \t and \n                                                                              
    data_reader = csv.reader(epiParaSeq, delimiter='\t') #transform lines into lists 
    sequences = []
    labels = []
    for line in data_reader:
        #print(line)
        if(not (line[0].startswith("#") or line[0].startswith("seqAGEpitope") or line[0].startswith("EpitopeEnc"))):
            if(removeStar):
                sequences.append(RemoveEndingStar(line[colI]))
                labels.append(RemoveEndingStar(line[colO]))
            else:
                sequences.append(line[colI])
                labels.append(line[colO])
    
    print("Got ", len(sequences), " sequences and ", len(labels), " from ", fileName)
    return (sequences, labels)


# In[12]:


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


# In[13]:


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


# In[14]:


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


# In[15]:


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


# In[16]:


def padNewSentences(sentences, maxLength, word2IdDictionary, paddingOption):
    input_tensor = [[word2IdDictionary[s] for s in sentence.split(' ')] for sentence in sentences]
    paddedSequences = tf.keras.preprocessing.sequence.pad_sequences(input_tensor,maxlen=maxLength,padding=paddingOption)
    return paddedSequences


# In[ ]:





# In[17]:


def get_angles(pos, i, d_model):
    angle_rates = 1 / np.power(10000, (2 * (i//2)) / np.float32(d_model))
    return pos * angle_rates


# In[18]:


def positional_encoding(position, d_model):
    angle_rads = get_angles(np.arange(position)[:, np.newaxis],
                          np.arange(d_model)[np.newaxis, :],
                          d_model)

    # apply sin to even indices in the array; 2i
    angle_rads[:, 0::2] = np.sin(angle_rads[:, 0::2])

    # apply cos to odd indices in the array; 2i+1
    angle_rads[:, 1::2] = np.cos(angle_rads[:, 1::2])

    pos_encoding = angle_rads[np.newaxis, ...]

    return tf.cast(pos_encoding, dtype=tf.float32)


# In[19]:


if(not(runningInCommandLine)):
    pos_encoding = positional_encoding(50, 512)
    print (pos_encoding.shape)

    plt.pcolormesh(pos_encoding[0], cmap='RdBu')
    plt.xlabel('Depth')
    plt.xlim((0, 512))
    plt.ylabel('Position')
    plt.colorbar()
    plt.show()


# In[20]:


def create_padding_mask(seq):
    seq = tf.cast(tf.math.equal(seq, 0), tf.float32)

    # add extra dimensions to add the padding
    # to the attention logits.
    return seq[:, tf.newaxis, tf.newaxis, :]  # (batch_size, 1, 1, seq_len)


# In[21]:


if(not(runningInCommandLine)):
    x = tf.constant([[7, 6, 0, 0, 1], [1, 2, 3, 0, 0], [0, 0, 0, 4, 5]])
    create_padding_mask(x)


# In[22]:


def create_look_ahead_mask(size):
    mask = 1 - tf.linalg.band_part(tf.ones((size, size)), -1, 0)
    return mask  # (seq_len, seq_len)


# In[23]:


if(not(runningInCommandLine)):
    x = tf.random.uniform((1, 3))
    temp = create_look_ahead_mask(x.shape[1])
    temp


# In[24]:


def scaled_dot_product_attention(q, k, v, mask):
    """Calculate the attention weights.
    q, k, v must have matching leading dimensions.
    k, v must have matching penultimate dimension, i.e.: seq_len_k = seq_len_v.
    The mask has different shapes depending on its type(padding or look ahead) 
    but it must be broadcastable for addition.

    Args:
    q: query shape == (..., seq_len_q, depth)
    k: key shape == (..., seq_len_k, depth)
    v: value shape == (..., seq_len_v, depth_v)
    mask: Float tensor with shape broadcastable 
          to (..., seq_len_q, seq_len_k). Defaults to None.

    Returns:
    output, attention_weights
    """

    matmul_qk = tf.matmul(q, k, transpose_b=True)  # (..., seq_len_q, seq_len_k)

    # scale matmul_qk
    dk = tf.cast(tf.shape(k)[-1], tf.float32)
    scaled_attention_logits = matmul_qk / tf.math.sqrt(dk)

    # add the mask to the scaled tensor.
    if mask is not None:
        scaled_attention_logits += (mask * -1e9)  

    # softmax is normalized on the last axis (seq_len_k) so that the scores
    # add up to 1.
    attention_weights = tf.nn.softmax(scaled_attention_logits, axis=-1)  # (..., seq_len_q, seq_len_k)

    output = tf.matmul(attention_weights, v)  # (..., seq_len_q, depth_v)

    return output, attention_weights


# In[25]:


def print_out(q, k, v):
    temp_out, temp_attn = scaled_dot_product_attention(
      q, k, v, None)
    print ('Attention weights are:')
    print (temp_attn)
    print ('Output is:')
    print (temp_out)


# In[26]:


if(not(runningInCommandLine)):
    np.set_printoptions(suppress=True)

    temp_k = tf.constant([[10,0,0],
                          [0,10,0],
                          [0,0,10],
                          [0,0,10]], dtype=tf.float32)  # (4, 3)

    temp_v = tf.constant([[   1,0],
                          [  10,0],
                          [ 100,5],
                          [1000,6]], dtype=tf.float32)  # (4, 2)

    # This `query` aligns with the second `key`,
    # so the second `value` is returned.
    temp_q = tf.constant([[0, 10, 0]], dtype=tf.float32)  # (1, 3)
    print_out(temp_q, temp_k, temp_v)


# In[27]:


if(not(runningInCommandLine)):
    # This query aligns with a repeated key (third and fourth), 
    # so all associated values get averaged.
    temp_q = tf.constant([[0, 0, 10]], dtype=tf.float32)  # (1, 3)
    print_out(temp_q, temp_k, temp_v)


# In[28]:


if(not(runningInCommandLine)):
    # This query aligns equally with the first and second key, 
    # so their values get averaged.
    temp_q = tf.constant([[10, 10, 0]], dtype=tf.float32)  # (1, 3)
    print_out(temp_q, temp_k, temp_v)


# In[29]:


if(not(runningInCommandLine)):
    temp_q = tf.constant([[0, 0, 10], [0, 10, 0], [10, 10, 0]], dtype=tf.float32)  # (3, 3)
    print_out(temp_q, temp_k, temp_v)


# In[30]:


class MultiHeadAttention(tf.keras.layers.Layer):
    def __init__(self, d_model, num_heads):
        super(MultiHeadAttention, self).__init__()
        self.num_heads = num_heads
        self.d_model = d_model

        assert d_model % self.num_heads == 0

        self.depth = d_model // self.num_heads

        self.wq = tf.keras.layers.Dense(d_model)
        self.wk = tf.keras.layers.Dense(d_model)
        self.wv = tf.keras.layers.Dense(d_model)

        self.dense = tf.keras.layers.Dense(d_model)

    def split_heads(self, x, batch_size):
        """Split the last dimension into (num_heads, depth).
        Transpose the result such that the shape is (batch_size, num_heads, seq_len, depth)
        """
        x = tf.reshape(x, (batch_size, -1, self.num_heads, self.depth))
        return tf.transpose(x, perm=[0, 2, 1, 3])

    def call(self, v, k, q, mask):
        batch_size = tf.shape(q)[0]

        q = self.wq(q)  # (batch_size, seq_len, d_model)
        k = self.wk(k)  # (batch_size, seq_len, d_model)
        v = self.wv(v)  # (batch_size, seq_len, d_model)

        q = self.split_heads(q, batch_size)  # (batch_size, num_heads, seq_len_q, depth)
        k = self.split_heads(k, batch_size)  # (batch_size, num_heads, seq_len_k, depth)
        v = self.split_heads(v, batch_size)  # (batch_size, num_heads, seq_len_v, depth)

        # scaled_attention.shape == (batch_size, num_heads, seq_len_q, depth)
        # attention_weights.shape == (batch_size, num_heads, seq_len_q, seq_len_k)
        scaled_attention, attention_weights = scaled_dot_product_attention(q, k, v, mask)

        scaled_attention = tf.transpose(scaled_attention, perm=[0, 2, 1, 3])  # (batch_size, seq_len_q, num_heads, depth)

        concat_attention = tf.reshape(scaled_attention, (batch_size, -1, self.d_model))  # (batch_size, seq_len_q, d_model)

        output = self.dense(concat_attention)  # (batch_size, seq_len_q, d_model)

        return output, attention_weights


# In[31]:


if(not(runningInCommandLine)):
    temp_mha = MultiHeadAttention(d_model=512, num_heads=8)
    y = tf.random.uniform((1, 60, 512))  # (batch_size, encoder_sequence, d_model)
    out, attn = temp_mha(y, k=y, q=y, mask=None)
    out.shape, attn.shape


# In[32]:


def point_wise_feed_forward_network(d_model, dff):
  return tf.keras.Sequential([
      tf.keras.layers.Dense(dff, activation='relu'),  # (batch_size, seq_len, dff)
      tf.keras.layers.Dense(d_model)  # (batch_size, seq_len, d_model)
  ])


# In[33]:


if(not(runningInCommandLine)):
    sample_ffn = point_wise_feed_forward_network(512, 2048)
    sample_ffn(tf.random.uniform((64, 50, 512))).shape


# In[34]:


class EncoderLayer(tf.keras.layers.Layer):
    def __init__(self, d_model, num_heads, dff, rate=0.1):
        super(EncoderLayer, self).__init__()

        self.mha = MultiHeadAttention(d_model, num_heads)
        self.ffn = point_wise_feed_forward_network(d_model, dff)

        self.layernorm1 = tf.keras.layers.LayerNormalization(epsilon=1e-6)
        self.layernorm2 = tf.keras.layers.LayerNormalization(epsilon=1e-6)

        self.dropout1 = tf.keras.layers.Dropout(rate)
        self.dropout2 = tf.keras.layers.Dropout(rate)

    def call(self, x, training, mask):

        attn_output, _ = self.mha(x, x, x, mask)  # (batch_size, input_seq_len, d_model)
        attn_output = self.dropout1(attn_output, training=training)
        out1 = self.layernorm1(x + attn_output)  # (batch_size, input_seq_len, d_model)

        ffn_output = self.ffn(out1)  # (batch_size, input_seq_len, d_model)
        ffn_output = self.dropout2(ffn_output, training=training)
        out2 = self.layernorm2(out1 + ffn_output)  # (batch_size, input_seq_len, d_model)

        return out2


# In[35]:


if(not(runningInCommandLine)):
    sample_encoder_layer = EncoderLayer(512, 8, 2048)

    sample_encoder_layer_output = sample_encoder_layer(
        tf.random.uniform((64, 43, 512)), False, None)

    sample_encoder_layer_output.shape  # (batch_size, input_seq_len, d_model)


# In[36]:


class DecoderLayer(tf.keras.layers.Layer):
    def __init__(self, d_model, num_heads, dff, rate=0.1):
        super(DecoderLayer, self).__init__()

        self.mha1 = MultiHeadAttention(d_model, num_heads)
        self.mha2 = MultiHeadAttention(d_model, num_heads)

        self.ffn = point_wise_feed_forward_network(d_model, dff)

        self.layernorm1 = tf.keras.layers.LayerNormalization(epsilon=1e-6)
        self.layernorm2 = tf.keras.layers.LayerNormalization(epsilon=1e-6)
        self.layernorm3 = tf.keras.layers.LayerNormalization(epsilon=1e-6)

        self.dropout1 = tf.keras.layers.Dropout(rate)
        self.dropout2 = tf.keras.layers.Dropout(rate)
        self.dropout3 = tf.keras.layers.Dropout(rate)


    def call(self, x, enc_output, training, 
               look_ahead_mask, padding_mask):
        # enc_output.shape == (batch_size, input_seq_len, d_model)

        attn1, attn_weights_block1 = self.mha1(x, x, x, look_ahead_mask)  # (batch_size, target_seq_len, d_model)
        attn1 = self.dropout1(attn1, training=training)
        out1 = self.layernorm1(attn1 + x)

        attn2, attn_weights_block2 = self.mha2(
            enc_output, enc_output, out1, padding_mask)  # (batch_size, target_seq_len, d_model)
        attn2 = self.dropout2(attn2, training=training)
        out2 = self.layernorm2(attn2 + out1)  # (batch_size, target_seq_len, d_model)

        ffn_output = self.ffn(out2)  # (batch_size, target_seq_len, d_model)
        ffn_output = self.dropout3(ffn_output, training=training)
        out3 = self.layernorm3(ffn_output + out2)  # (batch_size, target_seq_len, d_model)

        return out3, attn_weights_block1, attn_weights_block2


# In[37]:


#sample_decoder_layer = DecoderLayer(512, 8, 2048)

#sample_decoder_layer_output, _, _ = sample_decoder_layer(
#    tf.random.uniform((64, 50, 512)), sample_encoder_layer_output, 
#    False, None, None)

#sample_decoder_layer_output.shape  # (batch_size, target_seq_len, d_model)


# In[38]:


class Encoder(tf.keras.layers.Layer):
    def __init__(self, num_layers, d_model, num_heads, dff, input_vocab_size,
                   maximum_position_encoding, rate=0.1):
        super(Encoder, self).__init__()

        self.d_model = d_model
        self.num_layers = num_layers

        self.embedding = tf.keras.layers.Embedding(input_vocab_size, d_model)
        self.pos_encoding = positional_encoding(maximum_position_encoding, 
                                                self.d_model)


        self.enc_layers = [EncoderLayer(d_model, num_heads, dff, rate) 
                           for _ in range(num_layers)]

        self.dropout = tf.keras.layers.Dropout(rate)

    def call(self, x, training, mask):

        seq_len = tf.shape(x)[1]

        # adding embedding and position encoding.
        x = self.embedding(x)  # (batch_size, input_seq_len, d_model)
        x *= tf.math.sqrt(tf.cast(self.d_model, tf.float32))
        x += self.pos_encoding[:, :seq_len, :]

        x = self.dropout(x, training=training)

        for i in range(self.num_layers):
            x = self.enc_layers[i](x, training, mask)

        return x  # (batch_size, input_seq_len, d_model)


# In[39]:


if(not(runningInCommandLine)):
    sample_encoder = Encoder(num_layers=2, d_model=512, num_heads=8, 
                             dff=2048, input_vocab_size=8500,
                             maximum_position_encoding=10000)
    temp_input = tf.random.uniform((64, 62), dtype=tf.int64, minval=0, maxval=200)

    sample_encoder_output = sample_encoder(temp_input, training=False, mask=None)

    print (sample_encoder_output.shape)  # (batch_size, input_seq_len, d_model)


# In[40]:


class Decoder(tf.keras.layers.Layer):
    def __init__(self, num_layers, d_model, num_heads, dff, target_vocab_size,
                   maximum_position_encoding, rate=0.1):
        super(Decoder, self).__init__()

        self.d_model = d_model
        self.num_layers = num_layers

        self.embedding = tf.keras.layers.Embedding(target_vocab_size, d_model)
        self.pos_encoding = positional_encoding(maximum_position_encoding, d_model)

        self.dec_layers = [DecoderLayer(d_model, num_heads, dff, rate) 
                           for _ in range(num_layers)]
        self.dropout = tf.keras.layers.Dropout(rate)

    def call(self, x, enc_output, training, look_ahead_mask, padding_mask):

        seq_len = tf.shape(x)[1]
        attention_weights = {}

        x = self.embedding(x)  # (batch_size, target_seq_len, d_model)
        x *= tf.math.sqrt(tf.cast(self.d_model, tf.float32))
        x += self.pos_encoding[:, :seq_len, :]

        x = self.dropout(x, training=training)

        for i in range(self.num_layers):
            x, block1, block2 = self.dec_layers[i](x, enc_output, training, look_ahead_mask, padding_mask)

            attention_weights['decoder_layer{}_block1'.format(i+1)] = block1
            attention_weights['decoder_layer{}_block2'.format(i+1)] = block2

        # x.shape == (batch_size, target_seq_len, d_model)
        return x, attention_weights


# In[41]:


if(not(runningInCommandLine)):
    sample_decoder = Decoder(num_layers=2, d_model=512, num_heads=8, 
                             dff=2048, target_vocab_size=8000,
                             maximum_position_encoding=5000)
    temp_input = tf.random.uniform((64, 26), dtype=tf.int64, minval=0, maxval=200)

    output, attn = sample_decoder(temp_input, 
                                  enc_output=sample_encoder_output, 
                                  training=False,
                                  look_ahead_mask=None, 
                                  padding_mask=None)

    output.shape, attn['decoder_layer2_block2'].shape


# In[42]:


class Transformer(tf.keras.Model):
    def __init__(self, num_layers, d_model, num_heads, dff, input_vocab_size, 
                   target_vocab_size, pe_input, pe_target, rate=0.1):
        super(Transformer, self).__init__()

        self.encoder = Encoder(num_layers, d_model, num_heads, dff, 
                               input_vocab_size, pe_input, rate)

        self.decoder = Decoder(num_layers, d_model, num_heads, dff, 
                               target_vocab_size, pe_target, rate)

        self.final_layer = tf.keras.layers.Dense(target_vocab_size)

    def call(self, inp, tar, training, enc_padding_mask, 
               look_ahead_mask, dec_padding_mask):

        enc_output = self.encoder(inp, training, enc_padding_mask)  # (batch_size, inp_seq_len, d_model)

        # dec_output.shape == (batch_size, tar_seq_len, d_model)
        dec_output, attention_weights = self.decoder(
            tar, enc_output, training, look_ahead_mask, dec_padding_mask)

        final_output = self.final_layer(dec_output)  # (batch_size, tar_seq_len, target_vocab_size)

        return final_output, attention_weights


# In[43]:


if(not(runningInCommandLine)):
    sample_transformer = Transformer(
        num_layers=2, d_model=512, num_heads=8, dff=2048, 
        input_vocab_size=8500, target_vocab_size=8000, 
        pe_input=10000, pe_target=6000)

    temp_input = tf.random.uniform((64, 38), dtype=tf.int64, minval=0, maxval=200)
    temp_target = tf.random.uniform((64, 36), dtype=tf.int64, minval=0, maxval=200)

    fn_out, _ = sample_transformer(temp_input, temp_target, training=False, 
                                   enc_padding_mask=None, 
                                   look_ahead_mask=None,
                                   dec_padding_mask=None)

    fn_out.shape  # (batch_size, tar_seq_len, target_vocab_size)


# In[44]:


#Functions to transforme back a tokenized sentence into words (data processing) 


# In[45]:


def vocalise(tensor, dico, wordsize = 1):
    pulled = ''.join(dico[i] for i in tensor)
    if(wordsize < 0):
            pulled = ' '.join(dico[i] for i in tensor)
    return pulled.replace("<pad>", "").replace("<end>", "").replace("<start>", "")


# In[46]:


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





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[47]:


# All functions defined, now the script is starting


# In[48]:


# For testing, different commands:

if(not(runningInCommandLine)):
    if(True):
        d_model = 16
        dff = 510            #This is the hidden dimension, both encoder and decoder
        nSeqTot = 2000
        nRepeats = 10
        condition = 1
        wordSizeIn = 2
        wordSizeOut = 2
        FileInput = "Task4_D_EpiChem_ParaChemD2.txt"
        colInput = 0
        colOutput = 1
        simSuffix = "TestSim"
        #batch_size = 200 => so it's calcualted as 1/50 of the usable data size

        refMaxSeqFor20pctsTest = 10000
        paddingInputs = 'post'    #'pre' or 'post'
        paddingOutputs = 'post'
        EPOCHS = 25
        num_layers = 4      # Paper says 6, good range => we keep 4
        num_heads = 8


# In[49]:


# ============== Now opens the file and preprocess the data. ==================


#This step takes long time.
(paddedInputTensor, paddedOutputTensor, word2idxIN, idx2wordIN, word2idxOUT, idx2wordOUT, NI, NO) = processFile(FileInput, wordSizeIn, wordSizeOut, colInput, colOutput, paddingInputs, paddingOutputs)

#Note, we are actually working on np arrays. Tensorflow doesn't like the arrays from python.
type(paddedInputTensor)


# In[50]:


# ============== Calculates the size of the test and training dataset ==================
nAvailable = len(paddedInputTensor)
nSeqTot = min(nSeqTot,nAvailable)
train_size = int(0.8 * nSeqTot) 
test_size = int(0.2 * refMaxSeqFor20pctsTest)
if test_size + train_size > nAvailable:
   test_size = min(test_size, int(0.2*nAvailable))
print('Requested sequences numbers: dataset=', nSeqTot, 'train=', train_size, 'test=', test_size, 'totUsed=', train_size + test_size)


# In[51]:


#Calculates batch sizes.
#By default, we will take 50 batches, so the size of a batch is nr sequences / 50
batch_size = math.floor(nSeqTot / 50)
if(runningInCommandLine):
    if(len(sys.argv) > 13):
        batch_size = int(sys.argv[13])


# In[ ]:





# In[52]:


#Sum up of requested task
print('ML Task 4, suffix', simSuffix, 'd_model=', d_model, 'dff=', dff, 'nSeqTot=', nSeqTot, 'nRepeats=', nRepeats, 'condition=', condition, 'wordSizeIn=', wordSizeIn, 'wordSizeOut', wordSizeOut, 'FileInput', FileInput, 'colInput', colInput, 'colOutput', colOutput, 'batch_size=', batch_size, 'refMaxSeqFor20pctsTest', refMaxSeqFor20pctsTest)


# In[53]:


#This part can be repeated in the python script file

for repeat in range(nRepeats):
	#repeat = 1


	# In[54]:


	#Creates code and folder: ex: T4_E6_H500_S10.0k_C1_W1_1_CI0CO1_B200_ID1.0
	simCode = "TRA_" + simSuffix + "_E" + str(d_model) + "_Dff" + str(dff) + "_S" + str(nSeqTot/1000) + "k_C" + str(condition) + "_W" + str(wordSizeIn) + "_" + str(wordSizeOut) + "_CI" + str(colInput) + "CO" + str(colOutput) + "_B" + str(batch_size) + "_ID" + str(simID) + "." + str(repeat)   
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


	# In[55]:


	#Now randomly gets a subsampling 
	possibleDataIDs = [*range(0, nAvailable)]
	random.shuffle(possibleDataIDs)


	# In[56]:




	# In[57]:


	#Transform the lists into tensorflow datasets (not yet batched) - this step takes time so I do it only once (2 mins for 100 000 sequences)
	trainKeys = paddedInputTensor[possibleDataIDs[0:train_size]].astype(int);
	testKeys = paddedInputTensor[possibleDataIDs[train_size: train_size + test_size]].astype(int);

	if((condition == 3) or (condition == 13)):
		random.shuffle(possibleDataIDs)

	trainLabels = paddedOutputTensor[possibleDataIDs[0:train_size]].astype(int);
	testLabels = paddedOutputTensor[possibleDataIDs[train_size: train_size + test_size]].astype(int);


	# In[58]:


	#Transform into tensors
	TFkeysTrain = tf.data.Dataset.from_tensor_slices(tf.cast(trainKeys, tf.int64))
	TFvalTrain = tf.data.Dataset.from_tensor_slices(tf.cast(trainLabels, tf.int64))
	TFkeysTest = tf.data.Dataset.from_tensor_slices(tf.cast(testKeys, tf.int64))
	TFvalTest = tf.data.Dataset.from_tensor_slices(tf.cast(testLabels, tf.int64))


	# In[59]:


	#Transforms them into input and output elements (together) 
	tfTrainElements = tf.data.Dataset.zip((TFkeysTrain, TFvalTrain))
	tfTestElements = tf.data.Dataset.zip((TFkeysTest, TFvalTest))

	train_dataset = tfTrainElements.shuffle(train_size).batch(batch_size, drop_remainder=True)
	test_dataset = tfTestElements.shuffle(test_size).batch(batch_size, drop_remainder=True)


	# In[60]:


	#print(train_dataset)


	# In[61]:


	#print(train_dataset)


	# In[62]:


	#=============== Before was data prep, now is the architecture part =======


	# In[63]:


	# Each encoder layer has: EncoderLayer(d_model, num_heads, dff, rate)

	#num_layers = 4      # Paper says 6, good range => we keep 4
	#d_model = 128
	#dff = 512           # Number of neurons in the dense layers
	#num_heads = 8       # Number of parallel heads inside an encoder layer, I guess (multi-head attention).
	#dropout_rate = 0.1

	input_vocab_size = len(word2idxIN) 
	target_vocab_size = len(word2idxOUT) 


	# In[64]:


	class CustomSchedule(tf.keras.optimizers.schedules.LearningRateSchedule):
		def __init__(self, d_model, warmup_steps=4000):
			super(CustomSchedule, self).__init__()

			self.d_model = d_model
			self.d_model = tf.cast(self.d_model, tf.float32)

			self.warmup_steps = warmup_steps

		def __call__(self, step):
			arg1 = tf.math.rsqrt(step)
			arg2 = step * (self.warmup_steps ** -1.5)

			return tf.math.rsqrt(self.d_model) * tf.math.minimum(arg1, arg2)


	# In[65]:


	learning_rate = CustomSchedule(d_model)

	optimizer = tf.keras.optimizers.Adam(learning_rate, beta_1=0.9, beta_2=0.98, 
										 epsilon=1e-9)


	# In[66]:


	if(not(runningInCommandLine)):
		temp_learning_rate_schedule = CustomSchedule(d_model)

		plt.plot(temp_learning_rate_schedule(tf.range(40000, dtype=tf.float32)))
		plt.ylabel("Learning Rate")
		plt.xlabel("Train Step")


	# In[67]:


	loss_object = tf.keras.losses.SparseCategoricalCrossentropy(
		from_logits=True, reduction='none')


	# In[68]:


	def loss_function(real, pred):
		mask = tf.math.logical_not(tf.math.equal(real, 0))
		loss_ = loss_object(real, pred)

		mask = tf.cast(mask, dtype=loss_.dtype)
		loss_ *= mask

		return tf.reduce_sum(loss_)/tf.reduce_sum(mask)


	def accuracy_function(real, pred):
		accuracies = tf.equal(real, tf.argmax(pred, axis=2))

		mask = tf.math.logical_not(tf.math.equal(real, 0))
		accuracies = tf.math.logical_and(mask, accuracies)

		accuracies = tf.cast(accuracies, dtype=tf.float32)
		mask = tf.cast(mask, dtype=tf.float32)
		return tf.reduce_sum(accuracies)/tf.reduce_sum(mask)


	# In[69]:


	train_loss = tf.keras.metrics.Mean(name='train_loss')
	train_accuracy = tf.keras.metrics.Mean(name='train_accuracy')


	# In[ ]:





	# In[70]:


	transformer = Transformer(num_layers, d_model, num_heads, dff,
							  input_vocab_size, target_vocab_size, 
							  pe_input=trainKeys.shape[1], 
							  pe_target=trainLabels.shape[1],
							  rate=dropout_rate)


	# In[71]:


	def create_masks(inp, tar):
		# Encoder padding mask
		enc_padding_mask = create_padding_mask(inp)

		# Used in the 2nd attention block in the decoder.
		# This padding mask is used to mask the encoder outputs.
		dec_padding_mask = create_padding_mask(inp)

		# Used in the 1st attention block in the decoder.
		# It is used to pad and mask future tokens in the input received by 
		# the decoder.
		look_ahead_mask = create_look_ahead_mask(tf.shape(tar)[1])
		dec_target_padding_mask = create_padding_mask(tar)
		combined_mask = tf.maximum(dec_target_padding_mask, look_ahead_mask)

		return enc_padding_mask, combined_mask, dec_padding_mask


	# In[72]:


	#EPOCHS = 20


	# In[73]:


	# The @tf.function trace-compiles train_step into a TF graph for faster
	# execution. The function specializes to the precise shape of the argument
	# tensors. To avoid re-tracing due to the variable sequence lengths or variable
	# batch sizes (the last batch is smaller), use input_signature to specify
	# more generic shapes.

	train_step_signature = [
		tf.TensorSpec(shape=(None, None), dtype=tf.int64),
		tf.TensorSpec(shape=(None, None), dtype=tf.int64),
	]

	@tf.function(input_signature=train_step_signature)
	def train_step(inp, tar):
		tar_inp = tar[:, :-1]
		tar_real = tar[:, 1:]

		enc_padding_mask, combined_mask, dec_padding_mask = create_masks(inp, tar_inp)

		with tf.GradientTape() as tape:
			predictions, _ = transformer(inp, tar_inp, 
										 True, 
										 enc_padding_mask, 
										 combined_mask, 
										 dec_padding_mask)
			loss = loss_function(tar_real, predictions)

		gradients = tape.gradient(loss, transformer.trainable_variables)    
		optimizer.apply_gradients(zip(gradients, transformer.trainable_variables))

		train_loss(loss)
		train_accuracy(accuracy_function(tar_real, predictions))


	# In[74]:


	#================ Fitting is here! =====================


	# In[75]:


	for epoch in range(EPOCHS):
		start = time.time()

		train_loss.reset_states()
		train_accuracy.reset_states()

		# inp -> portuguese, tar -> english
		for (batch, (inp, tar)) in enumerate(train_dataset):
			train_step(inp, tar)

		if(batch % 50 == 0):
			print ('Epoch {} Batch {} Loss {:.4f} Accuracy {:.4f}'.format(epoch + 1, batch, train_loss.result(), train_accuracy.result()))

		#if (epoch + 1) % 5 == 0:
		#  ckpt_save_path = ckpt_manager.save()
		#  print ('Saving checkpoint for epoch {} at {}'.format(epoch+1,
		#                                                     ckpt_save_path))

		print ('Epoch {} Loss {:.4f} Accuracy {:.4f}'.format(epoch + 1, train_loss.result(), train_accuracy.result()))

		print ('Time taken for 1 epoch: {} secs\n'.format(time.time() - start))


	# In[76]:


	transformer.save_weights(simCode + "/" + 'transformer')


	# In[77]:


	#Manual example to run the transformer
	if(not(runningInCommandLine)):

		translated = []
		att_weights = []
		max_length_targ = trainLabels[0].shape[0]
		inp_sentence = tf.cast(trainKeys[0], tf.int64)
		encoder_input = tf.expand_dims(inp_sentence, 0)
		
		decoder_input = [word2idxOUT['<start>']]
		output = tf.expand_dims(decoder_input, 0)

		for i in range(max_length_targ):
			enc_padding_mask, combined_mask, dec_padding_mask = create_masks(
				encoder_input, output)

			# predictions.shape == (batch_size, seq_len, vocab_size)
			predictions, attention_weights = transformer(encoder_input, 
														 output,
														 False,
														 enc_padding_mask,
														 combined_mask,
														 dec_padding_mask)

			# select the last word from the seq_len dimension
			predictions = predictions[: ,-1:, :]  # (batch_size, 1, vocab_size)

			predicted_id = tf.cast(tf.argmax(predictions, axis=-1), tf.int32)


			# concatentate the predicted_id to the output which is given to the decoder
			# as its input.
			output = tf.concat([output, predicted_id], axis=-1)

			# return the result if the predicted_id is equal to the end token
			if predicted_id == word2idxOUT['<end>']:
				translated = tf.squeeze(output, axis=0)
				att_weights = attention_weights
				break


		translated = tf.squeeze(output, axis=0), 
		att_weights = attention_weights
		#return tf.squeeze(output, axis=0), attention_weights

		translated


	# In[78]:


	#Returns a list of 
	def evaluateBatch(input_tensor, output_tensor, transf, max_length_inp, max_length_targ, word2idIN, word2idOUT):
		input_batch_size = input_tensor.shape[0]
		#print(input_batch_size)

		encoder_input = input_tensor # tf.expand_dims(input_tensor, 0)
		#print(encoder_input)

		dec_input = tf.repeat(tf.expand_dims([word2idxOUT['<start>']], 0), input_batch_size, axis=0)
		output = dec_input #tf.expand_dims(dec_input, 0)
		
		#print(output)
		
		for t in range(max_length_targ):
			
			enc_padding_mask, combined_mask, dec_padding_mask = create_masks(encoder_input, output)
			#print(enc_padding_mask)
			#print(combined_mask)
			#print(dec_padding_mask)
		
			# predictions.shape == (batch_size, seq_len, vocab_size)
			predictions, attention_weights = transformer(encoder_input, 
													 output,
													 False,
													 enc_padding_mask,
													 combined_mask,
													 dec_padding_mask)

			# select the last word from the seq_len dimension
			predictions = predictions[: ,-1:, :]  # (batch_size, 1, vocab_size)

			predicted_id = tf.cast(tf.argmax(predictions, axis=-1), tf.int32)

		
			# concatentate the predicted_id to the output which is given to the decoder
			# as its input.
			output = tf.concat([output, predicted_id], axis=-1)

			# return the result if the predicted_id is equal to the end token
			#if predicted_id == word2idxOUT['<end>']:
			#    translated = tf.squeeze(output, axis=0)
			#    att_weights = attention_weights
			#    break
		
			#result = tf.concat([result, tf.expand_dims(predicted_id, 1)], axis = 1)

			dec_input = tf.expand_dims(predicted_id, 1)
	   
		return(output, attention_weights) #, np.array(predictionEachPos)) 


	# In[79]:


	if(not(runningInCommandLine)):
		input_tensor = tf.cast(trainKeys[0:10], tf.int64)
		output_tensor = tf.cast(trainLabels[0:10], tf.int64)
		max_length_inp = trainKeys[0].shape[0]
		max_length_targ = trainLabels[0].shape[0]
		output, attention_weights = evaluateBatch(input_tensor, output_tensor, transformer, max_length_inp, max_length_targ, word2idxIN, word2idxOUT)
		output


	# In[80]:


	if(not(runningInCommandLine)):

		def plot_attention_weights(attention, sentence, result, layer):
			fig = plt.figure(figsize=(16, 8))

			#sentence = tokenizer_pt.encode(sentence)

			attention = tf.squeeze(attention[layer], axis=0)

			for head in range(attention.shape[0]):
				ax = fig.add_subplot(2, 4, head+1)

				# plot the attention weights
				ax.matshow(attention[head][:-1, :], cmap='viridis')

				fontdict = {'fontsize': 10}

				ax.set_xticks(range(len(sentence)+2))
				ax.set_yticks(range(len(result)))

				ax.set_ylim(len(result)-1.5, -0.5)

				#ax.set_xticklabels(
				#    ['<start>']+[tokenizer_pt.decode([i]) for i in sentence]+['<end>'], 
				#    fontdict=fontdict, rotation=90)

				#ax.set_yticklabels([tokenizer_en.decode([i]) for i in result 
				#                    if i < tokenizer_en.vocab_size], 
				#                   fontdict=fontdict)

				ax.set_xlabel('Head {}'.format(head+1))

			plt.tight_layout()
			plt.show()


	# In[82]:


	#if(not(runningInCommandLine)):
	#    plot_attention_weights(attention_weights, inp_sentence, output[0], 'decoder_layer4_block2')


	# In[83]:


	def translate(sentence, plot=''):
		result, attention_weights = evaluate(sentence)

		predicted_sentence = tokenizer_en.decode([i for i in result 
												if i < tokenizer_en.vocab_size])  

		print('Input: {}'.format(sentence))
		print('Predicted translation: {}'.format(predicted_sentence))

		if plot:
			plot_attention_weights(attention_weights, sentence, result, plot)


	# In[ ]:





	# In[ ]:





	# In[ ]:





	# In[ ]:




	list_batched_inputs = np.array_split(trainKeys, 100, 0)
	list_batched_outputs =  np.array_split(trainLabels, 100, 0)
	
	for index in range(0,100):
		batched_inputs = list_batched_inputs[index]
		batched_outputs =  list_batched_outputs[index]
		
		(result, attention_plot) = evaluateBatch(batched_inputs, batched_outputs, transformer, NI, NO,  word2idxIN, word2idxOUT)
		batched_predictions = result.numpy()

		textIn = vocaliseBatch(batched_inputs, idx2wordIN, wordSizeIn)
		textPred = vocaliseBatch(result.numpy(),  idx2wordOUT, wordSizeOut)
		textExpected = vocaliseBatch(batched_outputs,  idx2wordOUT, wordSizeOut)


		# In[85]:


		AllResults = []
		f = open(simCode + "/" + "TrainResults.txt", "a")
		if(index == 0):
			f.write("lenout\tlenpred\tlenin\tld\tldnorm\tpred\tout\tin\n")


		# In[86]:


		for i in range(0,min(len(batched_predictions),maxOutput)):
			ld = jellyfish.levenshtein_distance(textPred[i],textExpected[i])
			ldnorm = ld/max(1,len(textPred[i]))

			lenpara = len(textExpected[i])
			lenpredpara = len(textPred[i])
			lenepi = len(textIn[i])

			#AllResults.append([lenpredpara, lenpara, lenepi, ld, ldnorm, textPred[i], textExpected[i], textIn[i]])
			f.write(str(lenpredpara) + "\t" + str(lenpara) + "\t" + str(lenepi) + "\t" + str(ld) + "\t" + str(ldnorm) + "\t" + str(textPred[i]) + "\t" + str(textExpected[i]) + "\t" + str(textIn[i]) + "\n")

		f.close()


		# In[88]:


	list_batched_inputs = np.array_split(testKeys, 100, 0)
	list_batched_outputs =  np.array_split(testLabels, 100, 0)
	
	if(len(testKeys) > 100000):
		list_batched_inputs = np.array_split(testKeys[0:100000], 100, 0)
		list_batched_outputs =  np.array_split(testLabels[0:100000], 100, 0)
	
	for index in range(0,100):
		batched_inputs = list_batched_inputs[index]
		batched_outputs =  list_batched_outputs[index]

		
		(result, attention_plot) = evaluateBatch(batched_inputs, batched_outputs, transformer, NI, NO,  word2idxIN, word2idxOUT)
		batched_predictions = result.numpy()


		# In[89]:


		textIn = vocaliseBatch(batched_inputs, idx2wordIN, wordSizeIn)
		textPred = vocaliseBatch(result.numpy(),  idx2wordOUT, wordSizeOut)
		textExpected = vocaliseBatch(batched_outputs,  idx2wordOUT, wordSizeOut)


		# In[90]:


		AllResults = []
		f = open(simCode + "/" + "TestResults.txt", "a")
		if(index == 0):
			f.write("lenout\tlenpred\tlenin\tld\tldnorm\tpred\tout\tin\n")


		# In[91]:


		for i in range(0,min(len(batched_predictions),maxOutput)):
			ld = jellyfish.levenshtein_distance(textPred[i],textExpected[i])
			ldnorm = ld/max(1,len(textPred[i]))

			lenpara = len(textExpected[i])
			lenpredpara = len(textPred[i])
			lenepi = len(textIn[i])

			#AllResults.append([lenpredpara, lenpara, lenepi, ld, ldnorm, textPred[i], textExpected[i], textIn[i]])
			f.write(str(lenpredpara) + "\t" + str(lenpara) + "\t" + str(lenepi) + "\t" + str(ld) + "\t" + str(ldnorm) + "\t" + str(textPred[i]) + "\t" + str(textExpected[i]) + "\t" + str(textIn[i]) + "\n")

		f.close()


		# In[ ]:





	# In[ ]:




