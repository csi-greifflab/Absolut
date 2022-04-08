# This script takes a two-columns file (epitope, paratope, more columns ignored) and performs umap
# between the paratopes using levenshtein distance and for each epitope value separately
# usage:
# thisScript.py tabSeparatedFile.txt outputFile=Same+Annotated.txt maxSizeBlocks=10000 colValues=0 colClasses=1
# Note: the first line of the input file is excluded (considered as header)

import numpy as np
import os
import csv
import sys  # for retrieving command line args
import Levenshtein as lv

#Local settings + default arguments
runAsCommandLine = True
LocationDataFiles = ""
maxLines = 100000000   #Biggest size (lines) of input file

# Decided from command line
fileIn = "D:/pprobert/Datasets/Task4/Task4_D_EpiChem_ParaChem.txt_Balanced_2000.txt"
fileOut = "TestProd.txt"

maxSizeBlocks = 10000
perEpi = True
colValues = 0
colClasses = 1

if(runAsCommandLine):
    if(len(sys.argv) > 1):
        fileIn = str(sys.argv[1])

    if(len(sys.argv) > 2):
        fileOut = str(sys.argv[2])

    if(len(sys.argv) > 3):
        maxSizeBlocks = int(sys.argv[3])

    if(len(sys.argv) > 4):
        colValues = int(sys.argv[4])

    if(len(sys.argv) > 5):
        colClasses = int(sys.argv[5])

        
print("Clusterize \nFileIn", fileIn, "\nFileOut", fileOut, "\nMaxSizeBlocks", maxSizeBlocks,  
      "colValues", colValues, "colClasses", colClasses)

        
print("Clusterize \nFileIn", fileIn, "\nFileOut", fileOut, "\nMaxSizeBlocks", maxSizeBlocks, "perEpi", perEpi, 
      "colValues", colValues, "colClasses", colClasses)
	  
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px

import scipy.spatial as sp, scipy.cluster.hierarchy as hc
import scipy
from scipy.cluster.hierarchy import ward, fcluster
#pip install umap.learn
import umap.umap_ as umap
import seaborn as sns; sns.set()
import hdbscan

from random import randint
np.random.seed(1984)

def readParaEpiFile(fileName):
    epiParaSeq = open(fileName, newline = '')   #one line is a text with \t and \n                                                                              
    data_reader = csv.reader(epiParaSeq, delimiter='\t') #transform lines into lists 
    paratopes = []
    epitopes = []
    paraepi = {}
    cpt = 0
    for line in data_reader:
        if((cpt > 0) and (cpt <= maxLines)): #Excludes the first line 
            paratopes.append(line[colClasses])
            epitopes.append(line[colValues])
            if line[colValues] in paraepi:
                paraepi[line[colValues]].append(line[colClasses])
            else:
                paraepi[line[colValues]] = [line[colClasses]] 
        cpt = cpt + 1

    print("Got ", len(paratopes), " paratopes and ", len(epitopes), " epitopes from ", fileName)
    return (paratopes, epitopes, paraepi)

[para, epi, paraepi] = readParaEpiFile(fileIn)
print("Got ", len(para), " unique paratopes and ", len(set(epi)), " unique epitopes")



def createClustersPerKey(mydict, maxdim):
    file_object = open(fileOut, 'w')

    id_cluster = 0;
    for k in mydict.keys():
        available = len(mydict[k]) 
        print(k, 'available', available, "start at cluster", id_cluster)
        ab_sequence = list(mydict[k])
        if(available > maxdim):
            ab_sequence = ab_sequence[0:maxdim]
            
        distmat = np.zeros((len(ab_sequence), len(ab_sequence)))
        distmat_norm = np.zeros((len(ab_sequence), len(ab_sequence)))

        for i in range(len(ab_sequence)):
            if i%500 == 0:
                print(i)
            for j in range(len(ab_sequence)):
                distmat[i,j] = lv.distance(ab_sequence[i],ab_sequence[j]) 
                distmat_norm[i,j] = distmat[i,j] / (max(1, max(len(ab_sequence[i]), len(ab_sequence[j]))))
                
        U = umap.UMAP(metric='precomputed', n_neighbors=30)
        XY = U.fit_transform(distmat_norm)

        labels2 = hdbscan.HDBSCAN(
            min_samples=50,
            min_cluster_size=50,
        ).fit_predict(XY)
        
        if(len(labels2) != len(ab_sequence)):
            print("ERR: labels have different sizes, key", k)
        
        for i in range(0,len(ab_sequence)):
            file_object.write(str(k) + "\t" + str(ab_sequence[i]) + "\t" + str(int(labels2[i]) + id_cluster) + "\n")            
        
        
        #Only one color since we are inside one 
        fig = px.scatter(XY, x=0, y=1, hover_data=[ab_sequence], color_continuous_scale=px.colors.sequential.Viridis)
        #fig.write_image(file=str(k) + "Supervised.png", format = '.png')
        #fig.show()
        
        fig2 = px.scatter(XY, x=0, y=1, color=labels2,hover_data=[ab_sequence], color_continuous_scale=px.colors.sequential.Viridis)
        #fig2.write_image(file=str(k) + "Automated.png", format = '.png')
    
        #fig, axs = plt.subplots(2)
        #fig.suptitle('Vertically stacked subplots')
        #axs[0].plot(x, y)
        #axs[1].plot(x, -y)


        newclusters = len(set(labels2))  
        id_cluster = id_cluster + newclusters
    
    file_object.close()

createClustersPerKey(paraepi, maxSizeBlocks)