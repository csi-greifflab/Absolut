# test svm

# import stuff
from sklearn.svm import SVC
from onehotencoder import batchhotEncodingAAStringflat as enc
from onehotencoder import label_binarizer as labin
import pandas as pd
from sklearn.model_selection import train_test_split
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
import sys
import random
from sklearn.decomposition import PCA

def do_pca(infile):
    '''
    performs pca, use 10k samples
    :return:
    '''
    fulldf = pd.read_csv(infile, sep='\t')
    df = fulldf.sample(10000)
    X = np.array(enc(df.Slide))
    y = np.array(labin(df.Label))
    pca = PCA(n_components=2)
    Xt = pca.fit_transform(X)
    print(pca.explained_variance_)
    # print(Xt[:,0])
    data = {'pc1':list(Xt[:,0]), 'pc2': list(Xt[:,1]), 'y':y}
    outdf = pd.DataFrame(data)
    exvar = '%s_%s' % (round(pca.explained_variance_[0],2), round(pca.explained_variance_[1],2))
    outdf['exvar'] = [exvar]*outdf.shape[0]
    outname = infile.split('/')[-1].split('.')[0]
    outfile = 'outfiles/' + outname + '_pc1_pc2.csv'
    outdf.to_csv(outfile, index=False)



def do_pca_aacomp(infile):
    '''
    performs pca, use 10k samples
    :return:
    '''
    fulldf = pd.read_csv(infile, sep='\t')
    df = fulldf.sample(10000)
    # X = np.array(enc(df.Slide))
    X = [[float(val) for val in row.split('_')] for row in df.AAcompoFullSlice]
    y = np.array(labin(df.Label))
    pca = PCA(n_components=2)
    Xt = pca.fit_transform(X)
    print(pca.explained_variance_)
    data = {'pc1':list(Xt[:,0]), 'pc2': list(Xt[:,1]), 'y':y}
    outdf = pd.DataFrame(data)
    exvar = '%s_%s' % (round(pca.explained_variance_[0],2), round(pca.explained_variance_[1],2))
    outdf['exvar'] = [exvar]*outdf.shape[0]
    outname = infile.split('/')[-1].split('.')[0]
    outfile = 'outfiles/' + outname + '_pc1_pc2_aacomp.csv'
    outdf.to_csv(outfile, index=False)


# run stuff
do_pca('dataset/1FBI_X_Task1_BalancedData.txt')
# do_pca_aacomp('dataset/1FBI_X_Task1_BalancedData.txt')
