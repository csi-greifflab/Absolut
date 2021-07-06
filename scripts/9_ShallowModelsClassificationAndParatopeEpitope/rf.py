# test svm

# import stuff
from sklearn.svm import SVC
from onehotencoder import batchhotEncodingAAStringflat as enc
from onehotencoder import label_binarizer as labin
import pandas as pd
from sklearn.model_selection import train_test_split
import numpy as np
from sklearn.ensemble import RandomForestClassifier
import sys
import random
from find_files import find_files as fifi
from data_splitter import split80_20 as ds
from sklearn.metrics import confusion_matrix

def svm_trainer(infile):
    '''
    train an svm model, uses phil's data, uses adapted encoder and binarizer from phil
    :param infile:
    :return:
    '''
    df = pd.read_csv(infile, sep='\t')
    print(df)
    df = df.sample(100)
    X = np.array(enc(df.Slide))
    print([len(item) for item in X[:3]])
    y = np.array(labin(df.Label))
    print(y[:13])
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.5, random_state = 0)
    print(X_train.shape)
    clf = SVC(kernel='linear')
    clf.fit(X_train,y_train)
    s = clf.score(X_test, y_test)
    print(s)


def rf_trainer_nsamples(infile):
    '''
    train an svm model, uses phil's data, uses adapted encoder and binarizer from phil
    :param infile:
    :return:
    '''
    fulldf = pd.read_csv(infile, sep='\t')
    nsamples = [100, 1000, 10000]
    outcontents = []
    for nsample in nsamples:
        df = fulldf.sample(nsample)
        X = np.array(enc(df.Slide))
        y = np.array(labin(df.Label))
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.5, random_state = 0)
        clf = RandomForestClassifier(max_depth=2, random_state=0)
        clf.fit(X_train,y_train)
        s = clf.score(X_test, y_test)
        outcontents.append([s, nsample, 'rf'])
    outcols = ['accuracy', 'nseqs', 'clf']
    outdf = pd.DataFrame(outcontents, columns=outcols)
    print(outdf)
    outname = 'outfiles/rf_acc_nsamples.csv'
    outdf.to_csv(outname, index=False)


def rf_trainer_nsamples_shuffled(infile):
    '''
    train an svm model, uses phil's data, uses adapted encoder and binarizer from phil
    :param infile:
    :return:
    '''
    fulldf = pd.read_csv(infile, sep='\t')
    nsamples = [100, 1000, 10000]
    outcontents = []
    for nsample in nsamples:
        df = fulldf.sample(nsample)
        X = np.array(enc(df.Slide))
        y = np.array(labin(df.Label))
        print(y)
        random.shuffle(y)
        print(y)
        # y = random.shuffle(y)
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.5, random_state = 0)
        clf = RandomForestClassifier(max_depth=2, random_state=0)
        clf.fit(X_train,y_train)
        s = clf.score(X_test, y_test)
        outcontents.append([s, nsample, 'rf'])
    outcols = ['accuracy', 'nseqs', 'clf']
    outdf = pd.DataFrame(outcontents, columns=outcols)
    print(outdf)
    outname = 'outfiles/rf_acc_nsamples_shuffled.csv'
    outdf.to_csv(outname, index=False)


def rf_trainer_nsamples_aacomp(infile):
    '''
    train an svm model, uses phil's data, uses adapted encoder and binarizer from phil
    :param infile:
    :return:
    '''
    fulldf = pd.read_csv(infile, sep='\t')
    nsamples = [100, 1000, 10000]
    outcontents = []
    for nsample in nsamples:
        df = fulldf.sample(nsample)
        # X = np.array(enc(df.Slide))
        X = [[float(val) for val in row.split('_')] for row in df.AAcompoFullSlice]
        y = np.array(labin(df.Label))
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.5, random_state = 0)
        clf = RandomForestClassifier(max_depth=2, random_state=0)
        clf.fit(X_train,y_train)
        s = clf.score(X_test, y_test)
        outcontents.append([s, nsample, 'rf'])
    outcols = ['accuracy', 'nseqs', 'clf']
    outdf = pd.DataFrame(outcontents, columns=outcols)
    print(outdf)
    outname = 'outfiles/rf_acc_nsamples_aacomp.csv'
    outdf.to_csv(outname, index=False)


def rf_trainer_nsamples_len(infile):
    '''
    train an svm model, uses phil's data, uses adapted encoder and binarizer from phil
    :param infile:
    :return:
    '''
    fulldf = pd.read_csv(infile, sep='\t')
    nsamples = [100, 1000, 10000]
    outcontents = []
    for nsample in nsamples:
        df = fulldf.sample(nsample)
        X = np.array(df.sizeCDR3).reshape(-1,1)
        y = np.array(labin(df.Label))
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.5, random_state = 0)
        clf = RandomForestClassifier(max_depth=2, random_state=0)
        clf.fit(X_train,y_train)
        s = clf.score(X_test, y_test)
        outcontents.append([s, nsample, 'rf'])
    outcols = ['accuracy', 'nseqs', 'clf']
    outdf = pd.DataFrame(outcontents, columns=outcols)
    print(outdf)
    outname = 'outfiles/rf_acc_nsamples_len.csv'
    outdf.to_csv(outname, index=False)


def rf_trainer_nsamples_ntimes(infile):
    '''
    train an svm model, uses phil's data, uses adapted encoder and binarizer from phil
    :param infile:
    :return:
    '''
    fulldf = pd.read_csv(infile, sep='\t')
    nsamples = [100, 1000, 10000]
    outcontents = []
    ntimes=10
    for nsample in nsamples:
        outns = []
        for i in range(ntimes):
            df = fulldf.sample(nsample)
            X = np.array(enc(df.Slide))
            y = np.array(labin(df.Label))
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.5, random_state = 0)
            # clf = SVC(kernel='linear')
            # clf = LogisticRegression(random_state=0)
            clf = RandomForestClassifier(max_depth=2, random_state=0)
            clf.fit(X_train,y_train)
            s = clf.score(X_test, y_test)
            outns.append(s)
        mean_s = np.mean(outns)
        sd_s = np.std(outns)
        print(mean_s, sd_s)
        outcontents.append([mean_s, sd_s, nsample, 'rf'])
    outcols = ['mean_accuracy', 'stdev', 'nseqs', 'clf']
    outdf = pd.DataFrame(outcontents, columns=outcols)
    outname = 'outfiles/rf_acc_nsamples_n%s.csv' % ntimes
    outdf.to_csv(outname, index=False)



def rf_trainer_nsamples_ntimes_antigen(infile):
    '''
    train an svm model, uses phil's data, uses adapted encoder and binarizer from phil
    :param infile:
    :return:
    '''
    fulldf = pd.read_csv(infile, sep='\t')
    nameparts = infile.split('/')
    filename = nameparts[-1]
    #antigen = filename.split('_')[0]
    antigen = '_'.join(filename.split('_')[:2])
    nsamples = [100, 1000, 10000]
    outcontents = []
    ntimes=3
    for nsample in nsamples:
        outns = []
        for i in range(ntimes):
            df = fulldf.sample(nsample)
            X = np.array(enc(df.Slide))
            y = np.array(labin(df.Label))
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.5, random_state = 0)
            # clf = SVC(kernel='linear')
            # clf = LogisticRegression(random_state=0)
            clf = RandomForestClassifier(max_depth=2, random_state=0)
            clf.fit(X_train,y_train)
            s = clf.score(X_test, y_test)
            outns.append(s)
        mean_s = np.mean(outns)
        sd_s = np.std(outns)
        print(mean_s, sd_s)
        outcontents.append([mean_s, sd_s, nsample, 'rf', antigen])
    outcols = ['mean_accuracy', 'stdev', 'nseqs', 'clf', 'antigen']
    outdf = pd.DataFrame(outcontents, columns=outcols)
    # outname = 'outfiles/rf_acc_nsamples_n%s.csv' % ntimes
    outname = 'clf_acc_antigens_outfiles/%s_rf_acc_nsamples_n%s.csv' % (antigen, ntimes)
    outdf.to_csv(outname, index=False)


def rf_trainer_nsamples_ntimes_antigen2(infile):
    '''
    train an svm model, uses phil's data, uses adapted encoder and binarizer from phil
    uses data splitter to 
    :param infile:
    :return:
    '''
    #fulldf = pd.read_csv(infile, sep='\t')
    traindf, testdf  = ds(infile)
    nameparts = infile.split('/')
    filename = nameparts[-1]
    #antigen = filename.split('_')[0]
    antigen = '_'.join(filename.split('_')[:2])
    nsamples = [100, 1000, 5000, 10000, 25000, 50000]
    outcontents = []
    ntimes=10
    for nsample in nsamples[:]:
        outns = []
        for i in range(ntimes):
            straindf = traindf.sample(nsample)
            X_train, y_train = np.array(enc(straindf.Slide)), np.array(labin(straindf.Label))
            X_test, y_test = np.array(enc(testdf.Slide)), np.array(labin(testdf.Label))
            ntrain  = len(X_train)
            ntest = len(X_test)
            #clf = LogisticRegression(random_state=0)
            clf = RandomForestClassifier(max_depth=2, random_state=0)
            clf.fit(X_train,y_train)
            y_pred = clf.predict(X_test)
            tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
            s = clf.score(X_test, y_test)
            #print(tn,fp,fn,tp, s)
            outns.append(s)
            nrep = 'rep%s' % (i + 1)
            outcontents.append([tn, fp, fn, tp, s, nrep, 'rf', antigen, ntrain,ntest])
        #mean_s = np.mean(outns)
        #sd_s = np.std(outns)
        #print(mean_s, sd_s)
        #outcontents.append([mean_s, sd_s, nsample, 'lr', antigen, ntrain,ntest])
    outcols = ['tn', 'fp', 'fn', 'tp', 'accuracy', 'nrep','clf', 'antigen', 'ntrain', 'ntest']
    outdf = pd.DataFrame(outcontents, columns=outcols)
    outname = 'clf_acc_antigens_outfiles/%s_rf_cm_acc_nsamples_nrep%s_pscheme.csv' % (antigen, ntimes)
    outdf.to_csv(outname, index=False)

# run stuff
# svm_trainer('dataset/1FBI_X_Task1_BalancedData.txt')
# rf_trainer_nsamples('dataset/1FBI_X_Task1_BalancedData.txt')
# rf_trainer_nsamples_shuffled('dataset/1FBI_X_Task1_BalancedData.txt')
# rf_trainer_nsamples_aacomp('dataset/1FBI_X_Task1_BalancedData.txt')
# rf_trainer_nsamples_len('dataset/1FBI_X_Task1_BalancedData.txt')
# rf_trainer_nsamples_ntimes('dataset/1FBI_X_Task1_BalancedData.txt')
infiles = fifi('/storage/pprobert/Task1', '.txt')
infiles = [item for item in infiles if 'MvsL_SlicesBalancedData.txt' in item]
print(infiles)
print(len(infiles))
for infile in infiles:
    rf_trainer_nsamples_ntimes_antigen2(infile)
