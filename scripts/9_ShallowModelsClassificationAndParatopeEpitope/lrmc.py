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
from prep_data_task2 import prep_data_task2 as pdt
from prep_data_task2 import prep_data_task2_nantigen as pdta
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


def lr_trainer_nsamples(infile):
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
        clf = LogisticRegression(random_state=0)
        clf.fit(X_train,y_train)
        s = clf.score(X_test, y_test)
        outcontents.append([s, nsample, 'lr'])
    outcols = ['accuracy', 'nseqs', 'clf']
    outdf = pd.DataFrame(outcontents, columns=outcols)
    print(outdf)
    outname = 'outfiles/lr_acc_nsamples.csv'
    outdf.to_csv(outname, index=False)


def lr_trainer_nsamples_shuffled(infile):
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
        random.shuffle(y)
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.5, random_state = 0)
        clf = LogisticRegression(random_state=0)
        clf.fit(X_train,y_train)
        s = clf.score(X_test, y_test)
        outcontents.append([s, nsample, 'lr'])
    outcols = ['accuracy', 'nseqs', 'clf']
    outdf = pd.DataFrame(outcontents, columns=outcols)
    print(outdf)
    outname = 'outfiles/lr_acc_nsamples_shuffled.csv'
    outdf.to_csv(outname, index=False)


def lr_trainer_nsamples_aacomp(infile):
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
        clf = LogisticRegression(random_state=0)
        clf.fit(X_train,y_train)
        s = clf.score(X_test, y_test)
        outcontents.append([s, nsample, 'lr'])
    outcols = ['accuracy', 'nseqs', 'clf']
    outdf = pd.DataFrame(outcontents, columns=outcols)
    print(outdf)
    outname = 'outfiles/lr_acc_nsamples_aacomp.csv'
    outdf.to_csv(outname, index=False)


def lr_trainer_nsamples_len(infile):
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
        X = np.array(df.sizeCDR3).reshape(-1,1)
        y = np.array(labin(df.Label))
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.5, random_state = 0)
        clf = LogisticRegression(random_state=0)
        clf.fit(X_train,y_train)
        s = clf.score(X_test, y_test)
        outcontents.append([s, nsample, 'lr'])
    outcols = ['accuracy', 'nseqs', 'clf']
    outdf = pd.DataFrame(outcontents, columns=outcols)
    print(outdf)
    outname = 'outfiles/lr_acc_nsamples_len.csv'
    outdf.to_csv(outname, index=False)


def lr_trainer_nsamples_ntimes(infile):
    '''
    train an svm model, uses phil's data, uses adapted encoder and binarizer from phil, uses phil's preproessing
    uses the multiclass variant of the classifier.
    :param infile:
    :return:
    '''
    fulldf = pd.read_csv(infile, sep='\t')
    nsamples = [250, 1000, 2500,10000,25000,100000,250000,1000000]
    outcontents = []
    ntimes=10
    outcontents = []
    for nsample in nsamples:
        outns = []
        for i in range(ntimes):
            X_train, y_train, X_test, y_test = pdt(infile, nsample)
            print(X_train.shape, y_train.shape)
            # clf = SVC(kernel='linear')
            clf = LogisticRegression(random_state=0)
            clf.fit(X_train,y_train)
            s = clf.score(X_test, y_test)
            y_pred = clf.predict(X_test)
            confs  = confusion_matrix(y_test, y_pred)
            confs_str = '\n'.join([' '.join([ str(val) for val in vals]) for vals in confs])
            print(s)
            print(confs)
            print(confs_str)
            print(confs.shape)
            nrep = 'rep%s' % (i + 1)
            ntrain  = len(X_train)
            ntest = len(X_test)
            outcontents.append([confs_str, s, nrep, 'lr', ntrain, ntest])
    outcols = ['conf_m', 'accuracy', 'nrep', 'clf', 'ntrain', 'ntest']
    outdf = pd.DataFrame(outcontents, columns=outcols)
    print(outdf)
    outname = 'clf_acc_antigens_outfiles/t2_lr_cm_acc_nsamples_nrep%s_pscheme.csv' % ntimes
    print(outname)
    outdf.to_csv(outname, index=False)


def lr_trainer_nsamples_ntimes_nantigen(infile):
    '''
    train an svm model, uses phil's data, uses adapted encoder and binarizer from phil, uses phil's preproessing
    uses the multiclass variant of the classifier.
    :param infile:
    :return:
    '''
    fulldf = pd.read_csv(infile, sep='\t')
    nsamples = [250, 1000, 2500,10000,25000,100000,250000,1000000]
    nantigens = [5, 10, 25, 50, 75, 100, 140]
    outcontents = []
    ntimes=10
    outcontents = []
    for nantigen in nantigens[:]:
        for nsample in nsamples[:]:
            outns = []
            for i in range(ntimes)[:]:
                X_train, y_train, X_test, y_test = pdta(infile, nsample, nantigen)
                print(X_train.shape, y_train.shape)
                # clf = SVC(kernel='linear')
                clf = LogisticRegression(random_state=0)
                clf.fit(X_train,y_train)
                s = clf.score(X_test, y_test)
                y_pred = clf.predict(X_test)
                confs  = confusion_matrix(y_test, y_pred)
                confs_str = '\n'.join([' '.join([ str(val) for val in vals]) for vals in confs])
                print(s)
                print(confs)
                print(confs_str)
                print(confs.shape)
                nrep = 'rep%s' % (i + 1)
                ntrain  = len(X_train)
                ntest = len(X_test)
                outcontents.append([confs_str, s, nrep, 'lr', ntrain, ntest, nantigen])
    outcols = ['conf_m', 'accuracy', 'nrep', 'clf', 'ntrain', 'ntest', 'nantigen']
    outdf = pd.DataFrame(outcontents, columns=outcols)
    print(outdf)
    outname = 'clf_acc_antigens_outfiles/t2_lr_cm_acc_nsamples_nantigens_nrep%s_pscheme.csv' % ntimes
    print(outname)
    outdf.to_csv(outname, index=False)


def lr_trainer_nsamples_ntimes_shuffled(infile):
    '''
    train an svm model, uses phil's data, uses adapted encoder and binarizer from phil, uses phil's preproessing
    uses the multiclass variant of the classifier.
    :param infile:
    :return:
    '''
    fulldf = pd.read_csv(infile, sep='\t')
    nsamples = [250, 1000, 2500,10000,25000,100000,250000,1000000]
    outcontents = []
    ntimes=10
    outcontents = []
    for nsample in nsamples:
        outns = []
        for i in range(ntimes):
            X_train, y_train, X_test, y_test = pdt(infile, nsample)
            np.random.shuffle(y_train)
            np.random.shuffle(y_test)
            print(X_train.shape, y_train.shape)
            # clf = SVC(kernel='linear')
            clf = LogisticRegression(random_state=0)
            clf.fit(X_train,y_train)
            s = clf.score(X_test, y_test)
            y_pred = clf.predict(X_test)
            confs  = confusion_matrix(y_test, y_pred)
            confs_str = '\n'.join([' '.join([ str(val) for val in vals]) for vals in confs])
            print(s)
            print(confs)
            print(confs_str)
            print(confs.shape)
            nrep = 'rep%s' % (i + 1)
            ntrain  = len(X_train)
            ntest = len(X_test)
            outcontents.append([confs_str, s, nrep, 'lr', ntrain, ntest])
    outcols = ['conf_m', 'accuracy', 'nrep', 'clf', 'ntrain', 'ntest']
    outdf = pd.DataFrame(outcontents, columns=outcols)
    print(outdf)
    outname = 'clf_acc_antigens_outfiles/t2_lr_cm_acc_nsamples_nrep%s_pscheme_shuffled.csv' % ntimes
    print(outname)
    outdf.to_csv(outname, index=False)


def lr_trainer_nsamples_ntimes_nantigen_shuffled(infile):
    '''
    train an svm model, uses phil's data, uses adapted encoder and binarizer from phil, uses phil's preproessing
    uses the multiclass variant of the classifier.
    :param infile:
    :return:
    '''
    fulldf = pd.read_csv(infile, sep='\t')
    nsamples = [250, 1000, 2500,10000,25000,100000,250000,1000000]
    nantigens = [5, 10, 25, 50, 75, 100, 140]
    outcontents = []
    ntimes=10
    outcontents = []
    for nantigen in nantigens[:]:
        for nsample in nsamples[:]:
            outns = []
            for i in range(ntimes)[:]:
                X_train, y_train, X_test, y_test = pdta(infile, nsample, nantigen)
                np.random.shuffle(y_train)
                np.random.shuffle(y_test)
                print(X_train.shape, y_train.shape)
                # clf = SVC(kernel='linear')
                clf = LogisticRegression(random_state=0)
                clf.fit(X_train,y_train)
                s = clf.score(X_test, y_test)
                y_pred = clf.predict(X_test)
                confs  = confusion_matrix(y_test, y_pred)
                confs_str = '\n'.join([' '.join([ str(val) for val in vals]) for vals in confs])
                print(s)
                print(confs)
                print(confs_str)
                print(confs.shape)
                nrep = 'rep%s' % (i + 1)
                ntrain  = len(X_train)
                ntest = len(X_test)
                outcontents.append([confs_str, s, nrep, 'lr', ntrain, ntest, nantigen])
    outcols = ['conf_m', 'accuracy', 'nrep', 'clf', 'ntrain', 'ntest', 'nantigen']
    outdf = pd.DataFrame(outcontents, columns=outcols)
    print(outdf)
    outname = 'clf_acc_antigens_outfiles/t2_lr_cm_acc_nsamples_nantigens_nrep%s_pscheme_shuffled.csv' % ntimes
    print(outname)
    outdf.to_csv(outname, index=False)

def get_aacomp(X):
    '''
    transforms onehot to aacomp
    '''
    nchar = 11
    X_out = []
    for x in X:
        xs = np.array_split(x,nchar)
        sumxs = np.sum(xs, axis=0)
        aafrac = sumxs/nchar
        X_out.append(aafrac)
    return X_out 

def lr_trainer_nsamples_ntimes_nantigen_shuffled_aacomp(infile):
    '''
    train an svm model, uses phil's data, uses adapted encoder and binarizer from phil, uses phil's preproessing
    uses the multiclass variant of the classifier.
    :param infile:
    :return:
    '''
    fulldf = pd.read_csv(infile, sep='\t')
    nsamples = [250, 1000, 2500,10000,25000,100000,250000,1000000]
    nantigens = [5, 10, 25, 50, 75, 100, 140]
    outcontents = []
    ntimes=10
    outcontents = []
    for nantigen in nantigens[:]:
        for nsample in nsamples[:]:
            outns = []
            for i in range(ntimes)[:]:
                X_train, y_train, X_test, y_test = pdta(infile, nsample, nantigen)
                print(X_train.shape)
                X_train = get_aacomp(X_train)
                X_test = get_aacomp(X_test)
                np.random.shuffle(y_train)
                np.random.shuffle(y_test)
                # clf = SVC(kernel='linear')
                clf = LogisticRegression(random_state=0)
                clf.fit(X_train,y_train)
                s = clf.score(X_test, y_test)
                y_pred = clf.predict(X_test)
                confs  = confusion_matrix(y_test, y_pred)
                confs_str = '\n'.join([' '.join([ str(val) for val in vals]) for vals in confs])
                print(s)
                print(confs)
                print(confs_str)
                print(confs.shape)
                nrep = 'rep%s' % (i + 1)
                ntrain  = len(X_train)
                ntest = len(X_test)
                outcontents.append([confs_str, s, nrep, 'lr', ntrain, ntest, nantigen])
    outcols = ['conf_m', 'accuracy', 'nrep', 'clf', 'ntrain', 'ntest', 'nantigen']
    outdf = pd.DataFrame(outcontents, columns=outcols)
    print(outdf)
    outname = 'clf_acc_antigens_outfiles/t2_lr_cm_acc_nsamples_nantigens_nrep%s_pscheme_aacomp_shuffled.csv' % ntimes
    print(outname)
    outdf.to_csv(outname, index=False)


def lr_trainer_nsamples_ntimes_nantigen_aacomp(infile):
    '''
    train an svm model, uses phil's data, uses adapted encoder and binarizer from phil, uses phil's preproessing
    uses the multiclass variant of the classifier.
    :param infile:
    :return:
    '''
    fulldf = pd.read_csv(infile, sep='\t')
    nsamples = [250, 1000, 2500,10000,25000,100000,250000,1000000]
    nantigens = [5, 10, 25, 50, 75, 100, 140]
    outcontents = []
    ntimes=10
    outcontents = []
    for nantigen in nantigens[:]:
        for nsample in nsamples[:]:
            outns = []
            for i in range(ntimes)[:]:
                X_train, y_train, X_test, y_test = pdta(infile, nsample, nantigen)
                print(X_train.shape)
                X_train = get_aacomp(X_train)
                X_test = get_aacomp(X_test)
                #np.random.shuffle(y_train)
                #np.random.shuffle(y_test)
                # clf = SVC(kernel='linear')
                clf = LogisticRegression(random_state=0)
                clf.fit(X_train,y_train)
                s = clf.score(X_test, y_test)
                y_pred = clf.predict(X_test)
                confs  = confusion_matrix(y_test, y_pred)
                confs_str = '\n'.join([' '.join([ str(val) for val in vals]) for vals in confs])
                print(s)
                print(confs)
                print(confs_str)
                print(confs.shape)
                nrep = 'rep%s' % (i + 1)
                ntrain  = len(X_train)
                ntest = len(X_test)
                outcontents.append([confs_str, s, nrep, 'lr', ntrain, ntest, nantigen])
    outcols = ['conf_m', 'accuracy', 'nrep', 'clf', 'ntrain', 'ntest', 'nantigen']
    outdf = pd.DataFrame(outcontents, columns=outcols)
    print(outdf)
    outname = 'clf_acc_antigens_outfiles/t2_lr_cm_acc_nsamples_nantigens_nrep%s_pscheme_aacomp.csv' % ntimes
    print(outname)
    outdf.to_csv(outname, index=False)

# run stuff
# svm_trainer('dataset/1FBI_X_Task1_BalancedData.txt')
# lr_trainer_nsamples('dataset/1FBI_X_Task1_BalancedData.txt')
# lr_trainer_nsamples_shuffled('dataset/1FBI_X_Task1_BalancedData.txt')
# lr_trainer_nsamples_aacomp('dataset/1FBI_X_Task1_BalancedData.txt')
# lr_trainer_nsamples_len('dataset/1FBI_X_Task1_BalancedData.txt')
#lr_trainer_nsamples_ntimes('dataset/Treated142.txt')
#lr_trainer_nsamples_ntimes_shuffled('dataset/Treated142.txt')
#lr_trainer_nsamples_ntimes_nantigen('dataset/Treated142.txt')
#lr_trainer_nsamples_ntimes_nantigen_shuffled('dataset/Treated142.txt')
lr_trainer_nsamples_ntimes_nantigen_shuffled_aacomp('dataset/Treated142.txt')
#lr_trainer_nsamples_ntimes_nantigen_aacomp('dataset/Treated142.txt')
