# test svm

# import stuff
from mpi4py import MPI
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
from prep_data_task4 import prep_data_task4 as pdt
from prep_data_task2 import prep_data_task2_nantigen as pdta
from sklearn.metrics import confusion_matrix
from find_files import find_files as fifi

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
    nsamples = [200, 1000, 5000,10000,20000,50000, 100000,200000,500000, 1000000, 2000000]
    outcontents = []
    ntimes=10
    outseqdf = pd.DataFrame()
    for nsample in nsamples[:]:
        outns = []
        for i in range(ntimes):
            X_train, y_train, X_test, y_test, x_train, x_test = pdt(infile, nsample)
            print(X_train.shape, y_train.shape)
            #print(X_train)
            # clf = SVC(kernel='linear')
            clf = LogisticRegression(random_state=0)
            clf.fit(X_train,y_train)
            s = clf.score(X_test, y_test)
            y_pred = clf.predict(X_test)
            print(y_pred, len(y_pred))
            print(y_test, len(y_test))
            print(x_train, len(x_train))
            print(x_test, len(x_test))
            ntest = len(y_test)
            rep = 'rep%s' % str(i+1)
            seqdict = {'test_seq':x_test, 'pred_label':y_pred, 'test_label':y_test, 'ntrain':[nsample]*ntest, 'rep':[rep]*ntest}
            seqdf = pd.DataFrame(data=seqdict)
            outseqdf = pd.concat([outseqdf, seqdf])
            #confs  = confusion_matrix(y_test, y_pred)
            #confs_str = '\n'.join([' '.join([ str(val) for val in vals]) for vals in confs])
            print(s)
            #print(confs)
            #print(confs_str)
            #print(confs.shape)
            nrep = 'rep%s' % (i + 1)
            ntrain  = len(X_train)
            ntest = len(X_test)
            #outcontents.append([confs_str, s, nrep, 'lr', ntrain, ntest])
            outcontents.append([s, nrep, 'lr', ntrain, ntest])
        outcols = ['accuracy', 'nrep', 'clf', 'ntrain', 'ntest']
        outdf = pd.DataFrame(outcontents, columns=outcols)
        print(outdf)
        outname = 'clf_acc_antigens_outfiles/t4_lr_acc_nsamples_nrep%s_pscheme.csv' % ntimes
        outname_seq = 'clf_acc_antigens_outfiles/t4_lr_seq_nsamples_nrep%s_pscheme.csv' % ntimes
        print(outname)
        outdf.to_csv(outname, index=False)
        outseqdf.to_csv(outname_seq, index=False)


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
                #print(confs)
                #print(confs_str)
                #print(confs.shape)
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



def rf_trainer_nsamples_ntimes_t4(infile):
    '''
    train an svm model, uses phil's data, uses adapted encoder and binarizer from phil, uses phil's preproessing
    uses the multiclass variant of the classifier.
    :param infile:
    :return:
    '''
    fulldf = pd.read_csv(infile, sep='\t')
    nrows = fulldf.shape[0]
    #nsamples = [2000, 10000,50000, 250000]
    nsamples = [100, 2000, 10000,50000]
    #nsamples = [item for item in nsamples if item < nrows]
    print(nsamples)
    outcontents = []
    ntimes=10
    outseqdf = pd.DataFrame()
    for nsample in nsamples[:]:
        outns = []
        for i in range(ntimes):
            X_train, y_train, X_test, y_test, x_train, x_test = pdt(infile, nsample)
            print(X_train.shape, y_train.shape)
            print(X_test.shape, y_test.shape)
            print('total rows %s; iteration %s' % (nrows, i))
            #print(X_train)
            # clf = SVC(kernel='linear')
            #clf = LogisticRegression(random_state=0)
            clf = RandomForestClassifier(max_depth=2, random_state=0)
            clf.fit(X_train,y_train)
            s = clf.score(X_test, y_test)
            y_pred = clf.predict(X_test)
            print(y_pred, len(y_pred))
            print(y_test, len(y_test))
            print(x_train, len(x_train))
            print(x_test, len(x_test))
            ntest = len(y_test)
            rep = 'rep%s' % str(i+1)
            seqdict = {'test_seq':x_test, 'pred_label':y_pred, 'test_label':y_test, 'nsample':[nsample]*ntest, 'rep':[rep]*ntest}
            seqdf = pd.DataFrame(data=seqdict)
            outseqdf = pd.concat([outseqdf, seqdf])
            #confs  = confusion_matrix(y_test, y_pred)
            #confs_str = '\n'.join([' '.join([ str(val) for val in vals]) for vals in confs])
            print(s)
            #print(confs)
            #print(confs_str)
            #print(confs.shape)
            nrep = 'rep%s' % (i + 1)
            ntrain  = len(X_train)
            ntest = len(X_test)
            #outcontents.append([confs_str, s, nrep, 'lr', ntrain, ntest])
            outcontents.append([s, nrep, 'lr', ntrain, ntest])
        outcols = ['accuracy', 'nrep', 'clf', 'ntrain', 'ntest']
        outdf = pd.DataFrame(outcontents, columns=outcols)
        print(outdf)
        outtag = infile.split('/')[-1].split('.')[0]
        print(outtag)
        outname = 'clf_acc_antigens_outfiles_t4/%s_rf_acc_nsamples_nrep%s_pscheme.csv' % (outtag, ntimes)
        outname_seq = 'clf_acc_antigens_outfiles_t4/%s_rf_rawoutput_nsamples_nrep%s_pscheme.csv' % (outtag,ntimes)
        outdf.to_csv(outname, index=False)
        outseqdf.to_csv(outname_seq, index=False)



def rf_trainer_nsamples_ntimes_t4_shuffled(infile):
    '''
    train an svm model, uses phil's data, uses adapted encoder and binarizer from phil, uses phil's preproessing
    uses the multiclass variant of the classifier.
    :param infile:
    :return:
    '''
    fulldf = pd.read_csv(infile, sep='\t')
    nrows = fulldf.shape[0]
    #nsamples = [2000, 10000,50000, 250000]
    nsamples = [100, 2000, 10000,50000]
    #nsamples = [item for item in nsamples if item < nrows]
    print(nsamples)
    outcontents = []
    ntimes=10
    outseqdf = pd.DataFrame()
    for nsample in nsamples[:]:
        outns = []
        for i in range(ntimes):
            X_train, y_train, X_test, y_test, x_train, x_test = pdt(infile, nsample)
            random.shuffle(y_train)
            print(X_train.shape, y_train.shape)
            print(X_test.shape, y_test.shape)
            print('total rows %s; iteration %s' % (nrows, i))
            #print(X_train)
            # clf = SVC(kernel='linear')
            #clf = LogisticRegression(random_state=0)
            clf = RandomForestClassifier(max_depth=2, random_state=0)
            clf.fit(X_train,y_train)
            s = clf.score(X_test, y_test)
            y_pred = clf.predict(X_test)
            print(y_pred, len(y_pred))
            print(y_test, len(y_test))
            print(x_train, len(x_train))
            print(x_test, len(x_test))
            ntest = len(y_test)
            rep = 'rep%s' % str(i+1)
            seqdict = {'test_seq':x_test, 'pred_label':y_pred, 'test_label':y_test, 'nsample':[nsample]*ntest, 'rep':[rep]*ntest}
            seqdf = pd.DataFrame(data=seqdict)
            outseqdf = pd.concat([outseqdf, seqdf])
            #confs  = confusion_matrix(y_test, y_pred)
            #confs_str = '\n'.join([' '.join([ str(val) for val in vals]) for vals in confs])
            print(s)
            #print(confs)
            #print(confs_str)
            #print(confs.shape)
            nrep = 'rep%s' % (i + 1)
            ntrain  = len(X_train)
            ntest = len(X_test)
            #outcontents.append([confs_str, s, nrep, 'lr', ntrain, ntest])
            outcontents.append([s, nrep, 'lr', ntrain, ntest])
        outcols = ['accuracy', 'nrep', 'clf', 'ntrain', 'ntest']
        outdf = pd.DataFrame(outcontents, columns=outcols)
        print(outdf)
        outtag = infile.split('/')[-1].split('.')[0]
        print(outtag)
        outname = 'clf_acc_antigens_outfiles_t4/%s_rf_acc_nsamples_nrep%s_pscheme_shuffled.csv' % (outtag, ntimes)
        outname_seq = 'clf_acc_antigens_outfiles_t4/%s_rf_rawoutput_nsamples_nrep%s_pscheme_shuffled.csv' % (outtag,ntimes)
        outdf.to_csv(outname, index=False)
        outseqdf.to_csv(outname_seq, index=False)



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
#lr_trainer_nsamples_ntimes_nantigen_shuffled_aacomp('dataset/Treated142.txt')
#lr_trainer_nsamples_ntimes_nantigen_aacomp('dataset/Treated142.txt')
#lr_trainer_nsamples_ntimes('dataset/Task4_E_EpiSeq_ParaSeq_NoDeg_NoStar_v2.txt')
#lr_trainer_nsamples_ntimes('dataset/Task4_E_EpiSeq_ParaSeq_NoDeg_NoStar_n10000.txt')
#lr_trainer_nsamples_ntimes_shuffled('dataset/Treated142.txt')


# run files in dataset_t4
infiles = fifi('dataset_t4', 'ns.txt')
infiles = [item for item in infiles if 'NoDeg' in item]
#print(infiles)
#for infile in infiles:
#    df = pd.read_csv(infile, sep='\t')
#    print(df.shape, infile)
#sys.exit()
seqfiles = [item for item in infiles if 'Seq' in item]
agrfiles = [item for item in infiles if 'Agr'in item]
chemfiles = [item for item in infiles if 'Chem'in item]
motiffiles = [item for item in infiles if 'Motif'in item]
sorted_infiles = motiffiles + chemfiles + agrfiles + seqfiles
print(sorted_infiles, len(sorted_infiles))
# filter for motif input files
sorted_infiles=  [item for item in sorted_infiles if 'Motif' in item]
print(sorted_infiles)
#sys.exit()
#for infile in sorted_infiles:
#    rf_trainer_nsamples_ntimes_t4(infile)
#    sys.exit()
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
if rank == 0:
    data = sorted_infiles
else:
    data = None
data = comm.scatter(data, root=0)
infile = data
print(data)
#rf_trainer_nsamples_ntimes_t4(infile)
rf_trainer_nsamples_ntimes_t4_shuffled(infile)
