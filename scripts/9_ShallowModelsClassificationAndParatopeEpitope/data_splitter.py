#splits data prior to subsampling to follow phil's scheme

# import stuff
import pandas as pd
import sys

def split80_20(infile):
    '''
    splits data to 80 and 20 percent
    '''
    df = pd.read_csv(infile, sep='\t')
    bdf = df[df.Label=='Binder']
    nbdf =df[df.Label!='Binder']
    trainfrac = 0.8
    testfrac = 1-trainfrac
    bdf = bdf.sample(frac=1).reset_index(drop=True)
    nbdf = nbdf.sample(frac=1).reset_index(drop=True)
    bdftrainindex = int(trainfrac*bdf.shape[0])
    nbdftrainindex = int(trainfrac*nbdf.shape[0])
    traindf = pd.concat([bdf.iloc[0:bdftrainindex,:], nbdf.iloc[0:nbdftrainindex,:]])
    testdf = pd.concat([bdf.iloc[bdftrainindex:,:], nbdf.iloc[nbdftrainindex:,:]])
    return traindf, testdf

#run stuff
#split80_20('/storage/pprobert/Task1/3HI6_A_Task1aMvsL_BalancedData.txt')
