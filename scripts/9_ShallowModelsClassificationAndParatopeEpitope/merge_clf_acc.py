# merges acc files in clf_acc_antigens_outfiles

# import stuff
from find_files import find_files as fifi
import pandas as pd


def merge_clf_acc_files():
    '''
    merges acc files
    '''
    infiles = fifi('clf_acc_antigens_outfiles', 'n3.csv')
    outdf = pd.concat([pd.read_csv(infile) for infile in infiles])
    outname = 'clf_acc_antigens_outfiles/clf_acc_antigens_merged.csv'
    outdf.to_csv(outname, index=False)





# run stuff
merge_clf_acc_files()
