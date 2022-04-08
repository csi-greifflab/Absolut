#!/usr/bin/env python3

#Script that performs DLAB-vs on Absolut! lattice pairs (of size 6x6)
#Usage in command line:
#python ThisFile.py nAntigens=100 nPairsTot=10000 startegyNegatives=1 condition=1 nRotations=20 groupSize=50 fnatBind=0.9 fnatDLABNegative=0.1 balancingStrategy=0 nRepeats=1 nEpochs=5 batch_size=2000

#New things, now it includes 1 more strategy for negatives (separating Others(100) versus low fnat (1000) 

from __future__ import division, print_function, absolute_import
import numpy as np
import tensorflow as tf
import csv
import sys
import math
import random
import os
import math
from collections import Counter
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
random.seed(a=None, version=2)



#Meaning of command-line parameters
# nAntigens = 100       - number of antigens we learn from [so we could use the others as transferrability]
# nPairsTot = 10000     - total number of pairs BEFORE data enhancement by rotation
# strategyNegatives=1   - 1 = DLAB negatives other pairs. add 10 = with low affinity. add 100 = with absolut non-binders. add 1000 = with low fnat of binding pairs 
#                         interesting ones: 1001 (DLAB), 1111(all), and 110 (no DLAB, only non binders), 1 (DLAB, only nonbinding pairs but no low fnat)
# strategyNegativeTesting=1   - 1 = DLAB negatives other pairs. add 10 = with low affinity. add 100 = with absolut non-binders. add 1000 = with low fnat of binding pairs 
# condition=1           - 1 one-hot, 3, shuffled, 11 = Chemical, 13 = Chemical shuffled
# nRotations            - 0 for no data enhancement
# groupSize = 1         - when pick a pair, how many poses of the same pair we consider
# fnatBind = 0.9
# fnatDLABNegative=0.1
# balancingStrategy = 0/1/2    # 0: unbalanced; 1:Includes weights during fitting; 2:Find a way to balance the inputs inside batches
# nRepeats                - will repeat the SAME selection of antigens multiple times. Restart multiple times with nRepeats=1 to have different antigens each time
# nEpochs
# batch_size

# Manually set options
runningInCommandLine = True     # with False it runs some examples / tests
LocationDataFiles = ""            # The preprocessed files are expected to be there
size_lattice = 6
stopDataLeakage = True           # Avoids data enhancement of the same pose to be in train and test
doAllAntigens = True              # Allows to use data on all antigens (if False, only considers 1ADQ_A, for debugging)
removeFirstDimention = True        # In one-hot encoding, removes the dimension meaning "empty grid position", True => 20 dims. False => 21 dims


# Default parameters similar to DLAB
nAntigens = 1
nPairsTot = 5000            #this is before data augmentation => the training size will be 0.8 * nPairsTot * nRotations
strategyNegatives = 1001       #only DLAB negatives, no other ones
strategyNegativeTesting = 1110 # We will show two test: itself (test=strategyNegatives) and external (ext=strategyNegativeTesting)
                               # Further, we will also balance the external dataset to have the same amount of positives as the test
condition = 1
nRotations = 1 #20
groupSize = 50                #this is already done during preprocessing, so we will just read preprocessed files with this value of groupSize
fnatBind = 0.9
fnatDLABNegative=0.1
balancingStrategy = 0
nRepeats = 1
nEpochs = 5
batch_size = 2000





#Filtering command line arguments
repeat = 1
listAntigens = ["1ADQ_A"]
if(doAllAntigens):
    listAntigens = ["1ADQ_A", "1FBI_X", "1FNS_A", "1FSK_A", "1H0D_C", "1JPS_T", "1KB5_AB", "1NCB_N", "1NSN_S", "1OAZ_A", "1OB1_C", "1OSP_O", "1PKQ_J", "1QFW_AB", "1RJL_C", "1S78_B", "1TQB_A", "1WEJ_F", "1YJD_C", "1ZTX_E", "2ARJ_RQ", "2B2X_A", "2FD6_AU", "2HFG_R", "2IH3_C", "2JEL_P", "2Q8A_A", "2R29_A", "2R4R_A", "2R56_A", "2UZI_R", "2VXQ_A", "2VXT_I", "2W9E_A", "2WUC_I", "2XQB_A", "2XWT_C", "2YC1_C", "2YPV_A", "2ZCH_P", "3BGF_S", "3BN9_A", "3CVH_ABC", "3DVG_XY", "3EFD_K", "3GI9_C", "3HI6_A", "3JBQ_B", "3KJ4_A", "3KR3_D", "3KS0_J", "3L5X_A", "3L95_X", "3MJ9_A", "3NCY_A", "3NFP_I", "3NH7_A", "3Q3G_E", "3R08_E", "3R1G_B", "3RAJ_A", "3RKD_A", "3RVV_A", "3SKJ_E", "3SQO_A", "3TT1_A", "3U9P_C", "3UBX_A", "3V6O_A", "3VG9_A", "3VRL_C", "3WD5_A", "4AEI_A", "4CAD_C", "4DKE_A", "4H88_A", "4HC1_B", "4HJ0_B", "4I18_R", "4I77_Z", "4K24_A", "4K3J_A", "4K9E_C", "4KI5_M", "4KXZ_A", "4LU5_B", "4MXV_B", "4N9G_C", "4NP4_A", "4OII_A", "4OKV_E", "4PP1_A", "4QCI_D", "4QEX_A", "4QNP_A", "4QWW_A", "4R9Y_D", "4RGM_S", "4U1G_A", "4U6V_A", "4WV1_F", "4Y5V_C", "4YPG_C", "4YUE_C", "4ZFG_A", "4ZFO_F", "4ZSO_E", "5B8C_C", "5BVP_I", "5C0N_A", "5C7X_A", "5CZV_A", "5D93_A", "5DFV_A", "5DHV_M", "5DMI_A", "5DO2_B", "5E8D_A", "5E8E_LH", "5E94_G", "5EII_G", "5EPM_C", "5EU7_A", "5EZO_A", "5F3B_C", "5FB8_C", "5H35_C", "5HDQ_A", "5HI4_B", "5IKC_M", "5J13_A", "5JW4_A", "5JZ7_A", "5KN5_C", "5KTE_A", "5L0Q_A", "5LQB_A", "5MES_A", "5T5F_A", "5TH9_A", "5TLJ_X", "5TZ2_C"]
    listAntigens = ["1ADQ_A", "1FBI_X", "1FNS_A", "1FSK_A", "1H0D_C", "1JPS_T", "1KB5_AB", "1NCB_N", "1NSN_S", "1OAZ_A", "1OB1_C", "1OSP_O", "1PKQ_J", "1QFW_AB", "1RJL_C", "1TQB_A", "1WEJ_F", "1YJD_C", "1ZTX_E", "2ARJ_RQ", "2B2X_A", "2FD6_AU", "2HFG_R", "2IH3_C", "2JEL_P", "2Q8A_A", "2R29_A", "2R4R_A", "2R56_A", "2UZI_R", "2VXQ_A", "2VXT_I", "2W9E_A", "2WUC_I", "2XQB_A", "2XWT_C", "2YC1_C", "2YPV_A", "2ZCH_P", "3BGF_S", "3BN9_A", "3DVG_XY", "3EFD_K", "3GI9_C", "3HI6_A", "3JBQ_B", "3KJ4_A", "3KR3_D", "3KS0_J", "3L5X_A", "3L95_X", "3MJ9_A", "3NFP_I", "3NH7_A", "3Q3G_E", "4R9Y_D", "4RGM_S", "4U1G_A", "4U6V_A", "4WV1_F", "4Y5V_C", "4YPG_C", "4YUE_C", "4ZFG_A", "4ZFO_F", "4ZSO_E", "5B8C_C", "5BVP_I", "5C0N_A", "5C7X_A", "5CZV_A", "5D93_A", "5DFV_A", "5DHV_M", "5DMI_A", "5DO2_B", "5E8D_A", "5E8E_LH", "5E94_G", "5JW4_A"]
    listAntigens = ["1ADQ_A", "1FBI_X", "1FNS_A", "1FSK_A", "1H0D_C", "1JPS_T", "1KB5_AB", "1NCB_N", "1NSN_S", "1OAZ_A", "1OB1_C"]
    #For 0.9 0.1 x 10
    listAntigens = ["1ADQ_A", "1FBI_X", "1FSK_A", "1H0D_C", "1KB5_AB", "1NCB_N", "1OAZ_A", "1OB1_C", "1QFW_AB", "1RJL_C", "1TQB_A", "1WEJ_F", "1YJD_C", "1ZTX_E", "2B2X_A"] 

  
 
 
 
if(not(runningInCommandLine)):
    import matplotlib.pyplot as plt
    
if(runningInCommandLine):
    if(len(sys.argv) > 1):
        nAntigens = int(sys.argv[1])
        if(doAllAntigens == False and nAntigens > 1):
            print("ERR: you are in debug mode (doAllAntigens == False), only one antigen allowed (will be 1ADQ_A)")
            sys.exit()

    if(len(sys.argv) > 2):
        nPairsTot = int(sys.argv[2])

    if(len(sys.argv) > 3):
        strategyNegatives = int(sys.argv[3])

    if(len(sys.argv) > 4):
        strategyNegativeTesting = int(sys.argv[4])

    if(len(sys.argv) > 5):
        condition = int(sys.argv[5])

    if(len(sys.argv) > 6):
        nRotations = int(sys.argv[6])

    if(len(sys.argv) > 7):
        groupSize = int(sys.argv[7])
        
    if(len(sys.argv) > 8):
        fnatBind = float(sys.argv[8])

    if(len(sys.argv) > 9):
        fnatDLABNegative = float(sys.argv[9])

    if(len(sys.argv) > 10):
        balancingStrategy = int(sys.argv[10])

    if(len(sys.argv) > 11):
        nRepeats = int(sys.argv[11])

    if(len(sys.argv) > 12):
        nEpochs = int(sys.argv[12])

    if(len(sys.argv) > 13):
        batch_size = int(sys.argv[13])

#sanity checks
if(nRepeats > 10):
    print("ERR: We don't allow more than 10 repeats because we have only preprocessed 10 different files per antigen. Comment this line if you have preprocessed more")
    sys.exit()
    
    
    

saveStrategyNegatives = strategyNegatives

#translating strategyNegatives into booleans
includeLowFnatAsNeg = False
if(strategyNegatives >= 1000):
    includeLowFnatAsNeg = True
    strategyNegatives = strategyNegatives - 1000

includeNonBinders = False
if(strategyNegatives >= 100):
    includeNonBinders = True
    strategyNegatives = strategyNegatives - 100

includeLowAff = False
if(strategyNegatives >= 10):
    includeLowAff = True
    strategyNegatives = strategyNegatives - 10
    
includeDLABneg = False
if(strategyNegatives == 1):
    includeDLABneg = True


saveStrategyNegativesTesting = strategyNegativeTesting
includeLowFnatAsNegExt = False
if(strategyNegativeTesting >= 1000):
    includeLowFnatAsNegExt = True
    strategyNegativeTesting = strategyNegativeTesting - 1000

includeNonBindersExt = False
if(strategyNegativeTesting >= 100):
    includeNonBindersExt = True
    strategyNegativeTesting = strategyNegativeTesting - 100

includeLowAffExt = False
if(strategyNegativeTesting >= 10):
    includeLowAffExt = True
    strategyNegativeTesting = strategyNegativeTesting - 10
    
includeDLABnegExt = False
if(strategyNegativeTesting == 1):
    includeDLABnegExt = True



saveCondition = condition
useAAchem = False
if(condition > 10):
    useAAchem = True
    condition = condition - 10

if removeFirstDimention:
    encodingDimension = 20
else:
    encodingDimension = 21
    
if useAAchem == True:
    encodingDimension = encodingDimension - 16
    print("ERR: the chemical encoding (4 characters, c n p r) is not yet implemented")
    
    

# ========= 1: scripts for reading data preprocessing , performed in advance ===========
# Since the raw poses take too much space (4 TB), it is not possible to regenerate a train/test dataset from them anymore.
# This is done by a separate script, that takes reads raw poses and returns only a number of poses filtered by the options:
#                                     Antigen nPosesPerLabelType   thresFnatBind  thresFnatNegative  PosesPerCDRH3
# example: python PreprocessScript.py 4R9Y_D  10000                0.9            0.1                50
# Of note, this script selects each label type (P, L, N, O, I, F but not ?), with 10 000 of each label. 
# Therefore, the antigen is processed once, and one can later filter which types of labels one consider.
# However, for different thresholds, one need to geenrate new files
# Raw input file expected:
# "Poses_" + antigenID + "_PosesSourceSlices" + antigenID + ".txt_1000_6.txt")
# Output filtered files (10 of them)
# "PrePro" + antigenID + "_" + str(nPairsTot) + "_" + str(fnatBind) + "_" + str(fnatDLABNegative) + "_" + str(groupSize) + "v" + str(repeat) + ".txt"

# For info, meaning of labels that are already annotated to the sequences
# L: Low affinity => includeLowAff
# N: Non binders => includeNonBinders
# O: binding to other ones
# F: FALSE negatives (binding to other ones but actually binding this one as well) 
# P: Positive binding pose (always taken)
# I: Incorrect pose (negative < fnat limit) => includeDLABneg
# ?: Pose that is not taken (neither positive nor negative)
#    O + F + I are the "DLAB negatives". They don't know a sequence is F, in DLAB it belongs to O

def preProFileToOpen(antigenID = "1ADQ_A", nPairsTot = 10000, repeat = 1): #PrePro1ADQ_A_10000_0.9_0.1_50v1
    # We have generated the files for 10 000 poses and less. For more, the exact file should be around
    if(nPairsTot < 10000):
        nPairsTot = 10000 #remember this is a local variable
    return(LocationDataFiles + "PrePro" + antigenID + "_" + str(nPairsTot) + "_" + str(fnatBind) + "_" + str(fnatDLABNegative) + "_" + str(groupSize) + "v" + str(repeat) + ".txt")

#repeat from 1 to 9
def openPreprocessedPosesOneAntigen(antigenID, repeat=1, includeNonBinders = True, includeLowAff = True, includeDLABneg = True, 
                        includeLowFnatAsNeg = True, thresholdFnatBinder = 0.9, thresholdFnatNegative = 0.1):
    
    print("Opening", preProFileToOpen(antigenID, repeat))
    prepro = open(preProFileToOpen(antigenID, repeat), newline = '')   #one line is a text with \t and \n                                                                              
    data_reader = csv.reader(prepro, delimiter='\t') #transform lines into lists 
    #    IDs    antigenLattice    antibodyLattice    labels    fnats    gid
    # plus le numero de ligne, damned
    
    IDs = []
    antigenLattice = []
    antibodyLattice = []
    labels = [] # P = Positive /L = Low affinity /N = Nonbinder /O = Other AGs (DLAB negative) - only P is positive.
    fnats = []  # FNAT score, or -1 if non-binder.
    
    for line in data_reader:
        
        takeThisLine = True
        
        if(data_reader.line_num == 1):
            takeThisLine = False #headers
        
        currentLabel = line[4]
        
        
        if(currentLabel == 'N' and includeNonBinders == False):
            takeThisLine = False
        if(currentLabel == 'L' and includeLowAff == False):
            takeThisLine = False
        if(currentLabel == 'I' and includeLowFnatAsNeg == False):
            takeThisLine = False
        if(currentLabel == 'F' and includeDLABneg == False):
            takeThisLine = False
        if(currentLabel == 'O' and includeDLABneg == False):
            takeThisLine = False

        if(takeThisLine):
            IDs.append(line[1])
            antigenLattice.append(line[2])
            antibodyLattice.append(line[3])
            labels.append(currentLabel)
            fnats.append(line[5])
        
    return([IDs, antigenLattice, antibodyLattice, labels, fnats])


    

# ========= 2: scripts for one-hot encoding of lattices ===========    
    
# The first dimension will represent '_', that is an empty position. We will encode it as zero.
# we will then remove this dimension because dummy (if all other ones are zero, this one is 1), so it remains only the 20 AAs
# see removeFirstDimention = True in manually set options

#Defining the one hot encoding
alphabet = '_ACDEFGHIKLMNPQRSTVWY'

if useAAchem:
    alphabet = '_cpnr'

# define a mapping of chars to integers
char_to_int = dict((c, i) for i, c in enumerate(alphabet))
int_to_char = dict((i, c) for i, c in enumerate(alphabet))

chemCode = {
    'D': 'c',
    'E': 'c',
    'H': 'c',
    'R': 'c',
    'K': 'c',
    'S': 'p',
    'T': 'p',
    'N': 'p',
    'Q': 'p',
    'G': 'n',
    'A': 'n',
    'V': 'n',
    'L': 'n',
    'I': 'n',
    'P': 'n',
    'M': 'n',
    'C': 'n',
    'F': 'r',
    'W': 'r',
    'Y': 'r',
    'B': '!',
    'J': '!',
    'O': '!',
    'U': '!',
    'X': '!',
    'Z': '!',
    '_': '_'
}

def stringToAAChemProp(s):
    return "".join(chemCode[l] for l in s)

def hotEncodingAAString(myString):
    onehot_encoded = list()
    integer_encoded = [char_to_int[char] for char in myString]
    for value in integer_encoded:
        letter = [0 for _ in range(len(alphabet))]
        letter[value] = 1        
        if(removeFirstDimention == True):
            letter = letter[1:len(alphabet)]
        
        onehot_encoded.append(letter)

    return onehot_encoded  #Note, in previous versions I was returning [], which has shape [1, N, Vocab], here without 1

#input: [0,0,...,1,0,0]
def retrieveAA(onehot_encodedAA):
    foundAA = '_'
    for i in range(len(onehot_encodedAA)):
        if(onehot_encodedAA[i] == 1):
            if(removeFirstDimention == True):
                foundAA = int_to_char[i+1]
            else:
                foundAA = int_to_char[i]
    return foundAA

#input: [[0,0,...,1,0,0] , ... , [0,1,0,0... 0]]
def retrieveString(onehot_encodedString):
    foundString = ""
    for k in range(len(onehot_encodedString)):
        foundString = foundString + retrieveAA(onehot_encodedString[k])
    return foundString

# Converts the text description of a lattice into a 2D list of string blocks, that can be directly one-hot encoded. 
def textTo2Dtext(lattice):
    lineX = lattice.split(',')[0:-1]    #the last character is a coma, so there is an empty block after
    res = [X.strip().split(' ') for X in lineX]
    print(res)
    
# one-hot a 1D of strings
# HERE we are going to change into chemical based on the global variable useAAchem
def batchOneHot(vecOfStrings):
    if useAAchem:
        return [hotEncodingAAString(stringToAAChemProp(s)) for s in vecOfStrings]
    else:
        return [hotEncodingAAString(s) for s in vecOfStrings]

# one-hot a 2D of strings
def textToTensor(lattice):
    lineX = lattice.split(',')[0:-1]  #the last character is a coma, so there is an empty block after
    return [batchOneHot(X.strip().split(' ')) for X in lineX]

#examples
if(not(runningInCommandLine)):
    test = hotEncodingAAString("A_CD")
    print(test, np.array(test).shape)
    
    found = retrieveAA(test[0])
    print(found)
    
    got = retrieveString(test)
    print(got)
    
    example = "______ ______ ______ ______ ______ ______ ,______ H_____ ______ ______ ______ ______ ,______ ______ __W___ __Y___ ______ ______ ,______ ______ __L___ __D___ ______ ______ ,______ ______ ______ ______ ______ ______ ,______ ______ ______ ______ ______ ______ ,"
    v1 = textTo2Dtext(example)
    print(v1, np.array(v1).shape)
    
    v2 = batchOneHot(['__P___', '______', '______', '______', '______', '______'])
    print(v2, np.array(v2).shape)
    
    v3 = textToTensor(example)
    print(np.array(v3).shape)    



# ========= 3 Functions for data enhancement by taking a 6x6x6x, lattice and rotating "randomly". ==========
PI = 3.1415927535

# This function returns a list of [Phi, theta] where phi goes -Pi to Pi and 
# theta goes -Pi/2 to Pi/2 and starts from the XY plane (the Z axis is theta = pi/2)
def getDirections(resolution):
    res = []
    if (resolution < 1 ):
        return(res)

    # choice of latitude
    for n in range(1, resolution):
        thetaJ = PI * n / resolution - PI / 2.0 ;
    
        #choice of points
        nj = math.floor( 0.5 + math.sqrt( 3 ) * resolution * math.cos(thetaJ))

        shift = 0
        if (n % 2 ) != 0 : 
            shift += 0.5 ;
            for j in range(0, nj):
        
                # their phi is from equator (X,Y plane)
                res.append([-PI + 2.0 * PI * ( j + shift) / nj, thetaJ])
        
    return res;
    
# Generates a pool of uniform directions that we are going as a library to pick randomly
possibleRotations = getDirections(1000)
random.shuffle(possibleRotations)
totalRotations = len(possibleRotations)
print("Generated ", totalRotations, "possible uniformly distributed rotations")

def getRandomRotation():
    return possibleRotations[random.randint(0, totalRotations-1)]

# Example
# To get the idea of how many uniformly distributed directions we have depending on the resolution
if(not(runningInCommandLine)):
    for i in range(1, 51, 5):
        print(i, len(getDirections(i)))        
    print(getRandomRotation())

# Rotation function for [6x6x6] lattices
# First centers to (2.5, 2.5, 2.5) then rotates according to angles phi and theta, then translate backs 
# Theta is the angle that preserves the Y axis (latitude, -90 to 90), Phi preserves the Z azis (longitude, 0 to 360).
     # if rotates by theta = 90 deg (pi / 4) then X1 = -Z and Z1 = X, 
     # therefore, if the vector was (1, 0, 0) it becomes (0, 0, 1): Theta starts from flat (X,Y plane) towards Z:

def rotate(pos3D, theta, phi):
    x = pos3D[0] - 2.5;
    y = pos3D[1] - 2.5;
    z = pos3D[2] - 2.5;

    X1 =  x*math.cos(theta)-z*math.sin(theta);
    Y1 =  y;
    Z1 =  x*math.sin(theta)+z*math.cos(theta);

    X2 =  X1*math.cos(phi)+Y1*math.sin(phi);
    Y2 =  -X1*math.sin(phi)+Y1*math.cos(phi);
    Z2 =  Z1;

    return([X2 + 2.5, Y2 + 2.5, Z2 + 2.5])

# If a rotated position happens to be outside the lattice we bring it back in [we could also have decided to kick it out]
def boxize6(floatNr):
    rounded = round(floatNr)
    if(rounded < 0):
        return(0)
    if(rounded > 5):
        return(5)
    return(rounded)

# Of note, depending on removal of the dummy dimension '_' and on the chemical or AA encoding, the last dimension is 4, 5, 20 or 21
checkdimEncoding = len(alphabet)
startingRelevantDimension = 1
if(removeFirstDimention == True):
    checkdimEncoding = checkdimEncoding - 1
    startingRelevantDimension = 0

if checkdimEncoding != encodingDimension:
    print("ERR: the alphabet and the encoding dimension were different.")
    sys.exit()

# Main function to rotate a lattice 
def rotateLattice(latticeTensor4D, phi_theta = None):
    if(phi_theta == None):
        [phi, theta] = getRandomRotation()
    else:
        [phi, theta] = phi_theta
    
    res = np.zeros(shape = [6,6,6,encodingDimension])
    
    if(not removeFirstDimention):
        for i in range(0,6):
            for j in range(0,6):
                for k in range(0,6):
                    res[i][j][k][0] = 1 # by default
                
    for i in range(0,6):
        for j in range(0,6):
            for k in range(0,6):
                if(sum(latticeTensor4D[i][j][k][startingRelevantDimension:encodingDimension]) > 0):
                    [newX, newY, newZ] = rotate([i,j,k], theta, phi)
                    #print(i, j, k, " -> ", newX, newY, newZ)
                    res[boxize6(newX)][boxize6(newY)][boxize6(newZ)] = latticeTensor4D[i][j][k]
    return(res)


if(not(runningInCommandLine)):
    example = "______ ______ ______ ______ ______ ______ ,______ H_____ ______ ______ ______ ______ ,______ ______ __W___ __Y___ ______ ______ ,______ ______ __L___ __D___ ______ ______ ,______ ______ ______ ______ ______ ______ ,______ ______ ______ ______ ______ ______ ,"
    v = textToTensor(example)
    
    # Example of rotation where X becomes Z and Y stays the same.
    rotateLattice(v, [0, PI/2.])
    
    res = rotateLattice(v)
    #print(sum(res))
    
    
    
    
    
# ============ 4 - Definition of the DLAB architecture ================

class CNN3D(tf.keras.Model):
    def __init__(self, drop_rate = 0.2):
        super(CNN3D, self).__init__()
        
        #self.mp1 = tf.keras.layers.MaxPool3D(pool_size=(2, 2, 2), strides=(2,2,2))
                
        #initializer1 = tf.keras.initializers.GlorotUniform()
        self.conv1 = tf.keras.layers.Conv3D(
            filters = 32,
            kernel_size = 3,
            padding = 'same',
            strides = (1,1,1),
            kernel_initializer='glorot_uniform',
            activation='relu'
        )
        self.bn1 = tf.keras.layers.BatchNormalization()
        
        #self.mp2 = tf.keras.layers.MaxPool3D(pool_size=(2, 2, 2), strides=(2,2,2))

        #initializer2 = tf.keras.initializers.GlorotUniform()
        self.conv2 = tf.keras.layers.Conv3D(
            #in_channels=32, 
            filters = 64, 
            kernel_size=3, 
            padding='valid',
            strides = (1,1,1),
            kernel_initializer='glorot_uniform',
            activation='relu'
        )
        self.bn2 = tf.keras.layers.BatchNormalization()
        
        #self.mp3 = tf.keras.layers.MaxPool3D(pool_size=(2, 2, 2), strides=(2,2,2))
        
        #initializer3 = tf.keras.initializers.GlorotUniform()
        self.conv3 = tf.keras.layers.Conv3D(
            #in_channels=32, 
            filters = 128, 
            kernel_size=3, 
            padding='valid',
            strides = (1,1,1),
            kernel_initializer='glorot_uniform',
            activation='relu'
        )        
        self.bn3 = tf.keras.layers.BatchNormalization()
        

        # self.flat = tf.keras.layers.Flatten()
        
        #self.dropout = tf.keras.layers.Dropout(drop_rate)
        
        #self.fc1 = tf.keras.layers.Dense(1, activation='sigmoid', kernel_initializer='glorot_uniform')


        print("Initialized CNN3D")

    def call(self, x):

        #x = self.mp1(x)
        x = self.conv1(x)
        x = self.bn1(x)

        #x = self.mp2(x)
        x = self.conv2(x)
        x = self.bn2(x)

        #x = self.mp3(x)
        x = self.conv3(x)
        x = self.bn3(x)

        # flatten
        #x = self.flat(x)
         # x = x.view(-1, self.num_flat_features(x))

        #x = self.dropout(x)
        #x = self.fc1(x)

        return x


def createDLABmodel():
    inps1 = tf.keras.layers.Input(shape=(size_lattice, size_lattice, size_lattice, encodingDimension), name = "antibody")
    inps2 = tf.keras.layers.Input(shape=(size_lattice, size_lattice, size_lattice, encodingDimension), name = "antigen")

    CNN1 = CNN3D()
    features1 = CNN1(inps1)

    CNN2 = CNN3D()
    features2 = CNN1(inps2)

    flat1 = tf.keras.layers.Flatten()
    flat2 = tf.keras.layers.Flatten()

    concat = tf.keras.layers.Concatenate() #axis=1)

    dropout = tf.keras.layers.Dropout(0.2)

    fc1 = tf.keras.layers.Dense(1, activation='sigmoid', kernel_initializer='glorot_uniform')

    flatenned1 = flat1(features1)
    flatenned2 = flat2(features2)
    output = concat([flatenned1, flatenned2])

    output = dropout(output)
    output = fc1(output)

    model = tf.keras.models.Model(inputs=[inps1, inps2], outputs=output)
    return(model)    

# tool function to shuffle something and return it    
def internalShuffle(l):
    random.shuffle(l)
    return l

def f1(precision, recall):
    if((precision == 0) or (recall == 0)):
        return(0)
    
    return((2 * precision * recall) / (precision + recall))
# Functions defined!    
    
def printMetrics(evaluation, nEpochs):
    f1_train = f1(evaluation['Precision'][nEpochs-1], evaluation['Recall'][nEpochs-1])
    f1_val = f1(evaluation['val_Precision'][nEpochs-1], evaluation['val_Recall'][nEpochs-1])
    res = [evaluation['loss'][nEpochs-1],  # : [0.07811099290847778, 0.021633705124258995],
    evaluation['accuracy'][nEpochs-1],  #: [0.9834374785423279, 0.9943749904632568],
    evaluation['FalseNegatives'][nEpochs-1],  #: [120.0, 45.0],
    evaluation['FalsePositives'][nEpochs-1],  #: [145.0, 45.0],
    evaluation['TrueNegatives'][nEpochs-1],  #: [14315.0, 14415.0],
    evaluation['TruePositives'][nEpochs-1],  #: [1420.0, 1495.0],
    evaluation['Precision'][nEpochs-1],  #: [0.9073482155799866, 0.9707792401313782],
    evaluation['Recall'][nEpochs-1],  #: [0.9220778942108154, 0.9707792401313782],
    evaluation['AUC'][nEpochs-1],
    f1_train,
    evaluation['val_loss'][nEpochs-1],  #: [10.108085632324219, 6.631629943847656],
    evaluation['val_accuracy'][nEpochs-1],  #: [0.7537500262260437, 0.7742499709129333],
    evaluation['val_FalseNegatives'][nEpochs-1],  #: [983.0, 878.0],
    evaluation['val_FalsePositives'][nEpochs-1],  #: [2.0, 25.0],
    evaluation['val_TrueNegatives'][nEpochs-1],  #: [2538.0, 2515.0],
    evaluation['val_TruePositives'][nEpochs-1],  #: [477.0, 582.0],
    evaluation['val_Precision'][nEpochs-1],  #: [0.9958246350288391, 0.9588138461112976],
    evaluation['val_Recall'][nEpochs-1],  #: [0.32671234011650085, 0.39863014221191406],
    evaluation['val_AUC'][nEpochs-1],
    f1_val]

    return(res)

def getF1TrainVal(histDict):
    L = len(histDict['accuracy'])
    trainF1 = [f1(histDict['Precision'][i], histDict['Recall'][i]) for i in range(0, L)]
    valF1 = [f1(histDict['val_Precision'][i], histDict['val_Recall'][i]) for i in range(0, L)]
    return (trainF1, valF1)
    
    
    
    
# ========== 5 Now we read preprocessed files, make train / test and enhance by rotation ===========

# selecting nAntigens randomly
possibleIDs = [*range(0,len(listAntigens))]
random.shuffle(possibleIDs)
nAGsForTrain = min(len(listAntigens), nAntigens)

selectedAGs = possibleIDs[0:nAGsForTrain]
externalAGs = possibleIDs[nAGsForTrain:len(listAntigens)]
print("Selection of ", nAntigens, "/", len(listAntigens) , "antigens:", selectedAGs)
print("Other antigens are kept as external dataset", externalAGs)
selectedAGnames = [listAntigens[i] for i in selectedAGs]
print("Selected antigens: ", selectedAGnames)


# Summing up all options before we start
print("AllOptions\tstopDataLeakage\tsize_lattice\tbalancingStrategy\tnEpochs\tnAntigens\tnPairsTot\tstrategyNegatives\t...\tnRepeats\tcondition\tuseAAchem\tnRotations\tgroupSize\tfnatBind\tfnatDLABNegative\tselectedAGnames\n")
print("AllOptions\t", str(stopDataLeakage) + "\t" + str(size_lattice) + "\t" + str(balancingStrategy) + "\t" + str(nEpochs) + "\t" + str(nAntigens) + "\t" + str(nPairsTot) + "\t" + str(saveStrategyNegatives) + "\t" + "..." + "\t" + str(nRepeats)+ "\t" + str(saveCondition) + "\t" + str(useAAchem) + "\t" + str(nRotations) + "\t" + str(groupSize) + "\t" + str(fnatBind) + "\t" + str(fnatDLABNegative) + "\t" + str(selectedAGnames) + "\n")





for repeatLoop in range(1, nRepeats+1):


    # Reading the preprocessed file of each antigen, 
    pooledAGname = []
    pooledIDs = []
    pooledantigenLattice = []
    pooledantibodyLattice = []
    pooledlabels = []
    pooledfnats = []
    successedAntigens = []

    # When train and test have different composition of negatives, we will do as following:
    # - take all classes needed in either the train and/or test
    # - for the train, remove the non-wanted classes in train and separate into train and (pre-)test
    # - for the test, take the generated pre-test, remove the unwanted classes and add the possibly wanted
    #   that were not in the train.
    #   This will cause differences in the composition of the train and test. 

    for ag in selectedAGnames:
        try:
            #f = open(preProFileToOpen(antigenID = "1ADQ_A", repeat = repeatLoop))
            [IDs, antigenLattice, antibodyLattice, labels, fnats] = openPreprocessedPosesOneAntigen(ag, repeatLoop, includeNonBinders or includeNonBindersExt, includeLowAff or includeLowAffExt, includeDLABneg or includeDLABnegExt, includeLowFnatAsNeg or includeLowFnatAsNegExt, fnatBind, fnatDLABNegative)

            print("For antigen ", ag, ", got " , len(antigenLattice), " preprocessed antigen poses, strategy ", strategyNegatives)
            print(Counter(labels))

            pooledAGname = pooledAGname + ([ag] * len(IDs))
            pooledIDs = pooledIDs + IDs
            pooledantigenLattice = pooledantigenLattice + antigenLattice
            pooledantibodyLattice = pooledantibodyLattice + antibodyLattice
            pooledlabels = pooledlabels + labels
            pooledfnats = pooledfnats + fnats
            successedAntigens = successedAntigens + [ag]
        except IOError:
            print("ERR: couldn't find preprocessed data for antigen " + str(ag) + ", continuing")

    print("Amount of each elements read: AGname, ID, antigenLattice, antibodyLattice, labels, fnats.")
    print(len(pooledAGname), len(pooledIDs), len(pooledantigenLattice), len(pooledantibodyLattice), len(pooledlabels), len(pooledfnats))
    print("They should have the same value, it represents the number of different poses.")
    print("Note: the poses are already per group of " + str(groupSize) + " pose for a single CDRH3")

    # Now we will create 3 datasets: train (strategyNegatives, 80%), test (strategyNegatives, 20%), 
    #                                  and ext (test filtered for strategyNegativesExt + new negatives only in ext)


    
    
    # Now we will add a new status to each instance, whether it belongs to train, train+ext or only ext
    def selectLabel(label, includeNonBinders, includeLowAff, includeDLABneg, includeLowFnatAsNeg):
        if(includeNonBinders == True and label == 'N'): 
            return(1)
        if(includeLowAff == True and label == 'L'): 
            return(1)
        if(includeDLABneg == True and label == 'O'): 
            return(1)
        if(includeLowFnatAsNeg == True and label == 'I'): 
            return(1)  
        if(label == 'P'):
            return(1)
        return(0)

        
        
    pooledIsInTrain = [selectLabel(label, includeNonBinders, includeLowAff, includeDLABneg, includeLowFnatAsNeg) 
                       for label in pooledlabels]
    pooledIsInExt = [selectLabel(label, includeNonBindersExt, includeLowAffExt, includeDLABnegExt, includeLowFnatAsNegExt) 
                       for label in pooledlabels]
                       
                       
    def shuffleListPositionsPerID(listIDs):
        d = {}
        for i in range(0, len(listIDs)):
            if listIDs[i] in d.keys():
                d[listIDs[i]].append(i)
            else:
                d[listIDs[i]] = [i]

        # We will assume each key has the same amount of elements.
        dicoPoseNames = Counter(listIDs)
        possiblePosesIDs = list(dicoPoseNames.keys())
        random.shuffle(possiblePosesIDs)

        shuffledIDsgroupedPerPose = [ internalShuffle(d[i]) for i in possiblePosesIDs]
        ResultPossibleDataIDs = [item for sublist in shuffledIDsgroupedPerPose for item in sublist]
        return(ResultPossibleDataIDs)
      

      
      
    #Shuffling the data is done independently of separarting the external data

    # Before generating rotations, if stopDataLeakage == True, we group poses from the same CDRH3 and shuffle them per group 
    # it puts all generated rotations in a group that is going to be shuffled as a group
    # i.e. the rotations of the same pose stay together, so they don't end up in train AND test when cutting train = [indices 0... train size], etc

    # Makes dictionary: CDR3ID => list of IDs of the poses for it.
    if stopDataLeakage:
        possibleDataIDs = shuffleListPositionsPerID(pooledIDs)
        
    else:  # If false, We shuffle the rotations of the same pose, so they might be in train and test
        possibleDataIDs = [*range(0,len(pooledIDs))]
        random.shuffle(possibleDataIDs)

        
        
        
    #Now that the possibleDataIDs are shuffled, we need to 'Mask' those not for the training
    # complex written: possibleDataIDs = [possibleDataIDs[i] for i in range(0,len(possibleDataIDs)) if pooledIsInTrain[possibleDataIDs[i]] == 1]
    possibleDataIDs = [pos for pos in possibleDataIDs if pooledIsInTrain[pos] == 1]
    allDataExt = [pos for pos in possibleDataIDs if pooledIsInExt[pos] == 1]
    
    
    
    # selection of nPairsTot poses, and creation of nRotations copies of each pose. The first copy is not rotated.
    # Note, to make sure the same group of poses is not in train and test we will leave a gap of groupSize between them
    lenAvailable = len(possibleDataIDs)
    nPairsTot = min(nPairsTot,lenAvailable-groupSize)
    train_size = int(0.8 * nPairsTot)         # we already know how much will be the train and test
    test_size = int(0.2 * nPairsTot) 
    print("Wished dataset sizes:")
    print("Will keep only ", nPairsTot, " out of ", lenAvailable, "Available poses")
    print("BEFORE data enhancement per rotation, expect ", train_size , " training and ", test_size, "validation poses")


    # Converting all labels into 'P' => 1, all other ones => 0
    def binary(label):
        if(label == 'P'):
            return(1)

        return(0)

    # Block for just for visualizing the data and stopping if only one class. (nothing from this block will be used)
    pooledlabelsAll = [pooledlabels[i] for i in possibleDataIDs[0:train_size + test_size + groupSize]]
    print("Classes in training:", Counter(pooledlabelsAll[0:train_size]))
    print("Classes in validation:", Counter(pooledlabelsAll[train_size + groupSize: train_size + test_size + groupSize]))
    encodedLabels = np.array(list(map(binary, pooledlabelsAll)))
    print("Labels in training:", Counter(encodedLabels[0:train_size]))
    print("Labels in validation:", Counter(encodedLabels[train_size + groupSize: train_size + test_size + groupSize]))
    if(len(np.unique(encodedLabels[0:train_size])) < 2):
        print("ERR: The splitting of train and validation has raised only one class in the training")
        
    if(len(np.unique(encodedLabels[train_size + groupSize: train_size + test_size + groupSize])) < 2):
        print("ERR: The splitting of train and validation has raised only one class in the validation")

    # Now the "external test" will be: validation + (classes not in train) - (class in train but not in test) 
    
    
    
    
    
    
    # Now cutting the train and test according to block-shuffled data. We are still encoded as text, as in the input file
    PosElementsTrain = possibleDataIDs[0:train_size]
    PosElementsTest = possibleDataIDs[train_size + groupSize: train_size + test_size + groupSize]

    pooledAGnameShufTrain = [pooledAGname[i] for i in PosElementsTrain]
    pooledIDsShufTrain = [pooledIDs[i] for i in PosElementsTrain]
    pooledantigenLatticeShufTrain = [pooledantigenLattice[i] for i in PosElementsTrain]
    pooledantibodyLatticeShufTrain = [pooledantibodyLattice[i] for i in PosElementsTrain]
    pooledlabelsShufTrain = [pooledlabels[i] for i in PosElementsTrain]
    pooledfnatsShufTrain = [pooledfnats[i] for i in PosElementsTrain]

    pooledAGnameShufTest = [pooledAGname[i] for i in PosElementsTest]
    pooledIDsShufTest = [pooledIDs[i] for i in PosElementsTest]
    pooledantigenLatticeShufTest = [pooledantigenLattice[i] for i in PosElementsTest]
    pooledantibodyLatticeShufTest = [pooledantibodyLattice[i] for i in PosElementsTest]
    pooledlabelsShufTest = [pooledlabels[i] for i in PosElementsTest]
    pooledfnatsShufTest = [pooledfnats[i] for i in PosElementsTest]

    encodedLabelsTrain = np.array(list(map(binary, pooledlabelsShufTrain)))
    encodedLabelsTest = np.array(list(map(binary, pooledlabelsShufTest)))

    print("Train labels", Counter(encodedLabelsTrain))
    print("Test labels", Counter(encodedLabelsTest))


    # Data enhancement: repeating the positive poses as to have as many positive and negative
    if balancingStrategy == 2:
        # Now we look at the positions into the training only
        #List positive elements in the training
        listIDPos = [i for i in range(0,train_size) if encodedLabelsTrain[i] == 1]

        #List negative elements in the training
        listIDNeg = [i for i in range(0,train_size) if encodedLabelsTrain[i] == 0]

        Npos = len(listIDPos)
        Npos = len(listIDPos)
        Nneg = len(listIDNeg)
        NPosToAdd = Nneg - Npos
        SelectedRepeatedPositions = random.choices(listIDPos, k=NPosToAdd)

        # Now data enhancement
        pooledAGnameShufTrain += [pooledAGnameShufTrain[i] for i in SelectedRepeatedPositions]
        pooledIDsShufTrain += [pooledIDsShufTrain[i] for i in SelectedRepeatedPositions]
        pooledantigenLatticeShufTrain += [pooledantigenLatticeShufTrain[i] for i in SelectedRepeatedPositions]
        pooledantibodyLatticeShufTrain += [pooledantibodyLatticeShufTrain[i] for i in SelectedRepeatedPositions]
        pooledlabelsShufTrain += [pooledlabelsShufTrain[i] for i in SelectedRepeatedPositions]
        pooledfnatsShufTrain += [pooledfnatsShufTrain[i] for i in SelectedRepeatedPositions]

        encodedLabelsTrain = np.array(list(map(binary, pooledlabelsShufTrain)))

        print("The training data has been balanced (balancingStrategy == 2; not the test), new label distribution:")
        print("Train", Counter(encodedLabelsTrain))
        print("Test", Counter(encodedLabelsTest))    



    # Now the datasets are prepared, encoded labels are 0/1 but lattices are still as text. The generator function will create lattices as 3D tensors with AA or chem features
    # 1 one-hot, 3, shuffled, 11 = Chemical, 13 = Chemical shuffled
    if(condition == 3 or condition == 13):
        random.shuffle(pooledlabelsShufTrain)

    
    def myBalancedTrainGenerator():
        N = len(pooledantigenLatticeShufTrain)

        for k in range(0,nRotations):
            randomizedIDs = [*range(0,N)]
            random.shuffle(randomizedIDs)

            print("Train data pack", k, "/", nRotations, " ", k*N, "/", nRotations * N, " Processed")

            for i in range(0,N):

                if(i % 10000 == 0):
                    print("+10000 done")

                antigenLat = pooledantigenLatticeShufTrain[randomizedIDs[i]]
                antibodyLat = pooledantibodyLatticeShufTrain[randomizedIDs[i]]
                label = encodedLabelsTrain[randomizedIDs[i]]

                if((i + k) % nRotations != 0):  # one pose is not rotated, shift by 1 each block of data 
                    [phi, theta] = getRandomRotation()  # Same rotation for antigen and antibody
                    hotEncodedAntigen = rotateLattice(textToTensor(antigenLat), [phi, theta])
                    hotEncodedAntibody = rotateLattice(textToTensor(antibodyLat), [phi, theta])
                else: 
                    hotEncodedAntigen = textToTensor(antigenLat)
                    hotEncodedAntibody = textToTensor(antibodyLat)

                yield((hotEncodedAntigen, hotEncodedAntibody), label)

    class callableGenerator:
        def __init__(self):
            print("Init")

        def __call__(self):
            return myBalancedTrainGenerator()

    trainingDataset = tf.data.Dataset.from_generator(callableGenerator(), 
        output_types=((tf.float32,tf.float32), tf.float32), 
        output_shapes=(([6,6,6,encodingDimension], [6,6,6,encodingDimension]), [])).batch(batch_size)            

    def myBalancedTestGenerator():
        N = len(pooledantigenLatticeShufTest)

        for k in range(0,nRotations):

            randomizedIDs = [*range(0,N)]
            random.shuffle(randomizedIDs)

            print("Test data pack", k, "/", nRotations, " ", k*N, "/", nRotations * N, " Processed")

            for i in range(0,N):

                if(i % 10000 == 0):
                    print("+10000 done")

                antigenLat = pooledantigenLatticeShufTest[randomizedIDs[i]]
                antibodyLat = pooledantibodyLatticeShufTest[randomizedIDs[i]]
                label = encodedLabelsTest[randomizedIDs[i]]

                if((i + k) % nRotations != 0):  # one pose is not rotated, shift by 1 each block of data 
                    [phi, theta] = getRandomRotation()  # Same rotation for antigen and antibody
                    hotEncodedAntigen = rotateLattice(textToTensor(antigenLat), [phi, theta])
                    hotEncodedAntibody = rotateLattice(textToTensor(antibodyLat), [phi, theta])
                else: 
                    hotEncodedAntigen = textToTensor(antigenLat)
                    hotEncodedAntibody = textToTensor(antibodyLat)

                yield((hotEncodedAntigen, hotEncodedAntibody), label)

    class callableGeneratorTest:
        def __init__(self):
            print("Init")

        def __call__(self):
            return myBalancedTestGenerator()

    testDataset = tf.data.Dataset.from_generator(callableGeneratorTest(), 
        output_types=((tf.float32,tf.float32), tf.float32), 
        output_shapes=(([6,6,6,encodingDimension], [6,6,6,encodingDimension]), [])).batch(batch_size)




     # Instantiation of the ML architecture
    test = createDLABmodel()
    print(test.summary())

    optimizer = tf.keras.optimizers.Adam(learning_rate=0.01)
    loss = 'binary_crossentropy'
    metrics = ['accuracy', 'FalseNegatives', 'FalsePositives', 'TrueNegatives', 'TruePositives', 
               'Precision', 'Recall', 'AUC'] #, average_precision_score] #, tf.compat.v1.metrics.average_precision_at_k]
    test.compile(loss=loss, optimizer=optimizer, metrics=metrics)


    
    
    
    #createError
    # Now, fit with or without compensated loss for imbalance (strategy 1), or without weights (0 or 2, in 2 the data is pre-balanced)
    if balancingStrategy == 1:
        from sklearn.utils import class_weight
        class_weights = class_weight.compute_class_weight('balanced',np.unique(trainLabels), trainLabels)
        print("The loss function will be reweighted with weights", str(class_weights))

        history = test.fit(trainingDataset,validation_data=testDataset,
                            class_weight=class_weights,
                            epochs=nEpochs,
                            verbose = 1)
    elif balancingStrategy != 1: # for option 2, it was already balanced as dataset processing
        history = test.fit(trainingDataset,validation_data=testDataset, 
                            epochs=nEpochs,
                            verbose = 1)
        
            
    #evaluation = test.evaluate(({'antibody': testKeysAntigen, 'antigen': testKeysAntibody}, [testLabels]))
    print(history.history)
    f1hist = getF1TrainVal(history.history)
    # Writes the output in one line with all information
    file_object = open('HistoryDLABTrainVal.txt', 'a')
    file_object.write("successedAntigens\tnSuccessAg\tstopDataLeakage\t" + "size_lattice\t" + "balancingStrategy\t" + "nEpochs\t" + "nAntigens\t" + "nPairsTot\t" + "strategyNegatives\t" + "repeat\t" + "nRepeats" + "\t" + "condition\t" + "useAAchem\t" + "nRotations\t" + "groupSize\t" + "fnatBind\t" + "fnatDLABNegative\t" + "batch_size" + "\t" + "test.evaluate(([testKeysAntigen,testKeysAntibody], [testLabels]))" + "\n")
    file_object.write(str(len(successedAntigens)) + "\t" + str(stopDataLeakage) + 
                      "\t" + str(size_lattice) + "\t" + str(balancingStrategy) + "\t" + str(nEpochs) + "\t" + str(nAntigens) + 
                      "\t" + str(nPairsTot) + "\t" + str(saveStrategyNegatives) + "\t" + str(saveStrategyNegativesTesting) + "\t" + str(repeatLoop) + "\t" + str(nRepeats)+
                      "\t" + str(condition) + "\t" + str(useAAchem) + "\t" + str(nRotations) + "\t" + str(groupSize) + 
                      "\t" + str(fnatBind) + "\t" + str(fnatDLABNegative) + '\t' + str(batch_size) + 
                      '\t' + str(printMetrics(history.history, nEpochs)) + "\t" + str(selectedAGnames) + "\t" + str(f1hist) + 
                      "\t" + str(history.history['AUC']) + "\t" + str(history.history['val_AUC']) +               
                      "\t" + str(Counter(pooledlabelsAll[0:train_size])) + str(Counter(pooledlabelsAll[train_size + groupSize: train_size + test_size + groupSize])) + "\t" + str(successedAntigens) + "\n")
    file_object.close()    
            

            
            
    #PossibleDataID represents all positions that can be taken in the train/val. We have already taken

    #Elements from the test set that we reuse in the external
    PosInTestCanBeRetaken = [i for i in possibleDataIDs[train_size + groupSize: train_size + test_size + groupSize] if pooledIsInExt[i] == True]

    # Now, the instances that could have belonged in train/val but not taken
    PosInSharedClassesNotTakenTrainTest = [i for i in possibleDataIDs[train_size + test_size + groupSize:len(possibleDataIDs)] if (pooledIsInExt[i] == True and pooledIsInTrain[i] == True)]

    # Now, take all other positions that are in Ext but not in train (basically the positions that are inot in possibleDataIDs)
    # Note, this is not shuffled yet
    PosInExtButNotTrainTest = [i for i in range(len(pooledAGname)) if (pooledIsInExt[i] == True and pooledIsInTrain[i] == False)] 

    
    
    
    # This is only for printing
    pooledAGnameShufTestAndExt = [pooledAGname[i] for i in PosInTestCanBeRetaken]
    pooledIDsShufTestAndExt = [pooledIDs[i] for i in PosInTestCanBeRetaken]
    pooledantigenLatticeShufTestAndExt = [pooledantigenLattice[i] for i in PosInTestCanBeRetaken]
    pooledantibodyLatticeShufTestAndExt = [pooledantibodyLattice[i] for i in PosInTestCanBeRetaken]
    pooledlabelsShufTestAndExt = [pooledlabels[i] for i in PosInTestCanBeRetaken]
    pooledfnatsShufTestAndExt = [pooledfnats[i] for i in PosInTestCanBeRetaken]
    encodedLabelsTestAndExt = np.array(list(map(binary, pooledlabelsShufTestAndExt)))
    print("Labels in training+Ext:", Counter(pooledlabelsShufTestAndExt))
    print("Labels in training+Ext:", Counter(encodedLabelsTestAndExt))

    pooledAGnameShufTestAndExtAdd = [pooledAGname[i] for i in PosInSharedClassesNotTakenTrainTest]
    pooledIDsShufTestAndExtAdd = [pooledIDs[i] for i in PosInSharedClassesNotTakenTrainTest]
    pooledantigenLatticeShufTestAndExtAdd = [pooledantigenLattice[i] for i in PosInSharedClassesNotTakenTrainTest]
    pooledantibodyLatticeShufTestAndExtAdd = [pooledantibodyLattice[i] for i in PosInSharedClassesNotTakenTrainTest]
    pooledlabelsShufTestAndExtAdd = [pooledlabels[i] for i in PosInSharedClassesNotTakenTrainTest]
    pooledfnatsShufTestAndExtAdd = [pooledfnats[i] for i in PosInSharedClassesNotTakenTrainTest]
    encodedLabelsTestAndExtAdd = np.array(list(map(binary, pooledlabelsShufTestAndExtAdd)))
    print("Labels not yet taken by train+val:", Counter(pooledlabelsShufTestAndExtAdd))
    print("Labels not yet taken by train+val:", Counter(encodedLabelsTestAndExtAdd))
    # Now, only the classes that are new (not in the train)

    # Classes that were not in the train/eval but are desired in the external test
    # We will first take ALL data that has been ignored, then we are going to shuffle by blocks (IDs) 
    # Then we just take a 'fair' amount.
    # Anyways, the external is going to be quantified separately on each label separately, so one can recalculate metrics later.

    pooledAGnameShufExtOnly = [pooledAGname[i] for i in PosInExtButNotTrainTest]
    pooledIDsShufExtOnly = [pooledIDs[i] for i in PosInExtButNotTrainTest]
    pooledantigenLatticeShufExtOnly = [pooledantigenLattice[i] for i in PosInExtButNotTrainTest]
    pooledantibodyLatticeShufExtOnly = [pooledantibodyLattice[i] for i in PosInExtButNotTrainTest]
    pooledlabelsShufExtOnly = [pooledlabels[i] for i in PosInExtButNotTrainTest]
    pooledfnatsShufExtOnly = [pooledfnats[i] for i in PosInExtButNotTrainTest]
    encodedLabelsExtOnly = np.array(list(map(binary, pooledlabelsShufExtOnly)))
    print("Labels available in Ext only:", Counter(pooledlabelsShufExtOnly))
    print("Labels available in Ext only:", Counter(encodedLabelsExtOnly))
    
    
    
    
    PosForExternal = PosInTestCanBeRetaken + PosInSharedClassesNotTakenTrainTest + PosInExtButNotTrainTest 

    pooledAGnameExt = [pooledAGname[i] for i in PosForExternal]
    pooledIDsExt = [pooledIDs[i] for i in PosForExternal]
    pooledantigenLatticeExt = [pooledantigenLattice[i] for i in PosForExternal]
    pooledantibodyLatticeExt = [pooledantibodyLattice[i] for i in PosForExternal]
    pooledlabelsExt = [pooledlabels[i] for i in PosForExternal]
    pooledfnatsExt = [pooledfnats[i] for i in PosForExternal]
    encodedLabelsExt = np.array(list(map(binary, pooledlabelsExt)))
    print("Labels in Ext:", Counter(pooledlabelsExt))
    print("Labels in Ext:", Counter(encodedLabelsExt))

    shuffledPositionsPerGroup = shuffleListPositionsPerID(pooledIDsExt)

    pooledAGnameExt = [pooledAGnameExt[i] for i in shuffledPositionsPerGroup]
    pooledIDsExt = [pooledIDsExt[i] for i in shuffledPositionsPerGroup]
    pooledantigenLatticeExt = [pooledantigenLatticeExt[i] for i in shuffledPositionsPerGroup]
    pooledantibodyLatticeExt = [pooledantibodyLatticeExt[i] for i in shuffledPositionsPerGroup]
    pooledlabelsExt = [pooledlabelsExt[i] for i in shuffledPositionsPerGroup]
    pooledfnatsExt = [pooledfnatsExt[i] for i in shuffledPositionsPerGroup]
    encodedLabelsExt = np.array(list(map(binary, pooledlabelsExt)))
    print("Labels in Ext after shuffling:", Counter(pooledlabelsExt))
    print("Labels in Ext after shuffling:", Counter(encodedLabelsExt))
    
    
    
        
        # Now we are going to pick the same amount of each class
    nAvailable = len(pooledAGnameExt)
    listIDsO = [i for i in range(0, nAvailable) if pooledlabelsExt[i] == 'O']
    listIDsF = [i for i in range(0, nAvailable) if pooledlabelsExt[i] == 'F']
    listIDsI = [i for i in range(0, nAvailable) if pooledlabelsExt[i] == 'I']
    listIDsN = [i for i in range(0, nAvailable) if pooledlabelsExt[i] == 'N']
    listIDsL = [i for i in range(0, nAvailable) if pooledlabelsExt[i] == 'L']
    listIDsP = [i for i in range(0, nAvailable) if pooledlabelsExt[i] == 'P']


    # We just pick test_size of each class
    nPerClass = test_size
    balancedPositionsExt = []
    if(len(listIDsO)>0):
        balancedPositionsExt = balancedPositionsExt + list(np.random.choice(listIDsO, size=nPerClass, replace=True))
    if(len(listIDsF)>0):
        balancedPositionsExt = balancedPositionsExt + list(np.random.choice(listIDsF, size=nPerClass, replace=True))
    if(len(listIDsI)>0):
        balancedPositionsExt = balancedPositionsExt + list(np.random.choice(listIDsI, size=nPerClass, replace=True))
    if(len(listIDsN)>0):
        balancedPositionsExt = balancedPositionsExt + list(np.random.choice(listIDsN, size=nPerClass, replace=True))
    if(len(listIDsL)>0):
        balancedPositionsExt = balancedPositionsExt + list(np.random.choice(listIDsL, size=nPerClass, replace=True))
    if(len(listIDsP)>0):
        balancedPositionsExt = balancedPositionsExt + list(np.random.choice(listIDsP, size=nPerClass, replace=True))


        
    pooledAGnameExtBalanced = [pooledAGnameExt[i] for i in balancedPositionsExt]
    pooledIDsExtBalanced = [pooledIDsExt[i] for i in balancedPositionsExt]
    pooledantigenLatticeExtBalanced = [pooledantigenLatticeExt[i] for i in balancedPositionsExt]
    pooledantibodyLatticeExtBalanced = [pooledantibodyLatticeExt[i] for i in balancedPositionsExt]
    pooledlabelsExtBalanced = [pooledlabelsExt[i] for i in balancedPositionsExt]
    pooledfnatsExtBalanced = [pooledfnatsExt[i] for i in balancedPositionsExt]
    encodedLabelsExtBalanced = np.array(list(map(binary, pooledlabelsExtBalanced)))
    print("Labels in Ext after shuffling:", Counter(pooledlabelsExtBalanced))
    print("Labels in Ext after shuffling:", Counter(encodedLabelsExtBalanced))    
        
        
        
    def ExtGenerator():
        N = len(pooledantigenLatticeExtBalanced)

        for k in range(0,nRotations):

            randomizedIDs = [*range(0,N)]
            random.shuffle(randomizedIDs)

            print("External data pack", k, "/", nRotations, " ", k*N, "/", nRotations * N, " Processed")

            for i in range(0,N):

                if(i % 10000 == 0):
                    print("+10000 done")

                antigenLat = pooledantigenLatticeExtBalanced[randomizedIDs[i]]
                antibodyLat = pooledantibodyLatticeExtBalanced[randomizedIDs[i]]
                label = encodedLabelsExtBalanced[randomizedIDs[i]]

                if((i + k) % nRotations != 0):  # one pose is not rotated, shift by 1 each block of data 
                    [phi, theta] = getRandomRotation()  # Same rotation for antigen and antibody
                    hotEncodedAntigen = rotateLattice(textToTensor(antigenLat), [phi, theta])
                    hotEncodedAntibody = rotateLattice(textToTensor(antibodyLat), [phi, theta])
                else: 
                    hotEncodedAntigen = textToTensor(antigenLat)
                    hotEncodedAntibody = textToTensor(antibodyLat)

                yield((hotEncodedAntigen, hotEncodedAntibody), label)

    class callableGeneratorExt:
        def __init__(self):
            print("Init")

        def __call__(self):
            return ExtGenerator()

    extDataset = tf.data.Dataset.from_generator(callableGeneratorExt(), 
        output_types=((tf.float32,tf.float32), tf.float32), 
        output_shapes=(([6,6,6,encodingDimension], [6,6,6,encodingDimension]), [])).batch(batch_size)
            
            
            
    history2 = test.evaluate(extDataset)
    
    file_object = open('HistoryDLABExt.txt', 'a')
    file_object.write("successedAntigens\tnSuccessAg\tstopDataLeakage\t" + "size_lattice\t" + "balancingStrategy\t" + "nEpochs\t" + "nAntigens\t" + "nPairsTot\t" + "strategyNegatives\t" + "repeat\t" + "nRepeats" + "\t" + "condition\t" + "useAAchem\t" + "nRotations\t" + "groupSize\t" + "fnatBind\t" + "fnatDLABNegative\t" + "batch_size" + "\t" + "metrics" + "\n")
    file_object.write("Ext" + "\t" + str(len(successedAntigens)) + "\t" + str(stopDataLeakage) + 
                      "\t" + str(size_lattice) + "\t" + str(balancingStrategy) + "\t" + str(nEpochs) + "\t" + str(nAntigens) + 
                      "\t" + str(nPairsTot) + "\t" + str(saveStrategyNegatives) + "\t" + str(saveStrategyNegativesTesting) + "\t" + str(repeatLoop) + "\t" + str(nRepeats)+
                      "\t" + str(condition) + "\t" + str(useAAchem) + "\t" + str(nRotations) + "\t" + str(groupSize) + 
                      "\t" + str(fnatBind) + "\t" + str(fnatDLABNegative) + '\t' + str(batch_size) + 
                      '\t' + str(history2) + "\t" + str(selectedAGnames) +              
                      "\t" + str(Counter(pooledlabelsExtBalanced)) + "\t" +str(successedAntigens) + "\n")
    file_object.close()
