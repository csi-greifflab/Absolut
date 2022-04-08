#!/usr/bin/env python3

#python ThisFile.py antigenID nPerType fnatBind= 0.9 fnatDLABNegative=0.1 groupSize=50
# This script takes the output of ./Absolut poses and filters only for interesting poses, such that the new files are smaller
# input file name: LocationDataFiles + "Poses_" + antigenID + "_PSlices" + antigenID + ".txt_500_6.txt"
# generated output files: "PrePro" + antigenID + "_" + str(nPerType) + "_" + str(thresFnatBinder) + "_" + str(thresFnatNegative) + "_" + str(groupSize) + "v" + str(repeat) + ".txt"

# nPerType says how many poses of each type (P / I / N / L / O) are kept from the original large file
# groupSize says how many poses for the same antibody sequence (ab-ag pair) are kept as a group. (this will be important in train/test because a group is not separated between train and test)


from __future__ import division, print_function, absolute_import

LocationDataFiles = ""
runningInCommandLine = True
listAntigens = ["1ADQ_A", "1FBI_X", "1FNS_A", "1FSK_A", "1H0D_C", "1JPS_T", "1KB5_AB", "1NCB_N", "1NSN_S", 
                "1OAZ_A", "1OB1_C", "1OSP_O", "1PKQ_J", "1QFW_AB", "1RJL_C", "1S78_B", "1TQB_A", "1WEJ_F", 
                "1YJD_C", "1ZTX_E", "2ARJ_RQ", "2B2X_A", "2FD6_AU", "2HFG_R", "2IH3_C", "2JEL_P", "2Q8A_A", 
                "2R29_A", "2R4R_A", "2R56_A", "2UZI_R", "2VXQ_A", "2VXT_I", "2W9E_A", "2WUC_I", "2XQB_A", 
                "2XWT_C", "2YC1_C", "2YPV_A", "2ZCH_P", "3BGF_S", "3BN9_A", "3CVH_ABC", "3DVG_XY", "3EFD_K", 
                "3GI9_C", "3HI6_A", "3JBQ_B", "3KJ4_A", "3KR3_D", "3KS0_J", "3L5X_A", "3L95_X", "3MJ9_A", 
                "3NCY_A", "3NFP_I", "3NH7_A", "3Q3G_E", "3R08_E", "3R1G_B", "3RAJ_A", "3RKD_A", "3RVV_A", 
                "3SKJ_E", "3SQO_A", "3TT1_A", "3U9P_C", "3UBX_A", "3V6O_A", "3VG9_A", "3VRL_C", "3WD5_A", 
                "4AEI_A", "4CAD_C", "4DKE_A", "4H88_A", "4HC1_B", "4HJ0_B", "4I18_R", "4I77_Z", "4K24_A",
                "4K3J_A", "4K9E_C", "4KI5_M", "4KXZ_A", "4LU5_B", "4MXV_B", "4N9G_C", "4NP4_A", "4OII_A", 
                "4OKV_E", "4PP1_A", "4QCI_D", "4QEX_A", "4QNP_A", "4QWW_A", "4R9Y_D", "4RGM_S", "4U1G_A", 
                "4U6V_A", "4WV1_F", "4Y5V_C", "4YPG_C", "4YUE_C", "4ZFG_A", "4ZFO_F", "4ZSO_E", "5B8C_C", 
                "5BVP_I", "5C0N_A", "5C7X_A", "5CZV_A", "5D93_A", "5DFV_A", "5DHV_M", "5DMI_A", "5DO2_B", 
                "5E8D_A", "5E8E_LH", "5E94_G", "5EII_G", "5EPM_C", "5EU7_A", "5EZO_A", "5F3B_C", "5FB8_C", 
                "5H35_C", "5HDQ_A", "5HI4_B", "5IKC_M", "5J13_A", "5JW4_A", "5JZ7_A", "5KN5_C", "5KTE_A", 
                "5L0Q_A", "5LQB_A", "5MES_A", "5T5F_A", "5TH9_A", "5TLJ_X", "5TZ2_C"]

import numpy as np
import csv
import sys
import math
import random
import pandas as pd
from collections import Counter

random.seed(a=None, version=2)

def fileToOpen(antigenID = "1ADQ_A"):
    return(LocationDataFiles + "Poses_" + antigenID + "_PSlices" + antigenID + ".txt_500_6.txt")

antigenID = "5BVP_I"
nPerType = 1000
thresFnatBinder = 0.9
thresFnatNegative = 0.1
groupSize = 50

print(sys.argv)

if(runningInCommandLine):
    if(len(sys.argv) > 1):
        antigenID = sys.argv[1]

    if(len(sys.argv) > 2):
        nPerType = int(sys.argv[2])
        
    if(len(sys.argv) > 3):
        thresFnatBinder = float(sys.argv[3])

    if(len(sys.argv) > 4):
        thresFnatNegative = float(sys.argv[4])

    if(len(sys.argv) > 5):
        groupSize = int(sys.argv[5])


# DLAB nonbinders (binders to another antigen) might actually be binders, we will annotate on the fly when reading the poses
# It's interesting that we could also define DLAB binders based on the binding threshold, but this is not known when we get ZDOCK poses I guess
ThresholdsBinding = {
    "1ADQ_A" : -94.19,
    "1FBI_X" : -89.99,
    "1FNS_A" : -112.48,
    "1FSK_A" : -83.83,
    "1H0D_C" : -82.04,
    "1JPS_T" : -112.51,
    "1KB5_AB" : -109.21,
    "1NCB_N" : -113.37,
    "1NSN_S" : -89.99,
    "1OAZ_A" : -103.2,
    "1OB1_C" : -92.52,
    "1OSP_O" : -102.11,
    "1PKQ_J" : -96.51,
    "1QFW_AB" : -95.39,
    "1RJL_C" : -86.76,
    "1S78_B" : -117.33,
    "1TQB_A" : -106.71,
    "1WEJ_F" : -73.24,
    "1YJD_C" : -88.09,
    "1ZTX_E" : -96.7,
    "2ARJ_RQ" : -105.82,
    "2B2X_A" : -91.71,
    "2FD6_AU" : -114.58,
    "2HFG_R" : -68.76,
    "2IH3_C" : -97.78,
    "2JEL_P" : -71.41,
    "2Q8A_A" : -100.94,
    "2R29_A" : -90.43,
    "2R4R_A" : -110.01,
    "2R56_A" : -91.93,
    "2UZI_R" : -101.68,
    "2VXQ_A" : -87.38,
    "2VXT_I" : -94.52,
    "2W9E_A" : -94.88,
    "2WUC_I" : -95.7,
    "2XQB_A" : -88.99,
    "2XWT_C" : -114.22,
    "2YC1_C" : -91.66,
    "2YPV_A" : -95.18,
    "2ZCH_P" : -97.6,
    "3BGF_S" : -100.99,
    "3BN9_A" : -103.98,
    "3CVH_ABC" : -101.1,
    "3DVG_XY" : -105.79,
    "3EFD_K" : -59.58,
    "3GI9_C" : -125.24,
    "3HI6_A" : -107.21,
    "3JBQ_B" : -107.87,
    "3KJ4_A" : -104,
    "3KR3_D" : -74.4,
    "3KS0_J" : -96.59,
    "3L5X_A" : -85.18,
    "3L95_X" : -104.04,
    "3MJ9_A" : -97.76,
    "3NCY_A" : -138.29,
    "3NFP_I" : -104.87,
    "3NH7_A" : -86.82,
    "3Q3G_E" : -100.8,
    "3R08_E" : -86.5,
    "3R1G_B" : -121.65,
    "3RAJ_A" : -100.3,
    "3RKD_A" : -86.44,
    "3RVV_A" : -101,
    "3SKJ_E" : -99.42,
    "3SQO_A" : -106.73,
    "3TT1_A" : -133.1,
    "3U9P_C" : -112.6,
    "3UBX_A" : -107.35,
    "3V6O_A" : -106.76,
    "3VG9_A" : -126.11,
    "3VRL_C" : -96.46,
    "3WD5_A" : -105.5,
    "4AEI_A" : -74.47,
    "4CAD_C" : -132.847,
    "4DKE_A" : -118.19,
    "4H88_A" : -102.13,
    "4HC1_B" : -101.94,
    "4HJ0_B" : -91.07,
    "4I18_R" : -111.76,
    "4I77_Z" : -94.59,
    "4K24_A" : -84.3,
    "4K3J_A" : -102.7,
    "4K9E_C" : -113.93,
    "4KI5_M" : -99.25,
    "4KXZ_A" : -91.42,
    "4LU5_B" : -93.54,
    "4MXV_B" : -96.07,
    "4N9G_C" : -52.29,
    "4NP4_A" : -122.49,
    "4OII_A" : -118.82,
    "4OKV_E" : -67.08,
    "4PP1_A" : -82.58,
    "4QCI_D" : -99.1,
    "4QEX_A" : -102.65,
    "4QNP_A" : -104.86,
    "4QWW_A" : -117.6,
    "4R9Y_D" : -98.52,
    "4RGM_S" : -93.42,
    "4U1G_A" : -124.09,
    "4U6V_A" : -96.16,
    "4WV1_F" : -90.47,
    "4Y5V_C" : -117.79,
    "4YPG_C" : -88.74,
    "4YUE_C" : -89.98,
    "4ZFG_A" : -98.38,
    "4ZFO_F" : -80.57,
    "4ZSO_E" : -101.34,
    "5B8C_C" : -100.08,
    "5BVP_I" : -104.06,
    "5C0N_A" : -96.67,
    "5C7X_A" : -107.4,
    "5CZV_A" : -104.86,
    "5D93_A" : -105.93,
    "5DFV_A" : -88.85,
    "5DHV_M" : -63.29,
    "5DMI_A" : -87.93,
    "5DO2_B" : -97.05,
    "5E8D_A" : -79.39,
    "5E8E_LH" : -108.02,
    "5E94_G" : -96.31,
    "5EII_G" : -91.11,
    "5EPM_C" : -71.69,
    "5EU7_A" : -102.4,
    "5EZO_A" : -118.97,
    "5F3B_C" : -91.67,
    "5FB8_C" : -88.57,
    "5H35_C" : -116.39,
    "5HDQ_A" : -106.33,
    "5HI4_B" : -97.85,
    "5IKC_M" : -85.13,
    "5J13_A" : -87.35,
    "5JW4_A" : -117.49,
    "5JZ7_A" : -100.46,
    "5KN5_C" : -71.88,
    "5KTE_A" : -130.18,
    "5L0Q_A" : -101.66,
    "5LQB_A" : -100.33,
    "5MES_A" : -106.8,
    "5T5F_A" : -105.02,
    "5TH9_A" : -117.76,
    "5TLJ_X" : -91.28,
    "5TZ2_C" : -87.93
}

# Note, there are two types of DLABNeg, the low fnat and the others.
def openPosesOneAntigen(antigenID, includeNonBinders = True, includeLowAff = True, includeDLABneg = True, 
                        thresholdFnatBinder = 0.9, thresholdFnatNegative = 0.1):
    print(fileToOpen(antigenID))
    balancedDataset = open(fileToOpen(antigenID), newline = '')   #one line is a text with \t and \n                                                                              
    data_reader = csv.reader(balancedDataset, delimiter='\t') #transform lines into lists 
    IDs = []
    antigenLattice = []
    antibodyLattice = []
    labels = [] # P = Positive /L = Low affinity /N = Nonbinder /O = Other AGs (DLAB negative) - only P is positive.
    fnats = []  # FNAT score, or -1 if non-binder.
    
    #cpt = 0
    for line in data_reader:
        takeThisLine = True
        #cpt = cpt + 1
        label = line[1][0]
        
        if(line[0].startswith("#") or line[0].startswith("AGname")):
            takeThisLine = False
        
        elif(label == 'L' and includeLowAff == False):
            takeThisLine = False
        
        elif(label == 'N' and includeNonBinders == False):
            takeThisLine = False
        
        elif(label == 'O' and includeDLABneg == False):
            takeThisLine = False
            
        elif(label == 'O'):	
            if(float(line[6]) <= ThresholdsBinding[antigenID]):
                label = 'F'  #for False negative or false nonbinder
        
        elif(label == 'P'):
            if(float(line[10]) >= thresholdFnatBinder):
                label = 'P'
                # nothing to do
            elif (float(line[10]) <= thresholdFnatNegative):
                label = 'I'   #Incorrect pose
            else:
                takeThisLine = False
                label = '?'  # No mans land, we will anyways not take it
            
        if(takeThisLine == True):
            IDs.append(line[1])
            antigenLattice.append(line[12])
            antibodyLattice.append(line[13])
            labels.append(label) # The first letter of the slice ID is the label
            if(line[1][0] == 'P'):
                fnats.append(line[10]) # note, this is fracfnat, already normalized in [0:1] (we call it fnat)
            else:
                fnats.append(-1)
        
        #if(cpt == 250000):
        #    break
        
    return([IDs, antigenLattice, antibodyLattice, labels, fnats])

# Now loading the files
# Expect the following format:
#AGname	ID_CDR3	CDR3	slice	poseStr	poseCode	bindEnergy	totEnergy	code	fnatScore	fracFnat	centeringNr	antigenLattice	antibodyLattice
#1ADQ_A	P1	HYDYPLCLDYW			120675-SSSLLUSRDD	jWjViDgViVjQgVkVaNkKkVhLiNjYcfhk	-54.21	-63.75	j0039j0041i0042g0044i0044j0045g0046k0046a0048k0050k0070h0071i0074j0081	14	1	1	
#"______ ______ ______ ______ ______ ______ ,______ H_____ ______ ______ ______ ______ ,______ ______ __W___ __Y___ ______ ______ ,______ ______ __L___ __D___ ______ ______ ,______ ______ ______ ______ ______ ______ ,______ ______ ______ ______ ______ ______ ,"
#"______ ______ ______ ______ ______ ______ ,______ _N____ __K___ __W___ ______ ______ ,______ ______ _V_V__ _Q_Y__ __V___ ______ ,______ ______ ___L__ _V_N__ __D___ ______ ,______ ______ ______ ______ ______ ______ ,______ ______ ______ ______ ______ ______ ,"

[IDs, antigenLattice, antibodyLattice, labels, fnats] = openPosesOneAntigen(antigenID, True, True, True, thresFnatBinder, thresFnatNegative)

print(Counter(labels))   #=> returns Counter({'P': 8000, 'I': 32000, 'L': 199840, 'N': 199868, 'O': 359250})

# So, it's not possile to define a panda dataframe from columns, it seems, only by lines
# df = pd.concat(pd.DataFrame([1, 2, 3]) , pd.DataFrame([4, 5, 6]), axis = 0)
# print(df)
# So I first make a 2D numpy array that then will be converted into a dataframe. So much wasted memory and calculation for nothing
dataThisAntigen = np.column_stack((IDs, antigenLattice,antibodyLattice,labels,fnats))

for repeat in range(1,10):
	df = pd.DataFrame(dataThisAntigen, columns=['IDs', 'antigenLattice', 'antibodyLattice', 'labels', 'fnats'])

	df['gid'] = df.groupby(['IDs', 'labels']).ngroup()

	grouped = df.groupby(['IDs', 'labels'])
	grouped = grouped.apply(lambda x: x.sample(groupSize, replace=True))

	selectedGroups = grouped[grouped.gid.isin(random.sample(list(grouped.gid.unique()),nPerType))]

	selectedGroups = selectedGroups.reset_index( drop=True)
	selectedGroups.to_csv("PrePro" + antigenID + "_" + str(nPerType) + "_" + str(thresFnatBinder) + "_" + str(thresFnatNegative) + "_" + str(groupSize) + "v" + str(repeat) + ".txt", sep="\t")

for repeat in range(1,10):
	df = pd.DataFrame(dataThisAntigen, columns=['IDs', 'antigenLattice', 'antibodyLattice', 'labels', 'fnats'])

	df['gid'] = df.groupby(['IDs']).ngroup()

	grouped = df.groupby(['IDs'])
	grouped = grouped.apply(lambda x: x.sample(groupSize, replace=True))

	selectedGroups = grouped[grouped.gid.isin(random.sample(list(grouped.gid.unique()),nPerType))]

	selectedGroups = selectedGroups.reset_index( drop=True)
	selectedGroups.to_csv("PreProSameID" + antigenID + "_" + str(nPerType) + "_" + str(thresFnatBinder) + "_" + str(thresFnatNegative) + "_" + str(groupSize) + "v" + str(repeat) + ".txt", sep="\t")



