#!/usr/bin/env Rscript

# Pools raw binding files to the same antigen from a folder (called with the name of this antigen)

#Usage (for antigen 1ADQ_A) : ./thisFile.R 1ADQ_A
# This script should be called where the raw bindings folder is located, i.e. where 1ADQ_A is located as subfolder
# Output: creates a new folder 1ADQ_AAnalyses/ containing all CDR3 and slices selected

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
    stop("This script needs one argument (name of antigen)", call.=FALSE)
}

antigenID <- args[1]

#Define input and output folders
fold <- paste(antigenID,'/', sep="")
outFolder <- paste(antigenID,"Analyses/",sep="")
dir.create(outFolder)


#List raw binding files in a folder 
getFilesForAntigen <- function(antigenID, folder){
    pattern <- paste("*", antigenID, ".*FinalBindings.*txt$", sep="");
    if(nchar(antigenID) == 0) { #stupid R says length("") = 1
        pattern <- ".*FinalBindings.*txt$"
    }
    list <- list.files(path = folder, pattern = pattern, recursive = TRUE, ignore.case = FALSE, include.dirs = FALSE);
    return(list)
}

#listFiles <- getFilesForAntigen(antigen, fold)

# This function reads one raw bindings file and only keeps the "best" 11-mer per CDRH3 (or more than one if equally best)
getData <- function(folder, fname){
    #print(fname)
    a <- read.table(paste(folder,fname, sep=""), header=TRUE)
    b <- a[which(a$Best == 'true'), ]
   
    # Add an additional column with the name of the Antigen (take the file name)
    NL = length(b$Structure)
    nameCol = rep(fname, NL)
    b <- cbind(b, Antigen=nameCol)
    a <- NULL
    return(b)
}

# Option 1: Takes a folder and reads all raw bindings files inside, keeps only best 11-mers per CDRH3
poolDataOneAntigen <- function(antigenID, folder){
    nrFilesMax = 1000
    filesToProcess <- getFilesForAntigen(antigenID, folder)
    #print(paste("Got ", length(filesToProcess), "files in", folder))
    # creating a larger dataframe with the info of all antigen
    #https://stackoverflow.com/questions/10689055/create-an-empty-data-frame
    df <- data.frame(ID_slide_Variant=character(),
                     CDR3=character(),
                     Best.=factor(),
                     Slide=character(),
                     Energy=double(),
                     Structure=character(),
                     Antigen=character(),
                     stringsAsFactors=TRUE);
    for(i in 1:min(nrFilesMax,length(filesToProcess))){
        fname <- filesToProcess[i]  
        print(paste("   ...",fname))
        fullData <- getData(fold, fname) 
        df <- rbind(df, fullData)
        fullData <- NULL
    }
    
    return(df)
}

# Option 2: Keeps all 11-mers that are above an energy threshold (independently of best=true or false)
poolSlicesOverThreshold <- function(antigenID, folder, threshold, reverse = FALSE){
    nrFilesMax = 1000
    filesToProcess <- getFilesForAntigen(antigenID, folder)
    #print(paste("Got ", length(filesToProcess), "files"))
    df <- data.frame(ID_slide_Variant=character(),
                     CDR3=character(),
                     Best.=factor(),
                     Slide=character(),
                     Energy=double(),
                     Structure=character(),
                     Antigen=character(),
                     stringsAsFactors=TRUE);
   
    
    for(i in 1:min(nrFilesMax,length(filesToProcess))){
        fname <- filesToProcess[i]  
        print(paste("   ...",fname))
        a <- read.table(paste(folder,fname, sep=""), header=TRUE)
        #print(paste("Got", length(a$Energy),"lines before filtering"))
        if(reverse){
            fullData <- a[which(a$Energy > threshold), ]
        } else {
            fullData <- a[which(a$Energy <= threshold), ]
        }
        a <- NULL
        df <- rbind(df, fullData)
        fullData <- NULL
    }
    return(df)
}



# Examples
#bigDataThisAntigen <- poolSlicesOverThreshold("1FBI_X", "C:/Users/pprobert/Desktop/1FBI_Manual/25000_1FBI.txt", -85)
bigDataThisAntigen <- poolDataOneAntigen (antigenID, fold)
