
# This file is generating the data files for ML Task 1 and Task 2
# It needs to be called in the folder where ll the 1FBI_XAnalyses/ folders are
# It browses along all the XXXX_XAnalyses/ folders and creates:
#     - datasets for Task 1 for each antigen (inside their respective Analyses folder)
#     - merged datasets along all antigens for ML Task 2

library(ggplot2)
library(dplyr)

#setwd("D:/pprobert/AnalyzedV2")
#setwd("C:/Users/pprobert/Desktop/TestScripts/")

#List of antigens we did process
toParse <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
listAvailable <- c()
for (rep in toParse){
    #print(substring(rep, nchar(rep) - 7, nchar(rep)))
    if((nchar(rep) > 1) & (substring(rep, nchar(rep) - 7, nchar(rep)) == "Analyses")){
        #print(rep)
        #Remove ./ and Analyses to get back antigenID
        listAvailable <- c(listAvailable, substring(rep, 3, nchar(rep) - 8))
    }
}
listAvailable    

#Now we will have two variants: Either we generate datasets of CDR3s,
# or we generate from slices. The source files will be different:
# - antigenID, "_MascotteFeatures.txt" for the CDR3 datasets, then we output CDR3s and the AA usage of CDR3s
# - antigenID, "_MascotteSlicesFeatures.txt" for the CDR3 datasets, then we output SLICES and the AA usage of SLICES

generateMLTask1 <- function(antigenID, nameLabel, useSlices=TRUE, uniqueSlices=TRUE, binderRef="Mascotte", nonbinderRef="500kNonMascotte", nameOut="Task1"){
    
    #Step 1: read filtered files with either the top 1% CDR3s or (if useSlices=true) all the slices that are above this same threshold
    if(useSlices){
        fileFeaturesBinders = paste(antigenID, "Analyses/", antigenID, "_", binderRef, "SlicesFeatures.txt", sep="")
        fileFeaturesNonBinders = paste(antigenID, "Analyses/", antigenID, "_", nonbinderRef, "SlicesFeatures.txt", sep="")
    } else {
        fileFeaturesBinders = paste(antigenID, "Analyses/", antigenID, "_", binderRef, "Features.txt", sep="")
        fileFeaturesNonBinders = paste(antigenID, "Analyses/", antigenID, "_", nonbinderRef, "Features.txt", sep="")
    }
    print(fileFeaturesBinders)
    print(fileFeaturesNonBinders)
    
    bind <- read.table(fileFeaturesBinders, header=TRUE)
    head(bind)
    nonbind <- read.table(fileFeaturesNonBinders, header=TRUE)
    head(nonbind)
    
    
    print(paste(length(bind$ID_slide_Variant), nameLabel))
    
    if(useSlices){
        # A duplicate Slice should have
        if(uniqueSlices){
            bind <- bind[!duplicated(bind$Slide),]   
            nonbind <- nonbind[!duplicated(nonbind$Slide),]   
        }
    }
    print(paste(length(nonbind$ID_slide_Variant), "NonBinders"))
    
    #merged now contains the full dataset (note, I don't use it)
    #merged <- rbind(cbind(bind, Label=nameColBind), cbind(nonbind, Label=nameColNonBind))
    #AAsSlice <- merged$AAcompoFullSlice
    #AACDR <- merged$AAcompoFullCDR
    
    # Get the distribution of CDR3 size on both binders and nonbinders
    distribCDR <- bind %>% group_by(sizeCDR3) %>% mutate(number=n()) %>% slice(1) %>% select(sizeCDR3, number)
    distribCDRnonbind <- nonbind %>% group_by(sizeCDR3) %>% mutate(number=n()) %>% slice(1) %>% select(sizeCDR3, number)
    
    # Create a new column with either the label or "Nonbinder"
    nameColBind = rep(nameLabel, length(bind$ID_slide_Variant))
    nameColNonBind = rep("NonBinder", length(nonbind$ID_slide_Variant))

    # Now shorten the files with only the columns of interest (Slide or CDR3, with the respective AA composition)
    # We keep track of the CDR3 size to later sample from a wished CDR3 lentgth distribution among non-binders (we have a lot of them)
    if(useSlices){
        AllBinders <- cbind(bind, Label=nameColBind) %>% select(Slide, Label, AAcompoFullSlice, sizeCDR3)
        NonBindersSmaller <- cbind(nonbind, Label=nameColNonBind) %>% select(Slide, Label, AAcompoFullSlice, sizeCDR3)
    } else {
        AllBinders <- cbind(bind, Label=nameColBind) %>% select(CDR3, Label, AAcompoFullCDR, sizeCDR3)
        NonBindersSmaller <- cbind(nonbind, Label=nameColNonBind) %>% select(CDR3, Label, AAcompoFullCDR, sizeCDR3)
    }
    
    #Dataset 1: nonbalanced
    nControl = length(AllBinders$Slide)
    sampledNonBinders <- NonBindersSmaller[sample(length(NonBindersSmaller$Slide), nControl),]
    nonBalancedDataset <- rbind(AllBinders, sampledNonBinders)
    #This step is just shuffling the lines randomly
    shuffledNonBalancedDataset <- nonBalancedDataset[sample(length(nonBalancedDataset$Label), length(nonBalancedDataset$Label)),]
    
    #Dataset 2: balanced per CDR3 size distribution-. We build the control dataset by picking same amount of non-binders for each CDR3 length 
    df <- NULL
    for(i in 1:length(distribCDR$sizeCDR3)){
        # l = possible length of a CDR3
        l = distribCDR$sizeCDR3[i]
        # nb = amount of sequences with this CDR3 length in the binder dataset.
        nb = distribCDR$number[i]
        
        #Note: the selected class already has the good columns
        
        # 1- list of non-binders with this size
        thisSize <- NonBindersSmaller[nonbind$sizeCDR3 == l,]
        nbnonbind = length(thisSize$sizeCDR3)
        
        # 2- selects the good amount (nb) from the nonbinders 
        selected <- thisSize[sample(nrow(thisSize), min(nb, nbnonbind)),]
        nselected = length(selected$sizeCDR3)
        df <- rbind(df, selected)
        
        print(paste("For CDR3 L=", l, "nb binders=", nb, " got nonbinders:", nbnonbind, "and selected", nselected))
    }
    
    # merge the dataset
    balancedDataset <- rbind(AllBinders, df)
    # shuffle the lines
    shuffledBalancedDataset <- balancedDataset[sample(length(balancedDataset$Label), length(balancedDataset$Label)),]
    
    #Dataset 3: Mock shuffling of labels (will make a function to redo it multiple times from a file)
    PreviousLabels <- shuffledBalancedDataset$Label
    Label <- PreviousLabels[sample(length(PreviousLabels), length(PreviousLabels))]
    # This step removes the Label, creates a new Label column, and reselects the columns in the good order
    if(useSlices){
        balancedShuffledLabels <- cbind(shuffledBalancedDataset %>% select(Slide, AAcompoFullSlice, sizeCDR3), Label) %>% select(Slide, Label, AAcompoFullSlice, sizeCDR3) 
    } else {
        balancedShuffledLabels <- cbind(shuffledBalancedDataset %>% select(CDR3, AAcompoFullCDR, sizeCDR3), Label) %>% select(CDR3, Label, AAcompoFullCDR, sizeCDR3) 
    }
    
    if(useSlices){
        write.table(shuffledNonBalancedDataset, paste(antigenID, "_", nameOut, "_SlicesImbalancedData.txt", sep=""), row.names = FALSE, sep = "\t", quote = FALSE)
        write.table(shuffledBalancedDataset, paste(antigenID, "_", nameOut, "_SlicesBalancedData.txt", sep=""), row.names = FALSE, sep = "\t", quote = FALSE)
        write.table(balancedShuffledLabels, paste(antigenID, "_", nameOut, "_SlicesBalancedShuffledLabels.txt", sep=""), row.names = FALSE, sep = "\t", quote = FALSE)
    } else {
        write.table(shuffledNonBalancedDataset, paste(antigenID, "_", nameOut, "_ImbalancedData.txt", sep=""), row.names = FALSE, sep = "\t", quote = FALSE)
        write.table(shuffledBalancedDataset, paste(antigenID, "_", nameOut, "_BalancedData.txt", sep=""), row.names = FALSE, sep = "\t", quote = FALSE)
        write.table(balancedShuffledLabels, paste(antigenID, "_", nameOut, "_BalancedShuffledLabels.txt", sep=""), row.names = FALSE, sep = "\t", quote = FALSE)
    }
    return(balancedDataset)
}

#antigenID = "1ADQ_A"
#Note: now, memory will explode, hihi
pulledCDR3 <- NULL
pulledSlices <- NULL
for(antigenID in listAvailable){
    resCDR <- generateMLTask1(antigenID, "Binder", useSlices = FALSE)
    resSlice <- generateMLTask1(antigenID, "Binder", useSlices = TRUE, uniqueSlices = TRUE)
    
	resCDRHarderV1 <- generateMLTask1(antigenID, "Binder", useSlices = FALSE, uniqueSlices=FALSE, binderRef="Mascotte", nonbinderRef="LooserExclusive", nameOut="Task1aMvsL")
	resSliceHarderV1 <- generateMLTask1(antigenID, "Binder", useSlices = TRUE, uniqueSlices=TRUE, binderRef="Mascotte", nonbinderRef="LooserExclusive", nameOut="Task1aMvsL")
    resCDRHarderV2 <- generateMLTask1(antigenID, "Binder", useSlices = FALSE, uniqueSlices = FALSE, binderRef="Heroes", nonbinderRef="MascotteExclusive", nameOut="Task1bHvsM")
    resSliceHarderV2 <- generateMLTask1(antigenID, "Binder", useSlices = TRUE, uniqueSlices = TRUE, binderRef="Heroes", nonbinderRef="MascotteExclusive", nameOut="Task1bHvsM")
	
    #pulledCDR3 <- rbind(pulledCDR3, resCDR)
    #pulledSlices <- rbind(pulledSlices, resSlice)
}

#write.table(pulledCDR3, "Task2A_mergedBalanced.txt", row.names = FALSE, sep = "\t", quote = FALSE)
#write.table(pulledSlice, "Task2A_SlicesmergedBalanced.txt", row.names = FALSE, sep = "\t", quote = FALSE)

    
    
