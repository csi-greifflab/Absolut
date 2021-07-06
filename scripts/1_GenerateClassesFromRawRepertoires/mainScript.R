#!/usr/bin/env Rscript

#Usage: ./mainScript.sh 1ADQ_A
# should be called in the raw bindings folder, i.e. where 1ADQ_A is located as subfolder
# Output: creates a new folder 1ADQ_AAbalyses/ containing all CDR3 and slices selected

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 1) {
    stop("This script analy", call.=FALSE)
}

antigenID <- args[1]

#for testing
#antigenID <- "1FBI_X"
#fold <- paste("C:/Users/pprobert/Desktop/TestScripts/", antigenID,'/', sep="")
#outFolder <- paste("C:/Users/pprobert/Desktop/TestScripts/", antigenID,"Analyses/",sep="")

fold <- paste(antigenID,'/', sep="")
outFolder <- paste(antigenID,"Analyses/",sep="")
dir.create(outFolder)



getFilesForAntigen <- function(antigenID, folder){
    pattern <- paste("*", antigenID, ".*FinalBindings.*txt$", sep="");
    if(nchar(antigenID) == 0) { #stupid R says length("") = 1
        pattern <- ".*FinalBindings.*txt$"
    }
    list <- list.files(path = folder, pattern = pattern, recursive = TRUE, ignore.case = FALSE, include.dirs = FALSE);
    return(list)
}

#listFiles <- getFilesForAntigen(antigen, fold)

getData <- function(folder, fname){
    #print(fname)
    a <- read.table(paste(folder,fname, sep=""), header=TRUE)
    b <- a[which(a$Best == 'true'), ]
    #this command is also fine
    #c <- subset(a, a$Best == 'true')
    # for testing,
    # b <- b[1:50,]
    
    # Add an additional column with the name of the Antigen (take the file name)
    NL = length(b$Structure)
    nameCol = rep(fname, NL)
    b <- cbind(b, Antigen=nameCol)
    a <- NULL
    return(b)
}

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



# tests
#bigDataThisAntigen <- poolSlicesOverThreshold("1FBI_X", "C:/Users/pprobert/Desktop/1FBI_Manual/25000_1FBI.txt", -85)
bigDataThisAntigen <- poolDataOneAntigen (antigenID, fold)

require(ggplot2)
require(ggridges)
library(dplyr)
require(ggseqlogo)


writeForAbsolut <- function(dataframe, nameFile, antigenID){
    shorter <- dataframe %>% select(ID_slide_Variant, CDR3,	Best.,	Slide,	Energy,	Structure) %>% rename(Best = Best.)
    fileConn<-file(nameFile)
    writeLines(paste("#Antigen\t",antigenID, sep=""), fileConn)
    close(fileConn)
    write.table(shorter, nameFile, row.names = FALSE, sep = "\t",  quote = FALSE, append = TRUE) 
}

print("Quantiles for 0.00001, 0.0001, 0.001, 0.01, 0.02");
lq <- c(0.00001, 0.0001, 0.001, 0.01, 0.02)
for(i in lq){
    print(quantile(bigDataThisAntigen$Energy, i * 0.05))
}

print("Quantiles for every 5 percent");
for(i in 1:20){
    print(quantile(bigDataThisAntigen$Energy, i * 0.05))
}

print("Energy for sUper Hero (0.01%)");
thresholdSuperHero = quantile(bigDataThisAntigen$Energy, 0.0001)
print(thresholdSuperHero)
SuperHero <- bigDataThisAntigen %>% filter(bigDataThisAntigen$Energy <= thresholdSuperHero) 
print(paste("Super Heroes CDR3s are ", length(SuperHero$Energy)))
#write.table(SuperHero, paste(outFolder, antigenID, "_superHeroes.txt", sep=""), row.names = FALSE, sep = "\t", quote = FALSE)
writeForAbsolut(SuperHero, paste(outFolder, antigenID, "_superHeroes.txt", sep=""), antigenID)

print("Energy for sUper Hero (0.1%)");
thresholdHero = quantile(bigDataThisAntigen$Energy, 0.001)
print(thresholdHero)
Hero <- bigDataThisAntigen %>% filter(bigDataThisAntigen$Energy <= thresholdHero)
print(paste("Heroes CDR3s are ", length(Hero$Energy)))
#write.table(Hero, paste(outFolder, antigenID, "_Heroes.txt", sep=""), row.names = FALSE, sep = "\t", quote = FALSE)
writeForAbsolut(Hero, paste(outFolder, antigenID, "_Heroes.txt", sep=""), antigenID)

print("Energy for Mascotte (1%)");
thresholdMascotte = quantile(bigDataThisAntigen$Energy, 0.01)
print(thresholdMascotte)
Mascotte <- bigDataThisAntigen %>% filter(bigDataThisAntigen$Energy <= thresholdMascotte) 
print(paste("Mascotte CDR3s are ", length(Mascotte$Energy)))
#write.table(Mascotte, paste(outFolder, antigenID, "_Mascotte.txt", sep=""), row.names = FALSE, sep = "\t", quote = FALSE)
writeForAbsolut(Mascotte, paste(outFolder, antigenID, "_Mascotte.txt", sep=""), antigenID)

print("Energy for Looser (5%)");
thresholdLooser = quantile(bigDataThisAntigen$Energy, 0.05)
print(thresholdLooser)
Looser <- bigDataThisAntigen %>% filter(bigDataThisAntigen$Energy <= thresholdLooser) 
print(paste("Looser CDR3s are ", length(Looser$Energy)))
#write.table(Looser, paste(outFolder, antigenID, "_Looser.txt", sep=""), row.names = FALSE, sep = "\t", quote = FALSE)
writeForAbsolut(Looser, paste(outFolder, antigenID, "_Looser.txt", sep=""), antigenID)

print("Generating Non Mascottes (> 1%)");
thresholdMascotte = quantile(bigDataThisAntigen$Energy, 0.01)
print(thresholdMascotte)
NonMascotte <- bigDataThisAntigen %>% filter(bigDataThisAntigen$Energy > thresholdMascotte) 
print(paste("Non Mascotte CDR3s are ", length(NonMascotte$Energy)))
nRows = 500000  #I think the Nth row is included - no repeat - error if sample less than exist
if(nrow(NonMascotte) >= nRows){
    SampledNonMascotte <- NonMascotte[sample(nrow(NonMascotte), nRows),]
} else {
    SampledNonMascotte <- NonMascotte[sample(nrow(NonMascotte), nRows, replace=TRUE),]
}
#write.table(SampledNonMascotte, paste(outFolder, antigenID, "_500kNonMascotte.txt", sep=""), row.names = FALSE, sep = "\t", quote = FALSE)
writeForAbsolut(SampledNonMascotte, paste(outFolder, antigenID, "_500kNonMascotte.txt", sep=""), antigenID)


print("Generating LooserExclusive (5% to 1%)");
LooserExclusive <- Looser %>% filter(Looser$Energy <= thresholdLooser & Looser$Energy > thresholdMascotte) 
print(paste("Looser Exclusive CDR3s are ", length(LooserExclusive$Energy)))
writeForAbsolut(LooserExclusive, paste(outFolder, antigenID, "_LooserExclusive.txt", sep=""), antigenID)

print("Generating MascotteExclusive (5% to 1%)");
MascotteExclusive <- Mascotte %>% filter(Mascotte$Energy <= thresholdMascotte & Mascotte$Energy > thresholdHero) 
print(paste("Mascotte Exclusive CDR3s are ", length(MascotteExclusive$Energy)))
writeForAbsolut(MascotteExclusive, paste(outFolder, antigenID, "_MascotteExclusive.txt", sep=""), antigenID)

print("Generating HeroExclusive (5% to 1%)");
HeroesExclusive <- Hero %>% filter(Hero$Energy <= thresholdHero & Hero$Energy > thresholdSuperHero) 
print(paste("Heroes Exclusive CDR3s are ", length(HeroesExclusive$Energy)))
writeForAbsolut(HeroesExclusive, paste(outFolder, antigenID, "_HeroesExclusive.txt", sep=""), antigenID)

library(ggplot2)       

library(ggseqlogo)       


logoSuperHero <- ggplot() + geom_logo( as.vector(SuperHero$Slide)) + theme_logo()
ggsave(plot = logoSuperHero, width = 6, height = 4, dpi = 300, filename = paste(antigenID, "Analyses/", antigenID, "_LogoSuperHero.pdf", sep=""))                     

logoHero <- ggplot() + geom_logo( as.vector(Hero$Slide)) + theme_logo()
ggsave(plot = logoHero, width = 6, height = 4, dpi = 300, filename = paste(antigenID, "Analyses/",antigenID, "_LogoHero.pdf", sep=""))                     

logoMascotte <- ggplot() + geom_logo( as.vector(Mascotte$Slide)) + theme_logo()
ggsave(plot = logoMascotte, width = 6, height = 4, dpi = 300, filename = paste(antigenID, "Analyses/",antigenID, "_LogoMascotte.pdf", sep=""))     

logoLooser <- ggplot() + geom_logo( as.vector(Looser$Slide)) + theme_logo()
ggsave(plot = logoLooser, width = 6, height = 4, dpi = 300, filename = paste(antigenID, "Analyses/",antigenID, "_LogoLooser.pdf", sep=""))    

logoAll <- ggplot() + geom_logo( as.vector(bigDataThisAntigen$Slide)) + theme_logo()
ggsave(plot = logoAll, width = 6, height = 4, dpi = 300, filename = paste(antigenID, "Analyses/",antigenID, "_LogoAll.pdf", sep=""))    

logoHeroesExclusive <- ggplot() + geom_logo( as.vector(HeroesExclusive$Slide)) + theme_logo()
ggsave(plot = logoHeroesExclusive, width = 6, height = 4, dpi = 300, filename = paste(antigenID, "Analyses/",antigenID, "_LogoHeroesExclusive.pdf", sep=""))                     

logoMascotteExclusive <- ggplot() + geom_logo( as.vector(MascotteExclusive$Slide)) + theme_logo()
ggsave(plot = logoMascotteExclusive, width = 6, height = 4, dpi = 300, filename = paste(antigenID, "Analyses/",antigenID, "_LogoMascotteExclusive.pdf", sep=""))     

logoLooserExclusive <- ggplot() + geom_logo( as.vector(LooserExclusive$Slide)) + theme_logo()
ggsave(plot = logoLooserExclusive, width = 6, height = 4, dpi = 300, filename = paste(antigenID, "Analyses/",antigenID, "_LogoLooserExclusive.pdf", sep=""))    

listStructures <- levels(as.factor(bigDataThisAntigen$Structure))
print(listStructures)

Looser <- NULL
Mascotte <- NULL
NonMascotte <- NULL
Hero <- NULL
SuperHero <- NULL
bigDataThisAntigen <- NULL
SampledNonMascotte <- NULL



# Now generating the datasets based on the slice only:
# Note: we don't use the function getData because it filters best=true, instead we use pool function that re-reads the files
print("Generating Looser Slices (That fall within 5% top CDR3s)");
dataSlicesLooser <- poolSlicesOverThreshold(antigenID, fold, thresholdLooser)
print(paste("Looser Slices are ", length(dataSlicesLooser$Energy)))
writeForAbsolut(dataSlicesLooser, paste(outFolder, antigenID, "_LooserSlices.txt", sep=""), antigenID)

print("Generating Mascotte Slices (That fall within 1% top CDR3s)");
dataSlicesMascotte <- poolSlicesOverThreshold(antigenID, fold, thresholdMascotte)
print(paste("Mascotte Slices are ", length(dataSlicesMascotte$Energy)))
writeForAbsolut(dataSlicesMascotte, paste(outFolder, antigenID, "_MascotteSlices.txt", sep=""), antigenID)

print("Generating Heroes Slices (That fall within 0.1% top CDR3s)");
dataSlicesHeroes <- poolSlicesOverThreshold(antigenID, fold, thresholdHero)
print(paste("Heroes Slices are ", length(dataSlicesHeroes$Energy)))
writeForAbsolut(dataSlicesHeroes, paste(outFolder, antigenID, "_HeroesSlices.txt", sep=""), antigenID)

print("Generating Super Heroes Slices (That fall within 0.01% top CDR3s)");
dataSlicesSuperHeroes <- poolSlicesOverThreshold(antigenID, fold, thresholdSuperHero)
print(paste("Super Heroes Slices are ", length(dataSlicesSuperHeroes$Energy)))
writeForAbsolut(dataSlicesSuperHeroes, paste(outFolder, antigenID, "_SuperHeroesSlices.txt", sep=""), antigenID)




print("Generating Non Mascottes Slices (> 1% top CDR3s)");
NonMascotteSlices <- poolSlicesOverThreshold(antigenID, fold, thresholdMascotte, reverse=TRUE)
print(paste("Non Mascotte CDR3s are ", length(NonMascotteSlices$Energy)))
nRows = 500000  #I think the Nth row is included - no repeat - error if sample less than exist
if(nrow(NonMascotteSlices) >= nRows){
    SampledNonMascotteSlices <- NonMascotteSlices[sample(nrow(NonMascotteSlices), nRows),]
} else {
    SampledNonMascotteSlices <- NonMascotteSlices[sample(nrow(NonMascotteSlices), nRows, replace=TRUE),]
}
#write.table(SampledNonMascotte, paste(outFolder, antigenID, "_500kNonMascotte.txt", sep=""), row.names = FALSE, sep = "\t", quote = FALSE)
writeForAbsolut(SampledNonMascotteSlices, paste(outFolder, antigenID, "_500kNonMascotteSlices.txt", sep=""), antigenID)




print("Generating Looser Exclusive Slices (That fall within 5% top CDR3s)");

print("Generating Looser Exclusive Slices (That fall between 5% and 1% top CDR3s)");
dataSlicesLooserExclusive <- dataSlicesLooser %>% filter(dataSlicesLooser$Energy > thresholdMascotte)
print(paste("Looser Exclusive Slices are ", length(dataSlicesLooserExclusive$Energy)))
writeForAbsolut(dataSlicesLooserExclusive, paste(outFolder, antigenID, "_LooserExclusiveSlices.txt", sep=""), antigenID)

print("Generating Mascotte Exclusive Slices (That fall between 1% and 0.1% top CDR3s)");
dataSlicesMascotteExclusive <- dataSlicesMascotte %>% filter(dataSlicesMascotte$Energy > thresholdHero)
print(paste("Mascotte Exclusive Slices are ", length(dataSlicesMascotteExclusive$Energy)))
writeForAbsolut(dataSlicesMascotteExclusive, paste(outFolder, antigenID, "_MascotteExclusiveSlices.txt", sep=""), antigenID)

print("Generating Heroes Exclusive Slices (That fall within 0.1% and 0.01% top CDR3s)");
dataSlicesHeroesExclusive <- dataSlicesHeroes %>% filter(dataSlicesHeroes$Energy > thresholdSuperHero)
print(paste("Heroes Exclusive Slices are ", length(dataSlicesHeroesExclusive$Energy)))
writeForAbsolut(dataSlicesHeroesExclusive, paste(outFolder, antigenID, "_HeroesExclusiveSlices.txt", sep=""), antigenID)

logoSuperHeroUniqueSlices <- ggplot() + geom_logo( as.vector(unique(dataSlicesSuperHeroes$Slide))) + theme_logo()
ggsave(plot = logoSuperHeroUniqueSlices, width = 6, height = 4, dpi = 300, filename = paste(antigenID, "Analyses/",antigenID, "_LogoSuperHeroUniqueSlices.pdf", sep=""))                     

logoHeroUniqueSlices <- ggplot() + geom_logo( as.vector(unique(dataSlicesHeroes$Slide))) + theme_logo()
ggsave(plot = logoHeroUniqueSlices, width = 6, height = 4, dpi = 300, filename = paste(antigenID, "Analyses/",antigenID, "_LogoHeroUniqueSlices.pdf", sep=""))                     

logoMascotteUniqueSlices <- ggplot() + geom_logo( as.vector(unique(dataSlicesMascotte$Slide))) + theme_logo()
ggsave(plot = logoMascotteUniqueSlices, width = 6, height = 4, dpi = 300, filename = paste(antigenID, "Analyses/",antigenID, "_LogoMascotteUniqueSlices.pdf", sep=""))     

logoLooserUniqueSlices <- ggplot() + geom_logo( as.vector(unique(dataSlicesLooser$Slide))) + theme_logo()
ggsave(plot = logoLooserUniqueSlices, width = 6, height = 4, dpi = 300, filename = paste(antigenID, "Analyses/",antigenID, "_LogoLooserUniqueSlices.pdf", sep=""))    

logoNonMascotteUniqueSlices <- ggplot() + geom_logo( as.vector(unique(NonMascotteSlices$Slide))) + theme_logo()
ggsave(plot = logoNonMascotteUniqueSlices, width = 6, height = 4, dpi = 300, filename = paste(antigenID, "Analyses/",antigenID, "_LogoNonMascotteUniqueSlices.pdf", sep=""))    


stop(paste("Finished Processing antigen ", antigenID))






#plot of affinities within one structure
q <- ggplot(SuperHero, aes(x = Energy, y = factor(Structure), fill = factor(Structure))) + geom_density_ridges(alpha = 0.5) + theme_classic() + theme(legend.text=element_text(size=4), axis.text = element_text(size = 4), axis.title = element_text(size = 5), text = element_text(size = 5), legend.key.size = unit(0.2,"line"))
ggsave(plot = q, width = 8, height = 12, dpi = 300, filename = paste(antigenID, "Analyses/",antigenID, "_DistribPerStructureSuperHero.pdf", sep=""))

q <- ggplot(Hero, aes(x = Energy, y = factor(Structure), fill = factor(Structure))) + geom_density_ridges(alpha = 0.5) + theme_classic() + theme(legend.text=element_text(size=4), axis.text = element_text(size = 4), axis.title = element_text(size = 5), text = element_text(size = 5), legend.key.size = unit(0.2,"line"))
ggsave(plot = q, width = 8, height = 12, dpi = 300, filename = paste(antigenID, "Analyses/",antigenID, "_DistribPerStructureHero.pdf", sep=""))

q <- ggplot(Mascotte, aes(x = Energy, y = factor(Structure), fill = factor(Structure))) + geom_density_ridges(alpha = 0.5) + theme_classic() + theme(legend.text=element_text(size=4), axis.text = element_text(size = 4), axis.title = element_text(size = 5), text = element_text(size = 5), legend.key.size = unit(0.2,"line"))
ggsave(plot = q, width = 8, height = 12, dpi = 300, filename = paste(antigenID, "Analyses/",antigenID, "_DistribPerStructureMascotte.pdf", sep=""))

q <- ggplot(Looser, aes(x = Energy, y = factor(Structure), fill = factor(Structure))) + geom_density_ridges(alpha = 0.5) + theme_classic() + theme(legend.text=element_text(size=4), axis.text = element_text(size = 4), axis.title = element_text(size = 5), text = element_text(size = 5), legend.key.size = unit(0.2,"line"))
ggsave(plot = q, width = 8, height = 12, dpi = 300, filename = paste(antigenID, "Analyses/",antigenID, "_DistribPerStructureLooser.pdf", sep=""))

q <- ggplot(arraggedBigData, aes(x = Energy, y = factor(Structure), fill = factor(Structure))) + geom_density_ridges(alpha = 0.5) + theme_classic() + theme(legend.text=element_text(size=4), axis.text = element_text(size = 4), axis.title = element_text(size = 5), text = element_text(size = 5), legend.key.size = unit(0.2,"line"))
ggsave(plot = q, width = 8, height = 12, dpi = 300, filename = paste(antigenID, "Analyses/",antigenID, "_DistribPerStructureALL.pdf", sep=""))

print("Total number of lines")
nLinesTotal = length(bigDataThisAntigen$ID_slide_Variant)
print(nLinesTotal)

print("Frequency of structures in the full data")
FreqStruct <- bigDataThisAntigen %>% group_by(Structure) %>% mutate(freqStructure = n() / nLinesTotal) %>% slice(1)
arraggedFreqStruct <- as.data.frame(FreqStruct) %>% arrange(desc(freqStructure))
write.table(arraggedFreqStruct, paste(antigenID, "_FreqStructuresALL.txt", sep=""), row.names = FALSE, sep = "\t", quote = FALSE)
     

FreqStruct <- bigDataThisAntigen %>% group_by(Structure) %>% mutate(freqStructure = n() / nLinesTotal) %>% slice(1)
arraggedFreqStruct <- as.data.frame(FreqStruct) %>% arrange(desc(freqStructure))
write.table(arraggedFreqStruct, paste(antigenID, "_FreqStructuresALL.txt", sep=""), row.names = FALSE, sep = "\t", quote = FALSE)


FreqStruct <- bigDataThisAntigen %>% group_by(Structure) %>% mutate(freqStructure = n() / nLinesTotal) %>% slice(1)
arraggedFreqStruct <- as.data.frame(FreqStruct) %>% arrange(desc(freqStructure))
write.table(arraggedFreqStruct, paste(antigenID, "_FreqStructuresALL.txt", sep=""), row.names = FALSE, sep = "\t", quote = FALSE)


FreqStruct <- bigDataThisAntigen %>% group_by(Structure) %>% mutate(freqStructure = n() / nLinesTotal) %>% slice(1)
arraggedFreqStruct <- as.data.frame(FreqStruct) %>% arrange(desc(freqStructure))
write.table(arraggedFreqStruct, paste(antigenID, "_FreqStructuresALL.txt", sep=""), row.names = FALSE, sep = "\t", quote = FALSE)


FreqStruct <- bigDataThisAntigen %>% group_by(Structure) %>% mutate(freqStructure = n() / nLinesTotal) %>% slice(1)
arraggedFreqStruct <- as.data.frame(FreqStruct) %>% arrange(desc(freqStructure))
write.table(arraggedFreqStruct, paste(antigenID, "_FreqStructuresALL.txt", sep=""), row.names = FALSE, sep = "\t", quote = FALSE)


p2 <- ggplot(data=FreqStruct, aes(x=Structure, y=freqStructure, fill=factor(Antigen))) +
    geom_bar(stat="identity", position=position_dodge())+
    #geom_text(aes(label=freqStructure), vjust=1.6, color="orange",
    #          position = position_dodge(0.9), size=0.6, angle = 90)+
    scale_fill_brewer(palette="Paired")+ 
    theme_classic() + theme(legend.text=element_text(size=4), axis.text = element_text(size = 4), axis.title = element_text(size = 5), text = element_text(size = 5), legend.key.size = unit(0.2,"line"), axis.text.x = element_text(angle = 90, hjust = 1))
p2


nLinesMascotte = length(Mascotte$ID_slide_Variant)
nLinesMascotte
FreqStructMascotte <- Mascotte %>% group_by(Structure) %>% mutate(freqStructure = n() / nLinesMascotte) %>% slice(1)
#FreqStruct

p3 <- ggplot(data=FreqStructMascotte, aes(x=Structure, y=freqStructure, fill=factor(Antigen))) +
    geom_bar(stat="identity", position=position_dodge())+
    #geom_text(aes(label=freqStructure), vjust=1.6, color="orange",
    #          position = position_dodge(0.9), size=0.6, angle = 90)+
    scale_fill_brewer(palette="Paired")+ 
    theme_classic() + theme(legend.text=element_text(size=4), axis.text = element_text(size = 4), axis.title = element_text(size = 5), text = element_text(size = 5), legend.key.size = unit(0.2,"line"), axis.text.x = element_text(angle = 90, hjust = 1))
p3



nLinesHeroes = length(Hero$CDR3)
FreqStructHeroes <- as.data.frame(Hero %>% group_by(Structure) %>% mutate(freqStructure = n() / nLinesHeroes) %>% slice(1))
essai <- FreqStructHeroes %>% mutate(icosile = ntile(Energy, 20))
essai


FreqStructBig <- bigDataThisAntigen %>% group_by(Structure) %>% mutate(freqStructure = n() / nLinesTotal)
testGroupQuantAndFreq <- as.data.frame(FreqStructBig) %>% mutate(icosile = ntile(Energy, 20))
length(testGroupQuantAndFreq$CDR3)
head(testGroupQuantAndFreq)


p3 <- ggplot(data=testGroupQuantAndFreq, aes(x=Structure, y=freqStructure, fill=factor(icosile), color=factor(icosile))) +
    geom_bar(stat="identity", position=position_dodge())+
    #geom_text(aes(label=freqStructure), vjust=1.6, color="orange", position = position_dodge(0.9), size=0.6, angle = 90)+
    scale_fill_brewer(palette="Paired")+ 
    theme_classic() + theme(legend.text=element_text(size=4), axis.text = element_text(size = 4), axis.title = element_text(size = 5), text = element_text(size = 5), legend.key.size = unit(0.2,"line"), axis.text.x = element_text(angle = 90, hjust = 1)) + facet_grid(icosile ~ .)
p3


ggsave(plot = p3, width = 10, height = 10, dpi = 500, filename = "MatrixStructuresV2.pdf")


r <- ggplot(Mascotte, aes(x = Energy, y = factor(Structure), fill = factor(Structure))) + geom_density_ridges(alpha = 0.5) + theme_classic() + theme(legend.text=element_text(size=4), axis.text = element_text(size = 4), axis.title = element_text(size = 5), text = element_text(size = 5), legend.key.size = unit(0.2,"line"))
r

sortStruct <- FreqStruct %>% select(Structure, freqStructure) %>% arrange(desc(freqStructure)) #for descending order
sortStruct

StructuresFromBest <- as.vector(sortStruct$Structure)
StructuresFromBest

