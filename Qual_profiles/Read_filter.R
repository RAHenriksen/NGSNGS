rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(tidyverse)
setwd(getwd())

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "Accumulative_frequencies.txt"
}

# Reading in the error profile for read 1 however i have removed all the nucleotide 
  #lines so only the "." lines

# test example "HiSeq2500L150R1filter_full.txt"
ncol <- max(count.fields(args[1], sep = "\t"))
profile <- read.table(args[1], fill = TRUE, header = FALSE,col.names=paste0('V', seq_len(ncol)))

##### CREATE THE ACCUMULATIVE FREQUENCIES WE USE FOR THE SAMPLING PROCEDURE

#Extracts every second line with the occurences of the qual
Occurence <- as.data.frame(profile[3:ncol(profile)])
Occurence <- tail(Occurence,-(nrow(profile)/6)) # removes those positions starting with "." (average quality)
toDelete <- seq(2, nrow(Occurence), 2)
Occurence <- Occurence[ toDelete ,]
Occurence[is.na(Occurence)] = 0

#sum of all occurences for each position 
Occurence$Sum <- c(rowSums(Occurence))
#frequency of the chosen qualities of each position (150)
Frequencies <- as.data.frame(Occurence[,1:ncol(Occurence)-1] / Occurence[,ncol(Occurence)])

#create a new dataframe with the rowsums of the values
mydata <- list()

mydata[[1]] <- Frequencies[,1]
for(i in 2:ncol(Frequencies)){
  mydata[[i]] <- rowSums(cbind(Frequencies[,1:i-1],Frequencies[,i]))
}

#rename colnames for sake of simplicity
colnames <- c()

for(i in 1:ncol(Frequencies)){
  nam <- paste("Col",i,sep="")
  colnames <- c(colnames, nam)
}

mydf <- data.frame(mydata)
colnames(mydf) <- colnames
mydf <- mydf %>% mutate_all(~replace(., is.nan(.), 1))

### ASCII VALUES OF THE QUALITIES

read_length <- (nrow(profile)/6)

#col 1 is nt, col 2 is position 1-150 each repeated twice (1st for mapQ, 2nd for freq)
Qualities <- as.data.frame(profile[1:ncol]) #you could remove nt and pos with [3:ncol]
Qualities <- head(Qualities,-read_length) # removes those starting wiht "." (average)

toDelete <- seq(1, nrow(Qualities), 2)
Qualities <- Qualities[ toDelete ,]
Qualities <- Qualities[3:ncol]

Qualities <- Qualities-1 #maximum qual is 40
Qualities<- Qualities %>% 
  filter(if_all(everything(), ~ !is.na(.x)))

Qualities <- head(Qualities,1)

### FINAL DATAFRAME 
acc_freq <- rbind(Qualities, setNames(mydf, names(Qualities)))

write.table(acc_freq,file=args[2],sep="\t",row.names=FALSE,col.names=FALSE)
