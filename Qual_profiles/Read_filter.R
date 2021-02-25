rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(tidyverse)
setwd("/home/wql443/WP1/SimulAncient/Qual_profiles/")
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "Qual_Freq.txt"
}

# Reading in the error profile for read 1 however i have removed all the nucleotide 
  #lines so only the "." lines
ncol <- 10
#"HiSeq2500L150R2filter_full.txt"
profile <- read.table(args[1], fill = TRUE, header = FALSE,
         col.names = c("nucleotide","Position",
                            "Mapq2","Mapq6","Mapq15",
                            "Mapq22","Mapq27","Map33",
                            "Mapq37","Mapq40"))

read_length <- (nrow(profile)/6)

#col 1 is nt, col 2 is position 1-150 each repeated twice (1st for mapQ, 2nd for freq)
Qualitities <- as.data.frame(profile[1:ncol]) #you could remove nt and pos with [3:ncol]
Qualitities <- head(Qualitities,-read_length) # removes those starting wiht "." (average)
Qualitities <- tail(Qualitities,-read_length) # removes qual freq from N

#Extracts every second line with the qualities
toDelete <- seq(1, nrow(Qualitities), 2)
Qualitities <- Qualitities[ toDelete ,]
Qualitities <- Qualitities[3:ncol]

Qualitities <- Qualitities-1 #maximum qual is 40
Qualitities[is.na(Qualitities)] = 0

#Extracts every second line with the occurences of the qual
Occurence <- as.data.frame(profile[3:ncol])
Occurence <- head(Occurence,-read_length) # removes those starting wiht "." (average)
Occurence <- tail(Occurence,-read_length) # removes qual freq from N
toDelete <- seq(2, nrow(Occurence), 2)
Occurence <- Occurence[ toDelete ,]
Occurence[is.na(Occurence)] = 0

#sum of all occurences for each position 
Occurence$Sum <- c(rowSums(Occurence))

#frequency of the chosen qualities of each position (150)
Frequencies <- as.data.frame(Occurence[,1:8] / Occurence[,9])
write.table(Frequencies,file=args[2],sep="\t",row.names=FALSE,col.names=FALSE)
