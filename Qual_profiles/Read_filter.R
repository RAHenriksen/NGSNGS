rm(list=ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(tidyverse)
library(purrr)
setwd(getwd())

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "Accumulative_frequencies.txt"
}


# test example "HiSeq2500L150R1filter_full.txt"
ncol <- max(count.fields(args[1], sep = "\t"))
profile <- read.table(args[1], fill = TRUE, header = FALSE,col.names=paste0('V', seq_len(ncol)))

read_length <- length(unique(profile$V2))

##### CREATE THE ACCUMULATIVE FREQUENCIES WE USE FOR THE SAMPLING PROCEDURE

#Extracts every second line with the occurences of the qual
Occurence <- as.data.frame(profile[3:ncol(profile)])
Occurence <- tail(Occurence,-(nrow(profile)/6)) # removes those starting wiht "." (average)
Occurence <- Occurence %>% mutate(across(everything(), .fns = ~replace_na(.,0))) 

# Align all the occurences according to their nucleotide quality,as some positions have different
# nucleotide qualities
Occurence_qual_align <- Occurence %>%
  mutate(pairID = 3 + 2*as.integer((row_number() - 1) /2)) %>%
  {split(select(., -pairID), pull(., pairID))} %>%
  purrr::map_dfr(~{data.frame(Quality = unlist(.x[1,]),
                              Occurrence = unlist(.x[2,]))},
                 .id = "pairID") %>%
  na.omit() %>%
  arrange(as.integer(pairID)) %>%
  pivot_wider(values_from = Occurrence,
              names_from = pairID,
              names_prefix = "Row") %>%
  arrange(Quality) %>%
  as.data.frame() %>% slice(-1)  %>%
  rownames_to_column() %>% 
  gather(variable, value, -rowname) %>% 
  spread(rowname, value)

# first line is the nucleotide qualities
Qualities <- head(Occurence_qual_align,1)[,2:ncol(Occurence_qual_align)]


# the occurences
Occurence_filter <- tail(Occurence_qual_align,nrow(Occurence_qual_align)-1)

rm(profile,Occurence,Occurence_qual_align)

# Convert first column with row numbers to actual numerical values: ROW102 -> 102
x<-as.numeric(substring(Occurence_filter[,1],4,nchar(Occurence_filter[,1])))
Occurence_filter$variable <- x

colnames(Occurence_filter) <- c(paste(rep("V",ncol(Occurence_filter)),1:ncol(Occurence_filter),sep = ""))

Occurence_filter <- Occurence_filter[order(Occurence_filter$V1),]  %>% 
  remove_rownames %>% column_to_rownames(var="V1")
Occurence_filter[Occurence_filter == "NULL"] <- 0

# convert to numerical values
i <- c(seq(1,ncol(Occurence_filter))) 
Occurence_filter[ , i] <- apply(Occurence_filter[ , i], 2,function(x) as.numeric(as.character(x)))

#sum of all occurences for each position 
Occurence_filter$Sum <- c(rowSums(Occurence_filter))

#frequency of the chosen qualities of each position (150)
Frequencies <- as.data.frame(Occurence_filter[,1:ncol(Occurence_filter)-1] / Occurence_filter[,ncol(Occurence_filter)])

# Ensure the first column to contain values != 0, V2 is due to previous naming
Frequencies$V2[Frequencies$V2==0] <- 1e-40

# Change other cells columns with no actual occurences to be NA
Frequencies[Frequencies == 0] <- NA

# NaN can be an artifact from previous operations on the profile, due to the quality profile values
Frequencies$V2[Frequencies$V2=="NaN"] <- 1 # to ensure values in the first column
Frequencies[Frequencies=="NaN"] <- NA

FreqList <- list()

FreqList[[1]] <- Frequencies[,1]
# create accumulative sum
for(i in 2:ncol(Frequencies)){
  FreqList[[i]] <- rowSums(cbind(Frequencies[,1:i-1],Frequencies[,i]))
}

#rename colnames for sake of simplicity
colnames_var <- c()

for(i in 1:ncol(Frequencies)){
  nam <- paste("Col",i,sep="")
  colnames_var <- c(colnames_var, nam)
}

rm(Occurence_filter,Frequencies)

Frequency_dataframe <- data.frame(FreqList)

colnames(Frequency_dataframe) <- colnames_var
row.names(Frequency_dataframe) <- NULL

Frequency_dataframe <- Frequency_dataframe %>% mutate(across(everything(), .fns = ~replace_na(.,0))) 
#Frequency_dataframe <- Frequency_dataframe %>% mutate_all(~replace(., is.nan(.), 1))

# Create the corresponding quality values
# ensure numerical values
i <- c(seq(1,ncol(Qualities))) 
Qualities[ , i] <- apply(Qualities[ , i], 2,function(x) as.numeric(as.character(x)))
Qualities <- Qualities-1
ProbErr <- 10^(-Qualities/10) # Convert to Probability of sequencing error


Qualities<-rbind(Qualities,ProbErr)
colnames(Qualities) <- colnames_var

Final_df <- rbind(Qualities,Frequency_dataframe)

write.table(Final_df,file=args[2],sep="\t",row.names=FALSE,col.names=FALSE)
