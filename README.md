 # [![make](https://github.com/RAHenriksen/SimulAncient/actions/workflows/make.yml/badge.svg)](https://github.com/RAHenriksen/NGSNGS/actions/workflows/make.yml) NGSNGS

# NEXT GENERATION SIMULATOR FOR NEXT GENERATION SEQUENCING DATA
Rasmus Amund Henriksen, Lei Zhao, Thorfinn Sand Korneliussen 
## Installation
* make
* make HTSSRC=../hstlib

## Usage
Next Generation Simulator for Next Generator Sequencing Data version 1.0.0 

Usage
./ngsngs [options] -i \<input reference.fa\> -r \<Number of reads\> -l/-lf \<fixed length or length file\> -seq \<SE/PE\> -f \<Output format\> -o \<Output name\>

Required: \
-i | --input: 			 Reference file in .fasta format to sample reads from \
-r | --reads: 			 Number of reads to simulate \
-l | --length: 			 Fixed length of simulated fragments, conflicts with -lf option \
-lf | --lengthfile: 		 CDF of a length distribution, conflicts with -l option \
-seq | --sequencing: 		 Simulate single-end or paired-end reads \
&emsp; \<SE\> &emsp;	 single-end \
&emsp; \<PE\> &emsp;	 paired-end \
-f | --format: 			 File format of the simulated outpur reads \
&emsp;	 <fa||fasta>&emsp;		 nucletide sequence \
&emsp; 	 <fa.gz||fasta.gz>&emsp;	 compressed nucletide sequence \
&emsp; 	 <fq||fastq>&emsp;		nucletide sequence with corresponding quality score \ 
&emsp; 	 <fq.gz||fastq.gz>&emsp;	 compressed nucletide sequence with corresponding quality score \ 
&emsp; 	 <bam>&emsp;			 Sequence Alignment Map format \
-o | --output: 			 Prefix of output file name. 

See helppage for further options and descriptions (./ngsngs -h)

## ERROR PROFILES AND FRAGMENT LENGTH DISTRIBUTIONS

### GENERATE ERROR PROFILES BASED ON ART's PROFILES
cd Qual_Profiles

Rscript Read_filter.R HiSeq2500L150R1filter.txt AccFreqL150R1.txt

### GENERATE SIZE DISTRIBUTIONS

## USAGE EXAMPLE
./ngsngs -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -r 100 -s 1 -f bam -o chr22 \
