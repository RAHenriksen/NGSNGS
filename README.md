 # [![make](https://github.com/RAHenriksen/SimulAncient/actions/workflows/make.yml/badge.svg)](https://github.com/RAHenriksen/NGSNGS/actions/workflows/make.yml) NGSNGS

# NEXT GENERATION SIMULATOR FOR NEXT GENERATION SEQUENCING DATA
Rasmus Amund Henriksen, Lei Zhao, Thorfinn Sand Korneliussen 
## Installation
* make
* make HTSSRC=../hstlib

## USAGE
Next Generation Simulator for Next Generator Sequencing Data version 0.0.0 

Usage: ./ngsngs [input reference] [numer of reads] [output file]

Options: \
-h | --help: 		 Print help page\
-i | --input: 		 Reference file in .fasta format to sample reads from\
-r | --reads: 		 Number of reads to simulate from\
-o | --output: 		 Prefix of output file name, with default extension in fasta format (.fa)\
-f | --format: 		 Adapter sequence to add for simulated reads\
	 <.fa||.fasta>	 nucletide sequence \
 	 <.fq||.fastq>	 nucletide sequence with corresponding quality score\ 
 	 <.sam||bam>	 Sequence Alignment Map format\
-t | --threads: 	 Number of threads to use for simulation\
-s | --seed: 		 Random seed, default value being computer time\
-a | --adapter: 	 Adapter sequence to add for simulated reads\

## TO DO
1. add adapter, sam functionality 
2. add fixed read size or empirial distribution option
3. SE vs PE option
4. VCF and msms

## GARGAMMEL
./gargammel.pl -c 1 --comp 0,0,1 -f src/sizefreq.size.gz -matfile src/matrices/double- -o data/simulation data/

./gargammel.pl -c 1 --comp 0,0,1 -f src/sizefreq.size.gz -damage 0.03,0.4,0.01,0.3 -o data/simulation data/

## ART
