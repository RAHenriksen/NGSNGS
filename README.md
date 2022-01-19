 # [![make](https://github.com/RAHenriksen/SimulAncient/actions/workflows/make.yml/badge.svg)](https://github.com/RAHenriksen/NGSNGS/actions/workflows/make.yml) NGSNGS

# NEXT GENERATION SIMULATOR FOR NEXT GENERATION SEQUENCING DATA
Rasmus Amund Henriksen, Lei Zhao, Thorfinn Sand Korneliussen 
## Installation
* make
* make HTSSRC=../hstlib

## USAGE
Next Generation Simulator for Next Generator Sequencing Data version 1.0.0 

Default usage: ./ngsngs [input reference] [numer of reads] 	 Generates single-end fasta file named output

Usage: ./ngsngs \<options> [input reference] [numer of reads] [output file] 

Required options: \
-i | --input: 			 Reference file in .fasta format to sample reads from\
-r | --reads: 			 Number of reads to simulate\

Optional options: 
-h | --help: 			 Print help page
-v | --version: 		 Print help page
-o | --output: 			 Prefix of output file name. No input generates output file 'output.fa'
-f | --format: 			 File format of the simulated outpur reads
	 <.fa||.fasta>		 nucletide sequence 
 	 <.fa.gz||.fasta.gz>	 compressed nucletide sequence 
 	 <.fq||.fastq>		 nucletide sequence with corresponding quality score 
 	 <.fq.gz||.fastq.gz>	 compressed nucletide sequence with corresponding quality score 
 	 <.bam>			 Sequence Alignment Map format
-l | --length: 			 Fixed length of simulated fragments, conflicts with -lf option
-lf | --lengthfile: 		 CDF of a length distribution, conflicts with -l option
-t1 | --threads1: 		 Number of threads to use for sampling sequence reads
-t2 | --threads2: 		 Number of threads to use write down sampled reads, default = threads1
-s | --seed: 			 Random seed, default value being computer time
-a1 | --adapter1: 		 Adapter sequence to add for simulated reads (SE) or first read pair (PE)
	 e.g. Illumina TruSeq Adapter 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG 
-a2 | --adapter2: 		 Adapter sequence to add for second read pair (PE) 
	 e.g. Illumina TruSeq Adapter 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT 
-q1 | --quality1: 		 Read Quality profile for simulated reads (SE) or first read pair (PE)
-q2 | --quality2: 		 Read Quality profile for for second read pair (PE)
-seq | --sequencing: 		 Simulate single-end or paired-end reads
	 <SE>	 single-end 
 	 <PE>	 paired-end
-b | --briggs: 			 Parameters for the damage patterns using the Briggs model
	 <nv,Lambda,Delta_s,Delta_d> : 0.024,0.36,0.68,0.0097 (from Briggs et al., 2007)
		 nv: Nick rate pr site 
 		 Lambda: Geometric distribution parameter for overhang length
 		 Delta_s: PMD rate in single-strand regions
 		 Delta_s: PMD rate in double-strand regions
		 
## USAGE EXAMPLE
./ngsngs -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -r 100 -s 1 -f bam -o chr22

## TO DO
1. add adapter, sam functionality 
2. add fixed read size or empirial distribution option
3. SE vs PE option
4. VCF and msms

## GENERATE ERROR PROFILES BASED ON ART's PROFILES
cd Qual_Profiles

Rscript Read_filter.R HiSeq2500L150R1filter.txt AccFreqL150R1.txt

## GARGAMMEL
./gargammel.pl -c 1 --comp 0,0,1 -f src/sizefreq.size.gz -matfile src/matrices/double- -o data/simulation data/

./gargammel.pl -c 1 --comp 0,0,1 -f src/sizefreq.size.gz -damage 0.03,0.4,0.01,0.3 -o data/simulation data/

## ART
