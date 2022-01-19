 # [![make](https://github.com/RAHenriksen/SimulAncient/actions/workflows/make.yml/badge.svg)](https://github.com/RAHenriksen/NGSNGS/actions/workflows/make.yml) NGSNGS

# NEXT GENERATION SIMULATOR FOR NEXT GENERATION SEQUENCING DATA
Rasmus Amund Henriksen, Lei Zhao, Thorfinn Sand Korneliussen 
## Installation
* make
* make HTSSRC=../hstlib

## Usage
Next Generation Simulator for Next Generator Sequencing Data version 1.0.0 

Usage
./ngsngs [options] -i <input_reference.fa> -r <Number of reads> -l/-lf <fixed length or length file> -seq <SE/PE> -f <Output format> -o <Output name>

Required: \
-i | --input: 			 Reference file in .fasta format to sample reads from \
-r | --reads: 			 Number of reads to simulate \
-l | --length: 			 Fixed length of simulated fragments, conflicts with -lf option \
-lf | --lengthfile: 		 CDF of a length distribution, conflicts with -l option \
-seq | --sequencing: 		 Simulate single-end or paired-end reads \
	 <SE>	 single-end \
 	 <PE>	 paired-end \
-f | --format: 			 File format of the simulated outpur reads \
	 <fa||fasta>		 nucletide sequence \
 	 <fa.gz||fasta.gz>	 compressed nucletide sequence \
 	 <fq||fastq>		 nucletide sequence with corresponding quality score \ 
 	 <fq.gz||fastq.gz>	 compressed nucletide sequence with corresponding quality score \ 
 	 <bam>			 Sequence Alignment Map format \
-o | --output: 			 Prefix of output file name. \

Optional: \
-h | --help: 			 Print help page \
-v | --version: 		 Print help page \
-t1 | --threads1: 		 Number of threads to use for sampling sequence reads \
-t2 | --threads2: 		 Number of threads to use write down sampled reads, default = 1 \
-s | --seed: 			 Random seed, default value being computer time \
-a1 | --adapter1: 		 Adapter sequence to add for simulated reads (SE) or first read pair (PE) \
	 e.g. Illumina TruSeq Adapter 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG  \
-a2 | --adapter2: 		 Adapter sequence to add for second read pair (PE)  \
	 e.g. Illumina TruSeq Adapter 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT \ 
-q1 | --quality1: 		 Read Quality profile for simulated reads (SE) or first read pair (PE) \
-q2 | --quality2: 		 Read Quality profile for for second read pair (PE) \
-b | --briggs: 			 Parameters for the damage patterns using the Briggs model \
	 <nv,Lambda,Delta_s,Delta_d> : 0.024,0.36,0.68,0.0097 (from Briggs et al., 2007) \
	 nv: Nick rate pr site  \
 	 Lambda: Geometric distribution parameter for overhang length \
 	 Delta_s: PMD rate in single-strand regions \
 	 Delta_s: PMD rate in double-strand regions \ 

## ERROR PROFILES AND FRAGMENT LENGTH DISTRIBUTIONS

### GENERATE ERROR PROFILES BASED ON ART's PROFILES
cd Qual_Profiles

Rscript Read_filter.R HiSeq2500L150R1filter.txt AccFreqL150R1.txt

### GENERATE SIZE DISTRIBUTIONS

## USAGE EXAMPLE
./ngsngs -i /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/chr22.fa -r 100 -s 1 -f bam -o chr22 \
