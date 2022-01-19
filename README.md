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

## EXAMPLES
### Simulate simple NGS paired-end output in fasta format with fixed length
~~~~bash
./ngsngs -i chr22.fa -r 10000 -l 100 -seq PE -f fa -o chr22pe
~~~~
### Simulate single-end NGS aDNA reads with deamination 
~~~~bash
./ngsngs -i chr22.fa -r 100000 -t1 2 -s 1 -lf Size_dist/Size_dist_sampling.txt -seq SE -b 0.024,0.36,0.68,0.0097 -f fa -o chr22se
# -t1: Number of threads to use for sampling sequence reads
# -s: Random seed
# -b: Parameters for the damage patterns using the Briggs model <nv,Lambda,Delta_s,Delta_d>
~~~~
### Simulate paired-end reads with fixed sized, quality profiles and adapters
~~~~bash
./ngsngs -i chr22.fa -r 100000 -t1 2 -s 1 -l 125 -q1 Qual_profiles/AccFreqL150R1.txt -q2 Qual_profiles/AccFreqL150R2.txt -a1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT -seq PE -f fq -o chr22pe
~~~~
### Create bam file with deamination pattern
~~~~bash
./ngsngs -i chr22.fa -r 100000 -t1 2 -s 1 -lf Size_dist/Size_dist_sampling.txt -seq SE -b 0.024,0.36,0.68,0.0097 -q1 Qual_profiles/AccFreqL150R1.txt -f bam -o chr22se
~~~~
Verify using mapDamage
~~~~bash
mapDamage -i chr22se.bam -r chr22.fa --no-stats
~~~~
or 
~~~~bash
samtools sort chr22se.bam -o chr22sesort.bam
samtools index chr22sesort.bam
mapDamage -i chr22sesort.bam -r chr22.fa --no-stats
~~~~
