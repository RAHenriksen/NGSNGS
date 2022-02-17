 # [![make](https://github.com/RAHenriksen/SimulAncient/actions/workflows/make.yml/badge.svg)](https://github.com/RAHenriksen/NGSNGS/actions/workflows/make.yml) NGSNGS

# NEXT GENERATION SIMULATOR FOR NEXT GENERATION SEQUENCING DATA
Rasmus Amund Henriksen, Lei Zhao, Thorfinn Sand Korneliussen \
Contact: rasmus.henriksen@sund.ku.dk

## Installation & Requirements
* Use local installation of htslib

git clone https://github.com/RAHenriksen/NGSNGS.git

git clone https://github.com/samtools/htslib.git

cd htslib; make; cd ../NGSNGS; make HTSSRC=../htslib

* Use systemwide installation of htslib

git clone https://github.com/RAHenriksen/NGSNGS.git

cd NGSNGS; make

## General
Next Generation Simulator for Next Generator Sequencing Data version 1.0.0 

Usage
~~~~bash
./ngsngs [options] -i <Reference.fa> -r/-c <Number of reads or Depth of coverage> -l/-lf <Fixed length or Length file> -seq <SE/PE> -f <Output format> -o <Prefix output name>

Options: \
-i   | --input: 		 Reference file in fasta format (.fa or .fasta) to sample reads.
-r   | --reads: 		 Number of reads to simulate, conflicts with -c option.
-c   | --coverage: 		 Depth of Coverage to simulate, conflics with -r option.
-l   | --length: 		 Fixed length of simulated fragments, conflicts with -lf option.
-lf  | --lengthfile: 		 CDF of a length distribution, conflicts with -l option.
-seq | --sequencing: 		 Simulate single-end or paired-end reads.
	 <SE>	 single-end 
 	 <PE>	 paired-end.
-f   | --format: 		 File format of the simulated output reads.
	 <fa||fasta>		 Nucletide sequence. 
 	 <fa.gz||fasta.gz>	 Compressed nucletide sequence. 
 	 <fq||fastq>		 Nucletide sequence with corresponding quality score. 
 	 <fq.gz||fastq.gz>	 Compressed nucletide sequence with corresponding quality score. 
 	 <bam>			 Sequence Alignment Map format.
-o   | --output: 		 Prefix of output file name.
-t1  | --threads1: 		 Number of threads to use for sampling sequence reads.
-t2  | --threads2: 		 Number of threads to use write down sampled reads, default = 1.
-s   | --seed: 			 Random seed, default = current calendar time (s).
-a1  | --adapter1: 		 Adapter sequence to add for simulated reads (SE) or first read pair (PE).
	 e.g. Illumina TruSeq Adapter 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG 

-a2  | --adapter2: 		 Adapter sequence to add for second read pair (PE). 
	 e.g. Illumina TruSeq Adapter 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT 

-p   | --poly: 			 Create Poly(X) tails for reads containing adapters with lengths below the inferred readcycle length. 
 	 e.g -p G or -p A 
-q1  | --quality1: 		 Read Quality profile for single-end reads (SE) or first read pair (PE).
-q2  | --quality2: 		 Read Quality profile for second read pair (PE).
-b   | --briggs: 		 Parameters for the damage patterns using the Briggs model.
	 <nv,Lambda,Delta_s,Delta_d> : 0.024,0.36,0.68,0.0097 (from Briggs et al., 2007).
	 nv: Nick rate pr site. 
 	 Lambda: Geometric distribution parameter for overhang length.
 	 Delta_s: PMD rate in single-strand regions.
 	 Delta_s: PMD rate in double-strand regions.
~~~~

## Output format
* When quality profiles are provided, NGSNGS infers the read lengths based on the position specific quality scores and creates a length limitations independent of the fixed length or the length disribution. 
* To simulate Fastq the quality profiles needs to be provided. For Sequence-Alignment-Map formats if no quality profile have been provided, then the nucleotide quality string will be the lowest quality.

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

## PIPELINE
Generate NGS reads
~~~~bash
./ngsngs -i chr22.fa -r 1000000 -t1 1 -s 1 -lf Size_dist/Size_dist_sampling.txt -q1 Qual_profiles/AccFreqL150R1.txt -seq SE -a1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG -b 0.024,0.36,0.68,0.0097 -f fq -o chr22out
~~~~
Trim adapters
~~~~bash
fastp -in1 chr22out.fq --ou1 chr22trim.fq --length_required 30
~~~~
Alignment and filtering
~~~~bash
bwa mem chr22.fa chr22trim.fq > chr22trim_se.bam
samtools view -F 4 -q 30 chr22trim_se.bam -b | samtools sort -o chr22trim_se_sort.bam
samtools index chr22trim_se_sort.bam
~~~~
Identify deamination pattern
~~~~bash
mapDamage -i chr22trim_se_sort.bam -r chr22.fa --no-stats
~~~~
