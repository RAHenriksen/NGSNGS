 # [![make](https://github.com/RAHenriksen/SimulAncient/actions/workflows/make.yml/badge.svg)](https://github.com/RAHenriksen/NGSNGS/actions/workflows/make.yml) 

# NEXT GENERATION SIMULATOR FOR NEXT GENERATION SEQUENCING DATA
Rasmus Amund Henriksen, Lei Zhao, Thorfinn Sand Korneliussen \
Contact: rasmus.henriksen@sund.ku.dk

## INSTALLATION & REQUIREMENTS
* Use local installation of htslib

git clone https://github.com/RAHenriksen/NGSNGS.git

git clone https://github.com/samtools/htslib.git

cd htslib; make; cd ../NGSNGS; make HTSSRC=../htslib

* Use systemwide installation of htslib

git clone https://github.com/RAHenriksen/NGSNGS.git

cd NGSNGS; make

**NOTE:** Newer version of htslib which includes bam_set1 is required

## GENERAL
Next Generation Simulator for Next Generator Sequencing Data version 1.0.0 

~~~~bash
Usage
Next Generation Simulator for Next Generator Sequencing Data version 0.5.0 
Usage
./ngsngs [options] -i <input_reference.fa> -r/-c <Number of reads or depth of coverage> -l/-lf <fixed length or length file> -seq <SE/PE> -f <output format> -o <output name prefix>

Example 
./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -r 100000 -t1 2 -s 1 -lf Test_Examples/Size_dist/Size_dist_sampling.txt -seq SE -b 0.024,0.36,0.68,0.0097 -q1 Test_Examples/Qual_profiles/AccFreqL150R1.txt -f bam -o MycoBactBamSEOut

./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -c 3 -t1 2 -s 1 -l 100 -seq PE -ne -a1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT -q1 Test_Examples/Qual_profiles/AccFreqL150R1.txt -q2 Test_Examples/Qual_profiles/AccFreqL150R2.txt -f fq -o MycoBactFqPEOut

./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -r 100000 -t1 1 -s 1 -ld Pois,78 -seq SE -mf Test_Examples/MisincorpFile.txt -f fa -o MycoBactFaSEOut

Parameters 

-h   | --help: 				 Print help page.

Required: 

-i   | --input: 			 Reference file in fasta format (.fa,.fasta) to sample reads.
-r   | --reads: 			 Number of reads to simulate, conflicts with -c option.
-c   | --coverage: 			 Depth of Coverage to simulate, conflics with -r option.
-l   | --length: 			 Fixed length of simulated fragments, conflicts with -lf & -ld option.
-lf  | --lengthfile: 			 CDF of a length distribution, conflicts with -l & -ld option.
-ld  | --lengthdist: 			 Discrete or continuous probability distributions, conflicts with -l & -lf option.
	<Uni,40,180 || Norm,80,30 || LogNorm,4,1 || Pois,165 || Exp,0.025 || Gam,20,2>
-seq | --sequencing: 			 Simulate single-end <SE> or paired-end <PE> reads.
	 <SE||PE>
-f   | --format: 			 File format with different compression levels of the simulated output reads.
	 <fa||fasta||fa.gz||fasta.gz>
	 <fq||fastq||fq.gz||fastq.gz>		 
	 <sam||bam||cram>
-o   | --output: 			 Prefix of output file name.

Optional: 

Nucleotide Alterations: 
-bcf: 					 Binary Variant Calling Format (.bcf), only incorporates SNP's
-id  | --indiv: 			 Unique ID for individual present in the bcf header
-b   | --briggs: 			 Parameters for the damage patterns using the Briggs model.
	 <nv,Lambda,Delta_s,Delta_d> : 0.024,0.36,0.68,0.0097 (from Briggs et al., 2007).
	 nv: Nick rate pr site. 
 	 Lambda: Geometric distribution parameter for overhang length.
 	 Delta_s: PMD rate in single-strand regions.
 	 Delta_d: PMD rate in double-strand regions.
-mf  | --mismatch: 			 Nucleotide substitution frequency file.
-ne  | --noerror: 			 Disabling the nucleotide subsitutions based on nucleotide qualities.

Read Specific: 
-chr | --chromosomes: 			 Specific chromosomes from input reference file. NB! Used in conjunction with input reference genome.
	 e.g -i Hg19 -chr chr2 
-a1  | --adapter1: 			 Adapter sequence to add for simulated reads (SE) or first read pair (PE).
	 e.g Illumina TruSeq Adapter 1: <AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG> 
-a2  | --adapter2: 			 Adapter sequence to add for second read pair (PE). 
	 e.g. Illumina TruSeq Adapter 2: <AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT> 
-p   | --poly: 				 Create Poly(X) tails for reads, containing adapters (-a1, a2) with lengths below the inferred readcycle length. 
 	 <A||G||C||T||N>
-q1  | --quality1: 			 Read Quality profile for single-end reads (SE) or first read pair (PE). NB! only required for output format different from fasta
-q2  | --quality2: 			 Read Quality profile for for second read pair (PE). NB! only required for output format different from fasta

Other: 
-t1  | --threads1: 			 Number of sampling threads, default = 1.
-t2  | --threads2: 			 Number of compression threads, default = 0.
-na  | --noalign: 			 Disable the storage of alignment information for the sam||bam||cram formats.
-s   | --seed: 				 Random seed, default = current calendar time (s).
-bl  | --bufferlength: 			 Length of the nucleotide string returned to the buffer nedded to store in the output file, default = 30000000.
-rnd | --random: 			 Pseudo-random number generator, OS specific
	 e.g. linux || unix -> drand48_r (-rand = 0) || MacOS (-rand = 0) is altered to -rand 3.
~~~~

## Output format - FQ, SAM
* When quality profiles are provided, NGSNGS infers the read lengths based on the position specific quality scores and creates a readlength limitation independent of the fixed length or the length disribution. 
* To simulate .Fastq the quality profiles needs to be provided. For Sequence-Alignment-Map formats if no quality profile have been provided, then the nucleotide quality string will be the lowest quality.
* When providing the '-e' option the nucleotides will be equally substituted between on of the four nucleotide depending on the error probability of the nucleotide quality score for a given position.
* The read ID's for all the reads follows this structure:
~~~~bash
T<threadNumber>_RID<randomID>_S<strandinfo>_chr<chromosome>:<start>-<end>_length:<sequencelength>
~~~~
e.g. in .fq format
~~~~bash
@T1_RID25_S1_chr22:21441362-21441421_length:60

S0 is the forward strand and S1 is the reverse strand
~~~~
## ERROR PROFILES AND FRAGMENT LENGTH DISTRIBUTIONS
Both the nucleotide quality profile and the read length distributions are the cumulative relative frequency distribution.

### CONVERSION TO GENERATE NUCLEOTIDE QUALITY PROFILES
#### QUAL PROFILES
See more nucleotide quality profiles on https://github.com/scchess/Art. 

#### MAPDAMAGE LENGTH AND ERROR RATES

### STRUCTURE
#### STRUCTURE OF THE NUCLEOTIDE QUALITY PROFILES
1st line: Information regarding the nucleotide qualities \
2nd line: Nucleotide quality converted into the probability of the base call being wrong \
3rd - end: The CDF of the nucleotide quality distribution for a given nulceotide at a given position.  \
~~~~bash
Line1: 2	6	15	22	27	33	37	40
Line2: 0.6309573	0.2511886	0.03162278	0.006309573	0.001995262	0.0005011872	0.0001995262	0.0001
Line3: 4.99876280620546e-07	0.000500876033181788	0.0192102454642476	0.0418921316974049	0.145226556427284	1	1	1
Line4: 3.57155740260631e-07	0.000417515060364677	0.0117939968548866	0.0281474438899403	0.0982699733097515	1	1	1
Line5: 2.62498333135585e-07	0.000614246099537268	0.0446979536679942	0.105743516028673	0.239278780579743	1	1	1
.
.
.
Line750: 1	0	0	0	0	0	0	0
Line751: 1	1	0	0		0	0	0
Line752: 1	1	0	0	0	0	0	0
~~~~
Line 3 - 152 -> A nucleotide \ 
Line 153 - 302 -> T  \
Line 303 - 452 -> G  \
Line 453 - 602 -> C  \
Line 603 -> 752 -> N.
#### FRAGMENT LENGTH DISTRIBUTIONS
The CDF of the fragment lengths of Ancient DNA.
~~~~bash
35	0.00540914
36	0.01326621
37	0.02248544
38	0.03442894
39	0.04907704
.
.
.
187	0.9997630062
188	0.9998814542
189	0.999920937
190	0.9999406784
191	1
~~~~

## EXAMPLE OF USAGE
### Simulate Single-end reads, using a provided number of reads, deamination and fragment length distribution using several threads in a .bam format
~~~~bash
../ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -r 100000 -t1 2 -s 1 -lf Test_Examples/Size_dist/Size_dist_sampling.txt -seq SE -b 0.024,0.36,0.68,0.0097 -q1 Test_Examples/Qual_profiles/AccFreqL150R1.txt -f bam -o MycoBactBamSEOut
~~~~
### Simulate Paired-end reads, with a desired coverage, without sequencing errors, a fixed fragment length using several threads in a .fq format
~~~~bash
./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -c 3 -t1 2 -s 1 -l 100 -seq PE -ne -a1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT -q1 Test_Examples/Qual_profiles/AccFreqL150R1.txt -q2 Test_Examples/Qual_profiles/AccFreqL150R2.txt -f fq -o MycoBactFqPEOut
~~~~
### Simulate Single-end reads, using a provided number of reads, a fragment length poisson distribution generating nucleotide substitutions from a misincorporation file in a .fa format 
~~~~bash
./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -r 100000 -t1 1 -s 1 -ld Pois,78 -seq SE -mf Test_Examples/MisincorpFile.txt -f fa -o MycoBactFaSEOut
~~~~


## MISC
### Some analysis requires MD tag which can be added after simulations
~~~~bash
samtools sort -@ 10 -m 2G MycoBactBamSEOut.bam -o MycoBactBamSEOut_sort.bam; 
samtools index MycoBactBamSEOut_sort.bam; 
samtools calmd -@ 10 -e -r -b MycoBactBamSEOut_sort.bam Test_Examples/Mycobacterium_leprae.fa.gz > MycoBactBamSEOut_sort_MD.bam; 
rm MycoBactBamSEOut_sort.bam.bai;rm MycoBactBamSEOut_sort.bam
~~~~
