 # [![make](https://github.com/RAHenriksen/SimulAncient/actions/workflows/make.yml/badge.svg)](https://github.com/RAHenriksen/NGSNGS/actions/workflows/make.yml) 
[![DOI](https://zenodo.org/badge/292509601.svg)](https://zenodo.org/badge/latestdoi/292509601)

# NEXT GENERATION SIMULATOR FOR NEXT GENERATION SEQUENCING DATA
Rasmus Amund Henriksen, Lei Zhao, Thorfinn Sand Korneliussen \
Contact: rasmus.henriksen@sund.ku.dk

## VERSION
The master branch represent a stable version of NGSNGS, and the development branch represent a newer and faster branch which eventually will become the master branch.

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
Next Generation Simulator for Next Generator Sequencing Data version 0.5.0 

~~~~bash
Next Generation Simulator for Next Generator Sequencing Data version 0.5.0 

Usage
./ngsngs [options] -i <input_reference.fa> -r/-c <Number of reads or depth of coverage> -l/-lf/-ld <fixed length, length file or length distribution> -seq <SE/PE> -f <output format> -o <output name prefix>

Example 
./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -r 100000 -t 2 -s 1 -lf Test_Examples/Size_dist/Size_dist_sampling.txt -seq SE -b 0.024,0.36,0.68,0.0097 -q1 Test_Examples/Qual_profiles/AccFreqL150R1.txt -f bam -o MycoBactBamSEOut

./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -c 3 -t 2 -s 1 -l 100 -seq PE -ne -a1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT -q1 Test_Examples/Qual_profiles/AccFreqL150R1.txt -q2 Test_Examples/Qual_profiles/AccFreqL150R2.txt -f fq -o MycoBactFqPEOut

./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -r 100000 -t 1 -s 1 -ld Pois,78 -seq SE -mf Test_Examples/MisincorpFile.txt -f fa -o MycoBactFaSEOut

-h   | --help: 			 Print help page.

Required: 

-i   | --input: 		 Reference file in fasta format (.fa,.fasta) to sample reads.
-r   | --reads: 		 Number of reads to simulate, conflicts with -c option.
-c   | --coverage: 		 Depth of Coverage to simulate, conflics with -r option.

-l   | --length: 		 Fixed length of simulated fragments, conflicts with -lf & -ld option.
-lf  | --lengthfile: 		 CDF of a length distribution, conflicts with -l & -ld option.
-ld  | --lengthdist: 		 Discrete or continuous probability distributions, given their Probability density function, conflicts with -l & -lf option.
	<Uni,Min,Max> 		 Uniform distribution from a closed interval given a minimum and maximum positive integer, e.g. Uni,40,180.
	<Norm,Mean,Variance> 	 Normal Distribution, given a mean and variance, e.g. Norm,80,30.
	<LogNorm,Mean,Variance>  Log-normal Distribution, given a mean and variance, e.g. LogNorm,4,1.
	<Pois,Rate> 		 Poisson distribution, given a rate, e.g. Pois,165.
	<Exp,Rate> 		 Exponential distribution, given a rate, e.g. Exp,0.025.
	<Gamma,Shape,Scale> 	 Gamma distribution, given a shape and scale, e.g. Gam,20,2.

-seq | --sequencing: 		 Simulate single-end or paired-end reads.
	 <SE>	 single-end 
 	 <PE>	 paired-end.
-f   | --format: 		 File format of the simulated output reads.
	 Nucletide sequence w. different compression levels. 
	 <fa||fasta> 
	 <fa.gz||fasta.gz>	
	 Nucletide sequence with corresponding quality score w. different compression levels. 
	 <fq||fastq> 
	 <fq.gz||fastq.gz>	
	 Sequence Alignment Map w. different compression levels. 
	 <sam||bam||cram>
-o   | --output: 		 Prefix of output file name.

Optional: 

Nucleotide Alterations: 
-bcf:| -vcf 			 Variant Calling Format (.vcf) or binary format (.bcf)
-id: | --indiv: 		 Integer value for the number of a specific individual defined in bcf header from -vcf/-bcf input file, default = -1 (no individual selected).
-b   | --briggs: 		 Parameters for the damage patterns using the Briggs model.
	 <nv,Lambda,Delta_s,Delta_d> : 0.024,0.36,0.68,0.0097 (from Briggs et al., 2007).
	 nv: Nick rate pr site. 
 	 Lambda: Geometric distribution parameter for overhang length.
 	 Delta_s: PMD rate in single-strand regions.
 	 Delta_d: PMD rate in double-strand regions.
-mf  | --mismatch: 		 Nucleotide substitution frequency file.
-ne  | --noerror: 		 Disabling the nucleotide subsitutions based on nucleotide qualities.
-na  | --noalign: 		 Using the SAM output as a sequence containing without alignment information.

Read Specific: 
-chr | --chromosomes: 		 Specific chromosomes from input reference file.
-a1  | --adapter1: 		 Adapter sequence to add for simulated reads (SE) or first read pair (PE).
	 e.g. Illumina TruSeq Adapter 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG 

-a2  | --adapter2: 		 Adapter sequence to add for second read pair (PE). 
	 e.g. Illumina TruSeq Adapter 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT 

-p   | --poly: 			 Create Poly(X) tails for reads, containing adapters with lengths below the inferred readcycle length. 
 	 e.g -p G or -p A 
-q1  | --quality1: 		 Read Quality profile for single-end reads (SE) or first read pair (PE).
-q2  | --quality2: 		 Read Quality profile for for second read pair (PE).

Other: 
-t   | --threads: 		 Number of sampling threads, default = 1.
-t2  | --threads2: 		 Number of compression threads, default = 0.
-s   | --seed: 			 Random seed, default = current calendar time (s).
-rng | --rand: 			 Pseudo-random number generator, OS specific
	 <0,1,2,3> 
	 0 :  			 drand48_r, default for linux or unix, not available for MacOS.
	 1 :  			 std::uniform_int_distribution
	 2 :  			 rand_r
	 3 :  			 erand48, default for MacOS.
~~~~

## Read ID
All formats shares a similar read ID structure
~~~~bash
T<ThreadNumber>_RID<RandomID>_S<Read1StrandInfo>_<Chromosome>:<Start>-<End>_length:<Fragmentlength> R<PairNumber>
~~~~
e.g. in .fq format
~~~~bash
@T1_RID25_S1_chr22:21441362-21441421_length:60 R1
~~~~
S0 is the forward strand and S1 is the reverse strand

## Nucleotide substitution models
### Nucleotide quality scores (-q1 and -q2)
Simulating a .fq or .sam format requires a provided nucleotide quality profile (-q1, -q2) with one example (Test_Examples/AccFreqL150R1.txt) and its structure: 
~~~~bash
Line1: 3	7	16	23	28	34	38	41
Line2: 0.501187	0.199526	0.025119	0.005012	0.001585	0.000398	0.000158	0.000079
Line3: 5.848056e-07 	5.853904e-04 	2.188810e-02 	2.653555e-02 	1.208910e-01 	1.000000e+00 	0.000000e+00 	0.000000 	
Line4: 3.960783e-07 	4.626195e-04 	1.261628e-02 	1.813564e-02 	7.776444e-02 	1.000000e+00 	0.000000e+00 	0.000000 
Line5: 3.450651e-07 	8.071072e-04 	5.794989e-02 	8.024696e-02 	1.755377e-01 	1.000000e+00 	0.000000e+00 	0.000000 	
.
.
.
Line750: 1.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00 	
Line751: 1.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00 	
Line752: 1.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00 	0.000000e+00
~~~~
* The first two lines signifies the possible quality scores to be simulated and their corresponding error probability.
* The remaining lines are seperated into 5 equal regions, i.e 750/5 -> 150. Signifying the number of positions for which a fragment can have simulated nucleotide quality scores. The first 150 lines -> A, next T, then G, then C and finally N.
* When quality profiles are provided (-q1,-q2), NGSNGS infers the maximum read length given the number of position with corresponding quality scores (number of lines) and creates a readlength limitation (e.g. (752-2)/5 = 150) independent of the fragment length options (-l,-lf,-ld). 
* To simulate .Fastq the quality profiles needs to be provided. For Sequence-Alignment-Map formats if no quality profile have been provided, then all the quality scores in the nucleotide quality string will be the lowest quality.
* From the error probability depending on the simulated quality score, sequencing errors will be simulated by equally substituting between the remaining three nucleotides. When providing the '-ne' option sequencing error substitution is disabled.
### Misincorporation file (-mf)
This misincorporation file represent the type specific probabilities of any of the bases transitioning to any other nucleotide or not transitioning for every cycle or position of the sequence.
The structure (Test_Examples/MisincorpFile.txt) is similar to that of the nucleotide quality profiles, with the number of bases for which nucleotide substitutions can occur being inferred from the dimension.
~~~~bash
	A				T			G			C
Line1: 0.865434 	0.888339 	0.953086 	1.000000
Line2: 0.882001 	0.894563 	0.979380 	1.000000
Line3: 0.915281 	0.932903 	0.983310 	1.000000
.
.
.
Line118: 0.004576 	0.023394 	0.029278 	1.000000
Line119: 0.002286 	0.017742 	0.021349 	1.000000
Line120: 0.003241 	0.015273 	0.019467 	1.000000
~~~~
* The substitution pattern from the misincorporation file, represent the substitution of all four nucleotide with the first half being from from both the 5’ termini and the second half being the 3’ termini of a given fragment. 
* With a dimension of 120 lines, the first 60 represent substitution frequencies of the first 15 positions within the read given the nucleotide belonging to A,T,G or C and the latter half being the last 15 nucleotides.
### FRAGMENT LENGTH DISTRIBUTIONS
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
* Given the above nucleotide quality profile, those reads with a fragment length above the inferred upper limit of the read length, will have an discrepancy between the read id and the simulated output sequence length. 
## EXAMPLE OF USAGE
### Simulate Single-end reads, using a provided number of reads, deamination and fragment length distribution using several threads in a .bam format
~~~~bash
../ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -r 100000 -t1 2 -s 1 -lf Test_Examples/Size_dist/Size_dist_sampling.txt -seq SE -b 0.024,0.36,0.68,0.0097 -q1 Test_Examples/AccFreqL150R1.txt -f bam -o MycoBactBamSEOut
~~~~
### Simulate Paired-end reads, with a desired coverage, without sequencing errors, a fixed fragment length using several threads in a .fq format
~~~~bash
./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -c 3 -t1 2 -s 1 -l 100 -seq PE -ne -a1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT -q1 Test_Examples/AccFreqL150R1.txt -q2 Test_Examples/AccFreqL150R2.txt -f fq -o MycoBactFqPEOut
~~~~
### Simulate Single-end reads, using a provided number of reads, a fragment length poisson distribution generating nucleotide substitutions from a misincorporation file in a .fa format 
~~~~bash
./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -r 100000 -t1 1 -s 1 -ld Pois,78 -seq SE -mf Test_Examples/MisincorpFile.txt -f fa -o MycoBactFaSEOut
~~~~
## MISC
### Example of adding an MD tag directly to the simulated bam files which can be added after simulations
~~~~bash
samtools sort -@ 10 -m 2G MycoBactBamSEOut.bam -o MycoBactBamSEOut_sort.bam; 
samtools index MycoBactBamSEOut_sort.bam; 
samtools calmd -@ 10 -e -r -b MycoBactBamSEOut_sort.bam Test_Examples/Mycobacterium_leprae.fa.gz > MycoBactBamSEOut_sort_MD.bam; 
~~~~
