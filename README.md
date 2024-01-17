 # [![make](https://github.com/RAHenriksen/NGSNGS/actions/workflows/make.yml/badge.svg)](https://github.com/RAHenriksen/NGSNGS/actions/workflows/make.yml) [![DOI:<10.1093/bioinformatics/btad041>](http://img.shields.io/badge/DOI-<10.1093/bioinformatics/btad041>-B31B1B.svg)](https://doi.org/10.1093/bioinformatics/btad041)

# NEXT GENERATION SIMULATOR FOR NEXT GENERATION SEQUENCING DATA
Rasmus Amund Henriksen, Lei Zhao, Thorfinn Sand Korneliussen \
Contact: rasmus.henriksen@sund.ku.dk

## VERSION
Next Generation Simulator for Next Generator Sequencing Data version 0.9.0, 

This article was [published](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btad041/6994180) in Oxford Academics as an application note the 20th of January 2023.

NGSNGS is a new program, therefore we are very interested in feedback to solve potential problems, as well as ideas for improvements or additions to specific and relevant features.

## INSTALLATION & REQUIREMENTS
* Use local installation of htslib

git clone https://github.com/RAHenriksen/NGSNGS.git

git clone https://github.com/samtools/htslib.git

cd htslib; make; cd ../NGSNGS; make HTSSRC=../htslib

* Use systemwide installation of htslib

git clone https://github.com/RAHenriksen/NGSNGS.git

cd NGSNGS; make

**NOTE:** Newer version of htslib which includes bam_set1 is required

## QUICK TUTORIAL
Examples of which parameters to include depending on the desired simulations
### simulate 1000 reads (-r) from human hg19.fa (-i), generate compressed fq.gz (-f), single end (--seq), make program use one threads (-t)
~~~~bash
./ngsngs -i hg19.fa -r 1000 -f fq -l 100 --seq SE -t 1 -q1 Test_Examples/AccFreqL150R1.txt -o HgSim
~~~~
### simulate 1000 reads (-r) from human hg19.fa (-i), generate fq (-f), single end (--seq), with a fixed quality score and lower fragment limit of 50
~~~~bash
./ngsngs -i hg19.fa -r 1000 -f fq -lf Test_Examples/Size_dist_sampling.txt --seq SE -t 1 -qs 40 -o HgSim
~~~~
### generate bam (-f), paired end (--seq), variable fragment length (-ld norm,350,20) but fixed readlength (-cl 100)
~~~~bash
./ngsngs -i hg19.fa -r 1000 -f bam -ld norm,350,20 -cl 100 --seq PE -t 1 -q1 Test_Examples/AccFreqL150R1.txt -q2 Test_Examples/AccFreqL150R2.txt -o HgSim
~~~~
### Disable platform specific errors (-ne), adding deamination pattern with Briggs 2007 model (-m b7,...), with ancient fragment length distribution (-lf), using seed 4 (-s)
~~~~bash
./ngsngs -i hg19.fa -r 1000 -f fq -s 4 -ne -lf Test_Examples/Size_dist_sampling.txt --seq SE -m b7,0.024,0.36,0.68,0.0097 -q1 Test_Examples/AccFreqL150R1.txt -o HgSim
~~~~
### Paired end reads, inferred cycle length from (-q1) to be 150, fragment length (-l) 400, inner distance of 100 (400-150*2)  
~~~~bash
./ngsngs -i hg19.fa -r 1000 -f fq -l 400 --seq PE -q1 Test_Examples/AccFreqL150R1.txt -q2 Test_Examples/AccFreqL150R1.txt -a1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT -o HgSim
~~~~
### Single end reads, adapter sequence (-a1) and poly-G tail (-p)
~~~~bash
./ngsngs -i hg19.fa -r 1000 -f fq -l 80 --seq SE -q1 Test_Examples/AccFreqL150R1.txt -a1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG -p G -o HgSim
~~~~
NB! the adapter sequences are only concatenated to the reads, if the inferred cycle length from quality profiles is greater than the fragment length, the poly-X tail is only added if the sequence with adapter length is still below the cycle length (Cycle length - fragment length - adapter length = 150 - 80 - 65 = 5).

## GENERAL
Next Generation Simulator for Next Generator Sequencing Data version 0.9.0 

~~~~bash
Next Generation Simulator for Next Generator Sequencing Data version 0.9.0 

Usage
./ngsngs [options] -i <input_reference.fa> -r/-c <Number of reads or depth of coverage> -l/-lf/-ld <fixed length, length file or length distribution> --seq <SE/PE> -f <output format> -o <output name prefix>

Example 
./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -r 100000 -t 2 -s 1 -lf Test_Examples/Size_dist_sampling.txt --seq SE -m b,0.024,0.36,0.68,0.0097 -q1 Test_Examples/AccFreqL150R1.txt -f bam -o MycoBactBamSEOut

./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -c 3 -t 2 -s 1 -l 100 --seq PE -ne -a1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT -q1 Test_Examples/AccFreqL150R1.txt -q2 Test_Examples/AccFreqL150R2.txt -f fq -o MycoBactFqPEOut

./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -r 100000 -t 1 -s 1 -ld Pois,78 --seq SE -mf Test_Examples/MisincorpFile.txt -f fa -o MycoBactFaSEOut

./ngsngs -i Test_Examples/hg19MSub.fa -r 1000 -t 1 -s 100 -l 150 --seq SE -ne -vcf Test_Examples/ChrMtSubDeletionDiploid.vcf -id 0 -q1 Test_Examples/AccFreqL150R1.txt -chr MT --out_haplo DeltionInfo -f fq -o MtDeletionOut 

-h   | --help: 			 Print help page.

----- Required -----

-i   | --input: 		 Reference file in fasta format (.fa,.fasta) to sample reads.

Sequence reads: 
-r   | --reads: 		 Number of reads to simulate, conflicts with -c option.
-c   | --coverage: 		 Depth of Coverage to simulate, conflics with -r option.

Fragment Length:
-l   | --length: 		 Fixed length of simulated fragments, conflicts with -lf & -ld option.
-lf  | --lengthfile: 		 CDF of a length distribution, conflicts with -l & -ld option.
-ld  | --lengthdist: 		 Discrete or continuous probability distributions, given their Probability density function, conflicts with -l & -lf option.
	 <Uni,Min,Max> 		 Uniform distribution from a closed interval given a minimum and maximum positive integer, e.g. Uni,40,180.
	 <Norm,Mean,Variance> 	 Normal Distribution, given a mean and variance, e.g. Norm,80,30.
	 <LogNorm,Mean,Variance> Log-normal Distribution, given a mean and variance, e.g. LogNorm,4,1.
	 <Pois,Rate> 		 Poisson distribution, given a rate, e.g. Pois,165.
	 <Exp,Rate> 		 Exponential distribution, given a rate, e.g. Exp,0.025.
	 <Gamma,Shape,Scale> 	 Gamma distribution, given a shape and scale, e.g. Gam,20,2.

Output characteristics:
--seq | --sequencing: 		 Simulate single-end or paired-end reads.
	 <SE||se||single||single-end>	 single-end. 
 	 <PE||pe||paired||paired-end>	 paired-end.
-f   | --format: 		 File format of the simulated output reads.
	Nucletide sequence w. different compression levels. 
	 <fa||fasta> 
	 <fa.gz||fasta.gz>
  	Nucletide sequence with corresponding quality score w. different compression levels.
	 <fq||fastq>
	 <fq.gz||fastq.gz>
  	Sequence Alignment Map w. different compression levels.
	 <sam||bam||cram> 

-o   | --out | --output:	 Output ALN file name (SAM format)

Format specific:
-q   | --quality:		 Read Quality profile for all reads for fastq or sequence alignment map formats.
-q1  | --quality1:		 Read Quality profile for single-end reads (SE) or first read pair (PE) for fastq or sequence alignment map formats.
-q2  | --quality2:		 Read Quality profile for for second read pair (PE) for fastq or sequence alignment map formats.
-qs  | --qualityscore:		 Fixed quality score, for both read pairs in fastq or sequence alignment map formats. It overwrites the quality profiles.

----- Optional -----

Genetic Variations:

--bcf | --vcf: 			 Variant Calling Format (.vcf) or binary format (.bcf)
--id  | --indiv: 		 Index (0 - based) of a specific individual from --vcf/--bcf input file, default = -1 (no individual selected).
--out_haplo:			 Output FASTA file with internally generated haplotypes with the variations from the provided `vcf` file; for diploid individuals
				 the fasta file contains two copies of the reference genome with the each allelic genotype.

Stochastic Variations:

--indel:			 Input probabilities and lambda values for a geometric distribution randomly generating insertions and deletions of a random length.
	 <InsProb,DelProb,InsParam,DelParam>
	 Insertions and deletions --indel 0.05,0.1,0.1,0.2
	 Only Insertions          --indel 0.05,0.0,0.1,0.0
	 Only Deletions           --indel 0.0,0.5,0.0,0.9
--out_indel:			 Output text file, containing the the read id, number of indels, the number of indel operations saving the position
				 before and after and length of the indel, simulated read length before and after.

Postmortem damage (PMD) - Deamination:

-m | --model:			 Choice of deamination model.
	 <b,nv,Lambda,Delta_s,Delta_d>  || <briggs,nv,Lambda,Delta_s,Delta_d>	Parameters for the damage patterns using the Briggs model altered to suit modern day library preparation.
	 <b7,nv,Lambda,Delta_s,Delta_d> || <briggs07,nv,Lambda,Delta_s,Delta_d>  Parameters for the damage patterns using the Briggs model 2007.
	 nv: Nick rate pr site.
	 Lambda: Geometric distribution parameter for overhang length.
	 Delta_s: PMD rate in single-strand regions.
	 Delta_d: PMD rate in double-strand regions.
	 e.g -m b,0.024,0.36,0.68,0.0097
--dup | --duplicates:	 	 Number of PCR duplicates, used in conjunction with briggs modern library prep -m <b,nv,Lambda,Delta_s,Delta_d>
	 <1,2,4>, Default = 1.

Nucleotide Alterations:

-mf  | --mismatch: 		 Nucleotide substitution frequency file.
-ne  | --noerror: 		 Disabling the nucleotide substitutions based on nucleotide qualities.

Read Specific:
-na  | --noalign: 		 Using the SAM output as a sequence containing without alignment information.
-cl  | --cycle:			 Read cycle length, the maximum length of sequence reads, if not provided the cycle length will be inferred from quality profiles (q1,q2).
-ll  | --lowerlimit:	 Lower fragment length limit, default = 30. The minimum fragment length for deamination is 30, so simulated fragments below will be fixed at 30.
-bl  | --bufferlength:		 Buffer length for generated sequence reads stored in the output files, default = 30000000.
-chr | --chromosomes: 		 Specific chromosomes from input reference file.
-a1  | --adapter1: 		 Adapter sequence to add for simulated reads (SE) or first read pair (PE).
	 e.g. Illumina TruSeq Adapter 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG 

-a2  | --adapter2: 		 Adapter sequence to add for second read pair (PE). 
	 e.g. Illumina TruSeq Adapter 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT 

-p   | --poly: 			 Create Poly(X) tails for reads, containing adapters with lengths below the inferred readcycle length. 
 	 e.g -p G or -p A 

Simulation Specific: 
-t    | --threads: 		 Number of sampling threads, default = 1.
-t2   | --threads2: 		 Number of compression threads, default = 0.
-s    | --seed: 		 Random seed, default = current calendar time (s).
--rng | --rand: 		 Pseudo-random number generator, OS specific
	 <0,1,2,3> 
	 0 :  			 drand48_r, default for linux or unix, not available for MacOS.
	 1 :  			 std::uniform_int_distribution
	 2 :  			 rand_r
	 3 :  			 erand48, default for MacOS.
~~~~

## Read ID
All formats shares a similar read ID structure
~~~~bash
T<ThreadNumber>_RID<RandomID>_S<Read1StrandInfo>_<Chromosome>:<Start>-<End>
_length:<Fragmentlength>_<modVal1Val2Val3Val4> F<FragmentNumber> R<PairNumber>
~~~~
e.g.
@T0_RID49_S0_NZ_CP029543.1:2236795-2236942_length:148_mod1000 F0 R1

S0 is the forward strand and S1 is the reverse strand, mod1000 equals read with deamination, F0 signifies the first fragment out of 4 possible PCR duplicates, R1 indicate the sequence is read 1 (See supplementary material for detailed description).

### Modification vector modVal1Val2Val3Val4
* Val1 signifies deamination of the simulated fragment, 0 = no deamination, 1 = deamination.
* Val2 signifies mismatches from file (-mf) of the simulated fragment, 0 = no mismatch, 1 = 5' mismatch, 2 = 3' mismatch, 3 = mismatch in both ends.
* Val3 signifies stochastic structural variations in the sequenced reads, 0 = no variation, 1 = insertions, 2 = deletions, 3 = insertions and deletions. 
* Val4 signifies sequencing error (dependent on -q1,-q2 or -qs) in the sequenced reads, 0 = no sequencing error, 1 = sequencing error.

## Nucleotide substitution models
### Nucleotide quality scores (-q1 and -q2)
Simulating a .fq or .sam format requires a provided nucleotide quality profile (-q1, -q2) with one example (Test_Examples/AccFreqL150R1.txt) and its structure: 
~~~~bash
Line1: 3		7		16		23		28		34		38		41
Line2: 0.501187		0.199526	0.025119	0.005012	0.001585	0.000398	0.000158	0.000079
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
* The remaining lines are seperated into 5 equal regions, i.e 750/5 -> 150. Signifying the number of positions for which a read can have simulated nucleotide quality scores. The first 150 lines -> A, next T, then G, then C and finally N.
* When quality profiles are provided (-q1,-q2), NGSNGS infers the cycle length given the number of position with corresponding quality scores (number of lines) and creates a read length limitation (e.g. (752-2)/5 = 150) independent of the fragment length options (-l,-lf,-ld). 
* To simulate .fastq outputs, the quality profiles needs to be provided. For Sequence-Alignment-Map formats if no quality profile have been provided, then all the quality scores in the nucleotide quality string will be the lowest quality.
* From the error probability depending on the simulated quality score, sequencing errors will be simulated by equally substituting between the remaining three nucleotides. When providing the '-ne' option sequencing error substitution is disabled.
### Misincorporation file (-mf)
This misincorporation file represent the type specific probabilities of any of the bases transitioning to any other nucleotide or not transitioning for every cycle or position of the sequence.
The structure (Test_Examples/MisincorpFile.txt) is similar to that of the nucleotide quality profiles, with the number of bases for which nucleotide substitutions can occur being inferred from the dimension.
~~~~bash
	A 	 	T 	 	G  		C
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

### Fragment length distribution file (-lf)
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
* Given the above nucleotide quality profile, those reads with a fragment length above the inferred the cycle length (150 bp), will have an discrepancy between the read id "length:&lt;Fragmentlength&gt;" and the simulated output sequence length. 

## CITATION
Bibtex citation
~~~~bash
@article{10.1093/bioinformatics/btad041,
    author = {Henriksen, Rasmus Amund and Zhao, Lei and Korneliussen, Thorfinn Sand},
    title = "{NGSNGS: Next generation simulator for next generation sequencing data}",
    journal = {Bioinformatics},
    year = {2023},
    month = {01},
    issn = {1367-4811},
    doi = {10.1093/bioinformatics/btad041},
    url = {https://doi.org/10.1093/bioinformatics/btad041},
    note = {btad041},
    eprint = {https://academic.oup.com/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btad041/48800233/btad041.pdf},
}
~~~~

## MISC
### Example of adding an MD tag directly to the simulated bam files which can be added after simulations
~~~~bash
samtools sort -@ 10 -m 2G MycoBactBamSEOut.bam -o MycoBactBamSEOut_sort.bam; 
samtools index MycoBactBamSEOut_sort.bam; 
samtools calmd -@ 10 -r -b MycoBactBamSEOut_sort.bam Test_Examples/Mycobacterium_leprae.fa.gz > MycoBactBamSEOut_sort_MD.bam; 
~~~~
 
