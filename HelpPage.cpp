#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <math.h>
#include "HelpPage.h"

int HelpPage(FILE *fp){
  fprintf(fp,"Next Generation Simulator for Next Generator Sequencing Data version 0.5.0 \n\n");
  fprintf(fp,"Usage\n./ngsngs [options] -i <input_reference.fa> -r/-c <Number of reads or depth of coverage> -l/-lf/-ld <fixed length, length file or length distribution> -seq <SE/PE> -f <output format> -o <output name prefix>\n");
  fprintf(fp,"\nExample \n./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -r 100000 -t 2 -s 1 -lf Test_Examples/Size_dist_sampling.txt -seq SE -m b,0.024,0.36,0.68,0.0097 -q1 Test_Examples/AccFreqL150R1.txt -f bam -o MycoBactBamSEOut\n");
  fprintf(fp,"\n./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -c 3 -t 2 -s 1 -l 100 -seq PE -ne -a1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG -a2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT -q1 Test_Examples/AccFreqL150R1.txt -q2 Test_Examples/AccFreqL150R2.txt -f fq -o MycoBactFqPEOut\n");  
  fprintf(fp,"\n./ngsngs -i Test_Examples/Mycobacterium_leprae.fa.gz -r 100000 -t 1 -s 1 -ld Pois,78 -seq SE -mf Test_Examples/MisincorpFile.txt -f fa -o MycoBactFaSEOut\n");  
  fprintf(fp,"\n./ngsngs -i Test_Examples/hg19MSub.fa -r 1000 -t 1 -s 100 -l 150 -seq SE -ne -vcf Test_Examples/ChrMtSubDeletionDiploid.vcf -id 0 -q1 Test_Examples/AccFreqL150R1.txt -chr MT -DumpVCF DeltionInfo -f fq -o MtDeletionOut \n");  
  fprintf(fp,"\n-h   | --help: \t\t\t Print help page.\n");
  fprintf(fp,"\n----- Required ----- \n\n");
  fprintf(fp,"\t-i   | --input: \t\t Reference file in fasta format (.fa,.fasta) to sample reads.\n");
  fprintf(fp,"\nSequence reads: \n\n");
  fprintf(fp,"\t-r   | --reads: \t\t Number of reads to simulate, conflicts with -c option.\n");
  fprintf(fp,"\t-c   | --coverage: \t\t Depth of Coverage to simulate, conflics with -r option.\n");
  fprintf(fp,"\nFragment Length: \n\n");
  fprintf(fp,"\t-l   | --length: \t\t Fixed length of simulated fragments, conflicts with -lf & -ld option.\n");
  fprintf(fp,"\t-lf  | --lengthfile: \t\t CDF of a length distribution, conflicts with -l & -ld option.\n");
  fprintf(fp,"\t-ld  | --lengthdist: \t\t Discrete or continuous probability distributions, given their Probability density function, conflicts with -l & -lf option.\n");
  fprintf(fp,"\t\t<Uni,Min,Max> \t\t Uniform distribution from a closed interval given a minimum and maximum positive integer, e.g. Uni,40,180.\n");
  fprintf(fp,"\t\t<Norm,Mean,Variance> \t Normal Distribution, given a mean and variance, e.g. Norm,80,30.\n");
  fprintf(fp,"\t\t<LogNorm,Mean,Variance>  Log-normal Distribution, given a mean and variance, e.g. LogNorm,4,1.\n");
  fprintf(fp,"\t\t<Pois,Rate> \t\t Poisson distribution, given a rate, e.g. Pois,165.\n");
  fprintf(fp,"\t\t<Exp,Rate> \t\t Exponential distribution, given a rate, e.g. Exp,0.025.\n");
  fprintf(fp,"\t\t<Gam,Shape,Scale> \t Gamma distribution, given a shape and scale, e.g. Gam,20,2.\n");
  fprintf(fp,"\nOutput characteristics: \n\n");
  fprintf(fp,"\t-seq | --sequencing: \t\t Simulate single-end or paired-end reads.\n");
  fprintf(fp,"\t\t <SE>\t single-end \n \t\t <PE>\t paired-end.\n");
  fprintf(fp,"\t-f   | --format: \t\t File format of the simulated output reads.\n");
  fprintf(fp,"\t\t Nucletide sequence w. different compression levels. \n\t\t <fa||fasta> \n\t\t <fa.gz||fasta.gz>\t\n"); 
  fprintf(fp,"\t\t Nucletide sequence with corresponding quality score w. different compression levels. \n\t\t <fq||fastq> \n\t\t <fq.gz||fastq.gz>\t\n"); 
  fprintf(fp,"\t\t Sequence Alignment Map w. different compression levels. \n\t\t <sam||bam||cram>\n"); 
  fprintf(fp,"\t-o   | --output: \t\t Prefix of output file name.\n");

  fprintf(fp,"\n----- Optional ----- \n");
  fprintf(fp,"\nGenetic Variations: \n\n");
  fprintf(fp,"\t-bcf | -vcf: \t\t\t Variant Calling Format (.vcf) or binary format (.bcf)\n");
  fprintf(fp,"\t-id  | --indiv: \t\t Integer value (0 - index) for the number of a specific individual defined in bcf header from -vcf/-bcf input file, default = -1 (no individual selected).\n");
  fprintf(fp,"\t\t e.g -id 0\t\t First individual in the provided vcf file.\n");
  fprintf(fp,"\t-DumpVCF:	\t\t The prefix of an internally generated fasta file, containing the sequences representing the haplotypes with the variations from the provided vcf file (-vcf|-bcf), for diploid individuals the fasta file contains two copies of the reference genome with the each allelic genotype.\n");
  fprintf(fp,"\nStochastic Variations: \n\n");
  fprintf(fp,"\t-indel: \t\t\t Input probabilities and lambda values for a geometric distribution randomly generating insertions and deletions of a random length.\n");
  fprintf(fp,"\t\t <InsProb,DelProb,InsParam,DelParam> \t \n");
  fprintf(fp,"\t\t Insertions and deletions \t-indel 0.05,0.1,0.1,0.2 \n");
  fprintf(fp,"\t\t Only Insertions \t\t-indel 0.05,0.0,0.1,0.0 \n");
  fprintf(fp,"\t\t Only Deletions \t\t-indel 0.0,0.5,0.0,0.9 \n");
  fprintf(fp,"\t-DumpIndel: \t\t The prefix of an internally generated txt file, containing the the read id, number of indels, the number of indel operations saving the position before and after and length of the indel, simulated read length before and after, see supplementary material for detailed example and description.\n");
  fprintf(fp,"\nPostmortem damage (PMD) - Deamination: \n\n");
  fprintf(fp,"\t-m   | --model: \t\t Choice of deamination model.\n");
  fprintf(fp,"\t\t <b,nv,Lambda,Delta_s,Delta_d>  || <briggs,nv,Lambda,Delta_s,Delta_d> \t Parameters for the damage patterns using the Briggs model altered to suit modern day library preparation.\n");
  fprintf(fp,"\t\t <b7,nv,Lambda,Delta_s,Delta_d> || <briggs07,nv,Lambda,Delta_s,Delta_d>  Parameters for the damage patterns using the Briggs model 2007.\n");
  fprintf(fp,"\t\t nv: Nick rate per site. \n \t\t Lambda: Geometric distribution parameter for overhang length.\n \t\t Delta_s: PMD rate in single-strand regions.\n \t\t Delta_d: PMD rate in double-strand regions.\n");
  fprintf(fp,"\t\t e.g -m b,0.024,0.36,0.68,0.0097\n\n");
  fprintf(fp,"\t-dup | --duplicates: \t\t Number of PCR duplicates, used in conjunction with briggs modern library prep -m <b,nv,Lambda,Delta_s,Delta_d>\n");
  fprintf(fp,"\t\t <1,2,4> \t\t Default = 1\n"); 
  fprintf(fp,"\nNucleotide Alterations: \n\n");
  fprintf(fp,"\t-mf  | --mismatch: \t\t Nucleotide substitution frequency file.\n");
  fprintf(fp,"\t-ne  | --noerror: \t\t Disabling the sequencing error substitutions based on nucleotide qualities.\n");
  
  fprintf(fp,"\nRead Specific: \n\n");
  fprintf(fp,"\t-na  | --noalign: \t\t Using the SAM output as a sequence container without alignment information.\n");
  fprintf(fp,"\t-cl  | --cycle: \t\t Read cycle length, the maximum length of sequence reads, if not provided the cycle length will be inferred from quality profiles (q1,q2).\n");
  fprintf(fp,"\t-bl  | --bufferlength: \t\t Buffer length for generated sequence reads stored in the output files, default = 30000000.\n");
  fprintf(fp,"\t-chr | --chromosomes: \t\t Specific chromosomes from input reference file to simulate from.\n\n");
  fprintf(fp,"\t-a1  | --adapter1: \t\t Adapter sequence to add for simulated reads (SE) or first read pair (PE).\n");
  fprintf(fp,"\t\t e.g. Illumina TruSeq Adapter 1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG \n\n");
  fprintf(fp,"\t-a2  | --adapter2: \t\t Adapter sequence to add for second read pair (PE). \n");
  fprintf(fp,"\t\t e.g. Illumina TruSeq Adapter 2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT \n\n");
  fprintf(fp,"\t-p   | --poly: \t\t\t Create Poly(X) tails for reads, containing adapters with lengths below the inferred readcycle length. \n \t\t e.g -p G or -p A \n");
  fprintf(fp,"\t-q1  | --quality1: \t\t Read Quality profile for single-end reads (SE) or first read pair (PE).\n");
  fprintf(fp,"\t-q2  | --quality2: \t\t Read Quality profile for for second read pair (PE).\n");
  fprintf(fp,"\nSimulation Specific: \n\n");
  fprintf(fp,"\t-t   | --threads: \t\t Number of sampling threads, default = 1.\n");
  fprintf(fp,"\t-t2  | --threads2: \t\t Number of compression threads, default = 0.\n");
  fprintf(fp,"\t-s   | --seed: \t\t\t Random seed, default = current calendar time (s).\n");
  fprintf(fp,"\t-rng | --rand: \t\t\t Pseudo-random number generator, OS specific\n");
  fprintf(fp,"\t\t <0,1,2,3> \n"); 
  fprintf(fp,"\t\t 0 :  \t\t\t drand48_r, default for linux or unix, not available for MacOS.\n"); 
  fprintf(fp,"\t\t 1 :  \t\t\t std::uniform_int_distribution\n"); 
  fprintf(fp,"\t\t 2 :  \t\t\t rand_r\n"); 
  fprintf(fp,"\t\t 3 :  \t\t\t erand48, default for MacOS.\n"); 
  exit(1);
  return 0;
}

float myatof(char *str){
  if (str==NULL)
    fprintf(stderr,"Could not parse Briggs parameters, provide <nv,Lambda,Delta_s,Delta_d>");
  
  return atof(str);
}

void Sizebreak(char *str){
  if (str==NULL)
    fprintf(stderr,"Could not parse the length parameters, provide either fixed length size (-l) \n or parse length distribution file (-lf)");
}

void ErrMsg(double messageno){
  if(messageno == 1.0){fprintf(stderr,"\nInput reference file not recognized, provide -i | --input\n");}
  else if(messageno == 2.0){fprintf(stderr,"\nNumber of reads or depth to simulate not recognized, provide either number of reads (-r) or depth of coverage (-c).\n");}
  else if(messageno == 2.2){fprintf(stderr,"\nUnable to utilize the desired number for depth of coverage (-c), use values above 0.\n");}
  else if(messageno == 2.4){fprintf(stderr,"\nUnable to utilize the desired number of reads to simulate (-r), use integers above 0.\n");}
  else if(messageno == 2.6){fprintf(stderr,"\nUnable to utilize the desired number of reads to simulate for the provided number of threads, i.e threads > no. reads when no. reads equals 1.\n");}
  else if(messageno == 2.99){fprintf(stderr,"\nCould not parse both the number reads and depth to simulate, provide either number of reads (-r) or depth of coverage (-c).\n");}
  else if(messageno == 3.0){fprintf(stderr,"\nCould not parse the length parameters, provide either fixed length size (-l) or parse length distribution file (-lf).\n");}
  else if(messageno == 3.2){fprintf(stderr,"\nUnable to simulate reads of the desired fixed length (-l), use integers above 0.\n");}
  else if(messageno == 5.0){fprintf(stderr,"\nCould not parse both length parameters, provide either fixed length size (-l) or parse length distribution file (-lf).\n");}
  else if(messageno == 6.0){fprintf(stderr,"\nSequence type not provided. provide -seq || --sequence : SE (single-end) or PE (paired-end).\n");}
  else if(messageno == 6.5){fprintf(stderr,"\nSequence type not recognized. provide either SE (single-end) or PE (paired-end).\n");}
  else if(messageno == 7.0){fprintf(stderr,"\nOutput format not recognized, provide -f | --format : <fa, fa.gz, fq, fq.gz, sam, bam, cram>.\n");}
  else if(messageno == 8.0){fprintf(stderr,"\nOutput filename not provided, provide -o.\n");}
  else if(messageno == 9.0){fprintf(stderr,"\nUnable to utilize the provided number of threads, use integers above 0.\n");}
  else if(messageno == 10.0){fprintf(stderr,"\nNucleotide for poly(x) homopolymers not recognized, provide -p : <A,G,C,T,N>.\n");}
  else if(messageno == 11.0){fprintf(stderr,"\nCould not parse the Nucleotide Quality profile(s), for format <fq, fq.gz, sam, bam, cram> provide -q1 for SE and -q1, -q2 for PE.\n");}
  else if(messageno == 12.0){fprintf(stderr,"\nBoth a mismatch file and briggs deamination parameters has been provided. Provide either briggs (-p) og mismatch (-mf).\n");}
  else if(messageno == 13.0){fprintf(stderr,"\nOnly variantion type has been provided (-v). Provide variant calling format (-bcf).\n");}
  else if(messageno == 13.5){fprintf(stderr,"\nProvide variantion type not recognized, provide snp (default)\n");}
  else if(messageno == 14.0){fprintf(stderr,"Poly tail error: Missing adapter sequence, provide adapter sequence (-a1,-a2) as well\n");}
  else {fprintf(stderr,"\nError with input parameters, see helppage (-h)");}
  fprintf(stderr,"see helppage (-h)\n");
  exit(0);
}
//else if(messageno == 13.5){fprintf(stderr,"\nProvide variantion type not recognized, provide either snp, indel or all (default)\n");}

void WarMsg(double messageno){
  if(messageno == 1.0){fprintf(stderr,"\nWarning: for the output format <fa> the quality profiles (-q1, -q2) remains unused\n");}
  else if(messageno == 2.0){fprintf(stderr,"\nWarning: for the output format <fa, sam, bam, cram> the parameter (-p) is rendered moot without the nucleotide quality profiles (-q1, -q2), since poly(x) tails cannot be added to the reads since the read length cannot be inferred\n");}
  else if(messageno == 3.0){fprintf(stderr,"\nWarning: sequencing errors (-e) are not added to the sequence reads without the nucleotide quality profiles (-q1, -q2) \n");}
  else if(messageno == 4.0){fprintf(stderr,"\nWarning: for the output format <fa> the provided nucleotide qualities (-q1, -q2) will not be used \n");}
}