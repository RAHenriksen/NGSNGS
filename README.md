# ANCIENT DNA SIMULATIONS

g++ SimulAncient_hslib.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl

Initial reporsitory for ancient DNA simulator

a.	Choose random location on the genome 

b.	Choose length of DNA either randomly or selected from a normal distribution – check the length of aDNA from articles (Thorfinn mentioned 40bp?). The random length will be similar to the fragmentation of the DNA

c.	Introduce errors – Deamination follows an age model? + check errors introduced by sequence platform

d.	Add adapters to the sequences

e.	Use multiple threads to selects chromosomes and the number of reads for the different positions

f.	Check the overlap of the regions to ensure enough coverage for the positions.

g.	Save all the reads continuously in the SAM format

# DNA damage 
1) Fragmentation: Remember that some nucleotides are more prone to breakage upon post mortem fragmentation: "It has been shown that sequence alignments of ancient DNA molecules preferentially start and end next to purines in the reference genome (Briggs et al. 2007; Sawyer et al. 2012), suggesting that
depurination and subsequent breakage of the sugar-phosphate backbone is at least partially responsible for postmortem DNA fragmentation. This pattern is also present in the Vindija 33.19 data analyzed here (Fig. 3; Supplemental Fig. S1), in which guanine and, to a lesser extent, adenine are overrepresented in the flanking bases of the reference genome, irrespective of the inferred structure of the double-stranded DNA fragment"

2) Deamination: here it shoould be possible to add or change nucleotides?

3) Cross-linking: im not sure this is possible to simulate
## Idea for input
~~~bash
Simulation_tool -ie XXX.fq -ib XXX.sam -ih XXX.fasta -dm 40 -ds 10 -cov XXX -f XXX -o XXX.fa
-ie / --endo : sequence input for endogenous DNA (fq,fastq,fa,fasta,sam,bam,cram,VCF)
-ib / --bact : sequence input for bacterical contamination (fq,fastq,fa,fasta,sam,bam,cram,VCF)
-ie / --hom : sequence input for present human contamination (fq,fastq,fa,fasta,sam,bam,cram,VCF)
-dm / --mean : mean value when creating the gaussian distribution for which the fragment sizes are collected (default == 40)
-ds / --std : standard error value when creating the gaussian distribution for which the fragment sizes are collected (default == 10)
-c / --cov : a given coverage to be obtained with simulated fragments
-f / --frag : number of fragments to simulate
-o / --output : the name and file format of simulated reads (fq,fastq,fa,fasta,sam,bam,cram,VCF).
~~~
