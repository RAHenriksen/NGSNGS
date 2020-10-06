# SimulAncient
DNAncient,

Initial reporsitory for ancient DNA simulator

a.	Choose random location on the genome
b.	Choose length of DNA either randomly or selected from a normal distribution – check the length of aDNA from articles (Thorfinn mentioned 40bp?). The random length will be similar to the fragmentation of the DNA
c.	Introduce errors – Deamination follows an age model? + check errors introduced by sequence platform
d.	Add adapters to the sequences
e.	Use multiple threads to selects chromosomes and the number of reads for the different positions
f.	Check the overlap of the regions to ensure enough coverage for the positions.
g.	Save all the reads continuously in the SAM format

## STEP 7 - Merge Chimeric eccDNA configurations using the coverage file and classified reads
~~~bash
SimulAncient -ie XXX.fq -ib XXX.sam -ih XXX.fasta -dm 40 -ds 10 -cov XXX -f XXX -o XXX.fa
-ie / -endo : sequence input for endogenous DNA (fq,fastq,fa,fasta,sam,bam,cram,VCF)
-ib / -bact : sequence input for bacterical contamination (fq,fastq,fa,fasta,sam,bam,cram,VCF)
-ie / -hom : sequence input for present human contamination (fq,fastq,fa,fasta,sam,bam,cram,VCF)
-dm / -mean : mean value when creating the gaussian distribution for which the fragment sizes are collected (default == 40)
-ds / -std : standard error value when creating the gaussian distribution for which the fragment sizes are collected (default == 10)
-cov : a given coverage to be obtained with simulated fragments
-f / -frag : number of fragments to simulate
-o / -output : the name and file format of simulated reads (fq,fastq,fa,fasta,sam,bam,cram,VCF).
~~~
