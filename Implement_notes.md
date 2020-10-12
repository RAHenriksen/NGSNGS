# SimulAncient, DNAncient - Implementation Notes

# loading the files
1) Use if , else if , else to check the proper file extension
2) Otherwise try with enum + switch and then just define all extensions in the enum. So fq,fastq,fa,fasta,sam,bam,cram,VCF.
	- What about .gz / gzip of the different files? 
	- Does capital letters influence the effect

# Number of fragments or coverage
3) Use while loop to state while coverage is below XXX continue, this is because we cannot know beforehand how 
many reads are needed to achieve a given coverage
4) Similar it would also be possible to use a while loop if we just want to simulate a given number of fragments. We can also use for loop for this,
since we know the number of reads we want to simulate

# Sampling and Nucleotide editing Output files
5) 
