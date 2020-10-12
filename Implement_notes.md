# SimulAncient, DNAncient - Implementation Notes

# Loading the files
1) Use if , else if , else to check the proper file extension
2) Otherwise try with enum + switch and then just define all extensions in the enum. So fq,fastq,fa,fasta,sam,bam,cram,VCF.
	- What about .gz / gzip of the different files? 
	- Does capital letters influence the effect

# Number of fragments or coverage
3) Use while loop to state while coverage is below XXX continue, this is because we cannot know beforehand how 
many reads are needed to achieve a given coverage
4) Similar it would also be possible to use a while loop if we just want to simulate a given number of fragments. We can also use for loop for this,
since we know the number of reads we want to simulate

# Sampling and Nucleotide editing
5) Using #include <string> we should be able to use functions like .find or .replace to find the index of cytosine. Which is useful if we want to
change it to uracil and then thymin to simulate the C->T substitution
6) When we're sampling we have to think about the strand too, so if the files gives us +/- then for + we use C->T but for - we have the complementary
G->A.
	- We also have to consider the replacements should be closer to one of the specific ends.

# Output files
7) If we're sampling DNA from exo, endo and bacterial we have to first create an weighted average since i believe the endo is the least frequent
8) When we're saving the output, since its simulated shouldn't we know the actual truth? So create the reads with the following sturcute
	>Bact_XXX	ATCGTATATGC
	>Endo_XXX	AGTCAGTACAA
	>Exo_XXX	ATGTGCTACGT
