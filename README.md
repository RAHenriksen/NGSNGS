 ## [![make-tests](https://github.com/RAHenriksen/SimulAncient/actions/workflows/make.yml/badge.svg)](https://github.com/RAHenriksen/SimulAncient/actions/workflows/make.yml) ngsngs
# ANCIENT DNA SIMULATIONS

g++ SimulAncient_hslib.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl

## Ideas i need to incorporate
1) giving as input real life data
  a) Fa, Fq, Bam -> I cant make Fa and Fq into Bam since i dont know the coordinates, neither giving them VCF variations
  b) But i should be able to take a bam sequence and convert to fa with VCF

2) Arguments parsed

3) A platform type, meaning 125bp illumina, 150bp illumina and then based on that input decide on the error files

4) Lei's briggs model.

# THEORY
## DNA damage 
1) Fragmentation: Remember that some nucleotides are more prone to breakage upon post mortem fragmentation: "It has been shown that sequence alignments of ancient DNA molecules preferentially start and end next to purines in the reference genome (Briggs et al. 2007; Sawyer et al. 2012), suggesting that
depurination and subsequent breakage of the sugar-phosphate backbone is at least partially responsible for postmortem DNA fragmentation. This pattern is also present in the Vindija 33.19 data analyzed here (Fig. 3; Supplemental Fig. S1), in which guanine and, to a lesser extent, adenine are overrepresented in the flanking bases of the reference genome, irrespective of the inferred structure of the double-stranded DNA fragment"

2) Nicking and overhang -> Lei's model

3) Cross-linking: im not sure this is possible to simulate

# GARGAMMEL
./gargammel.pl -c 1 --comp 0,0,1 -f src/sizefreq.size.gz -matfile src/matrices/double- -o data/simulation data/

./gargammel.pl -c 1 --comp 0,0,1 -f src/sizefreq.size.gz -damage 0.03,0.4,0.01,0.3 -o data/simulation data/

