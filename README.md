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
