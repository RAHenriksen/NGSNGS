#!/bin/bash

PRG=../ngsngs
IN=../Test_Examples/Mycobacterium_leprae.fa.gz
LF=../Test_Examples/Size_dist_sampling.txt
Q1=../Test_Examples/AccFreqL150R1.txt
Q2=../Test_Examples/AccFreqL150R2.txt
A1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG
A2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT
MF=../Test_Examples/MisincorpFile.txt
VCFDIR=../Test_Examples
VCFIN=../Test_Examples/hg19MSub.fa

echo "---------------------------------------------------------------------------------------------------------------"
echo "---------------------------------------- Testing '${PRG}' examples -----------------------------------------"
echo "---------------------------------------------------------------------------------------------------------------"
echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "1) Testing Single-end, length file, reads, briggs, sampling threads, sequencing error, sam"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t 1 -s 1 -lf ${LF} --seq SE -m b7,0.024,0.36,0.68,0.0097 -q1 ${Q1} -o MycoBactBamSEOut.sam

samtools view MycoBactBamSEOut.sam|cut -f1|sort > READID.txt
samtools view MycoBactBamSEOut.sam|cut -f6|sort > CIGAR.txt
samtools view MycoBactBamSEOut.sam|cut -f10|sort > SEQ.txt

#echo "From md5 file:"
#grep 'MycoBactBamSEOut' output.md5
#echo "New simulations"
#md5sum MycoBactBamSEOut.sam
#samtools sort MycoBactBamSEOut.sam -o MycoBactBamSEOutSort.sam
#echo "Original file after sorting"
#md5sum MycoBactBamSEOut.sam
#echo "Sorted file"
#md5sum MycoBactBamSEOutSort.sam
#md5sum MycoBactBamSEOut.sam > output.md5
echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "2) Testing Paired-end, fixed length, coverage, no sequencing error, adapters, sampling threads, fq"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -c 3 -t 1 -s 1 -l 100 --seq PE -ne -a1 ${A1} -a2 ${A2} -q1 ${Q1} -q2 ${Q2} -o MycoBactFqPEOut.sam
samtools fastq -1 MycoBactFqPEOut_R1.fq -2 MycoBactFqPEOut_R2.fq MycoBactFqPEOut.sam
#md5sum MycoBactFqPEOut_R1.fq >> output.md5
#md5sum MycoBactFqPEOut_R2.fq >> output.md5

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "3) Testing Single-end, length distributions, reads, sequencing error, misincorporation file, fa"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t 1 -s 1 -ld Pois,78 --seq SE -qs 1 -mf ${MF} -o MycoBactFaSEOut.sam
samtools fasta -0 MycoBactFaSEOut.fa MycoBactFaSEOut.sam
#md5sum MycoBactFaSEOut.fa >> output.md5

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "4) Testing Single-end, fixed length, reads, sequencing error, adapters, poly-G tail"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t 1 -s 1 -l 70 --seq SE -q1 ${Q1} -a1 ${A1} -p G -o MycoBactFqSEPolyOut.sam
samtools fastq -0 MycoBactFqSEPolyOut.fq MycoBactFqSEPolyOut.sam
#md5sum MycoBactFqSEPolyOut.fq >> output.md5

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "5) Testing Single-end, fixed length greater than error profile, reads, cycle length"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t 1 -s 1 -l 200 --seq SE -q1 ${Q1} -cl 90 -o MycoBactFqSEClOut.sam
samtools fastq -0 MycoBactFqSEClOut.fq MycoBactFqSEClOut.sam

Mean_read_length=$(awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' MycoBactFqSEClOut.fq)
if [ $Mean_read_length -ne 90 ]; then 
echo "Warning the mean read length isn't in accorandance with the cycle length"; exit 1;
fi
#md5sum MycoBactFqSEClOut.fq >> output.md5

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "6) Testing Single-end, simulating genetic variations, SNP, Insertion and deletions from Haploid vcf files, 
         for three indivdiuals with different genotype (represented as 0 1 or 2), homozygous for reference 
         and alternative and heterozygous. The alterede reference genome containing the variants are saved in 
         the HaploidRef*fa file and the simulated data are stored in the HaploidRes*fq."
echo "---------------------------------------------------------------------------------------------------------------"

for file in $(ls ${VCFDIR}/*Haploid*); do 
    for indiv in 0 1 2; do 
        Type=$(ls ${file}| sed 's/.*Sub//'|sed 's/Haploid.vcf//');
        echo ${file} Type ${Type} indiv ${indiv}; 
        ${PRG} -i ${VCFIN} -r 100 -s 1 -l 80 --seq SE -ne --vcf ${file} --id ${indiv} --out_haplo HaploidRef_${Type}_${indiv} -q1 ${Q1} --chr MT -o HaploidRes_${Type}_${indiv}.sam
	samtools fastq -0 HaploidRes_${Type}_${indiv}.fq HaploidRes_${Type}_${indiv}.sam
        #md5sum HaploidRef_${Type}_${indiv}.fa >> output.md5;
        #md5sum HaploidRes_${Type}_${indiv}.fq >> output.md5;
    done; 
done

echo "---------------------------------------------------------------------------------------------------------------"
echo "7) Testing Single-end, simulating genetic variations, SNP, Insertion and deletions from Diploid vcf files, 
         for three indivdiuals with different genotype (represented as 0 1 or 2), homozygous for reference 
         and alternative and heterozygous The alterede reference genome containing the variants are saved in 
         the DiploidRef*fa file and the simulated data are stored in the DiploidRes*fq."
echo "---------------------------------------------------------------------------------------------------------------"

for file in $(ls ${VCFDIR}/*Diploid*); do 
    for indiv in 0 1 2; do 
        Type=$(ls ${file}| sed 's/.*Sub//'|sed 's/Diploid.vcf//');
        echo ${file} Type ${Type} indiv ${indiv}; 
        ${PRG} -i ${VCFIN} -r 100 -s 1 -l 80 --seq SE -ne --vcf ${file} --id ${indiv} --out_haplo DiploidRef_${Type}_${indiv} -q1 ${Q1} --chr MT -o DiploidRes_${Type}_${indiv}.sam
	samtools fastq -0 DiploidRes_${Type}_${indiv}.fq DiploidRes_${Type}_${indiv}.sam
        #md5sum DiploidRef_${Type}_${indiv}.fa >> output.md5;
        #md5sum DiploidRes_${Type}_${indiv}.fq >> output.md5;
    done;
done;

echo "---------------------------------------------------------------------------------------------------------------"
echo "8) Testing Single-end, simulating stochastic variations, insertion or deletions, no error, length file "
echo "---------------------------------------------------------------------------------------------------------------"

${PRG} -i ${IN} -r 1000 -t 1 -s 1 -lf ${LF} --seq SE -ne --indel 0.0,0.05,0.0,0.9 -q1 ${Q1} --out_indel DelTmp -o DelOut.sam
samtools fastq -0 DelOut.fq DelOut.sam
${PRG} -i ${IN} -r 1000 -t 1 -s 1 -lf ${LF} --seq SE -ne --indel 0.05,0.0,0.9,0.0 -q1 ${Q1} --out_indel InsTmp -o InsOut.sam
samtools fastq -0 InsOut.fq InsOut.sam

echo "---------------------------------------------------------------------------------------------------------------"
echo "9) Testing Single-end, simulating fixed quality score of 40, with fragment lower-limit"
echo "---------------------------------------------------------------------------------------------------------------"

${PRG} -i ${IN} -r 100000 -t 1 -s 1 -lf ${LF} --seq SE -ll 50 -qs 40 -o MycoBactQSLLSEOUT.sam
samtools fastq -0 MycoBactQSLLSEOUT.fq MycoBactQSLLSEOUT.sam

No_SeqErr=$(cat MycoBactQSLLSEOUT.fq|grep 'mod0001'|wc -l)
if [ $No_SeqErr -ne 745 ]; then 
    echo "Warning different number of reads containing sequencing error with a fixed quality score of 40"; exit 1;
fi

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "--------------------------------------------------- MD5SUM ----------------------------------------------------"
echo "---------------------------------------------------------------------------------------------------------------"
echo " "
# Remove SAM headers (they include commit sha as version)
sed -i '/^@/d' *.sam
# Compare md5sum
md5sum -c output.md5 || exit 2;
