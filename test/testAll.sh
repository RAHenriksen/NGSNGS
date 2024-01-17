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
${PRG} -i ${IN} -r 100000 -t 1 -s 1 -lf ${LF} --seq SE -m b7,0.024,0.36,0.68,0.0097 -q1 ${Q1} -f sam -o MycoBactBamSEOut

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "2) Testing Paired-end, fixed length, coverage, no sequencing error, adapters, sampling threads, fq"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -c 3 -t 1 -s 1 -l 100 --seq PE -ne -a1 ${A1} -a2 ${A2} -q1 ${Q1} -q2 ${Q2} -f fq -o MycoBactFqPEOut1
${PRG} -i ${IN} -c 3 -t 1 -s 1 -l 100 --seq PE -ne -a1 ${A1} -a2 ${A2} -q ${Q1} -f fq -o MycoBactFqPEOut2

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "3) Testing Single-end, length distributions, reads, sequencing error, misincorporation file, fa"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t 1 -s 1 -ld Pois,78 --seq SE -mf ${MF} -f fa -o MycoBactFaSEOut

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "4) Testing Single-end, fixed length, reads, sequencing error, adapters, poly-G tail"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t 1 -s 1 -l 70 --seq SE -q1 ${Q1} -a1 ${A1} -p G -f fq -o MycoBactFqSEPolyOut

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "5) Testing Single-end, fixed length greater than error profile, reads, cycle length"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t 1 -s 1 -l 200 --seq SE -q1 ${Q1} -cl 90 -f fq -o MycoBactFqSEClOut

Mean_read_length=$(awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' MycoBactFqSEClOut.fq)
if [ $Mean_read_length -ne 90 ]; then 
echo "Warning the mean read length isn't in accorandance with the cycle length"; exit 1;
fi

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
        ../ngsngs -i ${VCFIN} -r 100 -s 1 -l 80 --seq SE -ne --vcf ${file} --id ${indiv} --out_haplo HaploidRef_${Type}_${indiv} -q1 ${Q1} --chr MT -f fq -o HaploidRes_${Type}_${indiv}    
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
        ../ngsngs -i ${VCFIN} -r 100 -s 1 -l 80 --seq SE -ne --vcf ${file} --id ${indiv} --out_haplo DiploidRef_${Type}_${indiv} -q1 ${Q1} --chr MT -f fq -o DiploidRes_${Type}_${indiv}   
    done;
done;

echo "---------------------------------------------------------------------------------------------------------------"
echo "8) Testing Single-end, simulating stochastic variations, insertion or deletions, no error, length file "
echo "---------------------------------------------------------------------------------------------------------------"

${PRG} -i ${IN} -r 1000 -t 1 -s 1 -lf ${LF} --seq SE -ne --indel 0.0,0.05,0.0,0.9 -q1 ${Q1} --out_indel DelTmp.txt -f fq -o DelOut
${PRG} -i ${IN} -r 1000 -t 1 -s 1 -lf ${LF} --seq SE -ne --indel 0.05,0.0,0.9,0.0 -q1 ${Q1} --out_indel InsTmp.txt -f fq -o InsOut

echo "---------------------------------------------------------------------------------------------------------------"
echo "9) Testing Single-end, simulating fixed quality score of 40, with fragment lower-limit"
echo "---------------------------------------------------------------------------------------------------------------"

${PRG} -i ${IN} -r 100000 -t 1 -s 1 -lf ${LF} --seq SE -ll 50 -qs 40 -f fq -o MycoBactQSLLSEOUT

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
