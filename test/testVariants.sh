#!/bin/bash

PRG=../ngsngs
LF=../Test_Examples/Size_dist_sampling.txt
Q1=../Test_Examples/AccFreqL150R1.txt
Q2=../Test_Examples/AccFreqL150R2.txt
VCFDIR=../Test_Examples
VCFIN=../Test_Examples/hg19MSub.fa

# Function to handle errors
handle_error() {
    if [ $? -ne 0 ]; then
        echo "Error: A segmentation fault or other error occurred. Exiting."
        exit 1
    fi
}

echo "---------------------------------------------------------------------------------------------------------------"
echo "------------------------------------------- Testing NGSNGS examples -------------------------------------------"
echo "---------------------------------------------------------------------------------------------------------------"
echo " "

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "-------------------------------- I) Testing genetic- and reference variations ---------------------------------"
echo "---------------------------------------------------------------------------------------------------------------"

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "1) Testing Single-end, simulating genetic variations, SNP, Insertion and deletions from Haploid vcf files, 
         for three indivdiuals with different genotype (represented as 0 1 or 2), homozygous for reference 
         and alternative and heterozygous. The alterede reference genome containing the variants are saved in 
         the HaploidRef*fa file and the simulated data are stored in the HaploidRes*fq."
echo "---------------------------------------------------------------------------------------------------------------"
echo " "
for file in $(ls ${VCFDIR}/*Haploid*); do 
    for indiv in 0 1 2; do 
        Type=$(ls ${file}| sed 's/.*Sub//'|sed 's/Haploid.vcf//');
        echo " ";
        echo ${file} Type ${Type} indiv ${indiv}; 
        ../ngsngs -i ${VCFIN} -r 100 -s 1 -l 80 -seq SE -ne -vcf ${file} -id ${indiv} -DumpVCF HaploidRef_${Type}_${indiv} -q1 ${Q1} -chr MT -f fq -o HaploidRes_${Type}_${indiv}    
        handle_error
    done; 
done

echo "---------------------------------------------------------------------------------------------------------------"
echo "2) Testing Single-end, simulating genetic variations, SNP, Insertion and deletions from Diploid vcf files, 
         for three indivdiuals with different genotype (represented as 0 1 or 2), homozygous for reference 
         and alternative and heterozygous The alterede reference genome containing the variants are saved in 
         the DiploidRef*fa file and the simulated data are stored in the DiploidRes*fq."
echo "---------------------------------------------------------------------------------------------------------------"

for file in $(ls ${VCFDIR}/*Diploid*); do 
    for indiv in 0 1 2; do 
        Type=$(ls ${file}| sed 's/.*Sub//'|sed 's/Diploid.vcf//');
        echo " ";
        echo ${file} Type ${Type} indiv ${indiv}; 
        ../ngsngs -i ${VCFIN} -r 100 -s 1 -l 80 -seq SE -ne -vcf ${file} -id ${indiv} -DumpVCF DiploidRef_${Type}_${indiv} -q1 ${Q1} -chr MT -f fq -o DiploidRes_${Type}_${indiv}   
        handle_error
    done;
done;

echo "---------------------------------------------------------------------------------------------------------------"
echo "2) Testing Single-end, simulating genetic variations, SNP, Insertion and deletions from Diploid vcf files, 
         for three indivdiuals with different genotype (represented as 0 1 or 2), homozygous for reference 
         and alternative and heterozygous The alterede reference genome containing the variants are saved in 
         the DiploidRef*fa file and the simulated data are stored in the DiploidRes*fq."
echo "---------------------------------------------------------------------------------------------------------------"

for indiv in 0 1 2; do 
    Type=$(ls ${file}| sed 's/.*Sub//'|sed 's/Diploid.vcf//');
    echo " ";
    echo ${file} Type ${Type} indiv ${indiv}; 
    ../ngsngs -i ${VCFIN} -r 100 -s 1 -l 80 -seq SE -ne -vcf ${file} -id ${indiv} -DumpVCF DiploidRef_${Type}_${indiv} -q1 ${Q1} -chr MT -f fq -o DiploidRes_${Type}_${indiv}   
    handle_error
done;


#for i in {0..2}; do ./ngsngs -i Test_Examples/hg19MSub.fa -r 1000 -t 1 -s 100 -l 150 -ne -chr MT -vcf ChrMT.vcf -qs 30 -seq SE -f fq -o sample${i} -id ${i} -DumpVCF Abby${i}Test;done
#for i in {0..2}; do awk '/^>/ {print; next} {print substr($0, 1, 40)}' Abby${i}Test.fa;done
#md5sum -c MycoBactTest.md5 || exit 2;
