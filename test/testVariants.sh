#!/bin/bash

PRG=../ngsngs
LF=../Test_Examples/Size_dist_sampling.txt
Q1=../Test_Examples/AccFreqL150R1.txt
Q2=../Test_Examples/AccFreqL150R2.txt
VCFDIR=../Test_Examples
MTREF=../Test_Examples/hg19MSub.fa
MYCOREF=../Test_Examples/Mycobacterium_leprae.fa.gz
# Function to handle errors
handle_error() {
    if [ $? -ne 0 ]; then
        echo "Error: A segmentation fault or other error occurred. Exiting."
        exit 1
    fi
}

echo "---------------------------------------------------------------------------------------------------------------"
echo "-------------------------------- Testing NGSNGS examples simulating variations --------------------------------"
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
        ../ngsngs -i ${MTREF} -r 100 -s 1 -l 80 -seq SE -ne -vcf ${file} -id ${indiv} -DumpVCF HaploidRef_${Type}_${indiv} -q1 ${Q1} -chr MT -f fq -o HaploidRes_${Type}_${indiv}    
        handle_error
        #md5sum HaploidRef_${Type}_${indiv}.fa
        #md5sum HaploidRes_${Type}_${indiv}.fq
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
        ../ngsngs -i ${MTREF} -r 100 -s 1 -l 80 -seq SE -ne -vcf ${file} -id ${indiv} -DumpVCF DiploidRef_${Type}_${indiv} -q1 ${Q1} -chr MT -f fq -o DiploidRes_${Type}_${indiv}   
        handle_error
        #md5sum DiploidRef_${Type}_${indiv}.fa
        #md5sum DiploidRes_${Type}_${indiv}.fq
    done;
done;

echo "----------------------------------------------------------------------------------------------------------------------------"
echo "3) Testing Single-end, simulating STR insertions, with reference allele larger than 1 and alternative larger than reference."
echo "----------------------------------------------------------------------------------------------------------------------------"

for indiv in 0 1 2; do 
    echo " ";
    ../ngsngs -i ${MTREF} -r 1000 -t 1 -s 100 -l 150 -ne -chr MT -vcf ${VCFDIR}/ChrMTcomplex.vcf -qs 30 -seq SE -f fq -o STR_${indiv} -id ${indiv} -DumpVCF STR_internal_${indiv}
    handle_error
    #md5sum STR_internal_${indiv}.fa
    #md5sum STR_${indiv}.fq
done;

for indiv in 0 1 2; do 
    awk '/^>/ {print; next} {print substr($0, 1, 40)}' STR_internal_${indiv}.fa &>> GAT_matches.txt;
done

awk -v OFS="\t" '{count=gsub(/GAT/, ""); print count}' GAT_matches.txt |  paste GAT_matches.txt - > STR_GAT_Count.txt
rm GAT_matches.txt


echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "4) Testing Single-end, simulating stochastic variations, insertion or deletions, no error, length file "
echo "---------------------------------------------------------------------------------------------------------------"

${PRG} -i ${MYCOREF} -r 1000 -t 1 -s 1 -lf ${LF} -seq SE -ne -indel 0.0,0.05,0.0,0.9 -q1 ${Q1} -DumpIndel DelTmp -f fq -o DelOut
handle_error
${PRG} -i ${MYCOREF} -r 1000 -t 1 -s 1 -lf ${LF} -seq SE -ne -indel 0.05,0.0,0.9,0.0 -q1 ${Q1} -DumpIndel InsTmp -f fq -o InsOut
handle_error

# Function to handle errors
#for i in {0..2}; do ./ngsngs -i Test_Examples/hg19MSub.fa -r 1000 -t 1 -s 100 -l 150 -ne -chr MT -vcf ChrMT.vcf -qs 30 -seq SE -f fq -o sample${i} -id ${i} -DumpVCF Abby${i}Test;done
#for i in {0..2}; do awk '/^>/ {print; next} {print substr($0, 1, 40)}' Abby${i}Test.fa;done
md5sum -c TestVar.md5 || exit 2;
