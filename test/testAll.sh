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
echo "------------------------------------------- Testing NGSNGS examples -------------------------------------------"
echo "---------------------------------------------------------------------------------------------------------------"
echo " "

echo "---------------------------------------------------------------------------------------------------------------"
echo "--------------------------------- I) Testing file format (-f) and compression----------------------------------"
echo "---------------------------------------------------------------------------------------------------------------"
echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "1) Testing Single-end, length file, reads, briggs, sampling threads, sequencing error, sam,bam,cram"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 10000 -t 1 -s 5643 -lf ${LF} -seq SE -m b7,0.024,0.36,0.68,0.0097 -q1 ${Q1} -f sam -o MycoBactSamSEOut
${PRG} -i ${IN} -r 10000 -t 1 -s 5643 -lf ${LF} -seq SE -m b7,0.024,0.36,0.68,0.0097 -q1 ${Q1} -f bam -o MycoBactBamSEOut
${PRG} -i ${IN} -r 10000 -t 1 -s 5643 -lf ${LF} -seq SE -m b7,0.024,0.36,0.68,0.0097 -q1 ${Q1} -f cram -o MycoBactCramSEOut

# Check file sizes and ordering
sam_size=$(stat -c%s MycoBactSamSEOut.sam)
bam_size=$(stat -c%s MycoBactBamSEOut.bam)
cram_size=$(stat -c%s MycoBactCramSEOut.cram)
echo " "
if [ $sam_size -gt $bam_size ] && [ $bam_size -gt $cram_size ]; then
    echo "File sizes are in correct order: sam ("${sam_size}") > bam ("${bam_size}") > cram ("${cram_size}")"
else
    echo "File sizes are in correct order: sam ("${sam_size}") & bam ("${bam_size}") & cram ("${cram_size}")"; exit 1;
fi

samtools view MycoBactSamSEOut.sam > SamSE.txt
samtools view MycoBactBamSEOut.bam > BamSE.txt
samtools view MycoBactCramSEOut.cram > CramSE.txt

SamMD5=$(md5sum SamSE.txt|cut -f1 -d ' ')
BamMD5=$(md5sum BamSE.txt|cut -f1 -d ' ')

if [ "$SamMD5" != "$BamMD5" ]; then
    echo "The decompression, or mate information differs from sam to bam"; exit 1;
fi
echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "2) Testing sequence alignment map format specific information for paired-end information"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 10000 -t 1 -s 81 -lf ${LF} -seq PE -qs 40 -f bam -o MycoBactBamPEOut
${PRG} -i ${IN} -r 10000 -t 1 -s 314 -lf ${LF} -seq PE -qs 40 -na -f bam -o MycoBactBamPEOutNa

samtools view MycoBactBamPEOut.bam > BamPE.txt
samtools view MycoBactBamPEOutNa.bam > BamPEna.txt

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "3) Testing fastq compression"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 10000 -t 1 -s 22 -l 100 -seq SE -qs 20 -f fq -o MycoBactFQout
${PRG} -i ${IN} -r 10000 -t 1 -s 22 -l 100 -seq SE -qs 20 -f fq.gz -o MycoBactFQGZout
# Check file sizes and ordering
fq_size=$(stat -c%s MycoBactFQout.fq)
fqgz_size=$(stat -c%s MycoBactFQGZout.fq.gz)

echo " "
if [ $fq_size -gt $fqgz_size ]; then
    echo "Fastq file sizes are in correct order: fq ("${fq_size}") > fq.gz ("${fqgz_size}")"
else
    echo "Fastq file sizes are in correct order: fq ("${fq_size}") < fq.gz ("${fqgz_size}")"; exit 1;
fi

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "--------------------------------------- II) Testing simulation options ----------------------------------------"
echo "---------------------------------------------------------------------------------------------------------------"

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "1) Testing Paired-end, fixed length, coverage, no sequencing error, adapters, sampling threads, fq"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -c 3 -t 1 -s 1 -l 100 -seq PE -ne -a1 ${A1} -a2 ${A2} -q1 ${Q1} -q2 ${Q2} -f fq -o MycoBactFqPEOut

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "2) Testing Single-end, length distributions, reads, sequencing error, misincorporation file, fa"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t 1 -s 1 -ld Pois,78 -seq SE -mf ${MF} -f fa -o MycoBactFaSEOut

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "3) Testing Single-end, fixed length, reads, sequencing error, adapters, poly-G tail"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t 1 -s 1 -l 70 -seq SE -q1 ${Q1} -a1 ${A1} -p G -f fq -o MycoBactFqSEPolyOut

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "4) Testing Single-end, fixed length greater than error profile, reads, cycle length"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t 1 -s 1 -l 200 -seq SE -q1 ${Q1} -cl 90 -f fq -o MycoBactFqSEClOut

Mean_read_length=$(awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' MycoBactFqSEClOut.fq)
if [ $Mean_read_length -ne 90 ]; then 
echo "Warning the mean read length isn't in accorandance with the cycle length"; exit 1;
fi

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "5) Testing alignments with cycle length below fragment length with deamination and position in terms of MD:Z"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 1000 -t 1 -s 1 -l 1000 -seq SE -q1 ${Q1} -m b,0.024,0.36,0.68,0.0097 -f bam -o L1000SEDeamin

samtools sort -@ 10 -m 2G -n L1000SEDeamin.bam -o L1000SEDeaminSort.bam;samtools calmd -@ 10 -r -b L1000SEDeaminSort.bam ${IN} > L1000SEDeaminSortMD.bam

MDcheck=$(samtools view L1000SEDeaminSortMD.bam |cut -f14|awk '{ md_index = index($0, "MD:Z:"); md_substr = substr($0, md_index + 5); md_length = length(md_substr); if (md_length > 60) print 0; else print 1}'|awk '{sum += $1} END {print sum}')

if [ $MDcheck -ne 1000 ]; then 
    echo "Warning issue with MD tag above 60 characters for some reads, suggesting wrong positions stored within bam with MD tag sum value " $MDcheck; exit 1;
fi

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "------------------------- III) Testing Stochastic-, genetic- and reference variations -------------------------"
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
        echo ${file} Type ${Type} indiv ${indiv}; 
        ../ngsngs -i ${VCFIN} -r 100 -s 1 -l 80 -seq SE -ne -vcf ${file} -id ${indiv} -DumpVCF HaploidRef_${Type}_${indiv} -q1 ${Q1} -chr MT -f fq -o HaploidRes_${Type}_${indiv}    
        #md5sum HaploidRef_${Type}_${indiv}.fa >> MycoBactTest.md5;
        #md5sum HaploidRes_${Type}_${indiv}.fq >> MycoBactTest.md5;
        echo " ";
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
        #md5sum DiploidRef_${Type}_${indiv}.fa >> MycoBactTest.md5;
        #md5sum DiploidRes_${Type}_${indiv}.fq >> MycoBactTest.md5;
    done;
done;
echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "3) Testing Single-end, simulating stochastic variations, insertion or deletions, no error, length file "
echo "---------------------------------------------------------------------------------------------------------------"

${PRG} -i ${IN} -r 1000 -t 1 -s 1 -lf ${LF} -seq SE -ne -indel 0.0,0.05,0.0,0.9 -q1 ${Q1} -DumpIndel DelTmp -f fq -o DelOut
${PRG} -i ${IN} -r 1000 -t 1 -s 1 -lf ${LF} -seq SE -ne -indel 0.05,0.0,0.9,0.0 -q1 ${Q1} -DumpIndel InsTmp -f fq -o InsOut

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "4) Testing Single-end, simulating fixed quality score of 40, with fragment lower-limit"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t 1 -s 1 -lf ${LF} -seq SE -ll 50 -qs 40 -f fq -o MycoBactQSLLSEOUT

No_SeqErr=$(cat MycoBactQSLLSEOUT.fq|grep 'mod0001'|wc -l)
if [ $No_SeqErr -ne 745 ]; then 
    echo "Warning different number of reads containing sequencing error with a fixed quality score of 40"; exit 1;
fi

if [ 0 -eq 1 ]; then
    md5sum MycoBactQSLLSEOUT.fq >> MycoBactTest.md5
    md5sum DelOut.fq >> MycoBactTest.md5
    md5sum InsOut.fq >> MycoBactTest.md5
    md5sum DelTmp.txt >> MycoBactTest.md5
    md5sum InsTmp.txt >> MycoBactTest.md5
fi
echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "5) Testing mutation rate (reference genome 1.9%), sampling with replacement when below lower limit (30), 
        adapters and poly G for fixed sequencing error "
echo "---------------------------------------------------------------------------------------------------------------"

${PRG} -i ${IN} -r 100000 -f fq.gz -seq PE -s 34532 -t 1 -ld Gam,20,2 -mr 0.019 -cl 100 -p G -a1 ${A1} -a2 ${A2} -qs 30 -o MutationRateSeqErrAdapters

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "--------------------------------------------------- MD5SUM ----------------------------------------------------"
echo "---------------------------------------------------------------------------------------------------------------"
echo " "
md5sum -c MycoBactTest.md5 || exit 2;

