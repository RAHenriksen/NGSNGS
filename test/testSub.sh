#!/bin/bash

PRG=../ngsngs
IN=../Test_Examples/Mycobacterium_leprae.fa.gz
LF=../Test_Examples/Size_dist_sampling.txt
Q1=../Test_Examples/AccFreqL150R1.txt
Q2=../Test_Examples/AccFreqL150R2.txt
MF=../Test_Examples/MisincorpFile.txt
M3=../Test_Examples/GlobalM3.bdamage.gz
# Function to handle errors
handle_error() {
    if [ $? -ne 0 ]; then
        echo "Error: A segmentation fault or other error occurred. Exiting."
        exit 1
    fi
}

echo "--------------------------------------------------------------------------------------------------------------"
echo "----------------------------- I) Testing different nucleotide substiution models -----------------------------"
echo "--------------------------------------------------------------------------------------------------------------"
echo " "
echo " "


echo " "
echo "-----------------------------------------------------------------------------------------------------------------------"
echo "1) Testing different Post-Mortem-Damage models representing different library preparation methods"
echo "-----------------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 1000 -t 1 -s 1 -l 100 -seq SE -qs 40 -ne -m Illumina,0.024,0.36,0.68,0.0097 -f fq -o IlluminaPMD
#md5sum IlluminaPMD.fq
handle_error

${PRG} -i ${IN} -r 1000 -t 1 -s 1 -l 100 -seq SE -qs 40 -ne -m Roche454,0.024,0.36,0.68,0.0097 -f fq -o RochePMD
#md5sum RochePMD.fq
handle_error

#${PRG} -i ${IN} -r 1000 -t 1 -s 1 -l 100 -seq SE -qs 40 -ne -m Illumina,0.024,0.36,0.68,0.0097 -f fq -o L1000PEDeamin
#handle_error


echo "---------------------------------------------------------------------------------------------------------------"
echo "2) Testing Single-end misincorporation file, fa"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 1000 -t 1 -s 1 -l 100 -seq SE -ne -mf ${MF} -f fa -o MycoBactMFSE
#md5sum MycoBactMFSE.fa
handle_error

echo "---------------------------------------------------------------------------------------------------------------"
echo "2) Testing Single-end bdamage - cycle specific mismatch matrix from metagenomic dataset"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 10000 -t 1 -s 1 -l 100 -seq SE -ne -m3 ${M3} -qs 35 -f fq -o MycoBactM3SE
#md5sum MycoBactM3SE.fq
handle_error

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "------------------------- III) Testing Stochastic-, genetic- and reference variations -------------------------"
echo "---------------------------------------------------------------------------------------------------------------"

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "4) Testing Single-end, simulating fixed quality score of 30, with fragment lower-limit"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 1000 -t 1 -s 1 -lf ${LF} -seq SE -ll 50 -qs 30 -f fq -o MycoBactQSLLSEOUT
#md5sum MycoBactQSLLSEOUT.fq

handle_error

No_SeqErr=$(cat MycoBactQSLLSEOUT.fq|grep 'mod0001'|wc -l)
if [ $No_SeqErr -ne 62 ]; then 
    echo "Warning different number of reads containing sequencing error with a fixed quality score of 30"; exit 1;
fi

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "5) Testing mutation rate (reference genome 1.9%), sampling with replacement when below lower limit (30), 
        adapters and poly G for fixed sequencing error "
echo "---------------------------------------------------------------------------------------------------------------"

${PRG} -i ${IN} -r 1000 -f fq.gz -seq PE -s 34532 -t 1 -ld Gam,20,2 -mr 0.019 -cl 100 -qs 30 -o MutationRateSeqErr
handle_error
#md5sum MutationRateSeqErr_R1.fq.gz
#md5sum MutationRateSeqErr_R2.fq.gz

md5sum -c TestSub.md5 || exit 2;
