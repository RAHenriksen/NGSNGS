#rm *.bam *.sam *.fa.gz *.fq.gz *.fq *.fa 
PRG=../ngsngs
IN=../Test_Examples/Mycobacterium_leprae.fa.gz
LF=../Test_Examples/Size_dist/Size_dist_sampling.txt
Q1=../Test_Examples/Qual_profiles/AccFreqL150R1.txt
Q2=../Test_Examples/Qual_profiles/AccFreqL150R2.txt
A1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG
A2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT
MF=../Test_Examples/MisincorpFile.txt

echo "---------------------------------------------------------------------------------------------------------------"
echo "---------------------------------------- Testing '${PRG}' examples -----------------------------------------"
echo "---------------------------------------------------------------------------------------------------------------"
echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "1) Testing Single-end, length file, reads, briggs, sampling threads, bam"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t1 2 -s 1 -lf ${LF} -seq SE -b 0.024,0.36,0.68,0.0097 -q1 ${Q1} -f bam -o MycoBactBamSEOut
samtools sort MycoBactBamSEOut.bam -o MycoBactBamSEOutSort.bam
Unsort_line=$(samtools view MycoBactBamSEOut.bam|wc -l)
Sort_line=$(samtools view MycoBactBamSEOutSort.bam|wc -l)
if [ $Unsort_line -ne $Sort_line ]; then 
echo "Warning different number of lines before and after sorting"; exit 1;
fi
rm MycoBactBamSEOut.bam
md5sum MycoBactBamSEOutSort.bam > MycoBactTest.md5
echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "2) Testing Paired-end, fixed length, coverage, no sequencing error, adapters, sampling threads, fq"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -c 3 -t1 2 -s 1 -l 100 -seq PE -ne -a1 ${A1} -a2 ${A2} -q1 ${Q1} -q2 ${Q2} -f fq -o MycoBactFqPEOut
sort MycoBactFqPEOut_R1.fq > MycoBactFqPEOut_R1_sort.fq
sort MycoBactFqPEOut_R2.fq > MycoBactFqPEOut_R2_sort.fq
rm MycoBactFqPEOut_R1.fq
rm MycoBactFqPEOut_R2.fq

md5sum MycoBactFqPEOut_R1_sort.fq >> MycoBactTest.md5
md5sum MycoBactFqPEOut_R2_sort.fq >> MycoBactTest.md5

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "3) Testing Single-end, length distributions, reads, misincorporation file, fa"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t1 1 -s 1 -ld Pois,78 -seq SE -mf ${MF} -f fa -o MycoBactFaSEOut
md5sum MycoBactFaSEOut.fa >> MycoBactTest.md5
echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "----------------------------------------- Additional '${PRG}' Test -----------------------------------------"
echo "---------------------------------------------------------------------------------------------------------------"
echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "4) Testing Single-end, fixed length, reads, sequencing error, adapters, poly-G tail"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t1 1 -s 1 -l 70 -seq SE -q1 ${Q1} -a1 ${A1} -p G -f fq -o MycoBactFqSEPolyOut
md5sum MycoBactFqSEPolyOut.fq >> MycoBactTest.md5

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "--------------------------------------------------- MD5SUM ----------------------------------------------------"
echo "---------------------------------------------------------------------------------------------------------------"
echo " "
### SAM + BCF TEST
#md5sum MycoBactFqPEOut_R1.fq >> MycoBactTest.md5
#md5sum MycoBactFqPEOut_R2.fq >> MycoBactFQ.md5
md5sum -c MycoBactTest.md5 || exit 2;
