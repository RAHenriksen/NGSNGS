#rm *.bam *.sam *.fa.gz *.fq.gz *.fq *.fa 
PRG=../ngsngs
IN=../Test_Examples/Mycobacterium_leprae.fa.gz
LF=../Test_Examples/Size_dist/Size_dist_sampling.txt
Q1=../Test_Examples/Qual_profiles/AccFreqL150R1.txt
Q2=../Test_Examples/Qual_profiles/AccFreqL150R2.txt
A1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG
A2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT
echo "--------------------"
echo "RUNNING NGS: '${PRG}'"
echo "--------------------"

echo "Testing Single-end, length file, reads, briggs, bam"
${PRG} -i ${IN} -r 1000 -t1 2 -s 1 -lf ${LF} -seq SE -b 0.024,0.36,0.68,0.0097 -q1 ${Q1} -f bam -o MycoBactBamSEOut
samtools sort MycoBactBamSEOut.bam -o MycoBactBamSEOutSort.bam
Unsort_line=$(samtools view MycoBactBamSEOut.bam|wc -l)
Sort_line=$(samtools view MycoBactBamSEOutSort.bam|wc -l)
if [ $Unsort_line -ne $Sort_line ]; then 
echo "Warning different number of lines before and after sorting"; exit 1;
fi
rm MycoBactBamSEOut.bam
##only do this once
#md5sum MycoBactBamSEOutSort.bam > MycoBactBam.md5
#md5sum -c MycoBactBam.md5 || exit 2;

echo "Testing Paired-end, fixed length, coverage, no sequence errors, adapters, fastq"
${PRG} -i ${IN} -c 1 -t1 1 -s 1 -l 100 -seq PE -ne -a1 ${A1} -a2 ${A2} -q1 ${Q1} -q2 ${Q2} -f fq -o MycoBactFqPEOut
#md5sum MycoBactFqPEOut_r1.fq > MycoBactFQ.md5
#md5sum MycoBactFqPEOut_r2.fq >> MycoBactFQ.md5
md5sum -c MycoBactFQ.md5 || exit 2;
