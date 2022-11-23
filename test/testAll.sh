#rm *.bam *.sam *.fa.gz *.fq.gz *.fq *.fa 
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
echo "1) Testing Single-end, length file, reads, briggs, sampling threads, sequencing error, bam"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t 1 -s 1 -lf ${LF} -seq SE -m b7,0.024,0.36,0.68,0.0097 -q1 ${Q1} -f bam -o MycoBactBamSEOut

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
${PRG} -i ${IN} -c 3 -t 1 -s 1 -l 100 -seq PE -ne -a1 ${A1} -a2 ${A2} -q1 ${Q1} -q2 ${Q2} -f fq -o MycoBactFqPEOut
sort MycoBactFqPEOut_R1.fq > MycoBactFqPEOut_R1_sort.fq
sort MycoBactFqPEOut_R2.fq > MycoBactFqPEOut_R2_sort.fq
rm MycoBactFqPEOut_R1.fq
rm MycoBactFqPEOut_R2.fq

md5sum MycoBactFqPEOut_R1_sort.fq >> MycoBactTest.md5
md5sum MycoBactFqPEOut_R2_sort.fq >> MycoBactTest.md5

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "3) Testing Single-end, length distributions, reads, sequencing error, misincorporation file, fa"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t 1 -s 1 -ld Pois,78 -seq SE -mf ${MF} -f fa -o MycoBactFaSEOut
md5sum MycoBactFaSEOut.fa >> MycoBactTest.md5

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "4) Testing Single-end, fixed length, reads, sequencing error, adapters, poly-G tail"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t 1 -s 1 -l 70 -seq SE -q1 ${Q1} -a1 ${A1} -p G -f fq -o MycoBactFqSEPolyOut
md5sum MycoBactFqSEPolyOut.fq >> MycoBactTest.md5

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "5) Testing Single-end, fixed length greater than error profile, reads, cycle length"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 100000 -t 1 -s 1 -l 200 -seq SE -q1 ${Q1} -cl 90 -f fq -o MycoBactFqSEClOut

Mean_read_length=$(awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' MycoBactFqSEClOut.fq)
if [ $Mean_read_length -ne 90 ]; then 
echo "Warning the mean read length isn't in accorandance with the cycle length"; exit 1;
fi
md5sum MycoBactFqSEClOut.fq >> MycoBactTest.md5

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "6) Testing Single-end, simulating genetic variations, SNP, Insertion and deletions from Haploid vcf files, 
         for three indivdiuals with different genotype, homozygous for reference and alternative and heterozygous "
echo "---------------------------------------------------------------------------------------------------------------"

for file in $(ls ${VCFDIR}/*Haploid*); do for indiv in {0,1,2}; do 
    echo ${file} indiv ${indiv}; 
    Type=$(ls ${file}| sed 's/.*Sub//'|sed 's/Haploid.vcf//');
    ${PRG} -i ${VCFIN} -r 100 -s 1 -l 80 -seq SE -ne -vcf ${file} -id ${indiv} -q1 ${Q1} -chr MT -—dump-internal VCF_temp_haploid_${Type}_${indiv} -f fq -o VCF_${Type}_${indiv};
    md5sum VCF_${Type}_${indiv}.fq >> MycoBactTest.md5; done; done

echo "---------------------------------------------------------------------------------------------------------------"
echo "7) Testing Single-end, simulating genetic variations, SNP, Insertion and deletions from Diploid vcf files, 
         for three indivdiuals with different genotype, homozygous for reference and alternative and heterozygous "
echo "---------------------------------------------------------------------------------------------------------------"

for file in $(ls ${VCFDIR}/*Diploid*); do for indiv in {0,1,2}; do 
    echo ${file} indiv ${indiv}; 
    Type=$(ls ${file}| sed 's/.*Sub//'|sed 's/Diploid.vcf//');
    ${PRG} -i ${VCFIN} -r 100 -t 1 -s 1 -l 150 -seq SE -ne -vcf ${file} -id ${indiv} -q1 ${Q1} -chr MT -—dump-internal VCF_temp_diploid_${Type}_${indiv} -f fq -o VCF_diploid_${Type}_${indiv};
    md5sum VCF_diploid_${Type}_${indiv}.fq >> MycoBactTest.md5; done; done

echo "---------------------------------------------------------------------------------------------------------------"
echo "8) Testing Single-end, simulating stochastic variations, insertion or deletions, no error, length file "
echo "---------------------------------------------------------------------------------------------------------------"

${PRG} -i ${IN} -r 1000 -t 1 -s 1 -lf ${LF} -seq SE -ne -indel 0.0,0.05,0.0,0.9 -q1 ${Q1} -—dump-indel DelTmp -f fq -o DelOut
${PRG} -i ${IN} -r 1000 -t 1 -s 1 -lf ${LF} -seq SE -ne -indel 0.05,0.0,0.9,0.0 -q1 ${Q1} -—dump-indel InsTmp -f fq -o InsOut

md5sum DelOut.fq >> MycoBactTest.md5
md5sum InsOut.fq >> MycoBactTest.md5

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "--------------------------------------------------- MD5SUM ----------------------------------------------------"
echo "---------------------------------------------------------------------------------------------------------------"
echo " "
md5sum -c MycoBactTest.md5 || exit 2;
rm DelTmp.txt
rm InsTmp.txt
rm *.fq
rm *.fa
