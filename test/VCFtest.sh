#rm *.bam *.sam *.fa.gz *.fq.gz *.fq *.fa 
VCFDIR=../Test_Examples/VCF
PRG=../ngsngs
IN=../Test_Examples/VCF/MTSub.fa
Q1=../Test_Examples/AccFreqL150R1.txt
Q2=../Test_Examples/AccFreqL150R2.txt

echo "---------------------------------------------------------------------------------------------------------------"
echo "------------------------------- Testing Haploid vcf files, for three indivdiuals ------------------------------"
echo "---------------------------------------------------------------------------------------------------------------"

for file in $(ls ${VCFDIR}/*Haploid*); do for indiv in {0,1,2}; do 
    echo ${file} indiv ${indiv}; 
    Type=$(ls ${file}| sed 's/.*Sub//'|sed 's/Haploid.vcf//');
    ${PRG} -i ${IN} -r 100 -t 1 -s 1 -l 150 -seq SE -ne -vcf ${file} -id ${indiv} -q1 ${Q1} -chr MT -—dump-internal VCF_temp_haploid_${Type}_${indiv} -f fq -o VCFtest;
    md5sum VCF_temp_haploid_${Type}_${indiv}.fa >> VCF.md5; done; done

echo "---------------------------------------------------------------------------------------------------------------"
echo "------------------------------- Testing Diploid vcf files, for three indivdiuals ------------------------------"
echo "---------------------------------------------------------------------------------------------------------------"

for file in $(ls ${VCFDIR}/*Diploid*); do for indiv in {0,1,2}; do 
    echo ${file} indiv ${indiv}; 
    Type=$(ls ${file}| sed 's/.*Sub//'|sed 's/Diploid.vcf//');
    ${PRG} -i ${IN} -r 100 -t 1 -s 1 -l 150 -seq SE -ne -vcf ${file} -id ${indiv} -q1 ${Q1} -chr MT -—dump-internal VCF_temp_diploid_${Type}_${indiv} -f fq -o VCFtest;
    md5sum VCF_temp_diploid_${Type}_${indiv}.fa >> VCF.md5; done; done

md5sum -c VCF.md5 || exit 2;
rm *.fq
rm *.fa