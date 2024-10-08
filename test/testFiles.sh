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
BEDIN=../Test_Examples/hg19MSubCapture.bed
VCFFILE=../Test_Examples/ChrMtSubSNPDiploid.vcf 

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

echo "---------------------------------------------------------------------------------------------------------------"
echo "--------------------------------- I) Testing file format (-f) and compression----------------------------------"
echo "---------------------------------------------------------------------------------------------------------------"
echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "1) Testing Single-end, length file, reads, sampling threads, sequencing error, sam,bam,cram"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 1000 -t 1 -s 5643 -lf ${LF} -seq SE -q1 ${Q1} -f sam -o MycoBactSamSEOut
handle_error
${PRG} -i ${IN} -r 1000 -t 1 -s 5643 -lf ${LF} -seq SE -q1 ${Q1} -f bam -o MycoBactBamSEOut
handle_error
${PRG} -i ${IN} -r 1000 -t 1 -s 5643 -lf ${LF} -seq SE -q1 ${Q1} -f cram -o MycoBactCramSEOut
handle_error

# Check file sizes and ordering
sam_size=$(stat -c%s MycoBactSamSEOut.sam)
bam_size=$(stat -c%s MycoBactBamSEOut.bam)
cram_size=$(stat -c%s MycoBactCramSEOut.cram)
echo " "
if [ $sam_size -gt $bam_size ] && [ $bam_size -gt $cram_size ]; then
    echo "File sizes are in correct order: sam ("${sam_size}") > bam ("${bam_size}") > cram ("${cram_size}")"
else
    echo "File sizes are in incorrect order: sam ("${sam_size}") & bam ("${bam_size}") & cram ("${cram_size}")"; exit 1;
fi

samtools view MycoBactSamSEOut.sam > SamSE.txt
handle_error
samtools view MycoBactBamSEOut.bam > BamSE.txt
handle_error

SamMD5=$(md5sum SamSE.txt|cut -f1 -d ' ')
BamMD5=$(md5sum BamSE.txt|cut -f1 -d ' ')

if [ "$SamMD5" != "$BamMD5" ]; then
    echo "The decompression, or mate information differs from sam to bam"; exit 1;
fi


echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "2) Testing sequence alignment map format specific information for paired-end information"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 1000 -t 1 -s 81 -lf ${LF} -seq PE -qs 40 -f bam -o MycoBactBamPEOut
handle_error
${PRG} -i ${IN} -r 1000 -t 1 -s 314 -lf ${LF} -seq PE -qs 40 -na -f bam -o MycoBactBamPEOutNa
handle_error

samtools view MycoBactBamPEOut.bam > BamPE.txt
handle_error
samtools view MycoBactBamPEOutNa.bam > BamPEna.txt
handle_error

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "3) Testing fastq compression"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 1000 -t 1 -s 22 -l 100 -seq SE -qs 20 -f fq -o MycoBactFqout
handle_error
${PRG} -i ${IN} -r 1000 -t 1 -s 22 -l 100 -seq SE -qs 20 -f fq.gz -o MycoBactFqGzout
handle_error
# Check file sizes and ordering
fq_size=$(stat -c%s MycoBactFqout.fq)
fqgz_size=$(stat -c%s MycoBactFqGzout.fq.gz)

echo " "
if [ $fq_size -gt $fqgz_size ]; then
    echo "Fastq file sizes are in correct order: fq ("${fq_size}") > fq.gz ("${fqgz_size}")"
else
    echo "Fastq file sizes are in incorrect order: fq ("${fq_size}") < fq.gz ("${fqgz_size}")"; exit 1;
fi

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "--------------------------------------- II) Testing simulation options ----------------------------------------"
echo "---------------------------------------------------------------------------------------------------------------"

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "1) Testing Paired-end, fixed length, coverage, no sequencing error, adapters, sampling threads, fq"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -c 1 -t 1 -s 1 -l 100 -seq PE -ne -a1 ${A1} -a2 ${A2} -q1 ${Q1} -q2 ${Q2} -f fq -o MycoBactFqPEOut
handle_error

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "3) Testing Single-end, fixed length, reads, sequencing error, adapters, poly-G tail"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 1000 -t 1 -s 1 -l 70 -seq SE -q1 ${Q1} -a1 ${A1} -p G -f fq -o MycoBactFqSEPolyOut
handle_error

echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "4) Testing Single-end, fixed length greater than error profile, reads, cycle length"
echo "---------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 1000 -t 1 -s 1 -l 200 -seq SE -q1 ${Q1} -cl 90 -f fq -o MycoBactFqSEClOut
handle_error

Mean_read_length=$(awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' MycoBactFqSEClOut.fq)
if [ $Mean_read_length -ne 90 ]; then 
echo "Warning the mean read length isn't in accorandance with the cycle length"; exit 1;
fi

echo "---------------------------------------------------------------------------------------------------------------"
echo "------------ IV) Testing MD tags to ensure accurate flags and position with- or without deamination -----------"
echo "---------------------------------------------------------------------------------------------------------------"

echo " "
echo "-----------------------------------------------------------------------------------------------------------------------"
echo "1) Testing single-end alignments with cycle length below fragment length with deamination and position in terms of MD:Z"
echo "-----------------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 1000 -t 1 -s 1 -l 1000 -seq SE -q1 ${Q1} -m Illumina,0.024,0.36,0.68,0.0097 -f bam -o L1000SEDeamin
handle_error

samtools sort -n L1000SEDeamin.bam -o L1000SEDeaminSort.bam;samtools calmd -@ 10 -r -b L1000SEDeaminSort.bam ${IN} > L1000SEDeaminSortMD.bam
handle_error

MDcheck1=$(samtools view L1000SEDeaminSortMD.bam |cut -f14|awk '{ md_index = index($0, "MD:Z:"); md_substr = substr($0, md_index + 5); md_length = length(md_substr); if (md_length > 60) print 0; else print 1}'|awk '{sum += $1} END {print sum}')

if [ $MDcheck1 -ne 1000 ]; then 
    echo "Warning issue with MD tag above 60 characters for some reads, suggesting wrong positions stored within bam with MD tag sum value " $MDcheck1; exit 1;
fi

samtools view L1000SEDeamin.bam > L1000SEDeamin.txt
handle_error

echo " "
echo "-----------------------------------------------------------------------------------------------------------------------"
echo "2) Testing paired-end alignments with cycle length above fragment length with deamination and position in terms of MD:Z"
echo "-----------------------------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 1000 -t 1 -s 1 -l 1000 -seq PE -q1 ${Q1} -q2 ${Q2} -m Illumina,0.024,0.36,0.68,0.0097 -f bam -o L1000PEDeamin
handle_error

samtools sort -n L1000PEDeamin.bam -o L1000PEDeaminSort.bam;samtools calmd -@ 10 -r -b L1000PEDeaminSort.bam ${IN} > L1000PEDeaminSortMD.bam
handle_error

MDcheck2=$(samtools view L1000PEDeaminSortMD.bam |cut -f14|awk '{ md_index = index($0, "MD:Z:"); md_substr = substr($0, md_index + 5); md_length = length(md_substr); if (md_length > 60) print 0; else print 1}'|awk '{sum += $1} END {print sum}')

if [ $MDcheck2 -ne 2000 ]; then 
    echo "Warning issue with MD tag above 60 characters for some reads, suggesting wrong positions stored within bam with MD tag sum value " $MDcheck2; exit 1;
fi

samtools view L1000PEDeamin.bam > L1000PEDeamin.txt
handle_error
echo " "
echo "------------------------------------------------------------------------------------------------"
echo "3) Testing MD:Z for single- and paired-end alignments without deamination and sequencing errors "
echo "------------------------------------------------------------------------------------------------"
${PRG} -i ${IN} -r 1000 -t 1 -s 200 -l 1000 -ne -seq SE -q1 ${Q1} -f bam -o L1000SE
handle_error

samtools sort -n L1000SE.bam -o L1000SESort.bam;samtools calmd -@ 10 -r -b L1000SESort.bam ${IN} > L1000SESortMD.bam
handle_error

#should only have one md tag equal to cycle length if positions are correct as we're not altering any nucleotides
MDcheck3=$(samtools view L1000SESortMD.bam |cut -f14|sort|uniq -c|wc -l)
if [ $MDcheck3 -ne 1 ]; then 
    echo "Warning issue with too many MD tags, suggesting wrong positions stored within bam with no nucleotide alterations" $MDcheck3; exit 1;
fi

samtools view L1000SESortMD.bam > L1000SESortMD.txt
handle_error

${PRG} -i ${IN} -r 1000 -t 1 -s 300 -l 1000 -ne -seq PE -q1 ${Q1} -q2 ${Q2} -f bam -o L1000PE
handle_error

samtools sort -n L1000PE.bam -o L1000PESort.bam;samtools calmd -@ 10 -r -b L1000PESort.bam ${IN} > L1000PESortMD.bam
handle_error

MDcheck4=$(samtools view L1000PESortMD.bam |cut -f14|sort|uniq -c|wc -l)
if [ $MDcheck4 -ne 1 ]; then 
    echo "Warning issue with too many MD tags, suggesting wrong positions stored within bam with no nucleotide alterations" $MDcheck4; exit 1;
fi

samtools view L1000PESortMD.bam > L1000PESortMD.txt
handle_error

echo "---------------------------------------------------------------------------------------------------------------"
echo "----------------------- V) Testing different simulation types, circular and capture ----------------------"
echo "---------------------------------------------------------------------------------------------------------------"

echo " "
echo "------------------------------------------------------------------------"
echo "1) Testing single-end capture simulations from provided bed file regions"
echo "------------------------------------------------------------------------"
${PRG} -i ${VCFIN} -r 1000 -t 1 -s 1 -l 100 -seq SE --circular -qs 35 -f fq -o TestMSubCircular
handle_error

awk '/@.*MT:[0-9]+-[0-9]+_length:[0-9]+/ {
    match($0, /MT:[0-9]+-([0-9]+)_length:[0-9]+/, arr)
    if (arr[1] > 500) {
        count++
    }
}
END {print count}' TestMSubCircular.fq &> BreakPointReadsMT.txt

echo " "
echo "------------------------------------------------------------------------"
echo "2) Testing single-end capture simulations from provided bed file regions"
echo "------------------------------------------------------------------------"
${PRG} -i ${VCFIN} -r 1000 -t 1 -s 1 -l 100 -seq SE -incl ${BEDIN} -f fa -o TestMSubCapture
handle_error

ranges=$(grep '>' TestMSubCapture.fa | cut -f4 -d_ | cut -f2 -d:)

awk -v ranges="$ranges" '
BEGIN {
    # Split the input ranges into an array
    split(ranges, arr, " ");
    for (i in arr) {
        split(arr[i], range, "-");
        start[i] = range[1];
        end[i] = range[2];
    }
    total_overlap_count = 0;
    total_non_overlap_count = 0;
    overlaps = "";
    non_overlaps = "";
    seen_non_overlap_reads = "";
}
{
    # Read start and end positions from the BED file
    bed_start = $2;
    bed_end = $3;
    overlap_found = 0;
    overlap_not_found = 0;

    # Check for overlap with each range
    for (i in start) {
        if (start[i] < bed_end && end[i] > bed_start) {
            overlap_found++;
        }
        else{
            overlap_not_found++;
        }
    }
    
    # If overlap is found
    if (overlap_found > 0) {
        total_overlap_count += overlap_found;
        overlaps = overlaps $1 "\t" bed_start "\t" bed_end "\n";
    } else {
        # If no overlap is found
        total_non_overlap_count += overlap_not_found;
        non_overlaps = non_overlaps $1 "\t" bed_start "\t" bed_end "\n";
    }
}
END {
    # Print the total number of overlaps and non-overlapping reads
    print "Total overlapping reads found for capture:" total_overlap_count > "overlaps.txt";
    print "Total non-overlapping reads found for capture:" total_non_overlap_count > "nonoverlaps.txt";
}' ${BEDIN}

echo " "
echo "---------------------------------------------------------------------------"
echo "3) Testing single-end capture simulations masking provided bed file regions"
echo "---------------------------------------------------------------------------"
${PRG} -i ${VCFIN} -r 1000 -t 1 -s 1 -l 100 -seq SE -excl ${BEDIN} -f fa -o TestMSubMaskCapture
handle_error

ranges=$(grep '>' TestMSubMaskCapture.fa | cut -f4 -d_ | cut -f2 -d:)

awk -v ranges="$ranges" '
BEGIN {
    # Split the input ranges into an array
    split(ranges, arr, " ");
    for (i in arr) {
        split(arr[i], range, "-");
        start[i] = range[1];
        end[i] = range[2];
    }
    total_overlap_count = 0;
    total_non_overlap_count = 0;
    overlaps = "";
    non_overlaps = "";
    seen_non_overlap_reads = "";
}
{
    # Read start and end positions from the BED file
    bed_start = $2;
    bed_end = $3;
    overlap_found = 0;
    overlap_not_found = 0;

    # Check for overlap with each range
    for (i in start) {
        if (start[i] < bed_end && end[i] > bed_start) {
            overlap_found++;
        }
        else{
            overlap_not_found++;
        }
    }
    
    # If overlap is found
    if (overlap_found > 0) {
        total_overlap_count += overlap_found;
        overlaps = overlaps $1 "\t" bed_start "\t" bed_end "\n";
    } else {
        # If no overlap is found
        total_non_overlap_count += overlap_not_found;
        non_overlaps = non_overlaps $1 "\t" bed_start "\t" bed_end "\n";
    }
}
END {
    # Print the total number of overlaps and non-overlapping reads
    print "Total overlapping reads found for masked capture:" total_overlap_count >> "overlaps.txt";
    print "Total non-overlapping reads found for masked capture:" total_non_overlap_count >> "nonoverlaps.txt";
}' ${BEDIN}

echo " "
echo "---------------------------------------------------------------------"
echo "4) Testing single-end capture simulations of variations from vcf file"
echo "---------------------------------------------------------------------"
${PRG} -i ${VCFIN} -r 1000 -s 1 -l 100 -seq SE -vcf ${VCFFILE} -id 0 -f fq -qs 40 -ne -fl 200 --capture -o DiploidCapture
handle_error

echo " "
echo "------------------------------------------------------------------------------------"
echo "5) Testing single-end linkage disequilibrium simulations of variations from vcf file"
echo "------------------------------------------------------------------------------------"
${PRG} -i ${VCFIN} -r 1000 -s 1 -l 100 -seq SE -vcf ${VCFFILE} -t 1 -id 2 -f fa -fl 200 -linkage -o DiploidLinkage
handle_error
echo " "
echo "---------------------------------------------------------------------------------------------------------------"
echo "--------------------------------------------------- MD5SUM ----------------------------------------------------"
echo "---------------------------------------------------------------------------------------------------------------"
echo " "
md5sum -c TestFiles.md5 || exit 2;
