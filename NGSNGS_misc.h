#ifndef NGSNGS_MISC_H
#define NGSNGS_MISC_H
#include <htslib/kstring.h>
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>

#define MAX_CHROM_NAME_LEN 100

typedef struct {
    char chromosome[MAX_CHROM_NAME_LEN];
    int start;
    int end;
} BedEntry;

void ReversComplement(char* seq);

void Complement(char* seq);

void reverseChar(char* str,int length);

void Complement_k(kstring_t* seq);

void ReversComplement_k(kstring_t* seq);

void CreateSeqQualKString(bam1_t *aln, kstring_t *Sequence, kstring_t *Quality,int offset);

char* PrintCigarBamSet1(size_t n_cigar,const uint32_t *cigar);

// Function to read BED file and store entries in an array of BedEntry
BedEntry* readBedFile(const char* filename, int* entryCount);

// create internal structure of a bed file from the vcf to create region from which we sample
int VCFtoBED(const char* bcffilename, int id,int range,BedEntry** bedentries, int* entryCount);

#endif
