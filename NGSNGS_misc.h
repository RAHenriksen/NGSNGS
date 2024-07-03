#ifndef NGSNGS_MISC_H
#define NGSNGS_MISC_H
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>

void ReversComplement(char* seq);

void Complement(char* seq);

void reverseChar(char* str,int length);

void Complement_k(kstring_t* seq);

void ReversComplement_k(kstring_t* seq);

void CreateSeqQualKString(bam1_t *aln, kstring_t *Sequence, kstring_t *Quality,int offset);

char* PrintCigarBamSet1(size_t n_cigar,const uint32_t *cigar);

#endif
