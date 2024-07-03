#include <algorithm>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <iostream>

#include "NGSNGS_misc.h"
#include <htslib/kstring.h>

int refToInt[256] = {
  0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//15
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//31
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//47
  0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//63
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//79
  4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//95
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//111
  4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//127
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//143
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//159
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//175
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//191
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//207
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//223
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//239
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4//255
};

char intToRef[4] = {'A','C','G','T'};

char refToChar[256] = {
    0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//15
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//31
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//47
    0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//63
    4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//79
    4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//95
    4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//111
    4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//127
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//143
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//159
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//175
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//191
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//207
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//223
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//239
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4//255
};

char NtComp[5] = {'T', 'G', 'C', 'A','N'};

const char *bass = "ACGTN";

void ReversComplement(char* seq) {
  // Reverse and complement input sequence in place
  reverseChar(seq, strlen(seq));
  Complement(seq);
}

void Complement(char* seq) {
  // Complemen input sequence in place
  for (; *seq; ++seq) {
    *seq = NtComp[refToInt[(unsigned char)*seq]];
  }
}

void reverseChar(char* str,int length) {
    std::reverse(str, str + length);
}

void Complement_k(kstring_t* seq) {
  // Complement input sequence in place
  char *s = seq->s;
  for (; *s; ++s) {
    *s = NtComp[refToInt[(unsigned char)*s]];
  }
}

void ReversComplement_k(kstring_t* seq) {
  // Reverse and complement input sequence in place
  reverseChar(seq->s, seq->l);
  Complement_k(seq);
}

char* PrintCigarBamSet1(size_t n_cigar,const uint32_t *cigar){
   // Allocate a buffer to store the CIGAR string
  char *cigar_str = (char *)malloc(n_cigar * 2); // a safe overestimate
  char *ptr = cigar_str;

  /*size_t l_qname, const char *qname,
  uint16_t flag, int32_t tid, hts_pos_t pos, uint8_t mapq,
  size_t n_cigar, const uint32_t *cigar,
  int32_t mtid, hts_pos_t mpos, hts_pos_t isize,
  size_t l_seq, const char *seq, const char *qual,
  size_t l_aux)*/

  // Convert CIGAR to string
  for (int i = 0; i < n_cigar; ++i) {
    int len = sprintf(ptr, "%d%c", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i]));
    ptr += len;
  }

  return cigar_str;
}

void CreateSeqQualKString(bam1_t *aln, kstring_t *Sequence, kstring_t *Quality, int offset){
  int len = aln->core.l_qseq;

  Sequence->s = (char *) realloc(Sequence->s, len + 1); // allocate memory for the sequence string (+1 for null terminator)
  Sequence->l = len;
  Sequence->m = len + 1;

  Quality->s = (char *) realloc(Quality->s, len + 1); // allocate memory for the sequence string (+1 for null terminator)
  Quality->l = len;
  Quality->m = len + 1;

  uint8_t *seq_data = bam_get_seq(aln);
  uint8_t *qual_data = bam_get_qual(aln);

  for (int i = 0; i < len; i++) {
    Sequence->s[i] = seq_nt16_str[bam_seqi(seq_data, i)]; // convert nucleotide id to IUPAC id
    Quality->s[i] = qual_data[i]+offset;
  }
  
  Sequence->s[len] = '\0'; // null-terminate the string
  Quality->s[len] = '\0';
}



#ifdef __WITH_MAIN__
//g++ NGSNGS_misc.cpp -D__WITH_MAIN__
int main(){
  // Example usage of kstring_t and Complement_k
  kstring_t seq;
  char sequence[] = "ACGTN";
  seq.s = sequence;
  seq.l = strlen(sequence);
  seq.m = seq.l + 1; // Just for safety, usually seq.m is set to the allocated size

  std::cout << "Original sequence: " << seq.s << std::endl;

  Complement_k(&seq);

  std::cout << "Complemented sequence: " << seq.s << std::endl;
  
  ReversComplement_k(&seq);

  std::cout << "Reverse complemented sequence: " << seq.s << std::endl;

  reverseChar(seq.s, seq.l);

  std::cout << "Reverse sequence: " << seq.s << std::endl;

  return 0;
}

#endif