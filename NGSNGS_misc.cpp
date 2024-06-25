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