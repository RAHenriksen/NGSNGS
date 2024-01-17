#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "RandSampling.h"
#include "mrand.h"
#include "NGSNGS_misc.h"

#define LENS 10000
#define MAXBINS 100
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

void ReversComplement(char seq[]){
  // generates the reverse complementary sequence from an input sequence
  char seq_intermediate[LENS];
  strcpy(seq_intermediate,seq);
  int seqlen = strlen(seq);
  //Complementing sequence
  for(int i=0;i<seqlen;i++){
    seq_intermediate[i] = NtComp[refToInt[(unsigned char) seq_intermediate[i]]]; //warning: array subscript has type 'char' [-Wchar-subscripts]
  }
  //reverse complement
  for(int i=seqlen-1;i>-1;i--){
    seq[seqlen-i-1] = seq_intermediate[i];
  }
  //Clearing out the intermediate sequence
  memset(seq_intermediate, 0, sizeof seq_intermediate);
}

void ReversComplement2(char seq[],size_t fraglength){
  // generates the reverse complementary sequence from an input sequence
  char seq_intermediate[fraglength];
  strcpy(seq_intermediate,seq);
  int seqlen = strlen(seq);
  //Complementing sequence
  for(int i=0;i<seqlen;i++){
    seq_intermediate[i] = NtComp[refToInt[(unsigned char) seq_intermediate[i]]]; //warning: array subscript has type 'char' [-Wchar-subscripts]
  }
  //reverse complement
  for(int i=seqlen-1;i>-1;i--){
    seq[seqlen-i-1] = seq_intermediate[i];
  }
  //Clearing out the intermediate sequence
  memset(seq_intermediate, 0, sizeof seq_intermediate);
}

void Complement(char seq[],size_t fraglength){
  // generates the complementary sequence from an input sequence
  char seq_intermediate[fraglength];
  strcpy(seq_intermediate,seq);
  int seqlen = strlen(seq);
  //Complementing sequence
  for(int i=0;i<seqlen;i++){
    seq_intermediate[i] = NtComp[refToInt[(unsigned char) seq_intermediate[i]]]; //warning: array subscript has type 'char' [-Wchar-subscripts]
    seq[i] = seq_intermediate[i];
  }

  //Clearing out the intermediate sequence
  memset(seq_intermediate, 0, sizeof seq_intermediate);
}

void reverseChar(char* str,int length) {
  std::reverse(str, str + length);
}
