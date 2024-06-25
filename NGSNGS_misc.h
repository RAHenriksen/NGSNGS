#ifndef NGSNGS_MISC_H
#define NGSNGS_MISC_H
#include <htslib/kstring.h>

void ReversComplement(char* seq);

void Complement(char* seq);

void reverseChar(char* str,int length);

void Complement_k(kstring_t* seq);

void ReversComplement_k(kstring_t* seq);

#endif
