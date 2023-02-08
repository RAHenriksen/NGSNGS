#ifndef NGSNGS_MISC_H
#define NGSNGS_MISC_H
#include "mrand.h"
#include "RandSampling.h"

void ReversComplement(char seq[]);

void ReversComplement2(char seq[],size_t fraglength);

void Complement(char seq[],size_t fraglength);

void reverseChar(char* str,int length);

#endif
