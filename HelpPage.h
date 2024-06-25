#ifndef HELPPAGE_H
#define HELPPAGE_H
#include <cstdio>
// -h or -v
int HelpPage(FILE *fp);

int AmpliconHelpPage(FILE *fp);

void ErrMsg(double messageno);

// general warnings
void WarMsg(double messageno);

float myatof(char *str);

void Sizebreak(char *str);

#endif
