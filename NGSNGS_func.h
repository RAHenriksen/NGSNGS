#ifndef SIMULANCIENT_FUNC_H
#define SIMULANCIENT_FUNC_H

void DNA_complement(char seq[]);

void reverseChar(char* str);

int Random_geometric_k(unsigned int  seed, const double p);

double myrand(unsigned int persistent);

void SimBriggsModel(char* reffrag, char* frag, int L, double nv, double lambda, double delta_s, double delta, unsigned int seed);

const char* Error_lookup(double a,double err[6000],int nt_offset, int read_pos);

double* Qual_array(double* freqval,const char* filename);

void Read_Qual_new(char *seq,char *qual,unsigned int seed,double* freqval);

int BinarySearch_fraglength(double* SearchArray,int low, int high, double key);

void FragArray(int& number,int*& Length, double*& Frequency,const char* filename);



#endif