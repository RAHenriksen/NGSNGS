#ifndef SIM_FRAGMENT_H
#define SIM_FRAGMENT_H

#include "mrand.h"

typedef struct{
  int type;
  int FixLength;
  int LengthDist;
  int RandType;
  int Thread_Seed;
  mrand_t *rand_alloc;
  int noRow;
  double* Frequency;
  int* Frag_len;
  std::default_random_engine Gen;
  std::uniform_int_distribution<int> UniDist;
  std::normal_distribution<double> NormDist;
  std::lognormal_distribution<double> LogNormDist;
  std::poisson_distribution<int> PoisDist;
  std::exponential_distribution<double> ExpDist;
  std::gamma_distribution<double> GammaDist;
}sim_fragment;

sim_fragment *sim_fragment_alloc(int type,double par1, double par2,int no_row,double*& FreqArray,int*& FragArray,
int RandType,unsigned int Thread_Seed,std::default_random_engine& generator);
void ReadLengthFile(int& number,int*& Length, double*& Frequency,const char* filename);

int getFragmentLength(sim_fragment *sf);
int BinarySearch_fraglength(double* SearchArray,int low, int high, double key);

#endif