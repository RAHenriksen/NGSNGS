#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <zlib.h>
#include <random>
#include <iostream>
#include <cassert>
#include <stdlib.h>

#include "getFragmentLength.h"
#include "mrand.h"

#define LENS 4096

int BinarySearch_fraglength(double* SearchArray,int low, int high, double key){
  int bin_idx = 0; 
  while (low <= high) {
    int mid = low + (high - low + 1) / 2;
    double midVal = SearchArray[mid];
    if (midVal < key) {bin_idx = mid;low = mid + 1;}
    else if (midVal > key) {high = mid - 1;}
    else if (midVal == key){high = mid - 1;}
  }
  return bin_idx+1;
}

void ReadLengthFile(int& number,int*& Length, double*& Frequency,const char* filename){
  int n =1;

  gzFile gz = Z_NULL;
  char buf[LENS];
  gz = gzopen(filename,"r");
  assert(gz!=Z_NULL);
  Length[0] = 0; Frequency[0] = (float) 0;
  while(gzgets(gz,buf,LENS)){
    Length[n] = atoi(strtok(buf,"\n\t "));
    Frequency[n] = atof(strtok(NULL,"\n\t "));
    n++;
  }
  gzclose(gz); 
  number = n;
}

int getFragmentLength(sim_fragment *sf){
  int res = 0;
  double rand_val;

  if(sf->type == 0){
    res = sf->FixLength;
  }
  else if(sf->type == 1){
    //fprintf(stderr,"type 1 lengthf file\n");
    rand_val = mrand_pop(sf->rand_alloc);
    int lengthbin = BinarySearch_fraglength(sf->Frequency,0, sf->noRow - 1, rand_val);
    int fraglength =  sf->Frag_len[lengthbin];
    res = fraglength;
  }
  else if(sf->type == 2){
    //fprintf(stderr,"type 2 uni\n");
    res = sf->UniDist(sf->Gen);//sf->LengthDist;
  }
  else if(sf->type == 3){
    //fprintf(stderr,"type 3 norm\n");
    res = sf->NormDist(sf->Gen);//sf->LengthDist;
  }
  else if(sf->type == 4){
    //fprintf(stderr,"type 4 log\n");
    res = sf->LogNormDist(sf->Gen);//sf->LengthDist;
  }
  else if(sf->type == 5){
    //fprintf(stderr,"type 5 pois\n");
    res = sf->PoisDist(sf->Gen);//sf->LengthDist;
  }
  else if(sf->type == 6){
    //fprintf(stderr,"type 6 exp\n");
    res = sf->ExpDist(sf->Gen);//sf->LengthDist;
  }
  else if(sf->type == 7){
    //fprintf(stderr,"type 7 gamma\n");
    res = sf->GammaDist(sf->Gen);//sf->LengthDist;
  }
  if(res<0)
    return getFragmentLength(sf);
  return res;
}

sim_fragment *sim_fragment_alloc(int type,double par1, double par2,int no_row,double*& FreqArray,
				 int*& FragArray,int RandType,unsigned int Thread_Seed,std::default_random_engine& generator){
  sim_fragment *fp = new sim_fragment;

  fp->rand_alloc= mrand_alloc(RandType,Thread_Seed);
  fp->type = type;
  fp->Gen =generator;

  if(fp->type == 0){
    fp->FixLength = par1;
  }
  else if(fp->type == 1){
    fp->Frequency = FreqArray;
    fp->Frag_len = FragArray;
    fp->noRow = no_row;
  }
  else if(fp->type == 2){
    std::uniform_int_distribution<int> distribution(par1,par2);
    fp->UniDist = distribution;
  }
  else if(fp->type == 3){
    std::normal_distribution<double> distribution(par1,par2);
    fp->NormDist = distribution;
  }
  else if(fp->type == 4){
    std::lognormal_distribution<double> distribution(par1,par2);
    fp->LogNormDist = distribution;
  }
  else if(fp->type == 5){
    std::poisson_distribution<int> distribution(par1);
    fp->PoisDist = distribution;
  }
  else if(fp->type == 6){
    std::exponential_distribution<double> distribution(par1);
    fp->ExpDist = distribution;
  }
  else if(fp->type == 7){
    std::gamma_distribution<double> distribution(par1,par2);
    fp->GammaDist = distribution;
  }
  return fp;
}


#ifdef __WITH_MAIN__

int main(int argc,char **argv){
  if(argc==1)
    fprintf(stderr,"./a.out dist par1 par2 nrep");

  int type = atoi(argv[1]);
  double par1 = atof(argv[2]);
  double par2 = atof(argv[3]);
  int nrep = atoi(argv[4]);
  const char *fname = "Test_Examples/Anc_size.txt";

  char buf[LENS];
  double* Frag_freq;
  int* Frag_len;
  int no_elem;
  Frag_len = new int[LENS];Frag_freq = new double[LENS];
  ReadLengthFile(no_elem,Frag_len,Frag_freq,fname);
  
  int Thread_Seed = 2;
  std::default_random_engine RndGen(Thread_Seed);
  sim_fragment *sf = sim_fragment_alloc(type,par1,par2,0,Frag_freq,Frag_len,0,Thread_Seed,RndGen);
  
  for(int i=0;i<nrep;i++)
    fprintf(stdout,"%d\n",getFragmentLength(sf));
  
  return 0;
}
#endif

//g++ getFragmentLength.cpp mrand.o -std=c++11 -lm -lz -D__WITH_MAIN__ -O3 -o FragLen
