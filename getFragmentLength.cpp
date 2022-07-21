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
    //fprintf(stderr,"first element %lf\n",SearchArray[low]);
    int ans = 0; 
    while (low <= high) {
        int mid = low + (high - low + 1) / 2;
        //fprintf(stderr,"test %lf\n",SearchArray[mid]);
        double midVal = SearchArray[mid];
 
        if (midVal < key) {
            ans = mid;
            low = mid + 1;
        }
        else if (midVal > key) {

            high = mid - 1;
        }
        else if (midVal == key) {
 
            high = mid - 1;
        }
    }
 
    return ans+1;
}

void ReadLengthFile(int& number,int*& Length, double*& Frequency,const char* filename){
  int n =1;

  gzFile gz = Z_NULL;
  char buf[LENS];
  gz = gzopen(filename,"r");
  assert(gz!=Z_NULL);
  Length[0] = 0; Frequency[0] = (float) 0;
  while(gzgets(gz,buf,LENS)){
    Length[n] = atoi(strtok(buf,"\n\t ")); //before it was Frag_len[n]
    Frequency[n] = atof(strtok(NULL,"\n\t "));
    n++;
  }
  gzclose(gz); 
  number = n;
}

int getFragmentLength(sim_fragment *sf){
  int res;
  double rand_val;

  if(sf->type == 0){
    //std::cout << " before "<< sf->FixLength << std::endl;
    res = sf->FixLength;
    //std::cout << " after " << res << std::endl;
  }
  else if(sf->type == 1){
    rand_val = mrand_pop(sf->rand_alloc);
    int lengthbin = BinarySearch_fraglength(sf->Frequency,0, sf->noRow - 1, rand_val);
    int fraglength =  sf->Frag_len[lengthbin];
    std::cout << lengthbin << " " << fraglength << std::endl;
    res = fraglength;
    std::cout << rand_val << " " << res << std::endl;
  }
  else if(sf->type == 2){
    res = sf->UniDist(sf->Gen);//sf->LengthDist;
  }
  return res;
}

sim_fragment *sim_fragment_alloc(int type,double par1, double par2,int no_row,double*& FreqArray,
int*& FragArray,int RandType,int Thread_Seed,std::default_random_engine& generator){
  sim_fragment *fp = new sim_fragment;

  fp->rand_alloc= mrand_alloc(RandType,Thread_Seed);
  fp->type = type;
  fp->Gen =generator;
  
  int res;
  if(fp->type == 0){
    //std::cout << " fp type before " << fp->type << " " << par1 << std::endl;
    fp->FixLength = par1;
    //std::cout << "fp fixlength " << fp->FixLength << std::endl;
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
    fp->UniDist = distribution;
  }
  else if(fp->type == 4){
    std::uniform_int_distribution<int> distribution(par1,par2);
    fp->UniDist = distribution;
  }
  else if(fp->type == 5){
    std::uniform_int_distribution<int> distribution(par1,par2);
    fp->UniDist = distribution;
  }
  else if(fp->type == 6){
    std::uniform_int_distribution<int> distribution(par1,par2);
    fp->UniDist = distribution;
  }
  else if(fp->type == 7){
    std::uniform_int_distribution<int> distribution(par1,par2);
    fp->UniDist = distribution;
  }
  return fp;
}


#ifdef __WITH_MAIN__

int main(int argc,char **argv){
  if(argc==1)
  fprintf(stderr,"./a.out dist par1 par2 nrep")
  int nrep = 10000000;

  int type = argv[1];
  double par1 = argv[2];
  double par2 = argv[3];
  double nrep = argv[3];
  const char *fname = "Test_Examples/Anc_size.txt";

  char buf[LENS];
  double* Frag_freq;
  int* Frag_len;
  int no_elem;
  Frag_len = new int[LENS];Frag_freq = new double[LENS];
  ReadLengthFile(no_elem,Frag_len,Frag_freq,fname);
  

  int Thread_Seed = 2;
  std::default_random_engine RndGen(Thread_Seed);
  sim_fragment *sf = sim_fragment_alloc(type,par1,par2,no_elem,Frag_freq,Frag_len,0,Thread_Seed,RndGen);

  
  for(i=0;i<nrep;i++)
    fprintf(stdout,"%d\n",getFragmentLength(sf));
  
  return 0;
}
#endif

//g++ getFragmentLength.cpp mrand.o -std=c++11 -lm -lz -D__WITH_MAIN__ -O3