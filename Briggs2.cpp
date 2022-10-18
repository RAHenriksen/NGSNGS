#include "Briggs2.h"
#include "mrand.h"
#include "fasta_sampler.h"
#include "NGSNGS_misc.h"

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>

#define LENS 4096
#define MAXBINS 100

extern int refToInt[256];
extern char NtComp[5];
extern const char *bass;

int SimBriggsModel2(char *ori, int L, double nv, double lambda, double delta_s, double delta, mrand_t *mr,char **res) {
  int IsDeam = 0;
  assert(L<1024);

  int l = 0;
  int r = L-1;
  
  char *rasmus = res[0];
  char *thorfinn = res[2];
  
  while (l+r > L-2){
    l = 0;
    r = 0;
    double u_l = mrand_pop(mr); // between 0 and 1
    double u_r = mrand_pop(mr); // between 0 and 1
    
    if (u_l > 0.5){
      l = (int) Random_geometric_k(lambda,mr);
    }
    if (u_r > 0.5){
      r = (int) Random_geometric_k(lambda,mr);
    }
  }
  
  strncpy(rasmus,ori+l,L-l-r);
  //fprintf(stderr,"SEQUNEC %s\n",rasmus);

  strncpy(thorfinn,ori+l,L-l-r); //Thorfinn equals Rasmus
  //fprintf(stderr,"SEQUNEC %s\n",thorfinn);

  ReversComplement(thorfinn); // revers complement of rasmus;
  //fprintf(stderr,"SEQUNEC %s\n",thorfinn);
  
  /*for(int i = 0; i < strlen(rasmus);i++){
    fprintf(stderr,"i value %d\n",i);
  }*/

  for (int i = 0; i<l; i++){
      // l means left overhang (ss)
    rasmus[i] = ori[i];
    thorfinn[i] = NtComp[refToInt[(unsigned char) rasmus[i]]]; //Since its a single nucleotide then i can just create the complementary of rasmus
    
    if (rasmus[i] == 'C' || rasmus[i] == 'c' ){
      double u = mrand_pop(mr);
      if (u < delta_s){
        IsDeam = 1;
        rasmus[i] = 'T'; //X
        thorfinn[i] = 'A';
      }
    }
  }
  
  for (int i =L-r; i<L ; i++){
    // 1 pos 3' (so last position in fragment)
    // r means right overhan (ss)
    rasmus[i] = ori[i];
    thorfinn[i] = NtComp[refToInt[(unsigned char) rasmus[i]]];//REVCOMPL
    
    if (thorfinn[i] == 'C' || thorfinn[L-i-1] == 'c'){
      double u = mrand_pop(mr);
      
      if (u < delta_s){
        IsDeam = 1;
        thorfinn[i] = 'T'; //'Y';
        rasmus[i] = 'A';
      }
    }
    
    //assumptions, only one nick in double stranded part
    double u_nick_m = mrand_pop(mr);
    
    double P_m = nv/((L-l-r-1)*nv+1-nv);
    int p_nick_m = l;
    double CumPm = P_m;

    while ((u_nick_m > CumPm) && (p_nick_m < L-r-1)){
      CumPm += P_m;
      p_nick_m +=1;
    }
    
    /*
    double u_nick_n = mrand_pop(mr);
    
    int p_nick_n = l;
    double P_n = nv*pow((1-v), double p_nick_m-p_nick_n-L-r-1);
    double CumPn = P_n;

    while ((u_nick_n > CumPn) && (u_nick_n < L-r-1)){
      CumPm += P_n;
      p_nick_m +=1;
    }*/

    /*for (int i = l; i < L-r; i++){
      if ((ori[i] == 'C' || ori[i] == 'c') && i<=p_nick){
        double urasmus = mrand_pop(mr);
        double uthorfinn = mrand_pop(mr);
        if (urasmus < delta){
          IsDeam = 1;
          rasmus[i] = 'T';
        }
        if (uthorfinn < delta){
          IsDeam = 1;
          thorfinn[i] = 'T';
        }
      }
    }*/

    // Insert joint distribution for selecting a nick position in both strands simultaneously 

  }
  
  char seq_intermediate[1024] = {0};
  char seq_intermediate2[1024] = {0};
  memset(seq_intermediate, 0, sizeof seq_intermediate);
  memset(seq_intermediate2, 0, sizeof seq_intermediate2);
  strcpy(seq_intermediate,rasmus);
  strcpy(seq_intermediate2,thorfinn);

  res[0] = seq_intermediate;
  ReversComplement(rasmus);
  res[1] = rasmus;
  res[2] = seq_intermediate2;
  ReversComplement(thorfinn);
  res[3] = thorfinn;
  
  return IsDeam;
}
  
#ifdef __WITH_MAIN__
//g++ Briggs2.cpp -D__WITH_MAIN__ mrand.o
int main(){
  int maxfraglength = 100;
  mrand_t *mr = mrand_alloc(2,777);
  char original[1024];
  char **results = new char *[4];
  for(int i=0;i<4;i++){
    results[i] = new char[1024];
    memset(results[i],'\0',1024);
  }

  while(1){

    int flen = mrand_pop_long(mr) % maxfraglength;
    //fprintf(stderr,"FLEN IS %d\n",flen);
    for(int i=0;i<flen;i++)
      original[i] = bass[mrand_pop_long(mr) %4];

    fprintf(stderr,"ori: %s\n",original);
    SimBriggsModel2(original,flen,0.024,0.36,0.68,0.9,mr,results);

    fprintf(stderr,"ori: %s\n",original);
    for(int i=0;i<4;i++)
      fprintf(stderr,"res[%d]: %s\n",i,results[i]);

    break;
  }

  return 0;
}

#endif
//g++ Briggs2.cpp NGSNGS_misc.cpp -D__WITH_MAIN__ mrand.o fasta_sampler.o RandSampling.o ../htslib/libhts.a -std=c++11 -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto -ggdb