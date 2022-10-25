#include "Briggs2.h"
#include "mrand.h"
#include "fasta_sampler.h"
#include "NGSNGS_misc.h"

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <math.h>

#define LENS 4096
#define MAXBINS 100

extern int refToInt[256];
extern char NtComp[5];
extern const char *bass;

 
int SimBriggsModel2(char *ori, int L, double nv, double lambda, double delta_s, double delta, mrand_t *mr,char **res,int strandR1) {
  int IsDeam = 0;
  assert(L<1024);
 
  int l = 0;
  int r = L-1;
 
  char *rasmus = res[0];
  char *thorfinn = res[2];
 
  while (l+r > L-2){
    l = 0;
    r = 0;
    double u_l = mrand_pop(mr);
    double u_r = mrand_pop(mr);
   
    if (u_l > 0.5){
      l = (int) Random_geometric_k(lambda,mr);
    }
    if (u_r > 0.5){
      r = (int) Random_geometric_k(lambda,mr);
    }
  }

  //Please do a check in the following part, since it may shift by 1
 
  strncpy(rasmus,ori,L);
  //fprintf(stderr,"SEQUNEC %s\n",rasmus);
  
  //fprintf(stderr,"SEQUNEC %s\n",rasmus)
  strncpy(thorfinn,ori,L); //Thorfinn equals Rasmus
  //fprintf(stderr,"SEQUNEC %s\n",thorfinn);
  
  if (strandR1 == 0){
      ReversComplement(thorfinn);
  }
  else if (strandR1 == 1){
    //if input sequence are reverse comp, then rasmus needs to be corresponding to the input *
    // reference genome and its orientation, whereas thorfinn remains reverse complemented
    ReversComplement(rasmus);
  }
  //fprintf(stderr,"SEQUNEC %s\n",thorfinn);

  // Contain everything strncpy(rasmus,ori,L);

 
  /*for(int i = 0; i < strlen(rasmus);i++){
    fprintf(stderr,"i value %d\n",i);
  }*/
 
  for (int i = 0; i<l; i++){
    // left 5' overhangs, Thorfinn's DMG pattern is fully dependent on that of Rasmus.

    if (rasmus[i] == 'C' || rasmus[i] == 'c' ){
      double u = mrand_pop(mr);
      if (u < delta_s){
        IsDeam = 1;
        rasmus[i] = 'T'; 
        thorfinn[i] = 'A';
      }
    }
  }
   
  for (int i = 0; i < r; i++){
    // right 5' overhangs, Rasmus's DMG pattern is fully dependent on that of Thorfinn.
    if (thorfinn[L-i-1] == 'C' || thorfinn[L-i-1] == 'c'){
      double u = mrand_pop(mr);
      if (u < delta_s){
        IsDeam = 1;
        thorfinn[L-i-1] = 'T';
        rasmus[L-i-1] = 'A';
      }
    }
  }
  
  
  // The nick positions on both strands are denoted as (m,n). m (The nick position on Rasmus) is sampled as the previous way, while n (The nick position on thorfinn) is sampled according to a 
  // conditional probability given m.
  double u_nick_m = mrand_pop(mr);
   
	// the counting starts from 0 rather than one so we shift
  double P_m = nv/((L-l-r-1)*nv+1-nv);
  int p_nick_m = l;
  double CumPm = P_m;
 
  while ((u_nick_m > CumPm) && (p_nick_m < L-r-1)){
    CumPm += P_m;
    p_nick_m +=1;
  }
  

  int p_nick_n;
  double u_nick_n = mrand_pop(mr);
  double CumPn;
 
  // Given m, sampling n
  if (p_nick_m < L-r-1){
      p_nick_n = L-p_nick_m-2; //we shift both n and m
      CumPn = nv;
      while((u_nick_n > CumPn) && (p_nick_n < L-l-1)){
         p_nick_n +=1;
         CumPn += nv*pow(1-nv,p_nick_m+p_nick_n-L+2);
      }
  }else if(p_nick_m == L-r-1){
      p_nick_n = r;
      CumPn = nv;
      while((u_nick_n > CumPn) && (p_nick_n < L-l-1)){
         p_nick_n +=1;
         CumPn += nv*pow(1-nv,p_nick_n-r);
      }
  }
  
  
  // Way 2 Complicated Way (should be a little bit faster)
  for (int i = l; i < L-r; i++){

    if (i<L-p_nick_n-1 && (rasmus[i] == 'C' || rasmus[i] == 'c')){
      //left of nick on thorfinn strand we change thorfinn according to rasmus
      double u = mrand_pop(mr);
      if (u < delta){
        IsDeam = 1;
        rasmus[i] = 'T';
        thorfinn[i] = 'A'; //Downstream nick one DMG pattern depends on the other strand
      }
    }
    else if (i>p_nick_m && (thorfinn[i] == 'C' || thorfinn[i] == 'c')){
      // right side of rasmus nick we change rasmus according to thorfinn
      double u = mrand_pop(mr);
      if (u < delta){
        IsDeam = 1;
        rasmus[i] = 'A';
        thorfinn[i] = 'T'; //Downstream nick one DMG pattern depends on the other strand
      }
    }
	
    // between the nick with rasmus showing DMG
    else if(i>=L-p_nick_n-1 && i<=p_nick_m && (rasmus[i] == 'C' || rasmus[i] == 'c')){
      double u = mrand_pop(mr);
      if (u < delta){
        IsDeam = 1;
        rasmus[i] = 'T'; //Upstream both nicks, DMG patterns are independent
      }
    }
    // between the nick with Thorfinn showing DMG
    else if(i>=L-p_nick_n-1 && i<=p_nick_m && (thorfinn[i] == 'C' || thorfinn[i] == 'c')){
      double u = mrand_pop(mr);
      if (u < delta){
        IsDeam = 1;
        thorfinn[i] = 'T'; //Upstream both nicks, DMG patterns are independent
      }
    }
  }
  
  char seq_intermediate[1024] = {0};
  char seq_intermediate2[1024] = {0};
  memset(seq_intermediate, 0, sizeof seq_intermediate);
  memset(seq_intermediate2, 0, sizeof seq_intermediate2);
  strcpy(seq_intermediate,rasmus);
  strcpy(seq_intermediate2,thorfinn);
 
  res[0] = seq_intermediate;
  ReversComplement(rasmus);
  res[3] = rasmus;
  res[1] = seq_intermediate2;
  ReversComplement(thorfinn);
  res[2] = thorfinn;
 
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