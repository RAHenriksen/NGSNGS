#include "Briggs2.h"
#include "mrand.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
const char *bass = "ACGTN";

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
  // thorfinn := revers complement of rasmus;

  
  
  for (int i = 0; i<l; i++){
      // l means left overhang (ss)
    rasmus[i] = ori[i];
    thorfinn[i] = rasmus[i];//REVCOMPL
    
    if (ori[i] == 'C' || ori[i] == 'c' ){
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
    thorfinn[i] = rasmus[i];//REVCOMPL
    
    if (ori[i] == 'C' || ori[L-i-1] == 'c'){
      
      double u = mrand_pop(mr);
      
      if (u < delta_s){
	IsDeam = 1;
	thorfinn[i] = 'T'; //'Y';
	rasmus[i] = 'A';
      }
      
    }
    
    //assumptions, only one nick?
    
    double u_nick = mrand_pop(mr);
    
    double d = nv/((L-l-r-1)*nv+1-nv);
    int p_nick = l;
    double cumd = d;
    while ((u_nick > cumd) && (p_nick < L-r-1)){
      cumd += d;
      p_nick +=1;
    }
    
    for (int i = l; i < L-r; i++){
      if ((ori[i] == 'C' || ori[i] == 'c') && i<=p_nick){
	double urasmus = mrand_pop(mr);
	double uthorfinn = mrand_pop(mr);
	if (urasmus < delta){
	  IsDeam = 1;
	  rasmus[i] = 'T'; //Q
	}
	if (uthorfinn < delta){
	  IsDeam = 1;
	  thorfinn[i] = 'T'; //Q
	}
      }
    }
    
    
    //result[0] = rasmus;
    res[1] = rasmus;//REVCOMPL
    // res[2] = thorfinn;
    res[3] = thorfinn; //REVCOMPL
    
    return IsDeam;
  }
  return 0;
}
  
#ifdef __WITH_MAIN__
//g++ Briggs2.cpp -D__WITH_MAIN__ mrand.o
int main(){
  
  int maxfraglength = 30;
  mrand_t *mr = mrand_alloc(2,777);
  char original[1024];
  char **results = new char *[4];
  for(int i=0;i<4;i++){
    results[i] = new char[1024];
    memset(results[i],'\0',1024);
  }

  while(1){

    int flen = mrand_pop_long(mr) % maxfraglength;
    for(int i=0;i<flen;i++)
      original[i] = bass[mrand_pop_long(mr) %4];

    SimBriggsModel2(original,flen,0.024,0.36,0.68,0.0097,mr,results);

    fprintf(stderr,"ori: %s\n",original);
    for(int i=0;i<4;i++)
      fprintf(stderr,"res[%d]: %s\n",i,results[i]);

    
    break;
  }
  return 0;

}

#endif
