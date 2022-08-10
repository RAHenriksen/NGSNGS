#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <cstdint>
#include <cmath>
#include <string>
#include <iostream>
#include <errno.h>
#include <random>

#include "Briggs.h"
#include "mrand.h"

int Random_geometric_k(unsigned int  seed, const double p)
{
  double u = ((double) rand_r(&seed)/ RAND_MAX); // this between 0 and 1
  int k;

  if (p == 1){k = 1;}
  else if(p == 0){k=0;}
  else{k = log (u) / log (1 - p);}

  return floor(k);
}

void SimBriggsModel(char* reffrag, char* frag, int L, double nv, double lambda, double delta_s, double delta, unsigned int seed,mrand_t *mr){

    double dtemp1;double dtemp2;
    dtemp1 = mrand_pop(mr); dtemp2 = mrand_pop(mr);
    int l = 0;
    int r = L-1;

    while (l+r > L-2){
      l = 0;
      r = 0;
      double u_l = dtemp1; // between 0 and 1
      double u_r = dtemp2; // between 0 and 1
      dtemp1 = mrand_pop(mr); dtemp2 = mrand_pop(mr);

      if (u_l > 0.5){
        l = Random_geometric_k((int) ((dtemp1*30000)+1),lambda); //Random_geometric_k(23424,lambda);//distribution1(generator1);
      }
      if (u_r > 0.5){
        r = Random_geometric_k((int) ((dtemp2*30000)+1),lambda); //Random_geometric_k(seed,lambda); //distribution1(generator2); //(int) ((rand_r(&seed)%30000)+1)
      }
    }
    for (int i = 0; i<l; i++){
      // l means left overhang (ss)
      if (reffrag[i] == 'C' || reffrag[i] == 'c' ){
        dtemp1 = mrand_pop(mr);//drand48_r(&buffer, &dtemp1);
        double u = dtemp1; //((double) rand_r(&seed)/ RAND_MAX);//
        if (u < delta_s){
          frag[i] = 'T'; //T
        }else{
          frag[i] = 'C'; //C
        }
      }else{
        frag[i] = reffrag[i];
      }
    }
    for (int i = 0; i < r; i++){
      // r means right overhan (ss)
      if (reffrag[L-i-1] == 'G' || reffrag[L-i-1] == 'g'){
        dtemp2 = mrand_pop(mr);//drand48_r(&buffer, &dtemp2);
        double u = dtemp2;//((double) rand_r(&seed)/ RAND_MAX);
        if (u < delta_s){
          frag[L-i-1] = 'A'; //A
        }
        else{
          frag[L-i-1] = 'G'; //G
        }
      }else{
        frag[L-i-1] = reffrag[L-i-1];
      }
    }
    dtemp1 = mrand_pop(mr);//drand48_r(&buffer, &dtemp1);
    double u_nick = dtemp1; //((double) rand_r(&seed)/ RAND_MAX);
    double d = nv/((L-l-r-1)*nv+1-nv);
    int p_nick = l;
    double cumd = d;
    while ((u_nick > cumd) && (p_nick < L-r-1)){
        cumd += d;
        p_nick +=1;
    }
    for (int i = l; i < L-r; i++){
      // The double strand part, the left and right hand overhang are probably cut, so only the midlle part of our DNA fragments (ds)
        if ((reffrag[i] == 'C' || reffrag[i] == 'c') && i<=p_nick){
          dtemp1 = mrand_pop(mr);//drand48_r(&buffer, &dtemp1);
          double u = dtemp1; //((double) rand_r(&seed)/ RAND_MAX);
          if (u < delta){
            frag[i] = 'T'; //T
          }
          else{
            frag[i] = 'C'; //C
          }
        }
        else if ((reffrag[i] == 'G' || reffrag[i] == 'g') && i>p_nick){
          dtemp2 = mrand_pop(mr);//drand48_r(&buffer, &dtemp2);
          double u = dtemp2; //((double) rand_r(&seed)/ RAND_MAX);
          if (u < delta){
            frag[i] = 'A'; //A
          }else{
            frag[i] = 'G'; //G
          }
        }else{
            frag[i] = reffrag[i];
        }
    }
}
