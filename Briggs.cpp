#include "Briggs.h"
#include "mrand.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>

int Random_geometric_k(const double p,mrand_t *mr)
{
  double u = mrand_pop(mr);
  int k;

  if (p == 1){k = 1;}
  else if(p == 0){k=0;}
  else{k = log (u) / log (1 - p);}

  return floor(k);
}

int SimBriggsModel(char seq[], int L, double nv, double lambda, double delta_s, double delta, mrand_t *mr){
    int IsDeam = 0;
    //fprintf(stderr,"INSIDE THE BRIGGS MODEL \n");
    double dtemp1;double dtemp2;
    dtemp1 = mrand_pop(mr);
    dtemp2 = mrand_pop(mr);
    int l = 0;
    int r = L-1;
    //fprintf(stderr,"SEQUENCE \t %s\n",seq);
    char seq_intermediate[1024];
    //just to ensure no issues arise in case of not clearing out the intermediate sequence
    memset(seq_intermediate, 0, sizeof seq_intermediate);
    strcpy(seq_intermediate,seq);
    //fprintf(stderr,"----------------------\n");
    //fprintf(stderr,"THE SEED IS %u and l : %d and r : %d\n",&seed,l,r);
    while (l+r > L-2){
      l = 0;
      r = 0;
      double u_l = dtemp1; // between 0 and 1
      double u_r = dtemp2; // between 0 and 1
      dtemp1 = mrand_pop(mr);//why ?
      dtemp2 = mrand_pop(mr);//why ?

      if (u_l > 0.5){
        l = (int) Random_geometric_k(lambda,mr); //Random_geometric_k(23424,lambda);//distribution1(generator1);
      }
      if (u_r > 0.5){
        r = (int) Random_geometric_k(lambda,mr); //Random_geometric_k(seed,lambda); //distribution1(generator2); //(int) ((rand_r(&seed)%30000)+1)
      }
    }
    //fprintf(stderr,"R and L values %d \t %d\n",r,l);
    for (int i = 0; i<l; i++){
      // l means left overhang (ss)
      //fprintf(stderr,"FIRST FOR LOOP \n");
      if (seq_intermediate[i] == 'C' || seq_intermediate[i] == 'c' ){
        dtemp1 = mrand_pop(mr);//drand48_r(&buffer, &dtemp1);
        double u = dtemp1; //((double) rand_r(&seed)/ RAND_MAX);//
        //fprintf(stderr,"Double u C 1 %f\n",u);
        if (u < delta_s){
          IsDeam = 1;
          seq[i] = 'T'; //X
        }else{
          seq[i] = seq_intermediate[i]; //C
        }
      }else{
        seq[i] = seq_intermediate[i];
      }
    }
    for (int i = 0; i < r; i++){
      // r means right overhan (ss)
      //fprintf(stderr,"SECOND FOR LOOP \n");
      if (seq_intermediate[L-i-1] == 'G' || seq_intermediate[L-i-1] == 'g'){
        dtemp2 = mrand_pop(mr);//drand48_r(&buffer, &dtemp2);
        double u = dtemp2;//((double) rand_r(&seed)/ RAND_MAX);
        //fprintf(stderr,"Double u G 1 %f\n",u);
        if (u < delta_s){
          IsDeam = 1;
          seq[L-i-1] = 'A'; //'Y';
        }
        else{
          seq[L-i-1] = seq_intermediate[L-i-1];
        }
      }else{
        seq[L-i-1] = seq_intermediate[L-i-1];
      }
    }
    dtemp1 = mrand_pop(mr);//drand48_r(&buffer, &dtemp1);
    double u_nick = dtemp1; //((double) rand_r(&seed)/ RAND_MAX);
    //fprintf(stderr,"Double u_nick %f\n",u_nick);
    double d = nv/((L-l-r-1)*nv+1-nv);
    int p_nick = l;
    double cumd = d;
    while ((u_nick > cumd) && (p_nick < L-r-1)){
        cumd += d;
        p_nick +=1;
    }
    for (int i = l; i < L-r; i++){
      // The double strand part, the left and right hand overhang are probably cut, so only the midlle part of our DNA fragments (ds)
      //fprintf(stderr,"THIRD FOR LOOP \n");
        if ((seq_intermediate[i] == 'C' || seq_intermediate[i] == 'c') && i<=p_nick){
          dtemp1 = mrand_pop(mr);//drand48_r(&buffer, &dtemp1);
          double u = dtemp1; //((double) rand_r(&seed)/ RAND_MAX);
          //fprintf(stderr,"Double u C 2 %f\n",u);
          if (u < delta){
            IsDeam = 1;
            seq[i] = 'T'; //Q
          }
          else{
            seq[i] = seq_intermediate[i]; //C
          }
        }
        else if ((seq_intermediate[i] == 'G' || seq_intermediate[i] == 'g') && i>p_nick){
          dtemp2 = mrand_pop(mr);//drand48_r(&buffer, &dtemp2);
          double u = dtemp2; //((double) rand_r(&seed)/ RAND_MAX);
          if (u < delta){
            IsDeam = 1;
            seq[i] = 'A'; //W
          }else{
            seq[i] = seq_intermediate[i]; //G
          }
        }else{
            seq[i] = seq_intermediate[i];
        }
    }
  //fprintf(stderr,"SEQUENCE \t %s\n-------------\n",seq);
  return IsDeam;
}