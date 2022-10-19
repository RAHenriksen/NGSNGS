#include "Briggs.h"
#include "mrand.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>

int SimBriggsModel(char seq[], int L, double nv, double lambda, double delta_s, double delta, mrand_t *mr,int strand,int& C_to_T_counter,int& G_to_A_counter,int& C_to_T_counter_rev,int& G_to_A_counter_rev){
    int IsDeam = 0;
    assert(L<1024);
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
          //if (i == 0&& strand == 0){C_to_T_counter++;}
          //else if(i == 0 && strand==1){C_to_T_counter_rev++;}
        }
        else{
          seq[i] = seq_intermediate[i]; //C
        }
      }else{
        seq[i] = seq_intermediate[i];
      }
    }
    for (int i = 0; i < r; i++){
      // 1 pos 3' (so last position in fragment)
      // r means right overhan (ss)
      //fprintf(stderr,"SECOND FOR LOOP \n");
      if (seq_intermediate[L-i-1] == 'G' || seq_intermediate[L-i-1] == 'g'){
        dtemp2 = mrand_pop(mr);//drand48_r(&buffer, &dtemp2);
        double u = dtemp2;//((double) rand_r(&seed)/ RAND_MAX);
        //fprintf(stderr,"Double u G 1 %f\n",u);
        if (u < delta_s){
          IsDeam = 1;
          seq[L-i-1] = 'A'; //'Y';
          //if (i == 0&& strand == 0){G_to_A_counter++;}
          //else if(i == 0 && strand==1){G_to_A_counter_rev++;}
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
    double d = nv/((L-l-r-1)*nv+1-nv); // P(m)
    int p_nick = l; //m
    double cumd = d; // P(m)

    while ((u_nick > cumd) && (p_nick < L-r-1)){
      // p_nick cannot be larger than L-r-1, so once p_nick is equal then we 
      // go out of this loop. So at that time p_nick = L-r-1. So positions within
      // the full fragment
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
              //if (i == 0 && strand == 0){C_to_T_counter++;}
              //else if(i == 0 && strand==1){C_to_T_counter_rev++;}
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
            if (i == 0){
              //if (i == 0 && strand==0){G_to_A_counter++;}
              //else if(i == 0 && strand==1){G_to_A_counter_rev++;}
            }
          }
          else{
            seq[i] = seq_intermediate[i]; //G
          }
        }else{
            seq[i] = seq_intermediate[i];
        }
    }
    
    if (seq_intermediate[0] == 'C' && seq[0] == 'T'){
      //fprintf(stderr," SUBSTITUTIONS %c > %c \n",seq_intermediate[0],seq[0]);
      if (strand == 0){C_to_T_counter++;}
      else if(strand==1){C_to_T_counter_rev++;}
    }
    if (seq_intermediate[L-1] == 'G' && seq[L-1] == 'A'){
      //fprintf(stderr," SUBSTITUTIONS %c > %c \n",seq_intermediate[L-1],seq[L-1]);
      if (strand == 0){G_to_A_counter++;}
      else if(strand==1){G_to_A_counter_rev++;}
    };
    
    // seq[0] = 
  //fprintf(stderr,"SEQUENCE \t %s\n-------------\n",seq);
  return IsDeam;
}
