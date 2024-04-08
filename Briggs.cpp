#include "Briggs.h"
#include "mrand.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <htslib/kstring.h>
#include <zlib.h>
#include <iostream>

/*
#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif
*/
int SimBriggsModel(char seq[], int L, double nv, double lambda, double delta_s, double delta, mrand_t *mr,int strand,int& C_to_T_counter,int& G_to_A_counter,int& C_to_T_counter_rev,int& G_to_A_counter_rev){
    int IsDeam = 0;
    assert(L<1024);
    //std::cout << "FRAGMENT LENGTH "<< L << std::endl;
    double dtemp1;double dtemp2;
    dtemp1 = mrand_pop(mr);
    dtemp2 = mrand_pop(mr);
    int l = 0;
    int r = L-1;
    char seq_intermediate[1024];
    //Clear out intermediate sequence
    memset(seq_intermediate, 0, sizeof seq_intermediate);
    strcpy(seq_intermediate,seq);

    while (l+r > L-2){
      l = 0;
      r = 0;
      double u_l = dtemp1; // between 0 and 1
      double u_r = dtemp2; // between 0 and 1
      dtemp1 = mrand_pop(mr);
      dtemp2 = mrand_pop(mr);

      if (u_l > 0.5){
        l = (int) Random_geometric_k(lambda,mr);
      }
      if (u_r > 0.5){
        r = (int) Random_geometric_k(lambda,mr);
      }
    }
    for (int i = 0; i<l; i++){
      // l means left overhang (ss)
      if (seq_intermediate[i] == 'C' || seq_intermediate[i] == 'c' ){
        dtemp1 = mrand_pop(mr);
        double u = dtemp1;
        if (u < delta_s){
          IsDeam = 1;
          seq[i] = 'T';
          if (i == 0&& strand == 0){C_to_T_counter++;}
          else if(i == 0 && strand==1){C_to_T_counter_rev++;}
        }
        else{
          seq[i] = seq_intermediate[i];
        }
      }else{
        seq[i] = seq_intermediate[i];
      }
    }
    for (int i = 0; i < r; i++){
      //std::cout << "deamin G>A" << std::endl;
      // 1 pos 3' (so last position in fragment)
      // r means right overhang (ss)
      if (seq_intermediate[L-i-1] == 'G' || seq_intermediate[L-i-1] == 'g'){
        dtemp2 = mrand_pop(mr);
        double u = dtemp2;
        if (u < delta_s){
          IsDeam = 1;
          seq[L-i-1] = 'A';
          if (i == 0&& strand == 0){G_to_A_counter++;}
          else if(i == 0 && strand==1){G_to_A_counter_rev++;}
        }
        else{
          seq[L-i-1] = seq_intermediate[L-i-1];
        }
      }else{
        seq[L-i-1] = seq_intermediate[L-i-1];
      }
    }
    dtemp1 = mrand_pop(mr);
    double u_nick = dtemp1;
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
        if ((seq_intermediate[i] == 'C' || seq_intermediate[i] == 'c') && i<=p_nick){
          dtemp1 = mrand_pop(mr);
          double u = dtemp1;
          if (u < delta){
              IsDeam = 1;
              seq[i] = 'T';
              if (i == 0 && strand == 0){C_to_T_counter++;}
              else if(i == 0 && strand==1){C_to_T_counter_rev++;}
          }
          else{
            seq[i] = seq_intermediate[i];
          }
        }
        else if ((seq_intermediate[i] == 'G' || seq_intermediate[i] == 'g') && i>p_nick){
          dtemp2 = mrand_pop(mr);
          double u = dtemp2;
          if (u < delta){
            IsDeam = 1;
            seq[i] = 'A';
            if (i == 0){
              if (i == 0 && strand==0){G_to_A_counter++;}
              else if(i == 0 && strand==1){G_to_A_counter_rev++;}
            }
          }
          else{
            seq[i] = seq_intermediate[i];
          }
        }else{
            seq[i] = seq_intermediate[i];
        }
    }
    
    if (seq_intermediate[0] == 'C' && seq[0] == 'T'){
      if (strand == 0){C_to_T_counter++;}
      else if(strand==1){C_to_T_counter_rev++;}
    }
    if (seq_intermediate[L-1] == 'G' && seq[L-1] == 'A'){
      if (strand == 0){G_to_A_counter++;}
      else if(strand==1){G_to_A_counter_rev++;}
    };

  return IsDeam;
}


int SimBriggsModel_kstring(kstring_t* seq, double nv, double lambda, double delta_s, double delta, mrand_t* mr, int strand, int& C_to_T_counter, int& G_to_A_counter, int& C_to_T_counter_rev, int& G_to_A_counter_rev) {
    int IsDeam = 0;
    fprintf(stderr,"kstring sequence lenght %zu\n",seq->l);

    double dtemp1, dtemp2;
    dtemp1 = mrand_pop(mr);
    dtemp2 = mrand_pop(mr);

    int l = 0;
    int r = seq->l - 1;

    while (l + r > seq->l - 2) {
        l = 0;
        r = 0;
        double u_l = dtemp1; // between 0 and 1
        double u_r = dtemp2; // between 0 and 1
        dtemp1 = mrand_pop(mr);
        dtemp2 = mrand_pop(mr);

        if (u_l > 0.5) {
            l = (int)Random_geometric_k(lambda, mr);
        }
        if (u_r > 0.5) {
            r = (int)Random_geometric_k(lambda, mr);
        }
    }

    kstring_t seq_intermediate;
    memset(&seq_intermediate, 0, sizeof(kstring_t)); // Initialize kstring_t

    kputs(seq->s, &seq_intermediate); // Copy seq to seq_intermediate

    for (int i = 0; i < l; i++) {
        if (seq_intermediate.s[i] == 'C' || seq_intermediate.s[i] == 'c') {
            dtemp1 = mrand_pop(mr);
            double u = dtemp1;
            if (u < delta_s) {
                IsDeam = 1;
                kputc('T', seq);
                if (i == 0 && strand == 0) { C_to_T_counter++; }
                else if (i == 0 && strand == 1) { C_to_T_counter_rev++; }
            }
            else {
                kputc(seq_intermediate.s[i], seq);
            }
        }
        else {
            kputc(seq_intermediate.s[i], seq);
        }
    }

    for (int i = 0; i < r; i++) {
        std::cout << "deamin G>A" << std::endl;
        if (seq_intermediate.s[seq->l - i - 1] == 'G' || seq_intermediate.s[seq->l - i - 1] == 'g') {
            dtemp2 = mrand_pop(mr);
            double u = dtemp2;
            if (u < delta_s) {
                IsDeam = 1;
                kputc('A', seq);
                if (i == 0 && strand == 0) { G_to_A_counter++; }
                else if (i == 0 && strand == 1) { G_to_A_counter_rev++; }
            }
            else {
                kputc(seq_intermediate.s[seq->l - i - 1], seq);
            }
        }
        else {
            kputc(seq_intermediate.s[seq->l - i - 1], seq);
        }
    }

    // Cleanup
    free(seq_intermediate.s);

    // Other parts of the function remain the same

    return IsDeam;
}

