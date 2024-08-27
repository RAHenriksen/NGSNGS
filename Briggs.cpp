#include "Briggs.h"
#include "mrand.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <htslib/kstring.h>
#include <zlib.h>
#include <iostream>

int Biotin_ds_454Roche(char seq[], int L, double nv, double lambda, double delta_s, double delta_d, mrand_t *mr,int strand,int& C_to_T_counter,int& G_to_A_counter,int& C_to_T_counter_rev,int& G_to_A_counter_rev){
  /*
  Biotin_ds_454Roche - Simulates deamination events in a DNA sequence based on given parameters and returns binary value whether any deamination occurred (1 occurred, 0 unaltered).

  @param seq: A character array representing the DNA sequence that will be modified by deamination. Alters sequence in place.
  @param L: An integer representing the length of the DNA sequence.
  @param nv: probability of introducing nicks in the DNA sequence.
  @param lambda: A double parameter used in the geometric distribution to determine the length of overhangs.
  @param delta_s: A double representing the probability threshold for a cytosine to be deaminated to uracil in the single-stranded regions (left overhang region).
  @param delta_d: A double representing the probability threshold for a cytosine to be deaminated to uracil or a guanine to adenine in the middle region of the sequence.
  @param mr: A pointer to an mrand_t struct used for generating random numbers.
  @param strand: An integer representing the DNA strand orientation (0 for forward, 1 for reverse).
  @param C_to_T_counter: integer tracking the number of C to T deamination events in the forward strand.
  @param G_to_A_counter: integer tracking the number of G to A deamination events in the forward strand.
  @param C_to_T_counter_rev: integer tracking the number of C to T deamination events in the reverse strand.
  @param G_to_A_counter_rev: integer tracking the number of G to A deamination events in the reverse strand.
  */
  
  double dtemp1;double dtemp2;
  dtemp1 = mrand_pop(mr);
  dtemp2 = mrand_pop(mr);

  int IsDeam = 0;
  
  //lengths of overhangs
  assert(L<1024);
  int left_overhang = 0;
  int right_overhang = L-1;
  char seq_intermediate[1024];
  
  //Clear out intermediate sequence
  memset(seq_intermediate, 0, sizeof seq_intermediate);
  strcpy(seq_intermediate,seq);

  // Determine left and right overhang lengths
  while (left_overhang+right_overhang > L-2){
    left_overhang = 0;
    right_overhang = 0;
    double uniform_prob_leftoverhang = dtemp1; // between 0 and 1
    double uniform_prob_rightoverhang = dtemp2; // between 0 and 1
    dtemp1 = mrand_pop(mr);
    dtemp2 = mrand_pop(mr);

    if (uniform_prob_leftoverhang > 0.5){
      left_overhang = (int) Random_geometric_k(lambda,mr);
    }
    if (uniform_prob_rightoverhang > 0.5){
      right_overhang = (int) Random_geometric_k(lambda,mr);
    }
  }
  
  for (int i = 0; i<left_overhang; i++){
    // left_overhang means left overhang (ss)
    if (seq_intermediate[i] == 'C' || seq_intermediate[i] == 'c' ){
      dtemp1 = mrand_pop(mr);
      double uniformprob_ss_left = dtemp1;
      if (uniformprob_ss_left < delta_s){
        // performs the deamination
        IsDeam = 1;
        seq[i] = 'T';
        if (i == 0&& strand == 0){C_to_T_counter++;}
        else if(i == 0 && strand==1){C_to_T_counter_rev++;}
      }
      else{
        // else copy the original sequence and continue
        seq[i] = seq_intermediate[i];
      }
    }
    else{
      // else copy the original sequence and continue
      seq[i] = seq_intermediate[i];
    }
  }
  
  for (int i = 0; i < right_overhang; i++){
    // 1 pos 3' (so last position in fragment)
    // r means right overhang (ss)
    if (seq_intermediate[L-i-1] == 'G' || seq_intermediate[L-i-1] == 'g'){
      dtemp2 = mrand_pop(mr);
      double uniformprob_ss_right = dtemp2;
      if (uniformprob_ss_right < delta_s){
        // performs the deamination from the 3' direction
        IsDeam = 1;
        seq[L-i-1] = 'A';
        if (i == 0&& strand == 0){G_to_A_counter++;}
        else if(i == 0 && strand==1){G_to_A_counter_rev++;}
      }
      else{
        // else copy the end of original sequence and continue
        seq[L-i-1] = seq_intermediate[L-i-1];
      }
    }
    else{
      // else copy the end of original sequence and continue
      seq[L-i-1] = seq_intermediate[L-i-1];
    }
  }

  dtemp1 = mrand_pop(mr);
  double uniformprob_nick = dtemp1;
  double prob_place_nick = nv/((L-left_overhang-right_overhang-1)*nv+1-nv); // P(m) probability of placing nick in the middle of the fragments
  int pos_nick = left_overhang; //m
  double cumd = prob_place_nick; // P(m)

  while ((uniformprob_nick > cumd) && (pos_nick < L-right_overhang-1)){
    // pos_nick cannot be larger than L-r-1, so once pos_nick is equal then we 
    // go out of this loop. So at that time p_nick = L-r-1. So positions within
    // the full fragment
    cumd += prob_place_nick; 
    pos_nick +=1; 
  }

  for (int i = left_overhang; i < L-right_overhang; i++){
    // The double strand part, the left and right hand overhang are probably cut, so only the midlle part of our DNA fragments (ds)
    if ((seq_intermediate[i] == 'C' || seq_intermediate[i] == 'c') && i<=pos_nick){
      dtemp1 = mrand_pop(mr);
      double uniformprob_ds_left = dtemp1;
      if (uniformprob_ds_left < delta_d){
        // perform deamination
        IsDeam = 1;
        seq[i] = 'T';
        if (i == 0 && strand == 0){
          C_to_T_counter++;
        }
        else if(i == 0 && strand==1){
          C_to_T_counter_rev++;
        }
      }
      else{
        // no alteration
        seq[i] = seq_intermediate[i];
      }
    }
    else if ((seq_intermediate[i] == 'G' || seq_intermediate[i] == 'g') && i>pos_nick){
      dtemp2 = mrand_pop(mr);
      double uniformprob_ds_right = dtemp2;
      if (uniformprob_ds_right < delta_d){
        // perform deamination
        IsDeam = 1;
        seq[i] = 'A';
        if (i == 0){
          if (i == 0 && strand==0){
            G_to_A_counter++;
          }
          else if(i == 0 && strand==1){
            G_to_A_counter_rev++;
          }
        }
      }
      else{
        // no alteration
        seq[i] = seq_intermediate[i];
      }
    }
    else{
      // no alteration
      seq[i] = seq_intermediate[i];
    }
  }
  
  // count number of deamination events
  if (seq_intermediate[0] == 'C' && seq[0] == 'T'){
    if (strand == 0){
      C_to_T_counter++;
    }
    else if(strand==1){
      C_to_T_counter_rev++;
    }
  }
  
  if (seq_intermediate[L-1] == 'G' && seq[L-1] == 'A'){
    if (strand == 0){
      G_to_A_counter++;
    }
    else if(strand==1){
      G_to_A_counter_rev++;
    }
  }

  return IsDeam;
}


int Biotin_ds_454Roche_kstring(kstring_t* seq, double nv, double lambda, double delta_s, double delta_d, mrand_t* mr, int strand, int& C_to_T_counter, int& G_to_A_counter, int& C_to_T_counter_rev, int& G_to_A_counter_rev) {
  /*
  Biotin_ds_454Roche - Simulates deamination events in a DNA sequence based on given parameters and returns binary value whether any deamination occurred (1 occurred, 0 unaltered).

  @param seq: A kstring representing the DNA sequence that will be modified by deamination.
  @param nv: probability of introducing nicks in the DNA sequence.
  @param lambda: A double parameter used in the geometric distribution to determine the length of overhangs.
  @param delta_s: A double representing the probability threshold for a cytosine to be deaminated to uracil in the single-stranded regions (left overhang region).
  @param delta_d: A double representing the probability threshold for a cytosine to be deaminated to uracil or a guanine to adenine in the middle region of the sequence.
  @param mr: A pointer to an mrand_t struct used for generating random numbers.
  @param strand: An integer representing the DNA strand orientation (0 for forward, 1 for reverse).
  @param C_to_T_counter: integer tracking the number of C to T deamination events in the forward strand.
  @param G_to_A_counter: integer tracking the number of G to A deamination events in the forward strand.
  @param C_to_T_counter_rev: integer tracking the number of C to T deamination events in the reverse strand.
  @param G_to_A_counter_rev: integer tracking the number of G to A deamination events in the reverse strand.
  */
  int IsDeam = 0;

  double dtemp1, dtemp2;
  dtemp1 = mrand_pop(mr);
  dtemp2 = mrand_pop(mr);

  int left_overhang = 0;
  int right_overhang = seq->l - 1;

  // Determine left and right overhang lengths
  while (left_overhang + right_overhang > seq->l - 2) {
    left_overhang = 0;
    right_overhang = 0;
    double uniform_prob_leftoverhang = dtemp1; // between 0 and 1
    double uniform_prob_rightoverhang = dtemp2; // between 0 and 1
    dtemp1 = mrand_pop(mr);
    dtemp2 = mrand_pop(mr);

    if (uniform_prob_leftoverhang > 0.5) {
      left_overhang = (int)Random_geometric_k(lambda, mr);
    }
    if (uniform_prob_rightoverhang > 0.5) {
      right_overhang = (int)Random_geometric_k(lambda, mr);
    }
  }

  // initialize intermediate kstrings to be altered
  kstring_t seq_intermediate;
  seq_intermediate.l = seq->l;
  seq_intermediate.m = seq->l;
  seq_intermediate.s = (char *)malloc((seq->l + 1) * sizeof(char)); // Allocate memory

  strcpy(seq_intermediate.s, seq->s); // Copy seq->s to seq_intermediate.s

  for (int i = 0; i < left_overhang; i++) {
    if (seq_intermediate.s[i] == 'C' || seq_intermediate.s[i] == 'c') {
      dtemp1 = mrand_pop(mr);
      double uniformprob_ss_left = dtemp1;
      if (uniformprob_ss_left < delta_s) {
        // performs the deamination
        IsDeam = 1;
        seq->s[i] = 'T'; //X
        if (i == 0 && strand == 0) {
          C_to_T_counter++;
        }
        else if (i == 0 && strand == 1){
          C_to_T_counter_rev++;
        }
      }
      else {
        // else copy the original sequence and continue
        seq->s[i] = seq_intermediate.s[i];
      }
    }
    else {
      // else copy the original sequence and continue
      seq->s[i] = seq_intermediate.s[i];
    }
  }

  for (int i = 0; i < right_overhang; i++) {
    if (seq_intermediate.s[seq->l - i - 1] == 'G' || seq_intermediate.s[seq->l - i - 1] == 'g') {
      dtemp2 = mrand_pop(mr);
      double uniformprob_ss_right = dtemp2;
      if (uniformprob_ss_right < delta_s) {
        // performs the deamination from the 3' direction
        IsDeam = 1;
        seq->s[seq->l-i-1] = 'A'; //Q
        if (i == 0 && strand == 0) {
          G_to_A_counter++;
        }
        else if (i == 0 && strand == 1) {
          G_to_A_counter_rev++;
        }
      }
      else{
      // else copy the end of original sequence and continue
        seq->s[seq->l-i-1] = seq_intermediate.s[seq->l-i-1];
      }
    }
    else {
      // else copy the end of original sequence and continue
      seq->s[seq->l-i-1] = seq_intermediate.s[seq->l-i-1];
    }
  }

  dtemp1 = mrand_pop(mr);
  double uniformprob_nick = dtemp1;
  double prob_place_nick = nv/((seq->l-left_overhang-right_overhang-1)*nv+1-nv); // P(m) probability of placing nick in the middle of the fragments
  int pos_nick = left_overhang; //m
  double cumd = prob_place_nick; // P(m)

  while ((uniformprob_nick > cumd) && (pos_nick < seq->l-right_overhang-1)){
    // pos_nick cannot be larger than L-r-1, so once pos_nick is equal then we 
    // go out of this loop. So at that time p_nick = L-r-1. So positions within
    // the full fragment
    cumd += prob_place_nick; 
    pos_nick +=1; 
  }

  for (int i = left_overhang; i < seq->l-right_overhang; i++){
    // The double strand part, the left and right hand overhang are probably cut, so only the midlle part of our DNA fragments (ds)
    if ((seq_intermediate.s[i] == 'C' || seq_intermediate.s[i] == 'c') && i<=pos_nick){
      dtemp1 = mrand_pop(mr);
      double uniformprob_ds_left = dtemp1;
      if (uniformprob_ds_left < delta_d){
        // perform deamination
        IsDeam = 1;
        seq->s[i] = 'T';
        if (i == 0 && strand == 0){
          C_to_T_counter++;
        }
        else if(i == 0 && strand==1){
          C_to_T_counter_rev++;
        }
      }
      else{
        // no alteration
        seq->s[i] = seq_intermediate.s[i];
      }
    }
    else if ((seq_intermediate.s[i] == 'G' || seq_intermediate.s[i] == 'g') && i>pos_nick){
      dtemp2 = mrand_pop(mr);
      double uniformprob_ds_right = dtemp2;
      if (uniformprob_ds_right < delta_d){
        // perform deamination
        IsDeam = 1;
        seq->s[i] = 'A';
        if (i == 0){
          if (i == 0 && strand==0){
            G_to_A_counter++;
          }
          else if(i == 0 && strand==1){
            G_to_A_counter_rev++;
          }
        }
      }
      else{
        // no alteration
        seq->s[i] = seq_intermediate.s[i];
      }
    }
    else{
      // no alteration
      seq->s[i] = seq_intermediate.s[i];
    }
  }

  // count number of deamination events
  if (seq_intermediate.s[0] == 'C' && seq->s[0] == 'T'){
    if (strand == 0){
      C_to_T_counter++;
    }
    else if(strand==1){
      C_to_T_counter_rev++;
    }
  }
  
  if (seq_intermediate.s[seq->l-1] == 'G' && seq->s[seq->l-1] == 'A'){
    if (strand == 0){
      G_to_A_counter++;
    }
    else if(strand==1){
      G_to_A_counter_rev++;
    }
  }

  free(seq_intermediate.s);
  seq_intermediate.s = NULL;
  seq_intermediate.l = seq_intermediate.m = 0;

  return IsDeam;
}

int PMD_Amplicon(kstring_t* seq, double nv, double lambda, double delta_s, double delta_d, mrand_t* mr){
  /*
  PMD_Amplicon - Simulates deamination events in a empirical DNA sequence based on given parameters and returns binary value whether any deamination occurred (1 occurred, 0 unaltered).

  @param seq: A kstring representing the sequence read from an input fasta,fastq,bam file that will be modified by deamination.
  @param nv: probability of introducing nicks in the DNA sequence.
  @param lambda: A double parameter used in the geometric distribution to determine the length of overhangs.
  @param delta_s: A double representing the probability threshold for a cytosine to be deaminated to uracil in the single-stranded regions (left overhang region).
  @param delta_d: A double representing the probability threshold for a cytosine to be deaminated to uracil or a guanine to adenine in the middle region of the sequence.
  @param mr: A pointer to an mrand_t struct used for generating random numbers.
  */

  /* 
  In order to insert deamination, we assume the sequence reads fully span the original deamination molecule, such that we can model what would have been the single-stranded and
  and double-stranded regions with the sequence read being created after blunt-end repair.
  */
    
  int IsDeam = 0;

  double dtemp1, dtemp2;
  dtemp1 = mrand_pop(mr);
  dtemp2 = mrand_pop(mr);

  int left_overhang = 0;
  int right_overhang = seq->l - 1;

  while (left_overhang + right_overhang > seq->l - 2) {
    left_overhang = 0;
    right_overhang = 0;
    double uniform_prob_leftoverhang = dtemp1; // between 0 and 1
    double uniform_prob_rightoverhang = dtemp2; // between 0 and 1
    dtemp1 = mrand_pop(mr);
    dtemp2 = mrand_pop(mr);

    if (uniform_prob_leftoverhang > 0.5) {
      left_overhang = (int)Random_geometric_k(lambda, mr);
    }
    if (uniform_prob_rightoverhang > 0.5) {
      right_overhang = (int)Random_geometric_k(lambda, mr);
    }
  }

  kstring_t seq_intermediate;
  seq_intermediate.l = seq->l;
  seq_intermediate.m = seq->l;
  seq_intermediate.s = (char *)malloc((seq->l + 1) * sizeof(char)); // Allocate memory

  strcpy(seq_intermediate.s, seq->s); // Copy seq->s to seq_intermediate.s

  for (int i = 0; i < left_overhang; i++) {
    if (seq_intermediate.s[i] == 'C' || seq_intermediate.s[i] == 'c') {
      dtemp1 = mrand_pop(mr);
      double uniformprob_ss_left = dtemp1;
      if (uniformprob_ss_left < delta_s) {
        // performs the deamination
        IsDeam = 1;
        seq->s[i] = 'T'; //X
      }
      else {
        // else copy the original sequence and continue
        seq->s[i] = seq_intermediate.s[i];
      }
    }
    else {
      // else copy the original sequence and continue
      seq->s[i] = seq_intermediate.s[i];
    }
  }

  for (int i = 0; i < right_overhang; i++) {
    if (seq_intermediate.s[seq->l - i - 1] == 'G' || seq_intermediate.s[seq->l - i - 1] == 'g') {
      dtemp2 = mrand_pop(mr);
      double uniformprob_ss_right = dtemp2;
      if (uniformprob_ss_right < delta_s) {
        // performs the deamination from the 3' direction
        IsDeam = 1;
        seq->s[seq->l-i-1] = 'A'; //Q
      }
      else {
        // else copy the end of original sequence and continue
        seq->s[seq->l-i-1] = seq_intermediate.s[seq->l-i-1];
      }
    }
    else {
      // else copy the end of original sequence and continue
      seq->s[seq->l-i-1] = seq_intermediate.s[seq->l-i-1];
    }
  }

  dtemp1 = mrand_pop(mr);
  double uniformprob_nick = dtemp1;
  double prob_place_nick = nv/((seq->l-left_overhang-right_overhang-1)*nv+1-nv); // P(m) probability of placing nick in the middle of the fragments
  int pos_nick = left_overhang; //m
  double cumd = prob_place_nick; // P(m)

  while ((uniformprob_nick > cumd) && (pos_nick < seq->l-right_overhang-1)){
    // pos_nick cannot be larger than L-r-1, so once pos_nick is equal then we 
    // go out of this loop. So at that time p_nick = L-r-1. So positions within
    // the full fragment
    cumd += prob_place_nick; 
    pos_nick +=1; 
  }

  for (int i = left_overhang; i < seq->l-right_overhang; i++){
    // The double strand part, the left and right hand overhang are probably cut, so only the midlle part of our DNA fragments (ds)
    if ((seq_intermediate.s[i] == 'C' || seq_intermediate.s[i] == 'c') && i<=pos_nick){
      dtemp1 = mrand_pop(mr);
      double uniformprob_ds_left = dtemp1;
      if (uniformprob_ds_left < delta_d){
        // perform deamination
        IsDeam = 1;
        seq->s[i] = 'T';
      }
      else{
        // no alteration
        seq->s[i] = seq_intermediate.s[i];
      }
    }
    else if ((seq_intermediate.s[i] == 'G' || seq_intermediate.s[i] == 'g') && i>pos_nick){
      dtemp2 = mrand_pop(mr);
      double uniformprob_ds_right = dtemp2;
      if (uniformprob_ds_right < delta_d){
        // perform deamination
        IsDeam = 1;
        seq->s[i] = 'A';
      }
      else{
        // no alteration
        seq->s[i] = seq_intermediate.s[i];
      }
    }
    else{
      // no alteration
      seq->s[i] = seq_intermediate.s[i];
    }
  }
  // Cleanup
  free(seq_intermediate.s);
  seq_intermediate.s = NULL;
  seq_intermediate.l = seq_intermediate.m = 0;

  return IsDeam;
}

