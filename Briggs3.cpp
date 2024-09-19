#include "Briggs3.h"
#include "mrand.h"
#include "fasta_sampler.h"
#include "NGSNGS_misc.h"

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <math.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>  // std::begin, std::end
#include <htslib/kstring.h>

extern int refToInt[256];
extern char NtComp[5];
extern const char *bass;

int SimBriggsModel3(char *ori, int L, double nv, double lambda, double delta_s, double delta_d,double end3primeoverhang, mrand_t *mr,char **&res,int* &frag_type,int strandR1,
                    int& C_to_T_counter,int& G_to_A_counter,int& C_total,int& G_total,
                    int& total_frag,int& fwd_fragment_count,int& rev_fragment_count){
  
  int IsDeam = 0;
  assert(L<1024);
  //fprintf(stderr,"\n------------------------\nINSIDE BRIGGS\n");

  //The input reference should always be equal to the 5' ---> fwrd ---> 3' orientation similar to the reference genome
  if (strandR1 == 1){
    //fprintf(stderr,"orig\t%s\n",ori);
    ReversComplement(ori);
  }

  // create out sequence containers
  char *rasmus = (char*)malloc((strlen(ori) + 1) * sizeof(char));
  char *thorfinn = (char*)malloc((strlen(ori) + 1) * sizeof(char));

  if (rasmus == nullptr || thorfinn == nullptr) {
    // Handle allocation failure
    if (rasmus) free(rasmus);
    if (thorfinn) free(thorfinn);
    fprintf(stderr,"fail 1\n");
    exit(1);
  }

  strncpy(rasmus,ori,L);
  rasmus[L] = '\0';

  strncpy(thorfinn,ori,L);
  thorfinn[L] = '\0';
  Complement(thorfinn);

  //fprintf(stderr,"rasmus and thorfinn succesfully created and allocated\n");
  // so now rasmus represent the forward strand of the original aDNA molecule
  // thorfinn represent the reverse strand of the original aDNA molecule but with the same orientation as rasmus


  int l5 = 0, r5 = L - 1, l3 = 0, r3 = 0;

  // Determine lengths of overhangs based on lambda
  while ((l5 + r5 > L - 2) || (l3 + r3 > L - 2) || (l3 + r5 > L - 2) || (r3 + l5 > L - 2)) {
    l5 = r5 = l3 = r3 = 0;
    double u_l = mrand_pop(mr);
    double u_r = mrand_pop(mr);

    // e.g. end3primeoverhang == 0.2 then 80% chance of 5' overhang and 20% 3' overhang
    if (u_l > end3primeoverhang){
      l5 = (int) Random_geometric_k(lambda,mr);
    }
    else{
      l3 = (int) Random_geometric_k(lambda,mr);
    } 

    if (u_r > end3primeoverhang){
      r5 = (int) Random_geometric_k(lambda,mr);
    }
    else{
      r3 = (int) Random_geometric_k(lambda,mr);
    }

  }

  //l3 = 20;
  //r5 = 20;
  std::vector<int> deamin_pos_vec;
  
  //fprintf(stderr,"rasmus sequence 1 %s \n",rasmus);
  //fprintf(stderr,"TK sequence 1 %s \n",thorfinn);

  // Deamination in single-stranded ends for 5' left overhang
  int tk_n_val = 0;
  for (int i = 0; i < l5; i++) {
    // Replace thorfinn sequence with N to represent non-existing overhang regions within thorfinn,
    // due to no blunt-end repair, while keeping the length of the fragment the same length.
    // Before storing these fragment sequences as reads we need to remove the N parts later on
    thorfinn[i] = 'N'; 
    tk_n_val = i;
    if (rasmus[i] == 'C' || rasmus[i] == 'c'){
      if (i == 0){C_total++;}
      double u = mrand_pop(mr);
      if (u < delta_s) {
        if (i == 0){C_to_T_counter++;G_to_A_counter++;}
        rasmus[i] = 'T';
        deamin_pos_vec.push_back(i);
      }
    }
  }
  //fprintf(stderr,"TK sequence 2 %s \n",thorfinn);

  // Deamination in single-stranded ends for 3' left overhang
  for (int i = 0; i < l3; i++){
    rasmus[i] = 'N';
    if (thorfinn[i] == 'C' || thorfinn[i] == 'c'){
      if (i == 0){C_total++;}
      double u = mrand_pop(mr);
      if (u < delta_s) {
        if (i == 0){C_to_T_counter++;G_to_A_counter++;}
        thorfinn[i] = 'T';
        deamin_pos_vec.push_back(i);
      }
    }
  }
  //fprintf(stderr,"rasmus sequence 2 %s \n",rasmus);

  // Deamination in single-stranded ends for 5' right overhang
  for (int j = 0; j < r5; j++){
    rasmus[L-j-1] = 'N';
    if (thorfinn[L-j-1] == 'C' || thorfinn[L-j-1] == 'c'){
      if (j == 0){C_total++;}
      double u = mrand_pop(mr);
      if (u < delta_s) {
        if (j == 0){C_to_T_counter++;G_to_A_counter++;}
        thorfinn[L-j-1] = 'T';
        deamin_pos_vec.push_back(L-j-1);
      }
    }
  }
  //fprintf(stderr,"rasmus sequence 3 %s \n",rasmus);

  // Deamination in single-stranded ends for 3' right overhang
  for (int j = 0; j < r3; j++){
    thorfinn[L-j-1] = 'N'; 
    if (rasmus[L-j-1] == 'C' || rasmus[L-j-1] == 'c'){
      if (j == 0){C_total++;}
      double u = mrand_pop(mr);
      if (u < delta_s) {
        if (j == 0){C_to_T_counter++;G_to_A_counter++;}
        rasmus[L-j-1] = 'T';
        deamin_pos_vec.push_back(L-j-1);
      }
    }
  }
  //fprintf(stderr,"TK sequence 3 %s \n",thorfinn);

  // Double stranded region deamination
  for (int i = std::max(l5, l3); i < L - std::min(r5, r3); i++) {
    if ((rasmus[i] == 'C' || rasmus[i] == 'c')){
      if (i == std::max(l5, l3)){C_total++;}
      double u = mrand_pop(mr);
      if (u < delta_d) {
        rasmus[i] = 'T';
        if (i == std::max(l5, l3)){C_to_T_counter++;G_to_A_counter++;}
        deamin_pos_vec.push_back(i);
      }
    }
     if ((thorfinn[i] == 'C' || thorfinn[i] == 'c')){
      if (i == std::max(l5, l3)){C_total++;}
      double u = mrand_pop(mr);
      if (u < delta_d) {
        if (i == std::max(l5, l3)){C_to_T_counter++;G_to_A_counter++;}

        thorfinn[i] = 'T';
        deamin_pos_vec.push_back(i);
      }
    }
  }
  //fprintf(stderr,"rasmus sequence 4 %s \n",rasmus);
  //fprintf(stderr,"TK sequence 4 %s \n",thorfinn);

  // placement of nick
  std::vector<int> nick_pos_rasmus;
  std::vector<int> nick_pos_thorfinn;

  // storing the index of potential nicks based on uniform random variable for the positions within L-1
  for (int i = 0; i < L - 1; i++) {
    if (mrand_pop(mr) < nv) {
      nick_pos_rasmus.push_back(i);
    }
    if (mrand_pop(mr) < nv) {
      nick_pos_thorfinn.push_back(i);
    }
  }

  // Add the starting and ending positions for both rasmus and thorfinn of the double stranded 
  // region when disregarding single-stranded overhangs in both 3' and 5'
  int l_rasmus = std::max(0, l3);
  int l_thorfinn = std::max(0, l5);
  int r_rasmus = std::min(L, L - r5 -1);
  int r_thorfinn = std::min(L, L - r3 -1);

  nick_pos_rasmus.insert(nick_pos_rasmus.begin(), l_rasmus);
  nick_pos_rasmus.push_back(r_rasmus);

  nick_pos_thorfinn.insert(nick_pos_thorfinn.begin(), l_thorfinn);
  nick_pos_thorfinn.push_back(r_thorfinn);

  // Now split the deamination aDNA strands into indepdenent fragments based on nick positions.
  std::vector<char*> fwd_fragments;
  std::vector<char*> fwd_revcomp_fragments;
  std::vector<char*> rev_fragments;
  std::vector<char*> rev_revcomp_fragments;

  // to store the potential number of fragments created for both rasmus and thorfinn, as they are combined into one fragments vector, and later one we still need to remember the orientation of the fragments
  // such that any potential sequencing reads can be generated from them.

  // for rasmus sequences store in the vector
  //fprintf(stderr,"rasmus sequence 5 %s \n",rasmus);
  //fprintf(stderr,"TK sequence 5 %s \n",thorfinn);

  //fprintf(stderr,"what is l3 and r5 and l5 and r3 \t %d & %d & %d & %d\n",l3,r5,l5,r3);
  for (int i = 0; i < (int)nick_pos_rasmus.size() - 1; i++) {
    int rassize = (int) strlen(rasmus);
    int start = nick_pos_rasmus[i];
    int end = nick_pos_rasmus[i + 1];
    int frag_length = end - start;
      
    //fprintf(stderr,"What is the positions %d & %d & %d \n",start,end,frag_length);
    // Only keep fragments larger than 30 bases
    if (start >= l3 && frag_length > 30 && frag_length <= (rassize-r5)) {
    //if (frag_length > 30) {
      // so now every first element is rasmus every second element is rasmus_rev_comp
      char* fwd_fragment = new char[frag_length + 1];  // +1 for null-terminator
      char* fwd_fragment_revcomp = new char[frag_length + 1];

      //&& start > l3 && end < (rassize-r5)
      strncpy(fwd_fragment, rasmus + start, frag_length);
      //fprintf(stderr,"rasmus sequence 6 %s \n",fwd_fragment);

      // store the deaminated sequences into the fragments
      strncpy(fwd_fragment_revcomp, rasmus + start, frag_length);
      ReversComplement(fwd_fragment_revcomp);
      
      fwd_fragment[strlen(fwd_fragment)] = '\0';  // Null-terminate the string
      fwd_fragment_revcomp[strlen(fwd_fragment_revcomp)] = '\0';  // Null-terminate the string

      fwd_fragments.push_back(fwd_fragment); // asign to vector
      //fprintf(stderr,"fwd_fragment sequence %s \n",fwd_fragment);
      fwd_revcomp_fragments.push_back(fwd_fragment_revcomp); // asign to vector

      fwd_fragment_count++;
    }
  }
  //fprintf(stderr,"before tk \n");
  // for thorfinn sequences store in the vector
  for (int i = 0; i < (int)nick_pos_thorfinn.size() - 1; i++) {
    int thorfinnsize = (int) strlen(thorfinn);
    int start = nick_pos_thorfinn[i];
    int end = nick_pos_thorfinn[i + 1];
    int frag_length = end - start;
    //fprintf(stderr,"What is the tk positions %d & %d & %d \n",start,end,frag_length);

    // Only keep fragments larger than 30 bases
    if (start >= l5 && frag_length > 30 && frag_length <= (thorfinnsize-r3)) {
    //if (frag_length > 30) {

      char* bwd_fragment = new char[frag_length + 1];  // +1 for null-terminator
      char* bwd_fragment_revcomp = new char[frag_length + 1];

      int TKstrlen = (int) strlen(thorfinn);
      reverseChar(thorfinn,TKstrlen); // reverse thorfinn to the orientation of reverse strand
      strncpy(bwd_fragment, thorfinn + start, frag_length);

      strncpy(bwd_fragment_revcomp, thorfinn + start, frag_length);
      ReversComplement(bwd_fragment_revcomp);

      bwd_fragment[strlen(bwd_fragment)] = '\0';  // Null-terminate the string
      bwd_fragment_revcomp[strlen(bwd_fragment_revcomp)] = '\0';  // Null-terminate the string

      rev_fragments.push_back(bwd_fragment); // asign to vector
      rev_revcomp_fragments.push_back(bwd_fragment_revcomp); // asign to vector

      rev_fragment_count++;
    }
  }
  //fprintf(stderr,"done\n");
  /*
  So the fragments vector contains

  [ rasmus fragment 1, reverse complement rasmus fragment 1, thorfinn fragment 1, reverse complement thorfinn fragment 1,
    rasmus fragment 2, reverse complement rasmus fragment 2, thorfinn fragment 2, reverse complement thorfinn fragment 2,
    .
    .
    .
    rasmus fragment n-1, reverse complement rasmus fragment n-1, thorfinn fragment n-1, reverse complement thorfinn fragment n-1,
    rasmus fragment n, reverse complement rasmus fragment n, thorfinn fragment n, reverse complement thorfinn fragment n,
  ]

  */

  // now just remember how many elements in the fragments are from rasmus and thorfinn, since these represent two different strands of the original double stranded molecule.

  // Allocate space in res for each valid fragment
  total_frag = (int) (fwd_fragments.size()+fwd_revcomp_fragments.size()+rev_fragments.size()+rev_revcomp_fragments.size());
  //fprintf(stderr,"fragment number %d \n",total_frag);

  res = new char*[total_frag];

  frag_type = new int[total_frag];
  
  //fprintf(stderr,"res continue\n");
  int tmp_count = 0;
  for (int i = 0; i < (int)fwd_fragments.size(); i++) {
    //fprintf(stderr,"i val %d \n",i);
    int len = (int) strlen(fwd_fragments[i]);
    //fprintf(stderr,"len %d \n",len);
    res[i] = new char[len];
    strncpy(res[i],fwd_fragments[i],len);
    frag_type[i] = 0;
    //fprintf(stderr,"res %s \n",res[i]);

    //fprintf(stderr,"fwd_fragments %d \t\t length %d \t %s \n",i,len,res[i]);
  }
  tmp_count = (int)fwd_fragments.size();
  //fprintf(stderr,"-----------\n");

  for (int i = tmp_count; i < (int)(tmp_count+fwd_revcomp_fragments.size()); i++) {
    //res[i] = fwd_revcomp_fragments[i-tmp_count];
    int len = (int) strlen(fwd_revcomp_fragments[i-tmp_count]);

    res[i] = new char[len];
    strncpy(res[i],fwd_revcomp_fragments[i-tmp_count],len);
    frag_type[i] = 1;

    //fprintf(stderr,"fwd_revcomp_fragments %d \t length %d \t %s \n",i,len,res[i]);
  }
  tmp_count += fwd_revcomp_fragments.size();

  //fprintf(stderr,"-----------\n");

  for (int i = tmp_count; i < (int)(tmp_count+rev_fragments.size()); i++) {
    //res[i] = rev_fragments[i-tmp_count];
    int len = (int) strlen(rev_fragments[i-tmp_count]);

    res[i] = new char[len];
    strncpy(res[i],rev_fragments[i-tmp_count],len);
    frag_type[i] = 2;

    //fprintf(stderr,"rev_fragments %d \t\t length %d \t %s \n",i,len,res[i]);
  }
  tmp_count += rev_fragments.size();

  //fprintf(stderr,"-----------\n");

  for (int i = tmp_count; i < (int)(tmp_count+rev_revcomp_fragments.size()); i++) {
    //res[i] = rev_revcomp_fragments[i-tmp_count];    
    int len = (int) strlen(rev_revcomp_fragments[i-tmp_count]);
    res[i] = new char[len];
    strncpy(res[i],rev_revcomp_fragments[i-tmp_count],len);
    frag_type[i] = 3;

    //fprintf(stderr,"rev_revcomp_fragments %d \t length %d \t %s \n",i,len,res[i]);
  }

  for (int i = 0;i<fwd_fragments.size();i++) {
    delete[] fwd_fragments[i];
  }
  for (int i = 0;i<fwd_revcomp_fragments.size();i++) {
    delete[] fwd_revcomp_fragments[i];
  }
  for (int i = 0;i<rev_fragments.size();i++) {
    delete[] rev_fragments[i];
  }
  for (int i = 0;i<rev_revcomp_fragments.size();i++) {
    delete[] rev_revcomp_fragments[i];
  }

  if(deamin_pos_vec.size()>1){
    IsDeam = 1;
  }

  /*for (int i = 0; i<nick_pos_rasmus.size();i++){
    fprintf(stderr,"Nick pos rasmus is %d \n",nick_pos_rasmus[i]);
  }
  for (int i = 0; i<nick_pos_thorfinn.size();i++){
    fprintf(stderr,"Nick pos thorfinn is %d \n",nick_pos_thorfinn[i]);
  }*/
  
  //fprintf(stderr,"WHAT IS DEAMIN %d \n",IsDeam);
  
  free(rasmus);
  free(thorfinn);

  deamin_pos_vec.clear();
  nick_pos_rasmus.clear();
  nick_pos_thorfinn.clear();
  fwd_fragments.clear();
  fwd_revcomp_fragments.clear();
  rev_fragments.clear();
  rev_revcomp_fragments.clear();

  return IsDeam;
}


#ifdef __WITH_MAIN__
int main(){
  int C_total = 0;int C_to_T_counter = 0;int C_to_T_counter_rev = 0;int C_total_rev=0;
  int G_total = 0;int G_to_A_counter = 0;int G_to_A_counter_rev = 0;int G_total_rev=0;
  int RefCp1 = 0; int RefCTp1 = 0;int RefCp2 = 0; int RefCTp2 = 0;
  int maxfraglength = 150;
  int seed = 235;
  mrand_t *mr = mrand_alloc(0,seed);

  size_t reads = 100; // 1000000;
  int modulovalue = 10;
  size_t moduloread = reads/modulovalue;
  
  int fwdcount = 0;
  int revcount = 0;
  
  for (int i = 0; i < reads; i++){  
    char original[1024];

    int flen; 
    flen = 100 ; //mrand_pop_long(mr) % maxfraglength;
    if (flen < 30){flen = 35;}
    
    memset(original, 0, sizeof original);
    
    for(int i=0;i<flen;i++)
      original[i] = bass[mrand_pop_long(mr) %4];

    fprintf(stderr,"-------------\nSequence read %d \t length is %d \t \nsequence  %s \n",i,flen,original);

    int strand = 0;//mrand_pop(mr)>0.5?0:1;
    //0.024,0.36,0.68,0.0097

    //fprintf(stderr,"BEFORE SimBriggsModel3 \n");
    //fprintf(stderr,"Original sequence %s \n ------------- \n",original);

    char** results = nullptr;  // Pointer to hold results
    int* frag_type;
    int total_fragments;
    //fprintf(stderr,"before SimBriggsModel3\n");
    int deamn = SimBriggsModel3(original,flen,0.024,0.36,0.68,0.0097,mr,results,frag_type,strand,
      C_to_T_counter,G_to_A_counter,C_total,G_total,
      total_fragments,fwdcount,revcount,0.2); 
    //fprintf(stderr,"after SimBriggsModel3\n");

    //fprintf(stderr,"Fwd %d \t BWD %d \t  %d\n",fwdcount,revcount,total_fragments);
    //fprintf(stderr,"AFTER SimBriggsModel3 \n");
    //exit(1);
    for (int i = 0; i < total_fragments; i++) {
      fprintf(stderr,"Deam is %d \t Total fragments %d \t after sim %d \t %s \n",deamn,total_fragments,frag_type[i],results[i]);
      delete[] results[i];
    }
    delete[] results;
  }

  free(mr);

  double C_Deam = (double)C_to_T_counter/(double)C_total;
  double G_Deam = (double)G_to_A_counter/(double)C_total;
  double Pair1 = (double) RefCTp1/(double)RefCp1;
  double Pair2 = (double) RefCTp2/(double)RefCp1;
  /*
  fprintf(stderr,"RefCTp1 %d \t RefCp1 %d \t RefCTp2 %d \t RefCp2 %d\n",RefCTp1,RefCp1,RefCTp2,RefCp2);
  fprintf(stderr,"[R,T,T',R'] C>T freq %f and G>A freq %f\n",C_Deam,G_Deam);
  fprintf(stderr,"[R,T'] C>T freq %f and [T,R'] C>T freq %f\n",Pair1,Pair2);
  */
  return 0;
}

#endif
//g++ Briggs3.cpp -D__WITH_MAIN__ mrand.o fasta_sampler.o RandSampling.o NGSNGS_misc.o ../htslib/libhts.a -std=c++11 -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto -ggdb