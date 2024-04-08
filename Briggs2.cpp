#include "Briggs2.h"
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

#define LENS 4096
#define MAXBINS 100

extern int refToInt[256];
extern char NtComp[5];
extern const char *bass;

int SimBriggsModel2(char *ori, int L, double nv, double lambda, double delta_s, double delta, mrand_t *mr,char **res,int strandR1,
                    int& C_to_T_counter,int& G_to_A_counter,
                    int& C_total,int& G_total){
  int IsDeam = 0;
  assert(L<1024);
  //fprintf(stderr,"\n------------------------\nINSIDE BRIGGS\n");

  //The input reference should always be equal to the 5' ---> fwrd ---> 3' orientation similar to the reference genome
  if (strandR1 == 1){
    //fprintf(stderr,"orig\t%s\n",ori);
    ReversComplement(ori);
  }
  /*
  strand = 0 //which are always in 5' -> fwd -> 3' so ref
  ori pre 	0	GACAGTGGAACTGGCCCTCAACGTATAGTGTGTAAAA
  ori post	0	GACAGTGGAACTGGCCCTCAACGTATAGTGTGTAAAA
  
  strand = 1 //which has previously been altered for those cases without deamination
  ori pre 	1	AGCGTTACCTAGAACAATTAGATCTGCTATAGGTATCT
  ori post	1	AGATACCTATAGCAGATCTAATTGTTCTAGGTAACGCT
  */

  //overhang lengths
  int l = 0;
  int r = L-1;
 
  while (l+r > L-2){
    l = 0;
    r = 0;
    double u_l = mrand_pop(mr);
    double u_r = mrand_pop(mr);
   
    if (u_l > 0.5){
      // Mean of 10^7 1.7789089147008
      l = (int) Random_geometric_k(lambda,mr);
    }
    if (u_r > 0.5){
      // Mean of 10^7 1.7800180303421
      r = (int) Random_geometric_k(lambda,mr);
    }
  }


  char *rasmus = res[0];
  char *thorfinn = res[1];
  char *thorfinn_rev_comp = res[2];
  char *rasmus_rev_comp = res[3];

  strncpy(rasmus,ori,L);
  strncpy(thorfinn,ori,L);
  Complement(thorfinn);
  

  // left 5' overhangs, Thorfinn's DMG pattern is fully dependent on that of Rasmus.
  /*
    5' CGTATACATAGGCACTATATCGACCACACT 3' Rasmus
    3'        TATCCGTGATATAGCTGGTGTGA 5' Thorfinn
  */

  std::vector<int> deamin_pos_vec={-1};

  for (int i = 0; i<l; i++){
    if (rasmus[i] == 'C' || rasmus[i] == 'c' ){
      if (i == 0){C_total++;}
      double u = mrand_pop(mr);
      if (u < delta_s){
        if (i == 0){C_to_T_counter++;G_to_A_counter++;}
        rasmus[i] = 'T'; 
        thorfinn[i] = 'A';
        deamin_pos_vec.push_back(i);
      }
    }
  }
  
  // right 5' overhangs, Rasmus's DMG pattern is fully dependent on that of Thorfinn.
  /*
    5' CGTATACATAGGCACTATATCGACC       3' 
    3' GCATATGTATCCGTGATATAGCTGGTGTGAC 5' 
  */
  for (int j = 0; j < r; j++){
    if (thorfinn[L-j-1] == 'C' || thorfinn[L-j-1] == 'c'){
      if (j == 0){C_total++;}
      double u = mrand_pop(mr);
      if (u < delta_s){
        if (j == 0){C_to_T_counter++;G_to_A_counter++;}
        rasmus[L-j-1] = 'A';
        thorfinn[L-j-1] = 'T';
        deamin_pos_vec.push_back((L-j-1));

      }
    }
  }

  if (nv > 0){
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
    }
    else if(p_nick_m == L-r-1){
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
        if (i == 0){C_total++;}
        //left of nick on thorfinn strand we change thorfinn according to rasmus

        /*
          5' CGTATACATAGGCACTATATCGACCACACT 3'
          3' GCATATGTA CCGTGATATAGCTGGTGTGA 5'
                        |
                        v
          5' CGTATACATAGGCACTATATCGACCACACT 3'
          3'           CCGTGATATAGCTGGTGTGA 5'                
        */
      
        double u = mrand_pop(mr);
        if (u < delta){
          rasmus[i] = 'T';
          thorfinn[i] = 'A'; //Downstream nick one DMG pattern depends on the other strand
          if (i == 0){C_to_T_counter++;G_to_A_counter++;}
          deamin_pos_vec.push_back((i));
        }
      }
      else if (i>p_nick_m && (thorfinn[i] == 'C' || thorfinn[i] == 'c')){
        // right side of rasmus nick we change rasmus according to thorfinn
        if (i == (L-1)){C_total++;}
        /*
          5' CGTATACAT GGCACTATATCGACCACACT 3'
          3' GCATATGTATCCGTGATATAGCTGGTGTGA 5'
                        |
                        v
          5' CGTATACAT                       3'
          3' GCATATGTATCCGTGATATAGCTGGTGTGA  5'                
        */

        double u = mrand_pop(mr);
        if (u < delta){
          rasmus[i] = 'A';
          thorfinn[i] = 'T'; //Downstream nick one DMG pattern depends on the other strand
          if (i == (L-1)){C_to_T_counter++;G_to_A_counter++;}
          deamin_pos_vec.push_back((i));
        }
      }

      // between the nick with rasmus showing DMG
      else if(i>=L-p_nick_n-1 && i<=p_nick_m && (rasmus[i] == 'C' || rasmus[i] == 'c')){
        if (i == 0){C_total++;}
        double u = mrand_pop(mr);
        if (u < delta){
          rasmus[i] = 'T'; //Upstream both nicks, DMG patterns are independent
          if (i == 0){C_to_T_counter++;}
          deamin_pos_vec.push_back((i));
        }
      }
      // between the nick with Thorfinn showing DMG
      else if(i>=L-p_nick_n-1 && i<=p_nick_m && (thorfinn[i] == 'C' || thorfinn[i] == 'c')){
        if (i == 0){C_total++;}
        double u = mrand_pop(mr);
        if (u < delta){
          thorfinn[i] = 'T'; //Upstream both nicks, DMG patterns are independent
          if (i == 0){G_to_A_counter++;}
          deamin_pos_vec.push_back((i));
        }
      }
    }
  }

  //Change orientation of Thorfinn to reverse strand
  reverseChar(thorfinn,strlen(thorfinn));
  strcpy(rasmus_rev_comp,rasmus);
  ReversComplement(rasmus_rev_comp);
  strcpy(thorfinn_rev_comp,thorfinn);
  ReversComplement(thorfinn_rev_comp);

  res[0] = rasmus;
  res[1] = thorfinn;
  res[2] = thorfinn_rev_comp;
  res[3] = rasmus_rev_comp;

  //fprintf(stderr, "DEAMIN %d \t length %d \t fraglength is %d \n", IsDeam, strlen(ori), L);
  //fprintf(stderr,"orig\t%.5s\nres[0]\t%.5s\nres[1]\t%.5s\n",ori,res[0],res[1]);
  //fprintf(stderr,"Ras res[0]\t%.*s\nTK res[1]\t%.*s\nRAS2 res[2]\t%.*s\nTK2 res[3]\t%.*s\n",5,res[0]+IsDeam-2,5,res[1]+IsDeam-2,5,res[2]+IsDeam-2,5,res[3]+IsDeam-2);

  //for(int i=0; i < deamin_pos_vec.size(); i++)
  // fprintf(stderr,"deamin pos %d\n",deamin_pos_vec.at(i));
  
  int mindeaminpos;
  int maxdeaminpos;
  
  if(deamin_pos_vec.size()>1){
    IsDeam = 1;
    //*std::minmax_element(deamin_pos_vec.begin()+1, deamin_pos_vec.end()).first;
    mindeaminpos = *std::minmax_element(deamin_pos_vec.begin()+1, deamin_pos_vec.end()).first;
    maxdeaminpos = *std::minmax_element(deamin_pos_vec.begin()+1, deamin_pos_vec.end()).second;
    //fprintf(stderr,"%d \t %d\n",mindeaminpos,maxdeaminpos);
    //fprintf(stderr,"orig\t%.5s\nres[0]\t%.5s\nres[1]\t%.5s\n",ori+mindeaminpos-2,res[0]+mindeaminpos-2,res[1]+mindeaminpos-2);
  }
  /*else{
    fprintf(stderr,"non deaminated \n");
  }

  if (strandR1 == 1){
    fprintf(stderr,"orig\t%s\nres[0]\t%s\nres[1]\t%s\nres[2]\t%s\nres[3]\t%s\n",ori,res[0],res[1],res[2],res[3]);
  }*/
  deamin_pos_vec.clear();

  //change the IsDeam to the minimum and maximum position at a given point
  return IsDeam;
}
  
#ifdef __WITH_MAIN__
//g++ Briggs2.cpp -D__WITH_MAIN__ mrand.o
int main(){
  int C_total = 0;int C_to_T_counter = 0;int C_to_T_counter_rev = 0;int C_total_rev=0;
  int G_total = 0;int G_to_A_counter = 0;int G_to_A_counter_rev = 0;int G_total_rev=0;
  int RefCp1 = 0; int RefCTp1 = 0;int RefCp2 = 0; int RefCTp2 = 0;
  int maxfraglength = 1000;
  int seed = 235;
  mrand_t *mr = mrand_alloc(0,seed);

  size_t reads = 1000000;
  int modulovalue = 10;
  size_t moduloread = reads/modulovalue;

  for (int i = 0; i < reads; i++){
    if (i%moduloread == 0){
      fprintf(stderr,"sequence %d \n",i);
    }
    
    char original[1024];
    char **results = new char *[4];
    for(int i=0;i<4;i++){
      results[i] = new char[1024];
      memset(results[i],'\0',1024);
    }

    int flen; 
    flen = mrand_pop_long(mr) % maxfraglength;
    if (flen < 30){flen = 35;}
    
    memset(original, 0, sizeof original);
    
    for(int i=0;i<flen;i++)
      original[i] = bass[mrand_pop_long(mr) %4];

    int strand = mrand_pop(mr)>0.5?0:1;
    //0.024,0.36,0.68,0.0097
    SimBriggsModel2(original,flen,0.024,0.36,0.68,0.0097,mr,results,strand,C_to_T_counter,G_to_A_counter,C_total,G_total); 
  }
  double C_Deam = (double)C_to_T_counter/(double)C_total;
  double G_Deam = (double)G_to_A_counter/(double)C_total;
  double Pair1 = (double) RefCTp1/(double)RefCp1;
  double Pair2 = (double) RefCTp2/(double)RefCp1;
  fprintf(stderr,"RefCTp1 %d \t RefCp1 %d \t RefCTp2 %d \t RefCp2 %d\n",RefCTp1,RefCp1,RefCTp2,RefCp2);
  fprintf(stderr,"[R,T,T',R'] C>T freq %f and G>A freq %f\n",C_Deam,G_Deam);
  fprintf(stderr,"[R,T'] C>T freq %f and [T,R'] C>T freq %f\n",Pair1,Pair2);
  return 0;
}

#endif
//g++ Briggs2.cpp NGSNGS_misc.cpp -D__WITH_MAIN__ mrand.o fasta_sampler.o RandSampling.o ../htslib/libhts.a -std=c++11 -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto -ggdb