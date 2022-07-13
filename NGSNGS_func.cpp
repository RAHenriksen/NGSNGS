#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm> //std::reverse
#include <iostream>
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>
#include <errno.h>
#include <random>
#include <map>
#include <math.h>
#include <pthread.h>

#include "NGSNGS_func.h"
#include "mrand.h"

//#if defined(__APPLE__) && defined(__MACH__) 
//#include "NGSNGS_Random.h"
//#endif /* __APPLE__ */

#define LENS 4096
#define MAXBINS 100

void FragDistArray(int& number,int*& Length, double*& Frequency,int SizeDistType,int seed,int val1, int val2){
  std::default_random_engine generator(seed);
  const int nrolls=10000;
  int p[nrolls]={};
  
  if (SizeDistType==1){std::uniform_int_distribution<int> distribution(val1,val2);for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}}
  else if (SizeDistType==2){std::normal_distribution<double> distribution(val1,val2);for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}}
  else if (SizeDistType==3){std::lognormal_distribution<double> distribution(val1,val2);for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}}
  else if (SizeDistType==4){std::poisson_distribution<int> distribution(val1);for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}}
  else if (SizeDistType==5){std::exponential_distribution<double> distribution(val1);for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}}
  else if (SizeDistType==6){std::gamma_distribution<double> distribution(val1,val2);for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}}

  
  int len = sizeof(p)/sizeof(p[0]);
  sort(p,p+len,std::less<int>());

  std::map<int,int> counts;
  for (int i =0; i<nrolls;i++){counts[p[i]]++;}

  int n =1;

  Length[0] = 0; Frequency[0] = (float) 0;
  double CumulativeCount = 0;
  for (auto const &key : counts){
      CumulativeCount += (double) key.second/(double)nrolls;
      Length[n] = key.first;
      Frequency[n] = CumulativeCount;
      n++;
  }
  number = n;
}

void RandDistLength(int seed,double distval1, double distval2,int DistType,int& mean, int*& Length, double*& Frequency){
  std::default_random_engine generator(seed);
  const int nrolls=10000;
  int p[nrolls]={};

  //uniform
  if(DistType==1){
    std::uniform_int_distribution<int> distribution(distval1,distval2);
    for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}
    mean = (int) (0.5*(distval1+distval2));
  }
  //normal based
  else if(DistType==2){
    std::normal_distribution<double> distribution(distval1,distval2);
    for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}
    mean = (int) distval1;
  }
  else if(DistType==3){
    //lognormal(4,1)
    std::lognormal_distribution<double> distribution(distval1,distval2);
    for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}
    double param = (distval1+((distval2*distval2)/2));
    mean = (int) exp(param);
  }
  //rate based
  else if(DistType==4){
    std::poisson_distribution<int> distribution(distval1);
    for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}
    mean = (int) distval1;
  }
  else if(DistType==5){
    std::exponential_distribution<double> distribution(distval1);
    for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}
    mean = (int) 1/distval1;
  }
  else if(DistType==6){
    std::gamma_distribution<double> distribution(distval1,distval2);
    for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}
    mean = (int) (distval1/distval2);
  }

  int len = sizeof(p)/sizeof(p[0]);
  sort(p,p+len,std::less<int>());

  std::map<int,int> counts;
  for (int i =0; i<nrolls;i++){counts[p[i]]++;}

  double CumulativeCount = 0;
  int n = 0;
  for (auto const &key : counts){
      CumulativeCount += (double) key.second/(double)nrolls;
      Length[n] = key.first;
      Frequency[n] = CumulativeCount;
      n++;
  }
}

static void delete_seq(char *str, size_t seq_len, size_t del_len, size_t pos,int alt_len){
    //instert_seq(char *str, size_t len, char insert_seq[],size_t ins_len, size_t pos){
    for (int i = alt_len; i < del_len; i++)
    {
        //fprintf(stderr,"deletion pos 1 %d \t and pos %zu\n",i,pos-seq_len);
        //std::cout << str[pos] << " " << str[pos+1] << std::endl;
        //std::cout << seq_len << " " << pos << " " << pos - seq_len << std::endl;
        memmove(&str[pos], &str[pos+1], pos-seq_len);
        //memmove(&str[pos-1], &str[pos], pos-seq_len-1);
        //fprintf(stderr,"deletion pos 2 %d\n",seq_len - pos);
        seq_len--;
    }
}

static void delete_seq_ins(char *str, size_t seq_len, size_t del_len, size_t pos){
    // after insertion it deleted the reference allele, e.g. REF A -> ALT AGGGGGG, which creates a 6 bp insertion
    for (int i = 0; i < del_len; i++)
    {
        memmove(&str[pos], &str[pos+1], seq_len - pos);
        seq_len--;
    }
}

static void instert_seq(char *str, size_t len, char insert_seq[],size_t ins_len, size_t pos){
    for (int i = 0; i < ins_len; i++)
    {
        memmove(&str[pos+1], &str[pos], len - pos + 1);
        str[pos] = insert_seq[i];
        len++;
    }
    //the insertion doesn't overlap so i remove the nucleotide which the insertion replace
    delete_seq_ins(str, strlen(str),1,pos+ins_len);
}

void DNA_CAPITAL(char seq[]){
  while (*seq) {
    switch(*seq) {
      case 'A':
      case 'a':
        *seq = 'A';
        break;
      case 'G':
      case 'g':
        *seq = 'G';
        break;
      case 'C':
      case 'c':
        *seq = 'C';
        break;
      case 'T':
      case 't':
        *seq = 'T';
        break;
      case 'N':
      case 'n':
        *seq = 'N';
        break;  
    }
    ++seq;
  }
}

void DNA_complement(char seq[]){
  while (*seq) {
    switch(*seq) {
      case 'A':
      case 'a':
        *seq = 'T';
        break;
      case 'G':
      case 'g':
        *seq = 'C';
        break;
      case 'C':
      case 'c':
        *seq = 'G';
        break;
      case 'T':
      case 't':
        *seq = 'A';
        break;
      case 'N':
      case 'n':
        *seq = 'N';
        break;  
    }
    ++seq;
  }
}

double uniform(){
    /*
     U(0,1): AS 183: Appl. Stat. 31:188-190
     Wichmann BA & Hill ID.  1982.  An efficient and portable
     pseudo-random number generator.  Appl. Stat. 31:188-190
     x, y, z are any numbers in the range 1-30000.  Integer operation up
     to 30323 required.
     
     Suggested to me by Ziheng Yang who also provided me with
     the source code used here.  I use it because it is both fast and portable.
     */
    static int x_rndu=11, y_rndu=23,z_rndu=137;
    double r;
    
    x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
    y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
    z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
    if (x_rndu<0) x_rndu+=30269;
    if (y_rndu<0) y_rndu+=30307;
    if (z_rndu<0) z_rndu+=30323;
    r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
    return (r-(int)r);
}

void reverseChar(char* str,int length) {
    std::reverse(str, str + length);
}

void deletechar(char* array,int seq_len, size_t index_to_remove, int del_len){
  /*memmove(&array[index_to_remove],
        &array[index_to_remove + 5],
        seq_len - index_to_remove - 5);*/
  std::string s(array);
  s.erase(s.begin()+index_to_remove, s.end()-(seq_len-index_to_remove-del_len));
  strcpy(array, s.c_str());
}

void InsertChar(char* array,std::string ins,int index){
  std::string s(array);
  s.insert(index,ins);  
  strcpy(array, s.c_str());    
}

void ErrorSub(double randval,char seqchar[], int pos){
  //fprintf(stderr,"SUBERROR\n"); X, W, Z, Y
  if (seqchar[pos] == 'A' || seqchar[pos] == 'a'){
    if (0 < randval && randval <= 0.3333333333){seqchar[pos] = 'C';} //X 
    else if (0.3333333333 < randval && randval <= 0.6666666667){seqchar[pos] = 'G';} //T
    else if (0.6666666667 < randval && randval <= 1){seqchar[pos]  = 'T';} //C
  }
  else if (seqchar[pos] == 'C'|| seqchar[pos] == 'c'){ //'Z'
    if (0 < randval && randval <= 0.3333333333){seqchar[pos] = 'G';} //Z
    else if (0.3333333333 < randval && randval <= 0.6666666667){seqchar[pos]  = 'T';} //G
    else if (0.6666666667 < randval && randval <= 1){seqchar[pos]  = 'A';} //G
  }
  else if (seqchar[pos] == 'G'|| seqchar[pos] == 'g'){
    if (0 < randval && randval <= 0.3333333333){seqchar[pos] = 'T';} //A
    else if (0.3333333333 < randval && randval <= 0.6666666667){seqchar[pos]  = 'A';} //T
    else if (0.6666666667 < randval && randval <= 1){seqchar[pos]  = 'C';} //C
  }
  else if (seqchar[pos] == 'T'|| seqchar[pos] == 't'){
    if (0 < randval && randval <= 0.3333333333){seqchar[pos] = 'A';} //G
    else if (0.3333333333 < randval && randval <= 0.6666666667){seqchar[pos]  = 'C';} //A
    else if (0.6666666667 < randval && randval <= 1){seqchar[pos]  = 'G';} //C
  }
}

double* DeamFileArray(double* freqval,const char* filename,int &deamcyclelength){
    char buf[LENS];
    int i = 0;
    gzFile gz = Z_NULL;
    gz = gzopen(filename,"r");
    assert(gz!=Z_NULL);
    while(gzgets(gz,buf,LENS)){
        double val1;double val2;double val3;double val4;
        sscanf(buf,"%lf\t%lf\t%lf\t%lf\n",&val1,&val2,&val3,&val4);
        //std::cout << "iter " << i << std::endl;// " " << buf << std::endl;
        freqval[i*4] = val1; freqval[i*4+1] = val2; freqval[i*4+2] = val3; freqval[i*4+3] = val4;
        i++;
    }
    gzclose(gz);
    deamcyclelength = i/8;
  return freqval;
}

void Deam_File(char seq[],mrand_t *mr,double* freqval,int LEN){
  char ntdeam[4] = {'A', 'T', 'G', 'C'};//{'R', 'Q', 'S', 'U'};//{'X', 'Y', 'Z', 'W'}; //{'A', 'T', 'G', 'C'};
  double dtemp1;
  // 5' moving downwards from the 1 postion in the sequence 
  int Astart = 0;
  int Tstart = LEN*4; //4*15
  int Gstart = LEN*8; //15*8
  int Cstart = LEN*12; //15*12

  // moving upstream from the last position
  int Aend3 = LEN*20;
  int Tend3 = LEN*24; //4*15
  int Gend3 = LEN*28; //15*8
  int Cend3 = LEN*32; //15*12

  int seqlen = strlen(seq);
  //5'
  for (int row_idx = 0; row_idx < LEN;row_idx++){
    dtemp1 = mrand_pop(mr);//0.99;// mrand_pop(mr);
    //fprintf(stderr,"RANDOM VALUE %lf\n",dtemp1);
    if (seq[row_idx] == 'A' || seq[row_idx] == 'a'){
      if (dtemp1 <= freqval[Astart+(row_idx*4)]){seq[row_idx] = ntdeam[0];}
      else if (freqval[Astart+(row_idx*4)] < dtemp1 && dtemp1 <= freqval[Astart+(row_idx*4)+1]){seq[row_idx] = ntdeam[1];}
      else if (freqval[Astart+(row_idx*4)+1] < dtemp1 && dtemp1 <= freqval[Astart+(row_idx*4)+2]){seq[row_idx] = ntdeam[2];}
      else if (freqval[Astart+(row_idx*4)+2] < dtemp1 && dtemp1 <= freqval[Astart+(row_idx*4)+3]){seq[row_idx] = ntdeam[3];}
    }
    else if (seq[row_idx] == 'T' || seq[row_idx] == 't'){
      if (dtemp1 <= freqval[Tstart+(row_idx*4)]){seq[row_idx] = ntdeam[0];}
      else if (freqval[Tstart+(row_idx*4)] < dtemp1 && dtemp1 <= freqval[Tstart+(row_idx*4)+1]){seq[row_idx] = ntdeam[1];}
      else if (freqval[Tstart+(row_idx*4)+1] < dtemp1 && dtemp1 <= freqval[Tstart+(row_idx*4)+2]){seq[row_idx] = ntdeam[2];}
      else if (freqval[Tstart+(row_idx*4)+2] < dtemp1 && dtemp1 <= freqval[Tstart+(row_idx*4)+3]){seq[row_idx] = ntdeam[3];}
    }
    else if (seq[row_idx] == 'G' || seq[row_idx] == 'g'){
      if (dtemp1 <= freqval[Gstart+(row_idx*4)]){seq[row_idx] = ntdeam[0];}
      else if (freqval[Gstart+(row_idx*4)] < dtemp1 && dtemp1 <= freqval[Gstart+(row_idx*4)+1]){seq[row_idx] = ntdeam[1];}
      else if (freqval[Gstart+(row_idx*4)+1] < dtemp1 && dtemp1 <= freqval[Gstart+(row_idx*4)+2]){seq[row_idx] = ntdeam[2];}
      else if (freqval[Gstart+(row_idx*4)+2] < dtemp1 && dtemp1 <= freqval[Gstart+(row_idx*4)+3]){seq[row_idx] = ntdeam[3];}
    }
    else if (seq[row_idx] == 'C' || seq[row_idx] == 'c'){
      //fprintf(stderr,"FREQ VAL INDEX %d\n",Cstart+(row_idx));
      if (dtemp1 <= freqval[Cstart+(row_idx*4)]){seq[row_idx] = ntdeam[0];}
      else if (freqval[Cstart+(row_idx*4)] < dtemp1 && dtemp1 <= freqval[Cstart+(row_idx*4)+1]){seq[row_idx] = ntdeam[1];}
      else if (freqval[Cstart+(row_idx*4)+1] < dtemp1 && dtemp1 <= freqval[Cstart+(row_idx*4)+2]){seq[row_idx] = ntdeam[2];}
      else if (freqval[Cstart+(row_idx*4)+2] < dtemp1 && dtemp1 <= freqval[Cstart+(row_idx*4)+3]){seq[row_idx] = ntdeam[3];}
    }
  }
  //3'
  for (int row_idx = 0; row_idx < LEN;row_idx++){
    int row_idx_3p = seqlen-(LEN-row_idx);
    if (seq[row_idx_3p] == 'A' || seq[row_idx_3p] == 'a'){
      if (dtemp1 <= freqval[Aend3-((row_idx)*4)-4]){seq[row_idx_3p] = ntdeam[0];}
      else if (freqval[Aend3-((row_idx)*4)-4] < dtemp1 && dtemp1 <= freqval[Aend3-((row_idx)*4)-3]){seq[row_idx_3p] = ntdeam[1];}
      else if (freqval[Aend3-((row_idx)*4)-3] < dtemp1 && dtemp1 <= freqval[Aend3-((row_idx)*4)-2]){seq[row_idx_3p] = ntdeam[2];}
      else if (freqval[Aend3-((row_idx)*4)-2] < dtemp1 && dtemp1 <= freqval[Aend3-((row_idx)*4)-1]){seq[row_idx_3p] = ntdeam[3];}
    }
    else if (seq[row_idx_3p] == 'T' || seq[row_idx_3p] == 't'){
      if (dtemp1 <= freqval[Tend3-((row_idx)*4)-4]){seq[row_idx_3p] = ntdeam[0];}
      else if (freqval[Tend3-((row_idx)*4)-4] < dtemp1 && dtemp1 <= freqval[Tend3-((row_idx)*4)-3]){seq[row_idx_3p] = ntdeam[1];}
      else if (freqval[Tend3-((row_idx)*4)-3] < dtemp1 && dtemp1 <= freqval[Tend3-((row_idx)*4)-2]){seq[row_idx_3p] = ntdeam[2];}
      else if (freqval[Tend3-((row_idx)*4)-2] < dtemp1 && dtemp1 <= freqval[Tend3-((row_idx)*4)-1]){seq[row_idx_3p] = ntdeam[3];}
    }
    else if (seq[row_idx_3p] == 'G' || seq[row_idx_3p] == 'g'){
      if (dtemp1 <= freqval[Gend3-((row_idx)*4)-4]){seq[row_idx_3p] = ntdeam[0];}
      else if (freqval[Gend3-((row_idx)*4)-4] < dtemp1 && dtemp1 <= freqval[Gend3-((row_idx)*4)-3]){seq[row_idx_3p] = ntdeam[1];}
      else if (freqval[Gend3-((row_idx)*4)-3] < dtemp1 && dtemp1 <= freqval[Gend3-((row_idx)*4)-2]){seq[row_idx_3p] = ntdeam[2];}
      else if (freqval[Gend3-((row_idx)*4)-2] < dtemp1 && dtemp1 <= freqval[Gend3-((row_idx)*4)-1]){seq[row_idx_3p] = ntdeam[3];}
    }
    else if (seq[row_idx_3p] == 'C' || seq[row_idx_3p] == 'c'){
      if (dtemp1 <= freqval[Cend3-((row_idx)*4)-4]){seq[row_idx_3p] = ntdeam[0];}
      else if (freqval[Cend3-((row_idx)*4)-4] < dtemp1 && dtemp1 <= freqval[Cend3-((row_idx)*4)-3]){seq[row_idx_3p] = ntdeam[1];}
      else if (freqval[Cend3-((row_idx)*4)-3] < dtemp1 && dtemp1 <= freqval[Cend3-((row_idx)*4)-2]){seq[row_idx_3p] = ntdeam[2];}
      else if (freqval[Cend3-((row_idx)*4)-2] < dtemp1 && dtemp1 <= freqval[Cend3-((row_idx)*4)-1]){seq[row_idx_3p] = ntdeam[3];}
    }
  }
}

int Random_geometric_k(unsigned int  seed, const double p)
{
  double u = ((double) rand_r(&seed)/ RAND_MAX); // this between 0 and 1
  int k;

  if (p == 1){k = 1;}
  else if(p == 0){k=0;}
  else{k = log (u) / log (1 - p);}

  return floor(k);
}

double myrand(unsigned int persistent){
  return ((double) rand_r(&persistent)/ RAND_MAX);
}

int myrandgenmodulo(unsigned int seed, int modulo){
  return ((double) (rand_r(&seed)%modulo)+1);
}

void SimBriggsModel(char* reffrag, char* frag, int L, double nv, double lambda, double delta_s, double delta, unsigned int seed,mrand_t *mr){

    double dtemp1;double dtemp2;
    dtemp1 = mrand_pop(mr); dtemp2 = mrand_pop(mr);
    int l = 0;
    int r = L-1;
    //fprintf(stderr,"----------------------\n");
    //fprintf(stderr,"THE SEED IS %u and l : %d and r : %d\n",&seed,l,r);
    while (l+r > L-2){
      l = 0;
      r = 0;
      double u_l = dtemp1; // myrand((unsigned int) (seed)); //((double) rand_r(&seed)/ RAND_MAX);//uniform(); // between 0 and 1
      double u_r = dtemp2; // myrand((unsigned int) (seed+1)); //((double) rand_r(&seed)/ RAND_MAX);//uniform(); //between 0 and 1
      //fprintf(stderr,"uniform %f \t %f \n",u_l,u_r);
      dtemp1 = mrand_pop(mr); dtemp2 = mrand_pop(mr);

      if (u_l > 0.5){
        l = Random_geometric_k((int) ((dtemp1*30000)+1),lambda); //Random_geometric_k(23424,lambda);//distribution1(generator1);
      }
      if (u_r > 0.5){
        r = Random_geometric_k((int) ((dtemp2*30000)+1),lambda); //Random_geometric_k(seed,lambda); //distribution1(generator2); //(int) ((rand_r(&seed)%30000)+1)
      }
    }
    //fprintf(stderr,"R and L values %d \t %d\n",r,l);
    for (int i = 0; i<l; i++){
      // l means left overhang (ss)
      //fprintf(stderr,"FIRST FOR LOOP \n");
      if (reffrag[i] == 'C' || reffrag[i] == 'c' ){
        dtemp1 = mrand_pop(mr);//drand48_r(&buffer, &dtemp1);
        double u = dtemp1; //((double) rand_r(&seed)/ RAND_MAX);//uniform();
        //fprintf(stderr,"Double u C 1 %f\n",u);
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
      //fprintf(stderr,"SECOND FOR LOOP \n");
      if (reffrag[L-i-1] == 'G' || reffrag[L-i-1] == 'g'){
        dtemp2 = mrand_pop(mr);//drand48_r(&buffer, &dtemp2);
        double u = dtemp2;//((double) rand_r(&seed)/ RAND_MAX);//uniform();
        //fprintf(stderr,"Double u G 1 %f\n",u);
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
    double u_nick = dtemp1; //((double) rand_r(&seed)/ RAND_MAX);//uniform();
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
        if ((reffrag[i] == 'C' || reffrag[i] == 'c') && i<=p_nick){
          dtemp1 = mrand_pop(mr);//drand48_r(&buffer, &dtemp1);
          double u = dtemp1; //((double) rand_r(&seed)/ RAND_MAX);//uniform();
          //fprintf(stderr,"Double u C 2 %f\n",u);
          if (u < delta){
            frag[i] = 'T'; //T
          }
          else{
            frag[i] = 'C'; //C
          }
        }
        else if ((reffrag[i] == 'G' || reffrag[i] == 'g') && i>p_nick){
          dtemp2 = mrand_pop(mr);//drand48_r(&buffer, &dtemp2);
          double u = dtemp2; //((double) rand_r(&seed)/ RAND_MAX);//uniform();
          //fprintf(stderr,"Double u G 2 %f\n",u);
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

const char* Error_lookup(double a,double err[6000],int nt_offset, int read_pos,int outputoffset){

  //const char* nt_qual[8] = {"#", "\'", "0", "7" ,"<", "B", "F","I"}; // 35,39,48,55,66,70,73 //then withouth & in nt_qual[0]
  //const char* nt_qual[8] = {"#", "\'", "0", "7" ,"<", "B", "F","I"};//{"#", "\'", "0", "7" ,"<", "B", "F","I"};

  char nt_qual[8] = {(char) (2+outputoffset),
                     (char) (6+outputoffset),
                     (char) (15+outputoffset),
                     (char) (22+outputoffset),
                     (char) (27+outputoffset),
                     (char) (33+outputoffset),
                     (char) (37+outputoffset),
                     (char) (40+outputoffset)};
  //const char* nt_qual[8] = (const char*) nt_qual_int[9];
  int offset = ((nt_offset+read_pos)*8);
  //printf("offset %d \n", offset);
  const char* nt_out;
  if (a <= err[offset]){
    nt_out = &nt_qual[0];
  }
  else if (err[offset] < a && a <= err[offset+1]){
    nt_out = &nt_qual[1];
    }
  else if (err[offset+1] < a && a <= err[offset+2]){
    nt_out = &nt_qual[2];
  }
  else if (err[offset+2] < a && a <= err[offset+3]){
    nt_out = &nt_qual[3];
  }
  else if (err[offset+3] < a && a<= err[offset+4]){
    nt_out = &nt_qual[4];
  }
  else if (err[offset+4] < a && a<= err[offset+5]){
    nt_out = &nt_qual[5];
  }
  else if (err[offset+5] < a && a<= err[offset+6]){
    nt_out = &nt_qual[6];
  }
  else if (err[offset+6] < a && a <= err[offset+7]){
    nt_out = &nt_qual[7];
  }
  return nt_out;
}

double* Qual_array(double* freqval,const char* filename){
  int length = 6000;
  char buf[length];
  gzFile gz = Z_NULL;
  gz = gzopen(filename,"r");
  assert(gz!=Z_NULL);
  int i = 0;
  while(gzgets(gz,buf,length)){
    double val1;double val2;double val3;double val4;double val5;double val6;double val7;double val8;
    sscanf(buf,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&val1,&val2,&val3,&val4,&val5,&val6,&val7,&val8);
    //std::cout << "iter " << i << std::endl;// " " << buf << std::endl;
    freqval[i*8] = val1; freqval[i*8+1] = val2; freqval[i*8+2] = val3; freqval[i*8+3] = val4;
    freqval[i*8+4] = val5; freqval[i*8+5] = val6; freqval[i*8+6] = val7; freqval[i*8+7] = val8;
    i++;
  }
  gzclose(gz);
  return freqval;
}

int BinarySearch_fraglength(double* SearchArray,int low, int high, double key){
    //fprintf(stderr,"first element %lf\n",SearchArray[low]);
    int ans = 0; 
    while (low <= high) {
        int mid = low + (high - low + 1) / 2;
        //fprintf(stderr,"test %lf\n",SearchArray[mid]);
        double midVal = SearchArray[mid];
 
        if (midVal < key) {
            ans = mid;
            low = mid + 1;
        }
        else if (midVal > key) {

            high = mid - 1;
        }
        else if (midVal == key) {
 
            high = mid - 1;
        }
    }
 
    return ans+1;
}

void FragArray(int& number,int*& Length, double*& Frequency,const char* filename){
  //int LENS = 4096;
  //int* Frag_len = new int[LENS];
  //double* Frag_freq = new double[LENS];
  int n =1;

  gzFile gz = Z_NULL;
  char buf[LENS];
  
  gz = gzopen(filename,"r");
  assert(gz!=Z_NULL);
  Length[0] = 0; Frequency[0] = (float) 0;
  while(gzgets(gz,buf,LENS)){
    Length[n] = atoi(strtok(buf,"\n\t ")); //before it was Frag_len[n]
    Frequency[n] = atof(strtok(NULL,"\n\t "));
    n++;
  }
  gzclose(gz); 

  number = n;
  //Length = Frag_len;
  //Frequency = Frag_freq;
  //delete[] Frag_len;
  //delete[] Frag_freq;
}

void printTime(FILE *fp){
  time_t rawtime;
  struct tm * timeinfo; 
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  fprintf (fp, "\t-> %s", asctime (timeinfo) );
}

void Header_func(htsFormat *fmt_hts,const char *outfile_nam,samFile *outfile,sam_hdr_t *header,faidx_t *seq_ref,int chr_total,int chr_idx_arr[],size_t genome_len,char CommandArray[1024],const char* version){
  // Creates a header for the bamfile. The header is initialized before the function is called //

  if (header == NULL) { fprintf(stderr, "sam_hdr_init");}
    
  // Creating header information
  
  char genome_len_buf[1024];
  for(int i=0;i<chr_total;i++){
    const char *name = faidx_iseq(seq_ref,chr_idx_arr[i]);
    //fprintf(stderr,"chromosomes added to bam header %d\n",chr_idx_arr[i]);

    int name_len =  faidx_seq_len(seq_ref,name);
    snprintf(genome_len_buf,1024,"%d", name_len);
    
    // reference part of the header, int r variable ensures the header is added
    int r = sam_hdr_add_line(header, "SQ", "SN", name, "LN", genome_len_buf, NULL);
    if (r < 0) { fprintf(stderr,"sam_hdr_add_line");}
    memset(genome_len_buf,0, sizeof(genome_len_buf));
    //snprintf(genome_len_buf,1024,"COMMAND ./ngsngs");
    //int r = sam_hdr_add_line(header, "TEST COMMAND", genome_len_buf, NULL);
    memset(genome_len_buf,0, sizeof(genome_len_buf));
  }
  // Adding PG tag
  sam_hdr_add_pg(header,"NGSNGS","VN",version,"CL",CommandArray,NULL);
  // saving the header to the file
  if (sam_hdr_write(outfile, header) < 0) fprintf(stderr,"writing headers to %s", outfile_nam); //outfile
}

char* vcf_info(char genome_data[],int chr_sizes,const char* bcf_file,const char *chr_names[],const char* VarType,FILE *VarInfoFile){

  DNA_CAPITAL(genome_data);

  htsFile *bcf_obj = NULL;
  bcf_hdr_t *bcf_head = NULL;
  // "/home/wql443/WP1/data/Chr14_subset.vcf  "/home/wql443/WP1/data/test2.vcf"   "/home/wql443/WP1/data/chr14_full.vcf"
  bcf_obj = bcf_open(bcf_file, "r");
  hts_idx_t *idx = bcf_index_load(bcf_file);
  bcf_head = bcf_hdr_read(bcf_obj);
  fprintf(stderr, "File contains %i samples\n",bcf_hdr_nsamples(bcf_head)); //bcf.nsamples(fp)
  int nseq = 0;

  const char **seqnames = NULL;
  seqnames = bcf_hdr_seqnames(bcf_head, &nseq);

  // check if the headers match
  int chr_exit = 0;
  for (int vcf_i = 0; vcf_i < nseq; vcf_i++){
    if(strcasecmp(seqnames[vcf_i],chr_names[0])==0){continue;}//fprintf(stderr,"number of chromosomes %d and names %s\n",chr_i,chr_names[chr_i]);
    else{chr_exit++;}
  }
  if(chr_exit == nseq){
    fprintf(stderr,"Reference file chromosome %s is not identified in vcf/bcf header\n",chr_names[0]);
    exit(0);
  }
  
  // Set sample number 

  int indiv = bcf_hdr_set_samples(bcf_head, "HG00096", 0);
  fprintf(stderr, "File contains %i samples\n",bcf_hdr_nsamples(bcf_head)); //bcf.nsamples(fp)

  //initializing variation sampling 
  bcf1_t *bcf_records = bcf_init();
  if (bcf_records == NULL){fprintf(stderr,"WARNING NO VARIANTS IN VCF FILE"); exit(0);}
  else if(bcf_records != NULL && bcf_read(bcf_obj, bcf_head, bcf_records)==0){
    // else if{bcf_records != NULL && bcf_read(bcf_obj, bcf_head, bcf_records)==0; 
    //bcf_read(bcf_obj, bcf_head, bcf_records)==0; //ignoring return value (int)
    size_t chr_pos_before;
    size_t chr_pos_after;
    size_t insert_total = 0;
    size_t del_total = 0;
      
    hts_itr_t *itr = bcf_itr_querys(idx, bcf_head,chr_names[0]);
    int record_indiv;
    size_t old_pos = 0;
    while ((record_indiv = bcf_itr_next(bcf_obj, itr, bcf_records)) == 0){
      fprintf(stderr,"WHILE LOOP\n");
      if(strcasecmp(bcf_hdr_id2name(bcf_head, bcf_records->rid),chr_names[0])==0){
        bcf_unpack((bcf1_t*)bcf_records, BCF_UN_ALL);

        int i, j, ngt, nsmpl = bcf_hdr_nsamples(bcf_head);
        char **dst = NULL;
        fprintf(stderr,"number of samples %d\n",nsmpl);
        int32_t *gt_arr = NULL, ngt_arr = 0;
        //ngt = bcf_get_genotypes(bcf_head, bcf_records, &gt_arr, &ngt_arr);
        //ngt = bcf_get_format_int32(hdr, rec, "GT", &gt, &ngt_arr);
        bcf_get_format_string;
        std::cout << ngt << std::endl;
        exit(0);
        fprintf(stderr,"INSIDE WhILE IF STATEMENT\n");
        size_t pos = (int) bcf_records->pos + 1;
        std::cout << chr_names[0] << " " << pos << std::endl;
        int n_allele = (int) bcf_records->n_allele;
        std::cout << n_allele << std::endl;
        char* ref = bcf_records->d.allele[0];
        std::cout << ref << std::endl;
        char* alt = bcf_records->d.allele[1];
        fprintf(stderr,"SNP\t%s\t%zu\t%s\n",chr_names[0],pos,ref);

        char seqbefore[1024];
        char seqafter[1024];
        char seq_indel[1024];

        int rand_alt = rand() % (n_allele-1) +1;
        char* alt_select = bcf_records->d.allele[rand_alt];

        int ref_size = strlen(ref);
        int alt_size = strlen(alt_select);
        size_t new_pos = pos + insert_total + del_total;
        fprintf(stderr,"SNP\t%s\t%zu\t%s\t%zu\t%s\n",chr_names[0],pos,ref,new_pos,alt_select);
        if(alt[0] != '<'){ //not CNV
          if(new_pos<old_pos){continue;}//If SNP has been removed due to a previous deletion
          else{
            fprintf(stderr,"INSIDE NOT CNV STATEMENT\n");
            if (strcasecmp(VarType,"snp")==0||strcasecmp(VarType,"all")==0){
              fprintf(stderr,"SNP\t%s\t%zu\t%s\t%zu\t%s\n",chr_names[0],pos,ref,new_pos,alt_select);
              if(ref_size == 1 && alt_size == 1){
                fprintf(VarInfoFile,"SNP\t%s\t%zu\t%s\t%zu\t%s\n",chr_names[0],pos,ref,new_pos,alt_select);
                char ref_nt = genome_data[new_pos-1];
                if(ref_nt == *ref){genome_data[new_pos-1] = *alt_select;}
                else{fprintf(stderr,"VCF reference allele doesn't match with the position in the genome data\n");exit(0);}
              }
            }
            //insertion
            if (strcasecmp(VarType,"indel")==0||strcasecmp(VarType,"all")==0){
              if(alt_size > ref_size){
                fprintf(VarInfoFile,"Indel\t%s\t%zu\t%s\t%zu\t%s\n",chr_names[0],pos,ref,new_pos,alt_select);
                instert_seq(genome_data, strlen(genome_data), alt_select,strlen(alt_select), new_pos-1);
                insert_total = insert_total + strlen(alt_select)-1;
              }
              //deletion
              else if(ref_size > alt_size){
                fprintf(VarInfoFile,"Indel\t%s\t%zu\t%s\t%zu\t%s\n",chr_names[0],pos,ref,new_pos,alt_select);
                delete_seq(genome_data, alt_size, ref_size, new_pos,alt_size);
                del_total = del_total - (ref_size-alt_size);
              }
            }
          }
        }
        old_pos = new_pos;
      }
      else{
        fprintf(stderr,"Reference chromosome %s has no variations in vcf file at chromosomes %s\n",chr_names[0],bcf_hdr_id2name(bcf_head, bcf_records->rid));
      }
    }
  }

  bcf_hdr_destroy(bcf_head);
  bcf_destroy(bcf_records); 
  bcf_close(bcf_obj);

  return genome_data;
}

char* full_genome_create(faidx_t *seq_ref,int chr_total,int chr_sizes[],const char *chr_names[],size_t chr_size_cumm[]){
  size_t genome_size = 0;
  chr_size_cumm[0] = 0;
  for (int i = 0; i < chr_total; i++){
    const char *chr_name = faidx_iseq(seq_ref,i);
    int chr_len = faidx_seq_len(seq_ref,chr_name);
    chr_sizes[i] = chr_len;
    chr_names[i] = chr_name;
    genome_size += chr_len;
    chr_size_cumm[i+1] = genome_size;
  }

  char* genome = (char*) malloc(sizeof(char) * (genome_size+chr_total+1));
  genome[0] = 0; //Init to create proper C string before strcat
  //chr_total
  for (int i = 0; i < chr_total; i++){

    const char *data = fai_fetch(seq_ref,chr_names[i],&chr_sizes[i]);
    //sprintf(&genome[strlen(genome)],data);
    //strcat(genome,data);  //Both gives conditional jump or move error
    if (data != NULL){
      sprintf(genome+strlen(genome),"%s",data); 
    }
    // several of the build in functions allocates memory without freeing it again.
    free((char*)data); //Free works on const pointers, so we have to cast into a const char pointer
  }
  return genome;
}

char* partial_genome_create(faidx_t *seq_ref,int chr_total,int chr_sizes[],const char *chr_names[],size_t chr_size_cumm[]){
  
  size_t genome_size = 0;
  chr_size_cumm[0] = 0;
  for (int i = 0; i < chr_total; i++){
    int chr_len = faidx_seq_len(seq_ref,chr_names[i]);
    //fprintf(stderr,"chr len %d\n",chr_len);
    chr_sizes[i] = chr_len;
    genome_size += chr_len;
    chr_size_cumm[i+1] = genome_size;
  }

  char* genome = (char*) malloc(sizeof(char) * (genome_size+chr_total));
  genome[0] = 0; //Init to create proper C string before strcat
  //chr_total
  for (int i = 0; i < chr_total; i++){

    const char *data = fai_fetch(seq_ref,chr_names[i],&chr_sizes[i]);
    //sprintf(&genome[strlen(genome)],data);
    //strcat(genome,data);  //Both gives conditional jump or move error
    if (data != NULL){
      sprintf(genome+strlen(genome),"%s",data);  
    }
    // several of the build in functions allocates memory without freeing it again.
    free((char*)data); //Free works on const pointers, so we have to cast into a const char pointer
  }
  return genome;
}

char* full_vcf_genome_create(faidx_t *seq_ref,int chr_total,int chr_sizes[],const char *chr_names[],size_t chr_size_cumm[],const char* bcf_file,const char* VarType){
  fprintf(stderr,"BCF FILE %s\n",bcf_file);
  size_t genome_size = 0;
  chr_size_cumm[0] = 0;
  for (int i = 0; i < chr_total; i++){
    const char *chr_name = faidx_iseq(seq_ref,i);
    int chr_len = faidx_seq_len(seq_ref,chr_name);
    chr_sizes[i] = chr_len;
    chr_names[i] = chr_name;
    genome_size += chr_len;
    chr_size_cumm[i+1] = genome_size;
  }

  char* genome = (char*) malloc(sizeof(char) * (genome_size+chr_total));
  genome[0] = 0; //Init to create proper C string before strcat
  //chr_total
  FILE *VarInfoFile = fopen("Variant_pos_log.txt", "w");
  fprintf(VarInfoFile,"Type\tChr\tRef_Pos\tRef\tRead_Pos\tAlt\n");
  for (int i = 0; i < chr_total; i++){

    char *data = fai_fetch(seq_ref,chr_names[i],&chr_sizes[i]);

    if (data != NULL){
      fprintf(stderr,"INSIDE NULL DATA \t WITH VARTYPE %s and chr %s\n",VarType,chr_names[i]);
      //char* vcf_info(char genome_data[],int chr_sizes,const char* bcf_file,const char *chr_names[]){
      vcf_info(data, chr_sizes[i],bcf_file,&chr_names[i],VarType,VarInfoFile); //Chr14_subset4_bcf
      sprintf(genome+strlen(genome),"%s",data); 
    }
    // several of the build in functions allocates memory without freeing it again.
    free((char*)data); //Free works on const pointers, so we have to cast into a const char pointer
  }
  fprintf(stderr,"DONE WITH DATA FOR LOOP \n");
  fclose(VarInfoFile);
  return genome;
}

//! Allocate workspace for random-number sampling.
ransampl_ws* ransampl_alloc( int n )
{
    ransampl_ws *ws;
    if ( !(ws = (ransampl_ws*) malloc( sizeof(ransampl_ws) )) ||
         !(ws->alias = (int*) malloc( n*sizeof(int) )) ||
         !(ws->prob = (double*) malloc( n*sizeof(double) )) ) {
        fprintf( stderr, "ransampl: workspace allocation failed\n" );
        exit(ENOMEM); // ENOMEM
    }
    ws->n = n;
    return ws;
}

//! Initialize workspace by precompute alias tables from given probabilities.
void ransampl_set( ransampl_ws *ws, const double* p )
{
    const int n = ws->n;
    int i, a, g;

    // Local workspace:
    double *P;
    int *S, *L;
    if ( !(P = (double*) malloc( n*sizeof(double) ) ) ||
         !(S = (int*) malloc( n*sizeof(int) ) ) ||
         !(L = (int*) malloc( n*sizeof(int) ) ) ) {
        fprintf( stderr, "ransampl: temporary allocation failed\n" );
        exit(ENOMEM);
    }

    // Normalise given probabilities:
    double sum=0;
    for ( i=0; i<n; ++i ) {
        if( p[i]<0 ) {
            fprintf( stderr, "ransampl: invalid probability p[%i]<0\n", i );
            exit(EINVAL);
        }
        sum += p[i];
    }
    if ( !sum ) {
        fprintf( stderr, "ransampl: no nonzero probability\n" );
        exit(EINVAL);
    }
    for ( i=0; i<n; ++i )
        P[i] = p[i] * n / sum;

    // Set separate index lists for small and large probabilities:
    int nS = 0, nL = 0;
    for ( i=n-1; i>=0; --i ) {
        // at variance from Schwarz, we revert the index order
        if ( P[i]<1 )
            S[nS++] = i;
        else
            L[nL++] = i;
    }

    // Work through index lists
    while ( nS && nL ) {
        a = S[--nS]; // Schwarz's l
        g = L[--nL]; // Schwarz's g
        ws->prob[a] = P[a];
        ws->alias[a] = g;
        P[g] = P[g] + P[a] - 1;
        if ( P[g] < 1 )
            S[nS++] = g;
        else
            L[nL++] = g;
    }

    while ( nL )
        ws->prob[ L[--nL] ] = 1;

    while ( nS )
        // can only happen through numeric instability
        ws->prob[ S[--nS] ] = 1;

    // Cleanup:
    free( P );
    free( S );
    free( L );
}

int ransampl_draw2( ransampl_ws *ws,double r1, double r2)
{   
    //fprintf(stderr,"%lf\t%lf\n",r1,r2);
    const int i = (int) (ws->n * r1);
    return r2 < ws->prob[i] ? i : ws->alias[i];
}

//! Free the random-number sampling workspace.
void ransampl_free( ransampl_ws *ws )
{
    free( ws->alias );
    free( ws->prob );
    free( ws );
}

ransampl_ws ***ReadQuality(char *ntqual, double *ErrProb, int ntcharoffset,const char *freqfile,unsigned long &readcycle){
  ransampl_ws ***dists = new ransampl_ws**[5];

  std::vector<char *> all_lines;
  gzFile gz = Z_NULL;
  assert(((gz = gzopen(freqfile,"rb")))!=Z_NULL);
  char buf[LENS];
  while(gzgets(gz,buf,LENS))
    all_lines.push_back(strdup(buf));
  gzclose(gz);

  //fprintf(stderr,"All lines: %lu\n",all_lines.size());

  unsigned long readcyclelength = (all_lines.size()-2)/5; //all_lines.size()-1

  fprintf(stderr,"\t-> Inferred read cycle lengths: %lu\n",readcyclelength);
  readcycle = readcyclelength; 
  //loop over inputdata
  int nbins = -1;
  double probs[MAXBINS];
  for(int b=0;b<5;b++){
    dists[b] = new ransampl_ws *[readcyclelength];
    for(unsigned long pos = 0 ; pos<readcyclelength;pos++){
      int at = 0;
      probs[at++] = atof(strtok(all_lines[2+b*readcyclelength+pos],"\n\t ")); //1+b*readcyclelength+pos
      char *tok = NULL;
      while(((tok=strtok(NULL,"\n\t ")))){
	      probs[at++] = atof(tok);
	      assert(at<MAXBINS);
      }
      if(nbins==-1){
	      nbins = at;
	      //fprintf(stderr,"Number of qualities/bins in inputfile: %d\n",nbins);
      }
      if(nbins!=at){
	      fprintf(stderr,"Problems, number of columns is different nbins: %d at: %d\n",nbins,at);
	      exit(0);
      }
      dists[b][pos] =  ransampl_alloc( nbins );
      ransampl_set(dists[b][pos],probs);
    }
  }
  //fprintf(stderr,"\t-> ransampl_ws done\n");

  //printf(all_lines[0]);
  int qualidx = 1;
  //extract the first token
  ntqual[0] = (char) (atoi(strtok(all_lines[0],"\n\t "))+ntcharoffset);
  char *qualtok = NULL;
  //extract the next
  while(((qualtok=strtok(NULL,"\n\t ")))){ntqual[qualidx++] = (char) (atoi(qualtok)+ntcharoffset);}
  
  int Err_idx = 1;
  ErrProb[0] = atof(strtok(all_lines[1],"\n\t"));
  char *Errtok = NULL;
  while ((Errtok=strtok (NULL,"\n\t"))){ErrProb[Err_idx++] = (double) atof(Errtok);}
  
  //strdup function allocate necessary memory to store the sourcing string implicitly, i need to free the returned string
  for (unsigned long i = 0; i < all_lines.size(); i++){free(all_lines[i]);}
  //fprintf(stderr,"\t-> Before return in READQUAL FUNC\n");
  return dists;
}