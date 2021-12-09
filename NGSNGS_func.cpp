#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <cmath>

#include <algorithm> //std::reverse

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>

#include <pthread.h>

#include "NGSNGS_func.h"

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

void reverseChar(char* str) {
    std::reverse(str, str + strlen(str));
}

int Random_geometric_k(unsigned int  seed, const double p)
{
  double u = ((double) rand_r(&seed)/ RAND_MAX);
  double k;

  if (p == 1){k = 1;}
  else if(p == 0){k=0;}
  else{k = log (u) / log (1 - p);}

  return floor(k);
}

double myrand(unsigned int persistent){
  return ((double) rand_r(&persistent)/ RAND_MAX);
}

void SimBriggsModel(char* reffrag, char* frag, int L, double nv, double lambda, double delta_s, double delta, unsigned int seed){
    int l = 0;
    int r = L-1;
    while (l+r > L-2){
        l = 0;
        r = 0;
        double u_l = ((double) rand_r(&seed)/ RAND_MAX);//uniform();
        double u_r = ((double) rand_r(&seed)/ RAND_MAX);//uniform();
        //std::cout << u_l << " " << u_r <<"\n";
        /*std::default_random_engine generator1(rand()%30000+1);
        std::geometric_distribution<int> distribution1(lambda);
        std::default_random_engine generator2(rand()%30000+1);
        std::geometric_distribution<int> distribution2(lambda);*/
        
        if (u_l > 0.5){
            l = Random_geometric_k(seed+0,0.36);//distribution1(generator1);
            //fprintf(stderr,"U_L LEI %d\n",l);
            //cout << l << "\n";
        }
        if (u_r > 0.5){
          //fprintf(stderr,"U_R");
          r = Random_geometric_k(seed+1,0.36); //distribution1(generator2);
          //fprintf(stderr,"U_R LEI %d\n",r);
          //fprintf(stderr,"U_R %lf", r);
          //cout << r << "\n";
        }
    }
    //cout << l << " " << r <<"\n";
    for (int i = 0; i<l; i++){
        if (reffrag[i] == 'C' || reffrag[i] == 'c' ){
            double u = ((double) rand_r(&seed)/ RAND_MAX);//uniform();
            if (u < delta_s){
                frag[i] = 'T';
            }else{
                frag[i] = 'C';
            }
        }else{
            frag[i] = reffrag[i];
        }
    }
    for (int i = 0; i < r; i++){
        if (reffrag[L-i-1] == 'G' || reffrag[L-i-1] == 'g'){
            double u = ((double) rand_r(&seed)/ RAND_MAX);//uniform();
            if (u < delta_s){
                frag[L-i-1] = 'A';
            }else{
                frag[L-i-1] = 'G';
            }
        }else{
            frag[L-i-1] = reffrag[L-i-1];
        }
    }
    double u_nick = ((double) rand_r(&seed)/ RAND_MAX);//uniform();
    double d = nv/((L-l-r-1)*nv+1-nv);
    int p_nick = l;
    double cumd = 0;
    while (u_nick > cumd & p_nick < L-r-1){
        cumd += d;
        p_nick +=1;
    }
    for (int i = l; i < L-r; i++){
        if ((reffrag[i] == 'C' || reffrag[i] == 'c') && i<=p_nick){
            double u = ((double) rand_r(&seed)/ RAND_MAX);//uniform();
            if (u < delta){
                frag[i] = 'T';
            }else{
                frag[i] = 'C';
            }
        }else if ((reffrag[i] == 'G' || reffrag[i] == 'g') && i>p_nick){
            double u = ((double) rand_r(&seed)/ RAND_MAX);//uniform();
            if (u < delta){
                frag[i] = 'A';
            }else{
                frag[i] = 'G';
            }
        }else{
            frag[i] = reffrag[i];
        }
    }
}

const char* Error_lookup(double a,double err[6000],int nt_offset, int read_pos,const char* OutputFormat){

  //const char* nt_qual[8] = {"#", "\'", "0", "7" ,"<", "B", "F","I"}; // 35,39,48,55,66,70,73 //then withouth & in nt_qual[0]
  //const char* nt_qual[8] = {"#", "\'", "0", "7" ,"<", "B", "F","I"};//{"#", "\'", "0", "7" ,"<", "B", "F","I"};

  if (OutputFormat == "bam"){char nt_qual[8] = {2,6,15,22,27,33,37,40};}
  if (OutputFormat == "fq"){char nt_qual[8] = {35,39,48,55,66,70,73};}
  
  char nt_qual[8] = {2,6,15,22,27,33,37,40};
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
  int LENS = 6000;
  char buf[LENS];
  gzFile gz = Z_NULL;
  gz = gzopen(filename,"r");
  assert(gz!=Z_NULL);
  int i = 0;
  while(gzgets(gz,buf,LENS)){
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

void Read_Qual_new(char *seq,char *qual,unsigned int seed,double* freqval,const char* OutputFormat){
  //fprintf(stderr,"INSIDE FUNCITON NOW \n");
  int Tstart = 150;
  int Gstart = 300;
  int Cstart = 450;
  int Nstart = 600;

  //srand(time( NULL ));
  //srand(seed);

  //srand48(seed);
  int seqlen = strlen(seq);

  for (int row_idx = 0; row_idx < seqlen; row_idx++){
    //fprintf(stderr,"index %d \n",row_idx);
    double r = ((double) rand_r(&seed)/ RAND_MAX);
    //double r = 0.43; // drand48();//0.43;//(double) rand()/RAND_MAX;
    //fprintf(stderr,"random value %lf \n",r);
    switch(seq[row_idx]){
      case 'A':
      case 'a':
        strncat(qual, Error_lookup(r,freqval,0,row_idx,OutputFormat), 1);
        break;
      case 'T':
      case 't':
        strncat(qual, Error_lookup(r,freqval,Tstart,row_idx,OutputFormat), 1);
        break;  
      case 'G':
      case 'g':
        strncat(qual, Error_lookup(r,freqval,Gstart,row_idx,OutputFormat), 1);
        break;
      case 'C':
      case 'c':
        strncat(qual, Error_lookup(r,freqval,Cstart,row_idx,OutputFormat), 1);
        break;
      case 'N':
      case 'n':
        strncat(qual, Error_lookup(r,freqval,Nstart,row_idx,OutputFormat), 1);
        break;
    }
  }
  //fprintf(stderr,"EXITING FUNCITON NOW \n");
}

int BinarySearch_fraglength(double* SearchArray,int low, int high, double key)
{
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
  int LENS = 4096;
  int* Frag_len = new int[LENS];
  double* Frag_freq = new double[LENS];
  int n =0;

  gzFile gz = Z_NULL;
  char buf[LENS];
  
  gz = gzopen(filename,"r");
  assert(gz!=Z_NULL);
  while(gzgets(gz,buf,LENS)){
    Frag_len[n] = atoi(strtok(buf,"\n\t "));
    Frag_freq[n] = atof(strtok(NULL,"\n\t "));
    n++;
  }
  gzclose(gz); 

  number = n;
  Length = Frag_len;
  Frequency = Frag_freq;
}

void printTime(FILE *fp){
  time_t rawtime;
  struct tm * timeinfo; 
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  fprintf (fp, "\t-> %s", asctime (timeinfo) );
}

