#include <algorithm>
#include <cstdio>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <zlib.h>

#include <cstdlib>
#include <ctime>

#include <cstdio>
#include <cassert>
#include <cstdint>

#include <random>
#include <iterator>
#include <cmath>

#include <thread>         // std::thread
#include <mutex>          // std::mutex mtx;

#include "SimulAncient_func.h"

void Deamin_char(char* str,char nt[],int seed,double alpha,double beta,int start,int end){   
  // Deamination of nucleotides
  std::vector<int> Index_vec;
  std::srand(seed+std::time(nullptr));
  std::random_device rd;
  std::default_random_engine generator(rd());
  std::gamma_distribution<double> distr(alpha,beta);

  int i = strcspn(str,nt);
  Index_vec.push_back(i);
  unsigned int vector_size;
  //consider just using signed int, just because they would not be negative..
  while(i < end) {
    int tmp = strcspn(str+i+1,nt);
    i += tmp + 1;
      Index_vec.push_back(i);
  }
  vector_size = Index_vec.size(); //Index_vec
    
  for (unsigned int i = 0; i < vector_size; i++){
    if (Index_vec.at(i) == int(distr(generator))) {
      //remember to create an input for the nt so it works for both 5' and 3' 
      str[Index_vec.at(i)]='T';
      }
		else {
      continue;
		}
  }
}

double** create2DArray(const char* filename,int width,int height){
  /* create 2d objects for a given error profile, with rows being position in read for specific nt
  and the cells being the frequency values */

  std::ifstream infile(filename);
  double** array2D = 0;
  array2D = new double*[height];
  for (int h = 0; h < height; h++){
    array2D[h] = new double[width];
    for (int w = 0; w < width; w++){
      infile >> array2D[h][w];
    }
  }
  infile.close();
  return array2D;
}

void Read_Qual2(char *seq,char *qual,std::discrete_distribution<>*Dist,std::default_random_engine &gen){
  /* creating the nucleotide quality string for fastq format and bam format using char array*/

  const char* nt_qual[8] = {"#", "\'", "0", "7" ,"<", "B", "F","I"};

  int read_length = strlen(seq);

  //the line offset for the distribution *Dist created from the 600*8 2Darray
  int Tstart = 150;
  int Gstart = 300;
  int Cstart = 450;
  
  // Dist[row_idx](gen) returns values between 0-7 (columns) sampling one out of 8 nt qual from the 
  // for each line (read positions) in the read using the error profile created using the 2d array 
  for (int row_idx = 0; row_idx < read_length; row_idx++){
    switch(seq[row_idx]){
      case 'A':
      case 'a':
        strncat(qual, nt_qual[Dist[row_idx](gen)], 1);
        break;
      case 'T':
      case 't':
        strncat(qual, nt_qual[Dist[row_idx + Tstart](gen)], 1);
        break;  
      case 'G':
      case 'g':
        strncat(qual, nt_qual[Dist[row_idx + Gstart](gen)], 1);
        break;
      case 'C':
      case 'c':
        strncat(qual, nt_qual[Dist[row_idx + Cstart](gen)], 1);
        break;
      case 'N':
      case 'n':
        strncat(qual, nt_qual[0], 1);
        break;
    }
  }
}

void Bam_baseQ(char *seq,char *qual,std::discrete_distribution<>*Dist,std::default_random_engine &gen){
  /* creating the nucleotide quality string for fastq format and bam format using char array*/

  //const char* nt_qual[8] = {"#", "\'", "0", "7" ,"<", "!", "%","("};

  // creating a char array and immediately casting the ASCII values into char
  char nt_qual[8] = {2,6,15,22,27,33,37,40};
  int read_length = strlen(seq);

  //the line offset for the distribution *Dist created from the 600*8 2Darray
  int Tstart = 150;
  int Gstart = 300;
  int Cstart = 450;
  
  // Dist[row_idx](gen) returns values between 0-7 (columns) sampling one out of 8 nt qual from the 
  // for each line (read positions) in the read using the error profile created using the 2d array 
  for (int row_idx = 0; row_idx < read_length; row_idx++){
    switch(seq[row_idx]){
      case 'A':
      case 'a':
        // const char* is a pointer to a char, so i am using reference & for the strcat
        strncat(qual, &nt_qual[Dist[row_idx](gen)], 1);
        break;
      case 'T':
      case 't':
        strncat(qual, &nt_qual[Dist[row_idx + Tstart](gen)], 1);
        break;  
      case 'G':
      case 'g':
        strncat(qual, &nt_qual[Dist[row_idx + Gstart](gen)], 1);
        break;
      case 'C':
      case 'c':
        strncat(qual, &nt_qual[Dist[row_idx + Cstart](gen)], 1);
        break;
      case 'N':
      case 'n':
        strncat(qual, &nt_qual[Dist[row_idx + Cstart](gen)], 1);
        break;
    }
  }
}


void Read_Qual_final(char *seq,char *qual,std::discrete_distribution<>*Dist,std::default_random_engine &gen,int qual_offset){
  /* creating the nucleotide quality string for fastq format and bam format using char array*/

  const char* nt_qual[8] = {"#", "\'", "0", "7" ,"<", "B", "F","I"};

  int read_length = strlen(seq);

  //the line offset for the distribution *Dist created from the 600*8 2Darray
  int Tstart = 150;
  int Gstart = 300;
  int Cstart = 450;
  //int Nstart = 600;
  
  for (int row_idx = 0; row_idx < read_length; row_idx++){
    switch(seq[row_idx]){
      case 'A':
      case 'a':
        strncat(qual, nt_qual[Dist[row_idx](gen)], 1);
        break;
      case 'T':
      case 't':
        strncat(qual, nt_qual[Dist[row_idx + Tstart](gen)], 1);
        break;  
      case 'G':
      case 'g':
        strncat(qual, nt_qual[Dist[row_idx + Gstart](gen)], 1);
        break;
      case 'C':
      case 'c':
        strncat(qual, nt_qual[Dist[row_idx + Cstart](gen)], 1);
        break;
      case 'N':
      case 'n':
        strncat(qual, nt_qual[0], 1);
        break;
    }
  }
}
/*
a=0, c=1, g=2, t=3, n=4

int revcom[5] = {3,2,1,0,4}

revcom[c]

char c= 'G'

revcom[c];

char revcom2[255];
memset(revcom2,'n',255);
revcom2['A'] = 't';
revcom2['C'] = 'g';*/


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

void Qual_dist(double** Array2d,std::discrete_distribution<> dist[],int size){
  /* creating a discrete distribution of nucleotide qualitie for each line in 2D array */
  for (int row_idx = 0; row_idx < size; row_idx++){
    std::discrete_distribution<> d({Array2d[row_idx][0],Array2d[row_idx][1],Array2d[row_idx][2],Array2d[row_idx][3],
                                    Array2d[row_idx][4],Array2d[row_idx][5],Array2d[row_idx][6],Array2d[row_idx][7]});  
    dist[row_idx] = d;
  }
}

void Seq_err(double** Array2d,std::discrete_distribution<> nt_sub[],int size){
  /* Similar to qual_dist creating a nucleotide distribution for each line in 2D array*/
  for (int row_idx = 0; row_idx < size; row_idx++){
    std::discrete_distribution<> d({Array2d[row_idx][0],
                                    Array2d[row_idx][1],
                                    Array2d[row_idx][2],
                                    Array2d[row_idx][3]});  
    nt_sub[row_idx] = d;
  }
}

void Ill_err(char *seq,std::discrete_distribution<>*Dist,std::default_random_engine &gen){
  /* Similar to Read_Qual2 but for creating substitutions in the string */
  int read_length = strlen(seq);
  int Tstart = 70;
  int Gstart = 140;
  int Cstart = 210;
  const char LookUp_nt[4] = {'A','T','G','C'};
  
  int top_idx;
  
  if (read_length <= 70){top_idx = read_length;} //fragments shorter than 70 can have potential sub in all pos
  else{top_idx = 70;} //for fragments longer than 70 then i'll just change substituions for first 70 nt
  
  for (int nt_idx = 0; nt_idx < top_idx; nt_idx++){
    switch(seq[nt_idx]){
      case 'A':
      case 'a':
        seq[nt_idx] = LookUp_nt[Dist[nt_idx](gen)];
        break;
      case 'T':
      case 't':
        seq[nt_idx] = LookUp_nt[Dist[nt_idx+Tstart](gen)];
        break;  
      case 'G':
      case 'g':
        seq[nt_idx] = LookUp_nt[Dist[nt_idx+Gstart](gen)];
        break;
      case 'C':
      case 'c':
        seq[nt_idx] = LookUp_nt[Dist[nt_idx+Cstart](gen)];
        break;
    }
  }
}

void Size_freq_dist(std::ifstream &infile, std::discrete_distribution<> dist[],int seed){
  // Creates a distribution with the frequencies for specific fragment lengths from the same file as used in Size_select_dist //
  int fraqment_length;
  double freq;
  
  std::vector<double> Freq_vec;

  // takes the freq from the file and creates a vector for that
  while (infile >> fraqment_length >> freq)
  {
    Freq_vec.push_back(freq);
  }

  // converts the vector of frequencies to a distribution
  //std::random_device rd;
  //std::default_random_engine gen(rd());
  std::default_random_engine gen(seed);
  std::discrete_distribution<> d(Freq_vec.begin(), Freq_vec.end());
  dist[1] = d;
}

int *Size_select_dist(std::ifstream &infile){
  // Creates an array of the fragments lengths from the same file as used in Size_freq_dist//
  int fraqment_length;
  double freq;
  //int sizearray[1024] = {};

  int* sizearray = new int[1024];
  int i = 0;

  // takes the fragment_length from the file and creates a array of that
  while (infile >> fraqment_length >> freq)
  {
    sizearray[i] = fraqment_length;
    i++;
  }
  return sizearray;
  
}

double uniform()
{
    /*
     U(0,1): AS 183: Appl. Stat. 31:188-190
     Wichmann BA & Hill ID.  1982.  An efficient and portable
     pseudo-random number generator.  Appl. Stat. 31:188-190
     x, y, z are any numbers in the range 1-30000.  Integer operation up
     to 30323 required.
     
     Suggested to me by Ziheng Yang who also provided me with
     the source code used here.  I use it because it is both fast and portable.
     */
    static int x_rndu=11, y_rndu=23;int static z_rndu=137;

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

int Random_geometric_k(unsigned int  seed, const double p)
{
  double u = ((double) rand_r(&seed)/ RAND_MAX);
  int k;

  if (p == 1){k = 1;}
  else if(p == 0){k=0;}
  else{k = log (u) / log (1 - p);}

  return k;
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

const char* Error_lookup(double a,double err[6000],int nt_offset, int read_pos){
  //const char* nt_qual[8] = {"#", "\'", "0", "7" ,"<", "!", "%","("};
  const char* nt_qual[8] = {"#", "\'", "0", "7" ,"<", "B", "F","I"};
  int offset = ((nt_offset+read_pos)*8);
  //printf("offset %d \n", offset);
  const char* nt_out;
  if (a <= err[offset]){
    nt_out = nt_qual[0];
  }
  else if (err[offset] < a && a <= err[offset+1]){
    nt_out = nt_qual[1];
    }
  else if (err[offset+1] < a && a <= err[offset+2]){
    nt_out = nt_qual[2];
  }
  else if (err[offset+2] < a && a <= err[offset+3]){
    nt_out = nt_qual[3];
  }
  else if (err[offset+3] < a && a<= err[offset+4]){
    nt_out = nt_qual[4];
  }
  else if (err[offset+4] < a && a<= err[offset+5]){
    nt_out = nt_qual[5];
  }
  else if (err[offset+5] < a && a<= err[offset+6]){
    nt_out = nt_qual[6];
  }
  else if (err[offset+6] < a && a <= err[offset+7]){
    nt_out = nt_qual[7];
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

void Read_Qual_new(char *seq,char *qual,unsigned int seed,double* freqval){
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
        strncat(qual, Error_lookup(r,freqval,0,row_idx), 1);
        break;
      case 'T':
      case 't':
        strncat(qual, Error_lookup(r,freqval,Tstart,row_idx), 1);
        break;  
      case 'G':
      case 'g':
        strncat(qual, Error_lookup(r,freqval,Gstart,row_idx), 1);
        break;
      case 'C':
      case 'c':
        strncat(qual, Error_lookup(r,freqval,Cstart,row_idx), 1);
        break;
      case 'N':
      case 'n':
        strncat(qual, Error_lookup(r,freqval,Nstart,row_idx), 1);
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