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
      infile >> array2D[h][w];}
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
        strncat(qual, nt_qual[0], 1);;
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
        strncat(qual, &nt_qual[Dist[row_idx + Cstart](gen)], 1);;
        break;
    }
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
    }
    ++seq;
  }
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
  
  for (int nt_idx = 0; nt_idx < read_length; nt_idx++){
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

void Size_freq_dist(std::ifstream &infile, std::discrete_distribution<> dist[]){
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
  std::random_device rd;
  std::default_random_engine gen(rd());
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