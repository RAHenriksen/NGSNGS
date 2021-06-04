#include <cstdio>
#include <cassert>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <typeinfo>
#include <random>
#include <iterator>
#include <cmath>
#include <chrono>
#include <time.h>
#include <algorithm>

void coverage(int cov){

  const int N = 10;
  int frag_len[N] = {0};

  float cov_current = 0;

  srand(time(NULL));
  int start;
  int length;
  
  int nread = 0;
  
  while (cov_current < cov) {
    int sum = 0;
    int count = N;
    start = rand() % N;
    length = 1 + rand() % 10;
    fprintf(stderr,"start %i and length %i \n",start,length);
    for (int j = start; j < start+length; j++){
      frag_len[j] += 1;
    }
    for (size_t i = 0; i < N; i++){std::cout << frag_len[i] << ' ';}
    std::cout << std::endl;
    for(int i = 0; i<N ; i++){
      if (frag_len[i] == 0){count--;}
      else{sum+=frag_len[i];}
    }
    nread++;
    fprintf(stderr,"count %i and sum %i \n",count,sum);
    cov_current = (float) sum / (float) count;
    std::cout << cov_current << std::endl;
    std::cout << "number of reads "<< nread << std::endl;
    //std::cout << "cov 2 "<< (float) sum / (float) N << std::endl;
    std::cout << "-----------" << std::endl;
  }
}

int Size_freq_dist2(const char* filename){
  std::ifstream infile(filename);
  int fraqment_length;
  std::string line;
  double freq;
  int sizearray[1024];
  std::vector<double> Freq_vec;

  int i = 0;
  while (infile >> fraqment_length >> freq)
  {
    //Freq_vec.push_back(b);
    //std::cout << freq << std::endl;
    //std::cout << "-----------" << std::endl;
    Freq_vec.push_back(freq);
    sizearray[i] = fraqment_length;
    i++;
    /*std::cout << a << "  " << b << std::endl;
    std::cout << "i " << i << std::endl;
    std::cout << "------------" << std::endl;*/
  }
  //std::cout << "loop done " << std::endl;
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::discrete_distribution<> d(Freq_vec.begin(), Freq_vec.end());
  return sizearray[d(gen)];
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

int main(int argc,char **argv){
  //coverage(5);

  /*for (size_t i = 0; i < 1000; i++)
  {
    int FL = Size_freq_dist("Size_freq.txt");
    int RL = 150;
    std::cout << "FL " << FL << std::endl;
    if (FL > 2*RL)
    {
      std::cout << "Case 1" << std::endl;
    }
    else if (FL > RL && FL < 2*RL)
    {
      std::cout << "Case 2 " << std::endl;
    }
    else if (FL < RL)
    {
      std::cout << "case 3 " << std::endl;
    }
    
  }*/
  std::ifstream infile("Size_freq.txt");
  int* sizearray = Size_select_dist(infile);
  infile.close();
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::discrete_distribution<> Qualdistr1[2]; // we only need 1 distribution but for some reason it fails when its 1.
  
  // Hvis jeg bare benytter samme infile så fucker den op og vælger kun index 0.
  std::ifstream infile2("Size_freq.txt");
  Size_freq_dist(infile2,Qualdistr1); //creates the distribution of all the frequencies
  // Qualdistr1[1](gen) // udvælger index af min den frekvens der er udtrukket 0-index
  // Qualdistr1[1] // print hele distributionen
  infile2.close();
  for (size_t i = 0; i < 10; i++)
  {
    int a = Qualdistr1[1](gen);
    std::cout << " a " << a  << std::endl;
    std::cout << sizearray[a] << std::endl;
    std::cout << "----------" << std::endl;
  }
  //std::cout << sizearray[0] << std::endl;
  return 0;
}

// g++ Cov_test.cpp FaFa_thread.cpp -std=c++11 

/*


void frag(int len){

  const int N = 10;
  int frag_len[N] = {0};

  float cov = 0;

  srand(time(NULL));
  int start;
  int length;
  
  int nread = 0;
  
  while (nread < 10) {
    int sum = 0;
    int count = N;
    start = rand() % N;
    length = 1 + rand() % 10;
    fprintf(stderr,"start %i and length %i \n",start,length);
    for (int j = start; j < start+length; j++){
      frag_len[j] += 1;
    }
    for (size_t i = 0; i < N; i++){std::cout << frag_len[i] << ' ';}
    std::cout << std::endl;
    for(int i = 0; i<N ; i++){
      if (frag_len[i] == 0){count--;}
      else{sum+=frag_len[i];}
    }
    nread++;
    fprintf(stderr,"count %i and sum %i \n",count,sum);
    cov = (float) sum / (float) count;
    std::cout << cov << std::endl;
    std::cout << "number of reads "<< nread << std::endl;
    
    std::cout << "-----------" << std::endl;
  }
}
*/