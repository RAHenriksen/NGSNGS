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

#include <random>
#include <iterator>
#include <cmath>

// g++ Char_array.cpp -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl -std=c++11
// ./a.out fafa fa 1607595921


double** create2DArray(int height, int width, const char* filename){
  // creates the 2d object with all the frequency values for a given positon of the nt
  // used to create the nt qual strings.
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

int Qual_random(double *a,std::default_random_engine &gen){
  //creates a random sequence of nt qualities based on a frequency distribution

  char Qualities[] = {'#', '\'', '0', '7' ,'<', 'B', 'F','I','\0'};
  int ASCII_qual[] = {35,39,48,55,60,66,70,73}; //Im choosing the ascii values, since im casting them into string later on

  std::discrete_distribution<> d({a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7]});
  return ASCII_qual[d(gen)];
  //Read_qual += ASCII_qual[d(gen)];
  //return Read_qual;
}

std::string Read_Qual(char *seq,double** Array2d,std::default_random_engine &gen){
  std::string qual;
  int read_length = strlen(seq);
  //std::cout << "read qual seq " << seq << std::endl;
  //std::cout << " length " << read_length << std::endl;
  int Tstart = 150;
  int Gstart = 300;
  int Cstart = 450;
  
  for (int row_idx = 0; row_idx < read_length; row_idx++){
    //std::cout << "row indx" << seq[row_idx] << std::endl;
    //std::cout << Qual_random(Array2d[row_idx],gen);
    switch(seq[row_idx]){
      // iterates through every element in seq and creates and extract qual from the given line 
        // of the sequence qual profile. based on the lines number given the Tstart etc..
      case 'A':
        qual += Qual_random(Array2d[row_idx],gen);
        break;
      case 'a':
        qual +=  Qual_random(Array2d[row_idx],gen);
        break;
      case 'T':
        qual +=  Qual_random(Array2d[row_idx + Tstart],gen);
        break;
      case 't':
        qual += Qual_random(Array2d[row_idx + Tstart],gen);
        break;  
      case 'G':
        qual +=  Qual_random(Array2d[row_idx + Gstart],gen);
        break;
      case 'g':
        qual +=  Qual_random(Array2d[row_idx + Gstart],gen);
        break;
      case 'C':
        qual +=  Qual_random(Array2d[row_idx + Cstart],gen);
        break;
      case 'c':
        qual +=  Qual_random(Array2d[row_idx + Cstart],gen);
        break;
    }
  }
  return qual;
}


// THIS IS NOT GOOD EITHER  -> THEN FOR EACH READ THE FILL DISTRIBUTION IS CREATED !!
std::string Read_Qual2(char *seq,double** Array2d,std::default_random_engine &gen){
  std::string qual;
  int read_length = strlen(seq);
  
  //Distribution for all 600 lines in 2
  int ASCII_qual[] = {35,39,48,55,60,66,70,73};
  std::discrete_distribution<> D[600];

  for (int row_idx = 0; row_idx < 600; row_idx++){
    std::discrete_distribution<> d({Array2d[row_idx][0],Array2d[row_idx][1],Array2d[row_idx][2],Array2d[row_idx][3],
                                    Array2d[row_idx][4],Array2d[row_idx][5],Array2d[row_idx][6],Array2d[row_idx][7]});  
    D[row_idx] = d;
  }
  std::cout << "size " << sizeof(D)/sizeof(D[0]) << std::endl;
  int Tstart = 150;
  int Gstart = 300;
  int Cstart = 450;
  
  for (int row_idx = 0; row_idx < read_length; row_idx++){
    //std::cout << "row indx" << seq[row_idx] << std::endl;
    //std::cout << Qual_random(Array2d[row_idx],gen);
    switch(seq[row_idx]){
      // iterates through every element in seq and creates and extract qual from the given line 
        // of the sequence qual profile. based on the lines number given the Tstart etc..
      case 'A':
        qual += ASCII_qual[D[row_idx](gen)];
        break;
      case 'a':
        qual += ASCII_qual[D[row_idx](gen)];
        break;
      case 'T':
        qual += ASCII_qual[D[row_idx + Tstart](gen)];
        break;
      case 't':
        qual += ASCII_qual[D[row_idx + Tstart](gen)];
        break;  
      case 'G':
        qual += ASCII_qual[D[row_idx + Gstart](gen)];
        break;
      case 'g':
        qual += ASCII_qual[D[row_idx + Gstart](gen)];
        break;
      case 'C':
        qual += ASCII_qual[D[row_idx + Cstart](gen)];
        break;
      case 'c':
        qual += ASCII_qual[D[row_idx + Cstart](gen)];
        break;
    }
  }
  return qual;
}

//-----------------------------------------//
// Here it works by seperating the distribution to have it once by Qual_dist, where the function is a pointer type,
// So i have to make the initialization outside...

std::discrete_distribution<> D[600];

std::discrete_distribution<>* Qual_dist(double** Array2d){
  // function returning a pointer to a discrete distribution. But the adress of the object dissappears after return statement.
  // as such the distribution we are returning are initialized outside of the function. So we dont have to create it multiple times
  
  // Distribution for all 600 lines in the 2D array with each line being a position in a read for a nt and each column being a
  // ascii value for an nt quality which should be cast into string

  for (int row_idx = 0; row_idx < 600; row_idx++){
    std::discrete_distribution<> d({Array2d[row_idx][0],Array2d[row_idx][1],Array2d[row_idx][2],Array2d[row_idx][3],
                                    Array2d[row_idx][4],Array2d[row_idx][5],Array2d[row_idx][6],Array2d[row_idx][7]});  
    D[row_idx] = d;
  }
  //std::cout << "size " << sizeof(D)/sizeof(D[0]) << std::endl;
  return D;
}

// Pass the distribution D as a pointer. but i am not actually dereferencing it. hmm..
std::string Read_Qual3(char *seq,std::discrete_distribution<>*Dist,std::default_random_engine &gen){
  std::string qual;
  int ASCII_qual[] = {35,39,48,55,60,66,70,73};
  int read_length = strlen(seq);

  int Tstart = 150;
  int Gstart = 300;
  int Cstart = 450;
  
  for (int row_idx = 0; row_idx < read_length; row_idx++){
    //std::cout << "row indx" << seq[row_idx] << std::endl;
    //std::cout << Qual_random(Array2d[row_idx],gen);
    switch(seq[row_idx]){
      // iterates through every element in seq and creates and extract qual from the given line 
        // of the sequence qual profile. based on the lines number given the Tstart etc..
      case 'A':
        qual += ASCII_qual[Dist[row_idx](gen)];
        break;
      case 'a':
        qual += ASCII_qual[Dist[row_idx](gen)];
        break;
      case 'T':
        qual += ASCII_qual[Dist[row_idx + Tstart](gen)];
        break;
      case 't':
        qual += ASCII_qual[Dist[row_idx + Tstart](gen)];
        break;  
      case 'G':
        qual += ASCII_qual[Dist[row_idx + Gstart](gen)];
        break;
      case 'g':
        qual += ASCII_qual[Dist[row_idx + Gstart](gen)];
        break;
      case 'C':
        qual += ASCII_qual[Dist[row_idx + Cstart](gen)];
        break;
      case 'c':
        qual += ASCII_qual[Dist[row_idx + Cstart](gen)];
        break;
    }
  }
  return qual;
}

std::discrete_distribution<> Distfunc(double** Array2d,int size){
  std::discrete_distribution<> dist[600];
  
  for (int row_idx = 0; row_idx < 600; row_idx++){
    std::discrete_distribution<> d({Array2d[row_idx][0],Array2d[row_idx][1],Array2d[row_idx][2],Array2d[row_idx][3],
                                    Array2d[row_idx][4],Array2d[row_idx][5],Array2d[row_idx][6],Array2d[row_idx][7]});  
    dist[row_idx] = d;
  }

  std::cout << "size " << sizeof(dist)/sizeof(dist[0]) << std::endl;
  std::cout << " test " << dist[0] << std::endl;
  return *dist;
}

void Distfunc2(double** Array2d,std::discrete_distribution<> dist[],int size){
  for (int row_idx = 0; row_idx < size; row_idx++){
    std::discrete_distribution<> d({Array2d[row_idx][0],Array2d[row_idx][1],Array2d[row_idx][2],Array2d[row_idx][3],
                                    Array2d[row_idx][4],Array2d[row_idx][5],Array2d[row_idx][6],Array2d[row_idx][7]});  
    dist[row_idx] = d;
  }
}

int main(){
  std::random_device rd;
  std::default_random_engine gen(rd()); 
  std::string qual;
  char seqtest[1024] = "AGCTaGCTAGTCAGCTaGCTAGTCAGCTaGCTAGTCAGCTaGCTAGTCAGCTaGCTAGTC";
  char seqtest2[1024] = "CAGCAGCAGAGCGATTAGA";

  double** a = create2DArray(600, 8,"/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt");
  // rækkefølgen er 0-150 "A" -> 150 - 300 "T" ->  300-450 "G" -> 450 - 600 "C"
  
  int distsize = 600;
  std::discrete_distribution<> test[distsize];
  Distfunc2(a,test,distsize);
  std::cout << "lol \n" << test[0] << std::endl;
  std::cout << "lol \n" << test[1] << std::endl;

  std::cout << Read_Qual3(seqtest,test,gen) << std::endl;
  std::cout << Read_Qual3(seqtest2,test,gen) << std::endl;

  //test = Distfunc(a,600);
  //std::cout << "test \n " << test << std::endl;


  std::cout << "lol" << std::endl;
  Qual_dist(a);
  //std::cout << D[1](gen) << std::endl;
  //std::cout << Read_Qual3(seqtest,D,gen) << std::endl;
  //std::cout << Read_Qual3(seqtest2,D,gen) << std::endl;

  /*
  int ASCII_qual[] = {35,39,48,55,60,66,70,73};

  for (int i = 0; i < 5; i++){
      qual += ASCII_qual[D[i](gen)];
      std::cout << "sampling " << D[i](gen) << std::endl;
      std::cout << " qual " << ASCII_qual[D[i](gen)] << std::endl;
      std::cout << " string "<< qual << std::endl;
    }
  
  std::string qual;

  char seqtest[1024] = "AGCTaGCTAGTCAGCTaGCTAGTCAGCTaGCTAGTCAGCTaGCTAGTCAGCTaGCTAGTC";

  int read_length = strlen(seqtest);
  std::random_device rd;
  std::default_random_engine gen(rd()); 
  
  std::cout << Read_Qual2(seqtest,a,gen);
  */
  return 0;
}

/*
int main(){
  //std::ifstream file("Freq.txt");
  double** a = create2DArray(600, 8,"/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt");
  // rækkefølgen er 0-150 "A" -> 150 - 300 "T" ->  300-450 "G" -> 450 - 600 "C"
  std::string Read_qual;
  std::string qual;

  char seqtest[1024] = "AGCTaGCTAGTCAGCTaGCTAGTCAGCTaGCTAGTCAGCTaGCTAGTCAGCTaGCTAGTC";

  int Tstart = 150;
  int Gstart = 300;
  int Cstart = 450;

  int read_length = strlen(seqtest);
  std::random_device rd;
  std::default_random_engine gen(rd()); 

  int ASCII_qual[] = {35,39,48,55,60,66,70,73};
  std::discrete_distribution<> D[600];

  for (int row_idx = 0; row_idx < 600; row_idx++){
    std::discrete_distribution<> d({a[row_idx][0],a[row_idx][1],a[row_idx][2],a[row_idx][3],a[row_idx][4],a[row_idx][5],a[row_idx][6],a[row_idx][7]});  
    D[row_idx] = d;
  }

  //std::cout << ASCII_qual[d(gen)] << std::endl;
  //D[1] = d;
  for (int i = 0; i < 5; i++){
    qual += ASCII_qual[D[i](gen)];
    std::cout << "sampling " << D[i](gen) << std::endl;
    std::cout << " qual " << ASCII_qual[D[i](gen)] << std::endl;
    std::cout << " string "<< qual << std::endl;
  }
  
  //std::cout << "dist " << D[1] << std::endl;
  //int b = D[1](gen);
  //std::cout << "sampling " << b << std::endl;
  //std::cout << "qual " << ASCII_qual[b] << std::endl;

  std::cout << "--------" << std::endl;
  for (int row_idx = 0; row_idx < 6; row_idx++){
    std::cout << Qual_random(a[row_idx],gen) << std::endl;
  }

  return 0;
}
*/