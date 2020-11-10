#include <cstdio>
#include <cassert>
#include <htslib/faidx.h>
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

// I would like to create a function with TK's code since its optimal in case we wish to 
//simulate a given number of fragments

void random_seq(faidx_t *seq_ref){
  // choose a random sequence -> still ned to change it so it saves the output to a single file.
  int readlength=35;
  int nreads = 8;
  for(int i=0;i<nreads;i++){
    char buf[96];//assume maxlength for readid is 96bytes
    int whichref = lrand48() % faidx_nseq(seq_ref);
    fprintf(stderr,"\t-> Whichref: %d\n",whichref);
    const char *name = faidx_iseq(seq_ref,whichref);
    int name_len =  faidx_seq_len(seq_ref,name);
    fprintf(stderr,"\t-> name: \'%s\' name_len: %d\n",name,name_len);

    int start = lrand48() % name_len;
    int stop = start+readlength;
    if(stop>name_len)
      stop = name_len;
    snprintf(buf,96,"%s:%d-%d",name,start,stop);
    fprintf(stderr,"buf: %s\n",buf);
    fprintf(stdout,"%s\n+\n",buf);
    char *data = fai_fetch(seq_ref,name,&name_len);

    for(int i=start;i<stop;i++)
      fprintf(stdout,"%c",data[i]);
    fprintf(stdout,"\n");

    for(int i=start;i<stop;i++)
      fprintf(stdout,"F");
    fprintf(stdout,"\n");
  }
}

void Deamin_char(char* str,char nt[],int seed,
              double alpha=1.0,double beta=2.0,int start=0,int end=25)
{   // use & to pass by reference to the vector
    // Get the first occurrence
  std::vector<int> Index_vec;

  int i = strcspn(str,nt);
  Index_vec.push_back(i);

  while(i < end) {
    int tmp = strcspn(str+i+1,nt);
    i += tmp + 1;
      //std::cout << "i" <<tmp <<std::endl;
      Index_vec.push_back(i);
      //std::cout << "i " << i << std::endl;
  }

  for (int i = 0; i < Index_vec.size(); i++){
    std::srand(std::time(nullptr)+i);
    std::default_random_engine generator(seed);
    std::gamma_distribution<double> distr(alpha,beta);

    if (Index_vec.at(i) == int(distr(generator))) {
      //std::cout << "INDEX " << Index_vec.at(i) << std::endl;
      //std::cout << "rand number " << int(distr(generator)) << std::endl;
      str[Index_vec.at(i)]='U';
      }
		else {
      continue;
		}
  }
}

std::string Qual_random(double *a){  
  char Qualities[] = {'#', '\'', '0', '7' ,'<', 'B', 'F','I','\0'};
  std::string Read_qual;
 
  //std::srand(std::time(nullptr)); //this interferes with the random number for the length.
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::discrete_distribution<> d({a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7]});
  Read_qual += Qualities[d(gen)];
  return Read_qual;
}

/*
double filename(const char* filename){
  std::ifstream infile(filename);
  int row = 150;
  int col = 8;
  double alldata[row][col];
  for(int i = 0; i < row; i++){
    for(int j = 0; j < col; j++){
      infile >> alldata[i][j];
    }
  }
  infile.close();
  return alldata;
}*/

int main(int argc,char **argv){
  
  std::ifstream infile("Freq.txt");
  int row = 150;
  int col = 8;
  double alldata[row][col];

  for(int i = 0; i < row; i++){
    for(int j = 0; j < col; j++){
      infile >> alldata[i][j];
    }
  }
  infile.close();
  std::string Read_qual;

  clock_t tStart = clock();
  const char *fastafile = "chr22.fa";
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *ref = NULL;
  ref  = fai_load(fastafile);
  assert(ref!=NULL);//check that we could load the file

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(ref));
  
  // choosing random sequences using -> random_seq(ref);

  // is lrand48() in order to pick a random sequence if containing more?
  int whichref = lrand48() % faidx_nseq(ref);
  std::cout << "reference number " << whichref << std::endl;
  const char *name = faidx_iseq(ref,whichref);
  int name_len =  faidx_seq_len(ref,name);
  std::cout << "chr name " << name << std::endl;
  std::cout << "size " << name_len << std::endl;
  
  int start_pos = 30000000;
  int end_pos = 30001000; //30001000
  
  std::ofstream outfa("output.fa");
  std::ofstream outfq("output.fq");
  while(start_pos <= end_pos){
    // creates random number in the range of the fragment size rand() % ( high - low + 1 ) + low
    int rand_len = (std::rand() % (80 - 30 + 1)) + 30;
    //std::cout << "random number " << rand_len << std::endl;
        
    char* sequence = faidx_fetch_seq(ref,name,start_pos,start_pos+rand_len,&name_len);
    //std::cout << "SEQUENCE \n" << sequence << std::endl;
    char * pch;
    pch = strchr(sequence,'N');
    if (pch != NULL){
      //Disregards any read with 'N' in.. change this to just change the reading position
      start_pos += rand_len;
      }
    else {
      //std::cout << pch << std::endl;
      char nt[] = "tT";
      Deamin_char(sequence,nt,rand_len);
      // std::cout << sequence ;
      // sequence.size(); //we can use .size if we did the std::string approach which we did with deamin_string calling it damage
      int length = strlen(sequence);
      if (length < 150){
        char adapter = 'X';
        double No = (150-length)/2.0;

        outfa << ">" << name << ":" << start_pos << "-" << start_pos+name_len << "_length:" << length << std::endl;
        outfa << std::string(floor(No),adapter) << sequence << std::string(ceil(No),adapter) << std::endl;

        outfq << "@" << name << ":" << start_pos << "-" << start_pos+name_len << "_length:" << length << std::endl;
        outfq << std::string(floor(No),adapter) << sequence << std::string(ceil(No),adapter) << std::endl;
        outfq << "+" << std::endl;
        for (int row_idx = 0; row_idx < row; row_idx++){Read_qual += Qual_random(alldata[row_idx]);}
        outfq << Read_qual << std::endl;
        Read_qual = "";
        }
        start_pos += rand_len;
      }
    }
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  return 0; 
}