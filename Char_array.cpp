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

void Deamin_char(char* str,char nt[],int seed,
              double alpha=1.0,double beta=2.0,int start=0,int end=25)
{   // use & to pass by reference to the vector
    // Get the first occurrence
  std::vector<int> Index_vec;
  std::srand(seed+std::time(nullptr));
  std::random_device rd;
  std::default_random_engine generator(rd());
  std::gamma_distribution<double> distr(alpha,beta);

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
    
    //std::cout << int(distr(generator)) << std::endl;
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

  //undersøg lige dette
  std::discrete_distribution<> d({a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7]});
  return ASCII_qual[d(gen)];

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


void DNA_complement(char seq[]){
  while (*seq) {
    switch(*seq) {
      case 'A':
        *seq = 'T';
        break;
      case 'a':
        *seq = 't';
        break;
      case 'G':
        *seq = 'C';
        break;
      case 'g':
        *seq = 'c';
        break;
      case 'C':
        *seq = 'G';
        break;
      case 'c':
        *seq = 'g';
        break;
      case 'T':
        *seq = 'A';
        break;
      case 't':
        *seq = 'a';
        break;  
    }
    ++seq;
  }
}

void fafq(const char* fastafile,const char* outname,const char* nt_profile,
          bool flag=false, double cov=1.0){
  
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *ref = NULL;
  ref  = fai_load(fastafile);
  assert(ref!=NULL);//check that we could load the file

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(ref));
  double init = 1.0;
  int chr_no = 0;
  //create the file for saving
  double** my2DArray = create2DArray(600, 8,"/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt");
  std::random_device rd;
  std::default_random_engine gen(rd()); 

  //std::string Qualities = Read_Qual(seqtest,my2DArray,gen);
  //std::cout << Qualities << std::endl;
  //std::cout << "---------------- " << std::endl;

  //std::ofstream outfqr1("output_r1.fq");
  //std::ofstream outfqr2("output_r2.fq");
  char *outname1 =(char*) malloc(10 + strlen(outname) + 1);    
  sprintf(outname1, "%s_r1.fq", outname);
  char *outname2 =(char*) malloc(10 + strlen(outname) + 1);    
  sprintf(outname2, "%s_r2.fq", outname);
  std::ofstream outfqr1(outname1);
  std::ofstream outfqr2(outname2);

  while (chr_no < faidx_nseq(ref)){
    const char *name = faidx_iseq(ref,chr_no);
    int name_len =  faidx_seq_len(ref,name);
    std::cout << "chr name " << name << std::endl;
    std::cout << "size " << name_len << std::endl;
    int start_pos = 1;
    int end_pos = name_len; //30001000

    while(start_pos <= end_pos){
      std::srand(start_pos+std::time(nullptr));
      // Seed random number generator
      int rand_len = (std::rand() % (80 - 30 + 1)) + 30;
      //std::cout << " random " << rand_len << std::endl;
      int dist = init/cov * rand_len; 
      char* sequence = faidx_fetch_seq(ref,name,start_pos,start_pos+rand_len,&name_len);
      // creates random number in the range of the fragment size rand() % ( high - low + 1 ) + low
      
      char * pch;
      pch = strchr(sequence,'N');
      if (pch != NULL){
        //Disregards any read with 'N' in.. change this to just change the reading position
        start_pos += dist + 1;
        }
      else {
        char nt[] = "tT";
        Deamin_char(sequence,nt,rand_len);
        //std::cout << "SEQUENCE \n" << sequence << std::endl;
        // sequence.size(); //we can use .size if we did the std::string approach which we did with deamin_string calling it damage
        int length = strlen(sequence);

        std::string Read_qual;
        if(flag==true){
          char adapter[] = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";

          char Ill_r1[75];
          char *read1 = (char*) malloc(1024);
          strcpy(read1, sequence);
          strcat(read1, adapter);
          strncpy(Ill_r1,read1, sizeof(Ill_r1));
          outfqr1 << "@" << name << ":" << start_pos << "-" << start_pos+name_len << "_length:" << length << std::endl;
          outfqr1 << Ill_r1 << std::endl;
          outfqr1 << "+" << std::endl;
          
          Read_qual = Read_Qual(Ill_r1,my2DArray,gen);
          outfqr1 << Read_qual << std::endl;
          Read_qual = "";

          char Ill_r2[75];
          char *read2 = (char*) malloc(1024);
          strcpy(read2, sequence);
          DNA_complement(read2);
          std::reverse(read2, read2 + strlen(read2));
          strcat(read2, adapter);
          strncpy(Ill_r2,read2, sizeof(Ill_r2));
          outfqr2 << "@" << name << ":" << start_pos << "-" << start_pos+name_len << "_length:" << length << std::endl;
          outfqr2 << Ill_r2 << std::endl;
          outfqr2 << "+" << std::endl;
          Read_qual = Read_Qual(Ill_r2,my2DArray,gen);
          outfqr2 << Read_qual << std::endl;
          Read_qual = "";
        }
        else{
          int length = strlen(sequence);
          char Ill_r1[1024];
          char *read1 = (char*) malloc(1024);
          strcpy(read1, sequence);
          strncpy(Ill_r1,sequence,sizeof(Ill_r1));
          //std::cout << "sequence " << Ill_r1 << std::endl;
          //std::cout << "lenght " << length << std::endl;
          outfqr1 << "@" << name << ":" << start_pos << "-" << start_pos+name_len << "_length:" << length << std::endl;
          outfqr1 << Ill_r1 << std::endl;
          outfqr1 << "+" << std::endl;
          //for (int row_idx = 0; row_idx < std::min(length,75); row_idx++){Read_qual += Qual_random(my2DArray[row_idx]);}
          //std::cout << "sequence2 " << Ill_r1 << std::endl;
          
          Read_qual = Read_Qual(Ill_r1,my2DArray,gen);

          outfqr1 << Read_qual << std::endl;
          Read_qual = "";

          char Ill_r2[75];
          char *read2 = (char*) malloc(1024);
          strcpy(read2, sequence);
          DNA_complement(read2);
          std::reverse(read2, read2 + strlen(read2));
          strncpy(Ill_r2,read2, sizeof(Ill_r2));
          outfqr2 << "@" << name << ":" << start_pos << "-" << start_pos+name_len << "_length:" << length << std::endl;
          outfqr2 << Ill_r2 << std::endl;
          outfqr2 << "+" << std::endl;
          //for (int row_idx = 150-std::min(length,75); row_idx < 150; row_idx++){Read_qual += Qual_random(my2DArray[row_idx]);}
          Read_qual = Read_Qual(Ill_r2,my2DArray,gen);
          outfqr2 << Read_qual << std::endl;
          Read_qual = "";
        }
      }
    start_pos += dist + 1;
    }
  chr_no++;
  }
  return; 
}
// ----------------------------
/*
int main(int argc,char **argv){
  clock_t tStart = clock();
  const char *fastafile = "/home/wql443/scratch/reference_genome/hg19/chr22.fa";
  const char *Profile = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt";
  fafq(fastafile,"test2",Profile,false);
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  return 0;
}
*/

void Distfunc2(double** Array2d,std::discrete_distribution<> dist[],int size){
  for (int row_idx = 0; row_idx < size; row_idx++){
    std::discrete_distribution<> d({Array2d[row_idx][0],Array2d[row_idx][1],Array2d[row_idx][2],Array2d[row_idx][3],
                                    Array2d[row_idx][4],Array2d[row_idx][5],Array2d[row_idx][6],Array2d[row_idx][7]});  
    dist[row_idx] = d;
  }
}

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


void fafq2(const char* fastafile,const char* outname,const char* nt_profile,
          bool flag=false, double cov=1.0){
  
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *ref = NULL;
  ref  = fai_load(fastafile);
  assert(ref!=NULL);//check that we could load the file

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(ref));
  double init = 1.0;
  int chr_no = 0;
  //create the file for saving
  double** my2DArray = create2DArray(600, 8,"/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt");
  std::random_device rd;
  std::default_random_engine gen(rd()); 

  int Qualdist_size = 600;
  std::discrete_distribution<> Qualdist[Qualdist_size];
  Distfunc2(my2DArray,Qualdist,Qualdist_size);

  //std::string Qualities = Read_Qual(seqtest,my2DArray,gen);
  //std::cout << Qualities << std::endl;
  //std::cout << "---------------- " << std::endl;

  //std::ofstream outfqr1("output_r1.fq");
  //std::ofstream outfqr2("output_r2.fq");
  char *outname1 =(char*) malloc(10 + strlen(outname) + 1);    
  sprintf(outname1, "%s_r1.fq", outname);
  char *outname2 =(char*) malloc(10 + strlen(outname) + 1);    
  sprintf(outname2, "%s_r2.fq", outname);
  std::ofstream outfqr1(outname1);
  std::ofstream outfqr2(outname2);

  while (chr_no < faidx_nseq(ref)){
    const char *name = faidx_iseq(ref,chr_no);
    int name_len =  faidx_seq_len(ref,name);
    std::cout << "chr name " << name << std::endl;
    std::cout << "size " << name_len << std::endl;
    int start_pos = 1;
    int end_pos = name_len; //30001000

    while(start_pos <= end_pos){
      std::srand(start_pos+std::time(nullptr));
      // Seed random number generator
      int rand_len = (std::rand() % (80 - 30 + 1)) + 30;
      //std::cout << " random " << rand_len << std::endl;
      int dist = init/cov * rand_len; 
      char* sequence = faidx_fetch_seq(ref,name,start_pos,start_pos+rand_len,&name_len);
      // creates random number in the range of the fragment size rand() % ( high - low + 1 ) + low
      
      char * pch;
      pch = strchr(sequence,'N');
      if (pch != NULL){
        //Disregards any read with 'N' in.. change this to just change the reading position
        start_pos += dist + 1;
        }
      else {
        char nt[] = "tT";
        Deamin_char(sequence,nt,rand_len);
        //std::cout << "SEQUENCE \n" << sequence << std::endl;
        // sequence.size(); //we can use .size if we did the std::string approach which we did with deamin_string calling it damage
        int length = strlen(sequence);

        std::string Read_qual;
        if(flag==true){
          char adapter[] = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";

          char Ill_r1[75];
          char *read1 = (char*) malloc(1024);
          strcpy(read1, sequence);
          strcat(read1, adapter);
          strncpy(Ill_r1,read1, sizeof(Ill_r1));
          outfqr1 << "@" << name << ":" << start_pos << "-" << start_pos+name_len << "_length:" << length << std::endl;
          outfqr1 << Ill_r1 << std::endl;
          outfqr1 << "+" << std::endl;
          
          Read_qual = Read_Qual3(Ill_r1,Qualdist,gen);
          outfqr1 << Read_qual << std::endl;
          Read_qual = "";

          char Ill_r2[75];
          char *read2 = (char*) malloc(1024);
          strcpy(read2, sequence);
          DNA_complement(read2);
          std::reverse(read2, read2 + strlen(read2));
          strcat(read2, adapter);
          strncpy(Ill_r2,read2, sizeof(Ill_r2));
          outfqr2 << "@" << name << ":" << start_pos << "-" << start_pos+name_len << "_length:" << length << std::endl;
          outfqr2 << Ill_r2 << std::endl;
          outfqr2 << "+" << std::endl;
          Read_qual = Read_Qual3(Ill_r2,Qualdist,gen);
          outfqr2 << Read_qual << std::endl;
          Read_qual = "";
        }
        else{
          int length = strlen(sequence);
          char Ill_r1[1024];
          char *read1 = (char*) malloc(1024);
          strcpy(read1, sequence);
          strncpy(Ill_r1,sequence,sizeof(Ill_r1));
          //std::cout << "sequence " << Ill_r1 << std::endl;
          //std::cout << "lenght " << length << std::endl;
          outfqr1 << "@" << name << ":" << start_pos << "-" << start_pos+name_len << "_length:" << length << std::endl;
          outfqr1 << Ill_r1 << std::endl;
          outfqr1 << "+" << std::endl;
          //for (int row_idx = 0; row_idx < std::min(length,75); row_idx++){Read_qual += Qual_random(my2DArray[row_idx]);}
          //std::cout << "sequence2 " << Ill_r1 << std::endl;
          
          Read_qual = Read_Qual3(Ill_r1,Qualdist,gen);

          outfqr1 << Read_qual << std::endl;
          Read_qual = "";

          char Ill_r2[75];
          char *read2 = (char*) malloc(1024);
          strcpy(read2, sequence);
          DNA_complement(read2);
          std::reverse(read2, read2 + strlen(read2));
          strncpy(Ill_r2,read2, sizeof(Ill_r2));
          outfqr2 << "@" << name << ":" << start_pos << "-" << start_pos+name_len << "_length:" << length << std::endl;
          outfqr2 << Ill_r2 << std::endl;
          outfqr2 << "+" << std::endl;
          //for (int row_idx = 150-std::min(length,75); row_idx < 150; row_idx++){Read_qual += Qual_random(my2DArray[row_idx]);}
          Read_qual = Read_Qual3(Ill_r2,Qualdist,gen);
          outfqr2 << Read_qual << std::endl;
          Read_qual = "";
        }
      }
    start_pos += dist + 1;
    }
  chr_no++;
  }
  return; 
}

void fafq3(const char* fastafile,FILE *fp1,FILE *fp2,const char* nt_profile,
          bool flag=false, double cov=1.0){
  
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);//check that we could load the file

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(seq_ref));
  double init = 1.0;
  int chr_no = 0;
  //create the file for saving
  double** my2DArray = create2DArray(600, 8,"/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt");
  std::random_device rd;
  std::default_random_engine gen(rd()); 

  int Qualdist_size = 600;
  std::discrete_distribution<> Qualdist[Qualdist_size];
  Distfunc2(my2DArray,Qualdist,Qualdist_size);

  while (chr_no < faidx_nseq(seq_ref)){
    const char *name = faidx_iseq(seq_ref,chr_no);
    int name_len =  faidx_seq_len(seq_ref,name);
    fprintf(stderr,"-> name: \'%s\' name_len: %d\n",name,name_len);
    char *data = fai_fetch(seq_ref,name,&name_len);

    int start_pos = 1;
    int end_pos = name_len; //30001000
    char seqmod[1024];

    while(start_pos <= end_pos){
      std::srand(start_pos+std::time(nullptr));
      // Seed random number generator
      int readlength = (std::rand() % (80 - 30 + 1)) + 30;
      //int readlength = drand48()*(80.0-30.0)+30.0;
      int stop = start_pos+(int) readlength;
      //int dist = init/cov * rand_len; 

      //extracts the sequence
      strncpy(seqmod,data+start_pos,readlength);

      char * pch; 
      pch = strchr(seqmod,'N');
      if (pch != NULL){start_pos += readlength + 1;}
      else {
        char nt[] = "tT";
        Deamin_char(seqmod,nt,readlength);
        std::string Read_qual;
        if(flag==true){
          char adapter[] = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
          char Ill_r1[75];
          char *read1 = (char*) malloc(1024);
          strcpy(read1, seqmod);
          strcat(read1, adapter);
          strncpy(Ill_r1,read1, sizeof(Ill_r1));

          fprintf(fp1,"@%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
          fprintf(fp1,"%s\n",seqmod);
          fprintf(fp1,"%s\n","+");
          Read_qual = Read_Qual3(Ill_r1,Qualdist,gen);
          fprintf(fp1,"%s\n",Read_qual.c_str());
          Read_qual = "";

          char Ill_r2[75];
          char *read2 = (char*) malloc(1024);
          strcpy(read2, seqmod);
          DNA_complement(read2);
          std::reverse(read2, read2 + strlen(read2));
          strcat(read2, adapter);
          strncpy(Ill_r2,read2, sizeof(Ill_r2));

          fprintf(fp2,"@%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
          fprintf(fp2,"%s\n",seqmod);
          fprintf(fp2,"%s\n","+");
          Read_qual = Read_Qual3(Ill_r2,Qualdist,gen);
          fprintf(fp2,"%s\n",Read_qual.c_str());
          Read_qual = "";
        }
        else{
          char Ill_r1[75];
          char *read1 = (char*) malloc(1024);
          strcpy(read1, seqmod);
          strncpy(Ill_r1,read1, sizeof(Ill_r1));

          fprintf(fp1,"@%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
          fprintf(fp1,"%s\n",seqmod);
          fprintf(fp1,"%s\n","+");
          Read_qual = Read_Qual3(Ill_r1,Qualdist,gen);
          fprintf(fp1,"%s\n",Read_qual.c_str());
          Read_qual = "";

          char Ill_r2[75];
          char *read2 = (char*) malloc(1024);
          strcpy(read2, seqmod);
          DNA_complement(read2);
          std::reverse(read2, read2 + strlen(read2));
          strncpy(Ill_r2,read2, sizeof(Ill_r2));

          fprintf(fp2,"@%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
          fprintf(fp2,"%s\n",seqmod);
          fprintf(fp2,"%s\n","+");
          Read_qual = Read_Qual3(Ill_r2,Qualdist,gen);
          fprintf(fp2,"%s\n",Read_qual.c_str());
          Read_qual = "";

        }
      }
    start_pos += readlength + 1;
    readlength = 0;
    memset(seqmod, 0, sizeof seqmod);
    }
  chr_no++;
  }
  return;
}


void fafq4(const char* fastafile,FILE *fp1,FILE *fp2,const char* r1_profile,const char* r2_profile,
          bool flag=false, double cov=1.0){
  
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);//check that we could load the file

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(seq_ref));
  double init = 1.0;
  int chr_no = 0;
  //create the file for saving
  double** R1_2Darray = create2DArray(600, 8,r1_profile);
  double** R2_2Darray = create2DArray(600, 8,r2_profile);

  std::random_device rd;
  std::default_random_engine gen(rd()); 

  int Qualdist_size = 600;
  std::discrete_distribution<> Qualdistr1[Qualdist_size];
  Distfunc2(R1_2Darray,Qualdistr1,Qualdist_size);
  std::discrete_distribution<> Qualdistr2[Qualdist_size];
  Distfunc2(R2_2Darray,Qualdistr2,Qualdist_size);

  while (chr_no < faidx_nseq(seq_ref)){
    const char *name = faidx_iseq(seq_ref,chr_no);
    int name_len =  faidx_seq_len(seq_ref,name);
    fprintf(stderr,"-> name: \'%s\' name_len: %d\n",name,name_len);
    char *data = fai_fetch(seq_ref,name,&name_len);

    int start_pos = 1;
    int end_pos = name_len; //30001000
    char seqmod[1024];

    while(start_pos <= end_pos){
      std::srand(start_pos+std::time(nullptr));
      // Seed random number generator
      int readlength = (std::rand() % (80 - 30 + 1)) + 30;
      //int readlength = drand48()*(80.0-30.0)+30.0;
      int stop = start_pos+(int) readlength;
      //int dist = init/cov * rand_len; 

      //extracts the sequence
      strncpy(seqmod,data+start_pos,readlength);

      char * pch; 
      pch = strchr(seqmod,'N');
      if (pch != NULL){start_pos += readlength + 1;}
      else {
        char nt[] = "tT";
        Deamin_char(seqmod,nt,readlength);
        std::string Read_qual;
        if(flag==true){
          char adapter[] = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
          char Ill_r1[150];
          char *read1 = (char*) malloc(1024);
          strcpy(read1, seqmod);
          strncpy(Ill_r1,read1, sizeof(Ill_r1));

          fprintf(fp1,"@%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
          fprintf(fp1,"%s\n",seqmod);
          fprintf(fp1,"%s\n","+");
          Read_qual = Read_Qual3(Ill_r1,Qualdistr1,gen);
          fprintf(fp1,"%s\n",Read_qual.c_str());
          Read_qual = "";

          char Ill_r2[150];
          char *read2 = (char*) malloc(1024);
          strcpy(read2, seqmod);
          DNA_complement(read2);
          std::reverse(read2, read2 + strlen(read2));
          strcat(read2, adapter);
          strncpy(Ill_r2,read2, sizeof(Ill_r2));

          fprintf(fp2,"@%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
          fprintf(fp2,"%s\n",seqmod);
          fprintf(fp2,"%s\n","+");
          Read_qual = Read_Qual3(Ill_r2,Qualdistr2,gen);
          fprintf(fp2,"%s\n",Read_qual.c_str());
          Read_qual = "";
        }
        else{
          char Ill_r1[150];
          strncpy(Ill_r1,seqmod, sizeof(Ill_r1));

          fprintf(fp1,"@%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
          fprintf(fp1,"%s\n",Ill_r1);
          fprintf(fp1,"%s\n","+");
          Read_qual = Read_Qual3(Ill_r1,Qualdistr1,gen);
          fprintf(fp1,"%s\n",Read_qual.c_str());
          Read_qual = "";

          char Ill_r2[150];
          DNA_complement(seqmod);
          std::reverse(seqmod, seqmod + strlen(seqmod));
          strncpy(Ill_r2,seqmod, sizeof(Ill_r2));

          fprintf(fp2,"@%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
          fprintf(fp2,"%s\n",seqmod);
          fprintf(fp2,"%s\n","+");
          Read_qual = Read_Qual3(seqmod,Qualdistr2,gen);
          fprintf(fp2,"%s\n",Read_qual.c_str());
          Read_qual = "";

        }
      }
    start_pos += readlength + 1;
    readlength = 0;
    memset(seqmod, 0, sizeof seqmod);
    }
  chr_no++;
  }
  return;
}

int main(int argc,char **argv){
  clock_t tStart = clock();
  const char *fastafile = "/home/wql443/scratch/reference_genome/hg19/chr22.fa";
  const char *Profile1 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt";
  const char *Profile2 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R2.txt";

  // fafq(fastafile,"test2",Profile,false);
  // fafq2(fastafile,"test2",Profile,false);
  
  FILE *fp1;
  FILE *fp2;
  fp1 = fopen("R1.fq","wb");
  fp2 = fopen("R2.fq","wb");
  fafq4(fastafile,fp1,fp2,Profile1,Profile2); 
  //fafq3(fastafile,fp1,fp2,Profile1);
  fclose(fp1);
  fclose(fp2);
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  return 0;
}





//const char* fastafile,FILE *fp1,FILE *fp2,const char* nt_profile

/*
int main(){
  char array1[] ="AAA";
  char array2[] = "XXXXXXXCCCCCCCC";
  char Ill_r1[150];
  std::strcat(array1,array2);
  std::strncpy(Ill_r1, array1, 14);
  std::cout << array1 << std::endl;
  std::cout << Ill_r1 << std::endl;
  return 0;
}*/