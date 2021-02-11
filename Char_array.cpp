#include <cstdio>
#include <stdlib.h>
#include <stdio.h>

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

void fafa(const char* fastafile, const char* outflag){
  //creates fasta structure
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  int whichref = lrand48() % faidx_nseq(seq_ref);
  const char *name = faidx_iseq(seq_ref,whichref);
  int name_len =  faidx_seq_len(seq_ref,name);
  fprintf(stderr,"-> name: \'%s\' name_len: %d\n",name,name_len);
  
  //creates randomized fragment length based on drand seed in main
  int readlength = drand48()*(80.0-30.0)+30.0;
  int start = lrand48() % name_len;
  int stop = start+(int) readlength;

  //read entire chromosome
  char *data = fai_fetch(seq_ref,name,&name_len);
  
  //creates sequence of random length for DNA damage
  char seqmod[1024];

  //adress of entire array, first element always the same but adding start for query source start, and stop-start : number of elements
  strncpy(seqmod,data+start,stop-start);

  char nt[] = "tT";
  printf(seqmod);
  std::cout << "----------" << std::endl;
  Deamin_char(seqmod,nt,readlength);
  printf(seqmod);

  //void pointer, casting to int and with size of chromosome
  int* cov2 = (int*) calloc(name_len,sizeof(int)); //some reason it didnt work with sizeof(int) * name_len
  for(int i = start; i <= stop;i++){cov2[i] += 1;}

  if (std::strcmp(outflag, "gz") == 0){
    kstring_t kstr;// kstring_t *kstr = new kstr;
    kstr.s = NULL; // kstr->s = NULL;
    kstr.l = kstr.m = 0; //kstr->l = kstr->m = 0;

    ksprintf(&kstr,">%s:%d-%d_length:%d\n",name,start,stop,readlength);
    ksprintf(&kstr,"%s\n",seqmod);

    gzFile gz=gzopen("test.fa.gz","wb");
    gzwrite(gz,kstr.s,kstr.l);kstr.l =0;
    gzclose(gz);
  }
  else
  {
    FILE *fp = fopen("test.fa","wb");
    fprintf(fp,">%s:%d-%d_length:%d\n",name,start,stop,readlength);
    fprintf(fp,"%s\n",seqmod);
    fclose(fp);
  }
}

/*
int main(int argc,char **argv){
  clock_t tStart = clock();
  const char *fastafile = "/home/wql443/scratch/reference_genome/hg19/chr22.fa";
  int seed = -1;
  if(argc>3)
    seed = atoi(argv[3]);
  else
    seed = time(NULL);

  fprintf(stderr,"seed is: %d\n",seed);
  srand48(seed);

  if (std::strcmp(argv[1], "fafa") == 0){
    if (std::strcmp(argv[2], "gz") == 0){
      fafa(fastafile,argv[2]); 
    }
    else if (std::strcmp(argv[2], "fa") == 0){
      fafa(fastafile,argv[2]); 
    }
  }
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
}
*/

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

  //std::string Read_qual;
 
  //std::srand(std::time(nullptr)); //this interferes with the random number for the length.

  //flyt udenfor funktionen, så du ikke lavert et nyt random objekt hver gang.
  std::random_device rd;
  std::default_random_engine gen(rd()); 

  //undersøg lige dette
  std::discrete_distribution<> d({a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7]});
  return ASCII_qual[d(gen)];
  //Read_qual += ASCII_qual[d(gen)];
  //return Read_qual;
}

int main(){
  //std::ifstream file("Freq.txt");
  double** my2DArray = create2DArray(600, 8,"Freq_R1.txt");
  // rækkefølgen er 0-150 "A" -> 150 - 300 "T" ->  300-450 "G" -> 450 - 600 "C"
  std::string Read_qual;
  
  char seqtest[1024] = "AGCTaGCTAGTCAGCTaGCTAGTCAGCTaGCTAGTCAGCTaGCTAGTCAGCTaGCTAGTC";

  int Tstart = 150;
  int Gstart = 300;
  int Cstart = 450;

  int read_length = strlen(seqtest);
  for (int row_idx = 0; row_idx < read_length; row_idx++){
    switch(seqtest[row_idx]){
      // iterates through every element in seq and creates and extract qual from the given line 
        // of the sequence qual profile. based on the lines number given the Tstart etc..
      case 'A':
        Read_qual += Qual_random(my2DArray[row_idx]);
        break;
      case 'a':
        Read_qual +=  Qual_random(my2DArray[row_idx]);
        break;
      case 'T':
        Read_qual +=  Qual_random(my2DArray[row_idx + Tstart]);
        break;
      case 't':
        Read_qual += Qual_random(my2DArray[row_idx + Tstart]);
        break;  
      case 'G':
        Read_qual +=  Qual_random(my2DArray[row_idx + Gstart]);
        break;
      case 'g':
        Read_qual +=  Qual_random(my2DArray[row_idx + Gstart]);
        break;
      case 'C':
        Read_qual +=  Qual_random(my2DArray[row_idx + Cstart]);
        break;
      case 'c':
        Read_qual +=  Qual_random(my2DArray[row_idx + Cstart]);
        break;
    }
  }
  std::cout << Read_qual << std::endl;  
  return 0;
}
