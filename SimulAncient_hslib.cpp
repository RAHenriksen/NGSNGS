#include <cstdio>
#include <cassert>
#include <htslib/faidx.h>
#include <htslib/sam.h>
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

int Qual_random(double *a){
  //creates a random sequence of nt qualities based on a frequency distribution

  char Qualities[] = {'#', '\'', '0', '7' ,'<', 'B', 'F','I','\0'};
  int ASCII_qual[] = {35,39,48,55,60,66,70,73}; //Im choosing the ascii values, since im casting them into string later on

  //std::string Read_qual;
 
  //std::srand(std::time(nullptr)); //this interferes with the random number for the length.
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::discrete_distribution<> d({a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7]});
  return ASCII_qual[d(gen)];
  //Read_qual += ASCII_qual[d(gen)];
  //return Read_qual;
}

void filename(const char* filename){
  //creates a nt qual string immediately from the frequency file
  std::ifstream infile(filename);
  int row = 150;
  int col = 8;
  double alldata[row][col];
  std::string Read_qual;
  for(int i = 0; i < row; i++){
    for(int j = 0; j < col; j++){
      infile >> alldata[i][j];
    }
  }
  for (int row = 0; row < 150; row++){ 
     Read_qual += Qual_random(alldata[row]);
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

/*
int main(int argc,char **argv){
  std::srand(std::time(nullptr));
  double** my2DArray = create2DArray(150, 8,"Freq.txt");
  clock_t tStart = clock();
  const char *fastafile = "chr22.fa";    
  faidx_t *ref = NULL;
  ref  = fai_load(fastafile);
  assert(ref!=NULL);//check that we could load the file
  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(ref));

  //Creates a pointer to allocated memomry for the format??
  htsFormat *fmt_hts =(htsFormat*) calloc(1,sizeof(htsFormat));
  //wb -> bam , wc -> cram
  char out_mode[5]="wb";
  const char *outfile_nam = "test.bam";
  samFile *outfile = NULL;
  // creates a pointer to generated header
  sam_hdr_t *header = sam_hdr_init();
  if (header == NULL) { fprintf(stderr, "sam_hdr_init"); return 0;}
  char *refName=NULL;

  if ((outfile = sam_open_format(outfile_nam, out_mode, fmt_hts)) == 0) {
    fprintf(stderr,"Error opening file for writing\n");
    exit(0);
  }

  // Creating sam_header
  char *name_len_char =(char*) malloc(1024);
  for(int i=0;i<faidx_nseq(ref);i++){
      const char *name = faidx_iseq(ref,i);
      int name_len =  faidx_seq_len(ref,name);
    // skal være c string i array så vi konvertere int om til char
    snprintf(name_len_char,1024,"%d",name_len);
    //    fprintf(stderr,"ref:%d %d %s\n",i,name,name_len_char);

    // reference part of the header, int r variable ensures the header is added
    int r = sam_hdr_add_line(header, "SQ", "SN", name, "LN", name_len_char, NULL);
    if (r < 0) { fprintf(stderr,"sam_hdr_add_line"); return 0; }
  }
  
  // saving the header to the file
  if (sam_hdr_write(outfile, header) < 0) fprintf(stderr,"writing headers to %s", outfile);
    
  // alignment delen, gemt i bam_1 type for at representere hver linje som 1 alignment
  bam1_t *b = bam_init1(); //initialisere bam1_t (type) til hukommelsen!

  char buf[96];
  int whichref = lrand48() % faidx_nseq(ref);
  const char *name = faidx_iseq(ref,whichref);
  int name_len =  faidx_seq_len(ref,name);
  
  int start_pos = 1;
  int end_pos = name_len; //30001000
  double cov = 1.0;
  double init = 1.0;

  while(start_pos <= end_pos){
    // creates random number in the range of the fragment size rand() % ( high - low + 1 ) + low
    int rand_len = (std::rand() % (80 - 30 + 1)) + 30;
    int dist = init/cov * rand_len;
    //we use structure faidx_t from htslib to load in a fasta
    int end_tmp = start_pos+rand_len;
    char* sequence = faidx_fetch_seq(ref,name,start_pos,end_tmp,&name_len);
    snprintf(buf,96,"%s:%d-%d",name,start_pos,end_tmp);
    
    char * pch;
    pch = strchr(sequence,'N');
    if (pch != NULL){
      //Disregards any read with 'N' in.. change this to just change the reading position
      start_pos += dist + 1;
      }
    else {
      std::string Read_qual;
      //std::cout << pch << std::endl;
      char nt[] = "tT";
      Deamin_char(sequence,nt,rand_len);
      //ensures proper bam format
      int unmap_flag = 4;
      int RNAME_chr_ID = 0; // using 0 takes info sam_hdr_t
      int mapQ = 255; // 255 if unavailable
      int no_cigar = 0; //number of cigar operations 
      const uint32_t *cigar = NULL; // cigar data, NULL if unavailable - But how do we then add actual cigar info?

      //casting the ascii values to strings, but im substracting -33 from the values, since bam used phred+33.
      for (int row_idx = 0; row_idx < strlen(sequence); row_idx++){Read_qual += Qual_random(my2DArray[row_idx])-33;}
      //std::cout << Read_qual << std::endl;
      const char* nt_qual = Read_qual.c_str();
      //std::cout << nt_qual << std::endl;
      //std::cout << (int) nt_qual[1] << (char) ((int) nt_qual[1] + 33) << std::endl;
      bam_set1(b,strlen(buf),buf,unmap_flag,RNAME_chr_ID,start_pos-1,mapQ,no_cigar,cigar,0,0,0,end_tmp-start_pos,sequence,nt_qual,0);
      sam_write1(outfile,header,b);
      Read_qual = "";

      start_pos += dist + 1;
    }
  }
  sam_hdr_destroy(header);
  sam_close(outfile);
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  return 0;
}*/

/*
int main(int argc,char **argv){

  double** my2DArray = create2DArray(150, 8,"Freq.txt");

  clock_t tStart = clock();
  const char *fastafile = "/home/wql443/scratch/reference_genome/hg19/chr22.fa";
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
  
  int start_pos = 1;
  int end_pos = name_len; //30001000
  double cov = 1.0;
  double init = 1.0;
  
  std::ofstream outfa("output.fa");
  std::ofstream outfq("output.fq");

  while(start_pos <= end_pos){
    // creates random number in the range of the fragment size rand() % ( high - low + 1 ) + low
    int rand_len = (std::rand() % (80 - 30 + 1)) + 30;
    int dist = init/cov * rand_len;
    
    char* sequence = faidx_fetch_seq(ref,name,start_pos,start_pos+rand_len,&name_len);
    //std::cout << "SEQUENCE \n" << sequence << std::endl;
    char * pch;
    pch = strchr(sequence,'N');
    if (pch != NULL){
      //Disregards any read with 'N' in.. change this to just change the reading position
      start_pos += dist + 1;
      }
    else {
      std::string Read_qual;
      //std::cout << pch << std::endl;
      char nt[] = "tT";
      Deamin_char(sequence,nt,rand_len);
      // std::cout << sequence ;
      // sequence.size(); //we can use .size if we did the std::string approach which we did with deamin_string calling it damage
      int length = strlen(sequence);
      if (length < 150){
        char adapter = 'X';
        double No = (150-length)/2.0;
        if (std::strcmp(argv[1], "fa") == 0){
          outfa << ">" << name << ":" << start_pos << "-" << start_pos+name_len << "_length:" << length << std::endl;
          outfa << std::string(floor(No),adapter) << sequence << std::string(ceil(No),adapter) << std::endl;
        }
        else if (std::strcmp(argv[1], "fq") == 0){
          outfq << "@" << name << ":" << start_pos << "-" << start_pos+name_len << "_length:" << length << std::endl;
          outfq << std::string(floor(No),adapter) << sequence << std::string(ceil(No),adapter) << std::endl;
          outfq << "+" << std::endl;
          for (int row_idx = 0; row_idx < 150; row_idx++){Read_qual += Qual_random(my2DArray[row_idx]);}
          outfq << Read_qual << std::endl;
          Read_qual = "";
        }
      }
      start_pos += dist + 1;
      //start_pos += rand_len;
    }
  }
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  return 0; 
}*/


int main(int argc,char **argv){

  double** my2DArray = create2DArray(150, 8,"Freq.txt");

  clock_t tStart = clock();
  const char *fastafile = "/home/wql443/scratch/reference_genome/hg19/chr22.fa";
  //const char *fastafile = "/home/wql443/scratch/reference_genome/hg19/chr22.fa";
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *ref = NULL;
  ref  = fai_load(fastafile);
  assert(ref!=NULL);//check that we could load the file

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(ref));
  std::cout << faidx_iseq(ref,1) << " lol" << std::endl;
  double cov = 1.0;
  double init = 1.0;
  
  int chr_no = 0;
  //if (std::strcmp(argv[1], "fa") == 0){std::ofstream outfa("output.fa");}
  //else if (std::strcmp(argv[1], "fq") == 0){std::ofstream outfq("output.fq");}
  std::ofstream outfa("output.fa");
  std::ofstream outfq("output.fq");
  std::ofstream outfqr1("output_r1.fq");
  std::ofstream outfqr2("output_r2.fq");

  std::srand(std::time(nullptr));
  while (chr_no < faidx_nseq(ref)){
    std::cout << "reference number " << chr_no << std::endl;
    const char *name = faidx_iseq(ref,chr_no);
    int name_len =  faidx_seq_len(ref,name);
    std::cout << "chr name " << name << std::endl;
    std::cout << "size " << name_len << std::endl;
    int start_pos = 1;
    int end_pos = name_len; //30001000

    while(start_pos <= end_pos){
      // Seed random number generator
      int rand_len = (std::rand() % (80 - 30 + 1)) + 30;
      //std::cout << " random " << rand_len << std::endl;
      int dist = init/cov * rand_len; 
      char* sequence = faidx_fetch_seq(ref,name,start_pos,start_pos+rand_len,&name_len);
      // creates random number in the range of the fragment size rand() % ( high - low + 1 ) + low
      
      //std::cout << "SEQUENCE \n" << sequence << std::endl;
      char * pch;
      pch = strchr(sequence,'N');
      if (pch != NULL){
        //Disregards any read with 'N' in.. change this to just change the reading position
        start_pos += dist + 1;
        }
      else {
        char nt[] = "tT";
        Deamin_char(sequence,nt,rand_len);
        // sequence.size(); //we can use .size if we did the std::string approach which we did with deamin_string calling it damage
        int length = strlen(sequence);

        if (std::strcmp(argv[1], "fa") == 0){
          outfa << ">" << name << ":" << start_pos << "-" << start_pos+name_len << "_length:" << length << std::endl;
          outfa << sequence << std::endl;
          }
        else if (std::strcmp(argv[1], "fq") == 0){
          std::string Read_qual;
          char adapter[] = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
          /*
          outfq << "@" << name << ":" << start_pos << "-" << start_pos+name_len << "_length:" << length << std::endl;
          outfq << sequence << std::endl;
          outfq << "+" << std::endl;
          for (int row_idx = 0; row_idx < 150; row_idx++){Read_qual += Qual_random(my2DArray[row_idx]);}
          outfq << Read_qual << std::endl;
          Read_qual = "";*/

          char Ill_r1[75];
          char *read1 = (char*) malloc(1024);
          strcpy(read1, sequence);
          strcat(read1, adapter);
          strncpy(Ill_r1,read1, sizeof(Ill_r1));
          outfqr1 << "@" << name << ":" << start_pos << "-" << start_pos+name_len << "_length:" << length << std::endl;
          outfqr1 << Ill_r1 << std::endl;
          outfqr1 << "+" << std::endl;
          for (int row_idx = 0; row_idx < 75; row_idx++){Read_qual += Qual_random(my2DArray[row_idx]);}
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
          for (int row_idx = 75; row_idx < 150; row_idx++){Read_qual += Qual_random(my2DArray[row_idx]);}
          outfqr2 << Read_qual << std::endl;
          Read_qual = "";
        }
      }
    start_pos += dist + 1;
    std::srand(start_pos);
    //std::cout << "start " << start_pos << std::endl;
    }
  chr_no++;
  std::cout << "---------------------" << std::endl;
  }
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  return 0; 
}