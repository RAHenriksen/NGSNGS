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


void fafq(const char* fastafile,const char* outflag,FILE *fp1,FILE *fp2,gzFile gz1,gzFile gz2,bool flag=false,
          const char* r1_profile = "Qual_profiles/Freq_R1.txt",
          const char* r2_profile = "Qual_profiles/Freq_R1.txt",
          char Adapter1[]=NULL,char Adapter2[]=NULL){
  
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);//check that we could load the file

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(seq_ref));
  double init = 1.0;
  int chr_no = 0;

  // more generalized way identify number of lines for read profile
  std::ifstream file(r1_profile);
  int Line_no = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
  file.close();

  // loading in error profiles for nt substitions and read qual
  double** Error_2darray = create2DArray("/home/wql443/WP1/Ill_err.txt",4,280);
  double** R1_2Darray = create2DArray(r1_profile,8,Line_no);
  double** R2_2Darray = create2DArray(r2_profile,8,Line_no);

  // creating random objects for all distributions.
  std::random_device rd;
  std::default_random_engine gen(rd()); 

  std::discrete_distribution<> Qualdistr1[Line_no];
  Qual_dist(R1_2Darray,Qualdistr1,Line_no);
  std::discrete_distribution<> Qualdistr2[Line_no];
  Qual_dist(R2_2Darray,Qualdistr2,Line_no);
  std::discrete_distribution<> Error[280];
  Seq_err(Error_2darray,Error,280);

  while (chr_no < faidx_nseq(seq_ref)){
    const char *name = faidx_iseq(seq_ref,chr_no);
    int name_len =  faidx_seq_len(seq_ref,name);
    fprintf(stderr,"-> name: \'%s\' name_len: %d\n",name,name_len);
    char *data = fai_fetch(seq_ref,name,&name_len);

    int start_pos = 1;
    int end_pos = name_len; //30001000
    char seqmod[1024];

    while(start_pos <= end_pos){
      //std::srand(start_pos+std::time(nullptr));
      //Seed random number generator 
      //int readlength = (std::rand() % (60 - 30 + 1)) + 30;
      int readlength = drand48()*(70.0-30.0)+30.0;
      int stop = start_pos+(int) readlength;
      //int dist = init/cov * rand_len; 

      //extracts the sequence
      strncpy(seqmod,data+start_pos,readlength);

      char * pch; 
      pch = strchr(seqmod,'N');
      if (pch != NULL){start_pos += readlength + 1;}
      else {
        char nt[] = "tT";
        char qual[1024] = "";

        //Deamination function is no longer needed due to the error profile
        Ill_err(seqmod,Error,gen);
        
        if(flag==true){
          size_t lib_size = Line_no/4; 
          char read1N[lib_size+1];
          char read2N[lib_size+1];

          memset(read1N,'N',lib_size);
          read1N[lib_size]='\0';
          memset(read2N,'N',lib_size);
          read2N[lib_size]='\0';

          char read1[lib_size + 1];
          char read2[lib_size + 1];
          //Copies sequence into both reads
          strcpy(read1, seqmod);
          strcpy(read2, seqmod);

          //creates read 1
          std::strcat(read1,Adapter1); //add adapter to read
          std::strncpy(read1N, read1, strlen(read1)); //copy read+adapter into N...

          //creates read 2
          DNA_complement(read2);
          std::reverse(read2, read2 + strlen(read2));
          std::strcat(read2,Adapter2);
          std::strncpy(read2N, read2, strlen(read2)); //copy read+adapter into N...

          if (std::strcmp(outflag, "gz") == 0){
            kstring_t kstr1;// kstring_t *kstr = new kstr;
            kstr1.s = NULL; // kstr->s = NULL;
            kstr1.l = kstr1.m = 0; //kstr->l = kstr->m = 0;

            ksprintf(&kstr1,">%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
            ksprintf(&kstr1,"%s\n",read1N);
            ksprintf(&kstr1,"%s\n","+");  
            Read_Qual2(read1N,qual,Qualdistr1,gen);
            ksprintf(&kstr1,"%s\n",qual);
            memset(qual, 0, sizeof(qual));    
            gzwrite(gz1,kstr1.s,kstr1.l);kstr1.l =0;

            kstring_t kstr2;// kstring_t *kstr = new kstr;
            kstr2.s = NULL; // kstr->s = NULL;
            kstr2.l = kstr2.m = 0; //kstr->l = kstr->m = 0;

            ksprintf(&kstr2,">%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
            ksprintf(&kstr2,"%s\n",read2N);
            ksprintf(&kstr2,"%s\n","+");  
            Read_Qual2(read2N,qual,Qualdistr2,gen);
            ksprintf(&kstr2,"%s\n",qual);    
            memset(qual, 0, sizeof(qual));
            gzwrite(gz2,kstr2.s,kstr2.l);kstr2.l =0;
          }
          else if (std::strcmp(outflag, "fq") == 0){
            kstring_t kstr1;// kstring_t *kstr = new kstr;
            kstr1.s = NULL; // kstr->s = NULL;
            kstr1.l = kstr1.m = 0; //kstr->l = kstr->m = 0;

            kstring_t kstr2;// kstring_t *kstr = new kstr;
            kstr2.s = NULL; // kstr->s = NULL;
            kstr2.l = kstr2.m = 0; //kstr->l = kstr->m = 0;

            ksprintf(&kstr1,">%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
            ksprintf(&kstr1,"%s\n",read1N);
            ksprintf(&kstr1,"%s\n","+");  
            Read_Qual2(read1N,qual,Qualdistr1,gen);
            ksprintf(&kstr1,"%s\n",qual);    
            fwrite(kstr1.s,sizeof(char),kstr1.l,fp1);kstr1.l =0;
            /*fprintf(fp1,"@%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
            fprintf(fp1,"%s\n",seqmod);
            fprintf(fp1,"%s\n","+");
            Read_Qual2(seqmod,qual,Qualdistr1,gen);
            fprintf(fp1,"%s\n",qual);*/
            memset(qual, 0, sizeof(qual));

            ksprintf(&kstr2,">%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
            ksprintf(&kstr2,"%s\n",read2N);
            ksprintf(&kstr2,"%s\n","+");  
            Read_Qual2(read2N,qual,Qualdistr2,gen);
            ksprintf(&kstr2,"%s\n",qual);  
            fwrite(kstr2.s,sizeof(char),kstr2.l,fp2);kstr1.l =0;

            /*fprintf(fp2,"@%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
            fprintf(fp2,"%s\n",read2N);
            fprintf(fp2,"%s\n","+");
            Read_Qual2(read2N,qual,Qualdistr2,gen);
            fprintf(fp2,"%s\n",qual);*/
            memset(qual, 0, sizeof(qual));
          }
        }
        else{
          if (std::strcmp(outflag, "gz") == 0){
            kstring_t kstr1;// kstring_t *kstr = new kstr;
            kstr1.s = NULL; // kstr->s = NULL;
            kstr1.l = kstr1.m = 0; //kstr->l = kstr->m = 0;

            ksprintf(&kstr1,">%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
            ksprintf(&kstr1,"%s\n",seqmod);
            ksprintf(&kstr1,"%s\n","+");  
            Read_Qual2(seqmod,qual,Qualdistr1,gen);
            ksprintf(&kstr1,"%s\n",qual);
            memset(qual, 0, sizeof(qual));    
            gzwrite(gz1,kstr1.s,kstr1.l);kstr1.l =0;

            kstring_t kstr2;// kstring_t *kstr = new kstr;
            kstr2.s = NULL; // kstr->s = NULL;
            kstr2.l = kstr2.m = 0; //kstr->l = kstr->m = 0;

            DNA_complement(seqmod);
            std::reverse(seqmod, seqmod + strlen(seqmod));
            Read_Qual2(seqmod,qual,Qualdistr2,gen);

            ksprintf(&kstr2,">%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
            ksprintf(&kstr2,"%s\n",seqmod);
            ksprintf(&kstr2,"%s\n","+");  
            ksprintf(&kstr2,"%s\n",qual);    
            memset(qual, 0, sizeof(qual));
            gzwrite(gz2,kstr2.s,kstr2.l);kstr2.l =0;
          }
          else if (std::strcmp(outflag, "fq") == 0){
            kstring_t kstr1;// kstring_t *kstr = new kstr;
            kstr1.s = NULL; // kstr->s = NULL;
            kstr1.l = kstr1.m = 0; //kstr->l = kstr->m = 0;

            kstring_t kstr2;// kstring_t *kstr = new kstr;
            kstr2.s = NULL; // kstr->s = NULL;
            kstr2.l = kstr2.m = 0; //kstr->l = kstr->m = 0;

            ksprintf(&kstr1,">%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
            ksprintf(&kstr1,"%s\n",seqmod);
            ksprintf(&kstr1,"%s\n","+");  
            Read_Qual2(seqmod,qual,Qualdistr1,gen);
            ksprintf(&kstr1,"%s\n",qual);    
            fwrite(kstr1.s,sizeof(char),kstr1.l,fp1);kstr1.l =0;
            /*fprintf(fp1,"@%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
            fprintf(fp1,"%s\n",seqmod);
            fprintf(fp1,"%s\n","+");
            Read_Qual2(seqmod,qual,Qualdistr1,gen);
            fprintf(fp1,"%s\n",qual);*/
            memset(qual, 0, sizeof(qual));
            
            DNA_complement(seqmod);
            std::reverse(seqmod, seqmod + strlen(seqmod));

            ksprintf(&kstr2,">%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
            ksprintf(&kstr2,"%s\n",seqmod);
            ksprintf(&kstr2,"%s\n","+");  
            Read_Qual2(seqmod,qual,Qualdistr2,gen);
            ksprintf(&kstr2,"%s\n",qual);  
            fwrite(kstr2.s,sizeof(char),kstr2.l,fp1);kstr1.l =0;
            /*fprintf(fp2,"@%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
            fprintf(fp2,"%s\n",seqmod);
            fprintf(fp2,"%s\n","+");
            Read_Qual2(seqmod,qual,Qualdistr2,gen);
            fprintf(fp2,"%s\n",qual);*/
            memset(qual, 0, sizeof(qual));
          }
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
  // /willerslev/users-shared/science-snm-willerslev-wql443/reference_files/Human/hg19canon.fa
  // "/home/wql443/scratch/reference_genome/hg19/chr2122.fa"
  const char *fastafile = "/home/wql443/scratch/reference_genome/hg19/chr2122.fa";
  int seed = -1;

  if (std::strcmp(argv[1], "fafq") == 0){
    // ./a.out fafq fq true /home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt /home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R2.txt
    if(argc>6){seed = atoi(argv[6]);}
    else{seed = time(NULL);} 
    fprintf(stderr,"seed is: %d\n",seed);
    srand48(seed);
    //const char *Profile1 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt";
    //const char *Profile2 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R2.txt";
    const char *Profile1 = argv[4];
    const char *Profile2 = argv[5];
    char Adapter1[] = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG";
    char Adapter2[] = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT";
    if (std::strcmp(argv[2], "gz") == 0){
      gzFile gz1;
      gz1 = gzopen("R1.fq.gz","wb");
      gzFile gz2;
      gz2 = gzopen("R2.fq.gz","wb");
      if (std::strcmp(argv[3], "False") == 0 || std::strcmp(argv[3], "false") == 0 || std::strcmp(argv[3], "F") == 0){
        fafq(fastafile,argv[2],NULL,NULL,gz1,gz2,false,Profile1,Profile2);
      }
      else if (std::strcmp(argv[3], "True") == 0 || std::strcmp(argv[3], "true") == 0 || std::strcmp(argv[3], "T") == 0){
        fafq(fastafile,argv[2],NULL,NULL,gz1,gz2,true,Profile1,Profile2,Adapter1,Adapter2);
      }
      gzclose(gz1);
      gzclose(gz2);
    }
    else if (std::strcmp(argv[2], "fq") == 0){
      FILE *fp1;
      FILE *fp2;
      fp1 = fopen("R1.fq","wb");
      fp2 = fopen("R2.fq","wb");
      //char Adapter1[] = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG";
      //char Adapter2[] = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT";
      //char Adapter1[] = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
      //char Adapter2[] = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";
      if (std::strcmp(argv[3], "False") == 0 ||
          std::strcmp(argv[3], "false") == 0 ||
          std::strcmp(argv[3], "F") == 0){fafq(fastafile,argv[2],fp1,fp2,NULL,NULL,false,Profile1,Profile2);}
      else if (std::strcmp(argv[3], "True") == 0 ||
              std::strcmp(argv[3], "true") == 0 ||
              std::strcmp(argv[3], "T") == 0){fafq(fastafile,argv[2],fp1,fp2,NULL,NULL,true,Profile1,Profile2,Adapter1,Adapter2);}
      fclose(fp1);
      fclose(fp2);
    }
  }
  else {
    printf("use mode fafq\n");
  }
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
}


/*
int main(){
  int sz_mem = 150;
  char *N = new char[sz_mem+1];
  strcpy(N, std::string(sz_mem, 'N').c_str());
  
  char *Frag = new char[1024];
  strcpy(Frag, "AAa");
  char Adapter1[1024] = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
  
  std::cout << N << std::endl;
  strcat(Frag,Adapter1);
  std::cout << Frag << std::endl;
  
  strncpy(N, Frag, strlen(Frag)); //copy read+adapter into N...
  std::cout << N << std::endl;

  return 0;
}
*/
// g++ XXX.cpp -std=c++11 -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lz -lbz2 -llzma -lcurl
//  ./a.out fafq fq T Qual_profiles/Freq_R1.txt Qual_profiles/Freq_R2.txt 1614849265