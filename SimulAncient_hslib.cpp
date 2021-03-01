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


// I would like to create a function with TK's code since its optimal in case we wish to 
//simulate a given number of fragments

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

std::string Read_Qual(char *seq,std::discrete_distribution<>*Dist,std::default_random_engine &gen){
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
      case 'a':
        qual += ASCII_qual[Dist[row_idx](gen)];
        break;
      case 'T':
      case 't':
        qual += ASCII_qual[Dist[row_idx + Tstart](gen)];
        break;  
      case 'G':
      case 'g':
        qual += ASCII_qual[Dist[row_idx + Gstart](gen)];
        break;
      case 'C':
      case 'c':
        qual += ASCII_qual[Dist[row_idx + Cstart](gen)];
        break;
      case 'N':
        qual += "!"; //character of 33 decimal or 21 hex
        break;
    }
  }
  return qual;
}


void Read_Qual2(char *seq,char *qual,std::discrete_distribution<>*Dist,std::default_random_engine &gen){
  //std::string qual;
  //char qual[1024] = "";
  const char* array[8] = {"#", "\'", "0", "7" ,"<", "B", "F","I"};

  int read_length = strlen(seq);

  int Tstart = 150;
  int Gstart = 300;
  int Cstart = 450;
  
  for (int row_idx = 0; row_idx < read_length; row_idx++){
    switch(seq[row_idx]){
      case 'A':
      case 'a':
        strncat(qual, array[Dist[row_idx](gen)], 1);
        break;
      case 'T':
      case 't':
        strncat(qual, array[Dist[row_idx + Tstart](gen)], 1);
        break;  
      case 'G':
      case 'g':
        strncat(qual, array[Dist[row_idx + Gstart](gen)], 1);
        break;
      case 'C':
      case 'c':
        strncat(qual, array[Dist[row_idx + Cstart](gen)], 1);
        break;
      case 'N':
        strncat(qual, array[0], 1);;
        break;
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

double** create2DArray(const char* filename,int width,int height){
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

void Distfunc(double** Array2d,std::discrete_distribution<> dist[],int size){
  for (int row_idx = 0; row_idx < size; row_idx++){
    std::discrete_distribution<> d({Array2d[row_idx][0],Array2d[row_idx][1],Array2d[row_idx][2],Array2d[row_idx][3],
                                    Array2d[row_idx][4],Array2d[row_idx][5],Array2d[row_idx][6],Array2d[row_idx][7]});  
    dist[row_idx] = d;
  }
}

void Ill_err(char *seq,std::discrete_distribution<>*Dist,std::default_random_engine &gen){
  int read_length = strlen(seq);
  int Tstart = 70;
  int Gstart = 140;
  int Cstart = 210;
  const char LookUp_nt[4] = {'A','T','G','C'};
  for (int nt_idx = 0; nt_idx < read_length; nt_idx++){
    //std::cout << "row indx" << seq[row_idx] << std::endl;
    //std::cout << Qual_random(Array2d[row_idx],gen);
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

void Seq_err(double** Array2d,std::discrete_distribution<> nt_sub[],int size){
  for (int row_idx = 0; row_idx < size; row_idx++){
    std::discrete_distribution<> d({Array2d[row_idx][0],
                                    Array2d[row_idx][1],
                                    Array2d[row_idx][2],
                                    Array2d[row_idx][3]});  
    nt_sub[row_idx] = d;
  }
}

//--------------------------FUNCTIONS FOR SEQUENCES----------------------

void fafa(const char* fastafile, const char* outflag,FILE *fp,gzFile gz){
  //creates fasta structure
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);
  
  double init = 1.0;
  int chr_no = 0;
  //creates randomized fragment length based on drand seed in main
  
  //read entire chromosome
  
  //creates sequence of random length for DNA damage

  while (chr_no < faidx_nseq(seq_ref)){
    //int whichref = lrand48() % faidx_nseq(seq_ref);
    //const char *name = faidx_iseq(seq_ref,whichref);
    const char *name = faidx_iseq(seq_ref,chr_no);
    int name_len =  faidx_seq_len(seq_ref,name);
    fprintf(stderr,"-> name: \'%s\' name_len: %d\n",name,name_len);
    char *data = fai_fetch(seq_ref,name,&name_len);
    
    int start_pos = 1;
    int end_pos = name_len; //30001000
    char seqmod[1024];
    
    //kstring_t kstr;// kstring_t *kstr = new kstr;
    //kstr.s = NULL; // kstr->s = NULL;
    //kstr.l = kstr.m = 0; //kstr->l = kstr->m = 0;

    while(start_pos <= end_pos){
      //std::cout << "------------ " << std::endl;
      //std::cout << "Start " <<start_pos << std::endl;
      //std::srand(start_pos+std::time(nullptr));
      int readlength = drand48()*(80.0-30.0)+30.0;
      //std::cout << " random " << readlength << std::endl;
      int stop = start_pos+(int) readlength;
      
      //adress of entire array, first element always the same but adding start for query source start, and stop-start : number of elements
      strncpy(seqmod,data+start_pos,readlength);

      char * pch;
      pch = strchr(seqmod,'N');
      if (pch != NULL){start_pos += readlength + 1;}
      else {
        char nt[] = "tT";
        Deamin_char(seqmod,nt,readlength);
        // sequence.size(); //we can use .size if we did the std::string approach which we did with deamin_string calling it damage
        if (std::strcmp(outflag, "gz") == 0){
          kstring_t kstr;// kstring_t *kstr = new kstr;
          kstr.s = NULL; // kstr->s = NULL;
          kstr.l = kstr.m = 0; //kstr->l = kstr->m = 0;

          ksprintf(&kstr,">%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
          ksprintf(&kstr,"%s\n",seqmod);
          
          //if (kstr.l > 4000){
          //  gzwrite(gz,kstr.s,kstr.l);kstr.l =0;
          //}

          gzwrite(gz,kstr.s,kstr.l);kstr.l =0;

        }
        else
        {
          fprintf(fp,">%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
          fprintf(fp,"%s\n",seqmod);
        }
        //fprintf(fp,">%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
        //fprintf(fp,"%s\n",seqmod);

        start_pos += readlength + 1;
        readlength = 0;
        memset(seqmod, 0, sizeof seqmod);
      }
    }
  chr_no++;
  }
  //gzclose(gz);
}

void fafq(const char* fastafile,FILE *fp1,FILE *fp2,const char* r1_profile,const char* r2_profile,
          bool flag=false, double cov=1.0){
  
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *seq_ref = NULL;
  seq_ref  = fai_load(fastafile);
  assert(seq_ref!=NULL);//check that we could load the file

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(seq_ref));
  double init = 1.0;
  int chr_no = 0;
  //create the file for savin
  std::ifstream file(r1_profile);
  // change Read_size to profile_line
  int Read_size = std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n');
  file.close();
  // int Qualdist_size = 600;

  double** Error_2darray = create2DArray("/home/wql443/WP1/Ill_err.txt",4,280);
  double** R1_2Darray = create2DArray(r1_profile,8,Read_size);
  double** R2_2Darray = create2DArray(r2_profile,8,Read_size);

  std::random_device rd;
  std::default_random_engine gen(rd()); 

  std::discrete_distribution<> Qualdistr1[Read_size];
  Distfunc(R1_2Darray,Qualdistr1,Read_size);
  std::discrete_distribution<> Qualdistr2[Read_size];
  Distfunc(R2_2Darray,Qualdistr2,Read_size);

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
      std::srand(start_pos+std::time(nullptr));
      // Seed random number generator 
      int readlength = (std::rand() % (60 - 30 + 1)) + 30;
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
        char qual[1024] = "";

        Ill_err(seqmod,Error,gen);
        //Deamin_char(seqmod,nt,readlength);
        std::string Read_qual;
        if(flag==true){
          char Adapter1[] = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG";
          char Adapter2[] = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT";
          char read1N[Read_size/4 + 1];
          char read2N[Read_size/4 + 1];
          strcpy(read1N, std::string(150, 'N').c_str()); // create read full of N
          strcpy(read2N, std::string(150, 'N').c_str()); // create read full of N
          
          char read1[Read_size/4 + 1];
          char read2[Read_size/4 + 1];
          //Copies sequence into both reads
          strcpy(read1, seqmod);
          strcpy(read2, seqmod);

          //creates read 1
          std::strcat(read1,Adapter1); //add adapter to read
          std::strncpy(read1N, read1, strlen(read1)); //copy read+adapter into N...
          fprintf(fp1,"@%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
          fprintf(fp1,"%s\n",read1N);
          fprintf(fp1,"%s\n","+");
          Read_Qual2(read1N,qual,Qualdistr1,gen);
          fprintf(fp1,"%s\n",qual);
          memset(qual, 0, sizeof(qual));

          //creates read 2
          DNA_complement(read2);
          std::reverse(read2, read2 + strlen(read2));
          std::strcat(read2,Adapter2);
          std::strncpy(read2N, read2, strlen(read2)); //copy read+adapter into N...

          fprintf(fp2,"@%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
          fprintf(fp2,"%s\n",read2N);
          fprintf(fp2,"%s\n","+");
          Read_Qual2(read2N,qual,Qualdistr2,gen);
          fprintf(fp2,"%s\n",qual);
          memset(qual, 0, sizeof(qual));
        }
        else{
          //char Ill_r1[150];
          //strncpy(Ill_r1,seqmod, sizeof(Ill_r1));

          fprintf(fp1,"@%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
          fprintf(fp1,"%s\n",seqmod);
          fprintf(fp1,"%s\n","+");
          Read_Qual2(seqmod,qual,Qualdistr1,gen);
          fprintf(fp1,"%s\n",qual);
          memset(qual, 0, sizeof(qual));

          //char Ill_r2[150];
          //strncpy(Ill_r2,seqmod, sizeof(Ill_r2));
          
          DNA_complement(seqmod);
          std::reverse(seqmod, seqmod + strlen(seqmod));

          fprintf(fp2,"@%s:%d-%d_length:%d\n",name,start_pos,start_pos+readlength,readlength);
          fprintf(fp2,"%s\n",seqmod);
          fprintf(fp2,"%s\n","+");
          Read_Qual2(seqmod,qual,Qualdistr2,gen);
          fprintf(fp2,"%s\n",qual);
          memset(qual, 0, sizeof(qual));
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


void faBam(const char* fastafile,const char* nt_profile,double cov=1.0){
  double** my2DArray = create2DArray(nt_profile,8,150);
  
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *ref = NULL;
  ref  = fai_load(fastafile);
  assert(ref!=NULL);//check that we could load the file

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(ref));
  double init = 1.0;
  int chr_no = 0;
  //create the file for saving
  
  //Creates a pointer to allocated memomry for the format??
  htsFormat *fmt_hts =(htsFormat*) calloc(1,sizeof(htsFormat));
  //wb -> bam , wc -> cram
  char out_mode[5]="wb";
    
  const char *outfile_nam = "test.bam";
  samFile *outfile = NULL;

  if ((outfile = sam_open_format(outfile_nam, out_mode, fmt_hts)) == 0) {
    fprintf(stderr,"Error opening file for writing\n");
    exit(0);
  }
  
  // creates a pointer to generated header
  sam_hdr_t *header = sam_hdr_init();
  if (header == NULL) { fprintf(stderr, "sam_hdr_init");}
  char *refName=NULL;
    
  // Creating sam_header
  char *name_len_char =(char*) malloc(1024);

  for(int i=0;i<faidx_nseq(ref);i++){
    const char *name = faidx_iseq(ref,i);
    int name_len =  faidx_seq_len(ref,name);
    // skal være c string i array så vi konvertere int om til char
    snprintf(name_len_char,1024,"%d",name_len);
    fprintf(stderr,"ref:%d %d %s\n",i,name,name_len_char);
    // reference part of the header, int r variable ensures the header is added
    int r = sam_hdr_add_line(header, "SQ", "SN", name, "LN", name_len_char, NULL);
    if (r < 0) { fprintf(stderr,"sam_hdr_add_line");}
  }
  
  // saving the header to the file
  if (sam_hdr_write(outfile, header) < 0) fprintf(stderr,"writing headers to %s", outfile);
  
  // alignment delen, gemt i bam_1 type for at representere hver linje som 1 alignment
  bam1_t *bam_file = bam_init1(); //initialisere bam1_t (type) til hukommelsen!

  while (chr_no < faidx_nseq(ref)){  
    const char *name = faidx_iseq(ref,chr_no);
    int name_len =  faidx_seq_len(ref,name);
    int start_pos = 1;
    int end_pos = name_len; //30000000
    
    while(start_pos <= end_pos){
      char Qname[96];
      std::srand(start_pos+(std::rand() % (80 - 30 + 1)) + 30 +std::time(nullptr));
      // Seed random number generator
      // creates random number in the range of the fragment size rand() % ( high - low + 1 ) + low
      int rand_len = (std::rand() % (80 - 30 + 1)) + 30;
      int dist = init/cov * rand_len;
      //we use structure faidx_t from htslib to load in a fasta
      char* sequence = faidx_fetch_seq(ref,name,start_pos,start_pos+rand_len,&name_len);
      snprintf(Qname,96,"%s:%d-%d",name,start_pos,start_pos+rand_len);

      char * pch;
      pch = strchr(sequence,'N');

      if (pch != NULL){
        //Disregards any read with 'N' in.. change this to just change the reading position
        start_pos += dist + 1;
      }
      else {
        std::string Read_qual;
        char nt[] = "tT";
        Deamin_char(sequence,nt,rand_len);

        //ensures proper bam format with 11 mandatory fields mandatory fields
        // QNAME = char Qname[96] 
        int Flag = 4; // 4 for unmapped
        int RNAME = chr_no; // Reference sequence name, chr_no takes chr from sam_hdr_t
        // POS = int start_pos -1 as it is already 1-based from faidx but bam_set1 converts it further to 1-based
        int mapQ = 255; // 255 if unavailable du to no mapping
        const uint32_t *cigar = NULL; // cigar string, NULL if unavailable - But how do we then add actual cigar info?
        // RNEXT = mtid -> chr_no
        int Pnext = 0; // position for next mate -> mpos
        int Tlen = 0; // template length -> isize
        // SEQ = sequence
        // QUAL = quality

        //casting the ascii values to strings, but im substracting -33 from the values, since bam used phred+33.
        for (int row_idx = 0; row_idx < strlen(sequence); row_idx++){Read_qual += Qual_random(my2DArray[row_idx])-33;}

        const char* nt_qual = Read_qual.c_str();
        int no_cigar = 0; //number of cigar operations 

        //bam, lqname, qname, flag, tid, pos, mapq, n_cigar, cigar, mtid, mpos, isize, l_seq, seq, qual, l_aux
        bam_set1(bam_file,strlen(Qname),Qname,Flag,RNAME,start_pos-1,mapQ,no_cigar,cigar,chr_no,Pnext-1,Tlen,rand_len,sequence,nt_qual,0);
        // -1 for the different positions is because they are made 1 - based with the bam_set
        
        sam_write1(outfile,header,bam_file);
        Read_qual = ""; 
      }
      start_pos += dist + 1;
    }
    chr_no++;
  }
  sam_hdr_destroy(header);
  sam_close(outfile);
  return; 
}

/*
int main(int argc,char **argv){
  clock_t tStart = clock();
  const char *fastafile = "/willerslev/users-shared/science-snm-willerslev-wql443/reference_files/Human/hg19canon.fa";

  if (std::strcmp(argv[1], "fa") == 0){
    fafa(fastafile,argv[2]);
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  }
  else if (std::strcmp(argv[1], "fq") == 0){
    const char *Profile = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq.txt";
    fafq(fastafile,argv[2],Profile,false);
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  }
  else if (std::strcmp(argv[1],"bam")==0){
    const char *Profile = "/home/wql443/WP1/hslib_test/Freq.txt";
    faBam(fastafile,Profile);
    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
  }
}*/


int main(int argc,char **argv){
  clock_t tStart = clock();
  // /willerslev/users-shared/science-snm-willerslev-wql443/reference_files/Human/hg19canon.fa
  // "/home/wql443/scratch/reference_genome/hg19/chr22.fa"
  const char *fastafile = "/home/wql443/scratch/reference_genome/hg19/chr2122.fa";
  int seed = -1;

  if (std::strcmp(argv[1], "fafa") == 0){
    // ./a.out fafa fa /home/wql443/scratch/reference_genome/hg19/chr14.fa test.fa
    if(argc>5){seed = atoi(argv[5]);}
    else{seed = time(NULL);} 
    fprintf(stderr,"seed is: %d\n",seed);
    srand48(seed);

    const char *fastafile = argv[3];
    if (std::strcmp(argv[2], "gz") == 0){
      gzFile gz;
      gz = gzopen(argv[4],"wb");
      fafa(fastafile,argv[2],NULL,gz); 
      gzclose(gz);
    }
    else if (std::strcmp(argv[2], "fa") == 0){
      FILE *fp;
      fp = fopen(argv[4],"wb");
      fafa(fastafile,argv[2],fp,NULL); 
      fclose(fp);
    }
  }
  else if (std::strcmp(argv[1], "fafq") == 0){
    // ./a.out fafq fq true /home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt /home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R2.txt
    if(argc>6){seed = atoi(argv[6]);}
    else{seed = time(NULL);} 
    fprintf(stderr,"seed is: %d\n",seed);
    srand48(seed);
    //const char *Profile1 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R1.txt";
    //const char *Profile2 = "/home/wql443/WP1/SimulAncient/Qual_profiles/Freq_R2.txt";
    if (std::strcmp(argv[2], "gz") == 0){
      std::cout << "YET TO BE IMPLEMENTED" << std::endl;
    }
    else if (std::strcmp(argv[2], "fq") == 0){
      FILE *fp1;
      FILE *fp2;
      fp1 = fopen("lol_R1.fq","wb");
      fp2 = fopen("lol_R2.fq","wb");
      const char *Profile1 = argv[4];
      const char *Profile2 = argv[5];
      char Adapter1[] = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG";
      char Adapter2[] = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTNNNNNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT";
      if (std::strcmp(argv[3], "false") == 0){fafq(fastafile,fp1,fp2,Profile1,Profile2,false);}
      else if (std::strcmp(argv[3], "true") == 0){fafq(fastafile,fp1,fp2,Profile1,Profile2,true);}
      fclose(fp1);
      fclose(fp2);
    }
  }
  else if (std::strcmp(argv[1],"bam")==0){
    if(argc>2){seed = atoi(argv[2]);}
    else{seed = time(NULL);} 
    fprintf(stderr,"seed is: %d\n",seed);
    srand48(seed);
    const char *Profile = "/home/wql443/WP1/hslib_test/Freq.txt";
    faBam(fastafile,Profile);
  }
  printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
}
