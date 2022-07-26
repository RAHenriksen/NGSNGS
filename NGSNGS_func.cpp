#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm> //std::reverse
#include <iostream>
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>
#include <errno.h>
#include <random>
#include <map>
#include <math.h>
#include <pthread.h>

#include "NGSNGS_func.h"
#include "NtSubModels.h"
#include "mrand.h"
#include "RandSampling.h"

//#if defined(__APPLE__) && defined(__MACH__) 
//#include "NGSNGS_Random.h"
//#endif /* __APPLE__ */

#define LENS 4096
#define MAXBINS 100

void FragDistArray(int& number,int*& Length, double*& Frequency,int SizeDistType,int seed,int val1, int val2){
  // Generate a similar structure as length distribution file, with fragment length and frequency using length distributions
  std::default_random_engine generator(seed);
  const int nrolls=10000;
  int p[nrolls]={};
  
  if (SizeDistType==1){std::uniform_int_distribution<int> distribution(val1,val2);for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}} //if(number > 30){p[i] = number;}}
  else if (SizeDistType==2){std::normal_distribution<double> distribution(val1,val2);for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}}
  else if (SizeDistType==3){std::lognormal_distribution<double> distribution(val1,val2);for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}}
  else if (SizeDistType==4){std::poisson_distribution<int> distribution(val1);for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}}
  else if (SizeDistType==5){std::exponential_distribution<double> distribution(val1);for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}}
  else if (SizeDistType==6){std::gamma_distribution<double> distribution(val1,val2);for (int i=0; i<nrolls; ++i) {int number = distribution(generator);p[i] = number;}}

  
  int len = sizeof(p)/sizeof(p[0]);
  sort(p,p+len,std::less<int>());

  std::map<int,int> counts;
  for (int i =0; i<nrolls;i++){counts[p[i]]++;}

  int n =1;

  Length[0] = 0; Frequency[0] = (float) 0;
  double CumulativeCount = 0;
  for (auto const &key : counts){
      CumulativeCount += (double) key.second/(double)nrolls;
      Length[n] = key.first;
      Frequency[n] = CumulativeCount;
      n++;
  }
  number = n;
}

void FragArray(int& number,int*& Length, double*& Frequency,const char* filename){
  //int LENS = 4096;
  //int* Frag_len = new int[LENS];
  //double* Frag_freq = new double[LENS];
  int n =1;

  gzFile gz = Z_NULL;
  char buf[LENS];
  
  gz = gzopen(filename,"r");
  assert(gz!=Z_NULL);
  Length[0] = 0; Frequency[0] = (float) 0;
  while(gzgets(gz,buf,LENS)){
    Length[n] = atoi(strtok(buf,"\n\t ")); //before it was Frag_len[n]
    Frequency[n] = atof(strtok(NULL,"\n\t "));
    n++;
  }
  gzclose(gz); 

  number = n;
  //Length = Frag_len;
  //Frequency = Frag_freq;
  //delete[] Frag_len;
  //delete[] Frag_freq;
}

void delete_seq(char *str, int seq_len, int del_len, size_t pos,int alt_len){
    //instert_seq(char *str, size_t len, char insert_seq[],size_t ins_len, size_t pos){
    for (int i = alt_len; i < del_len; i++)
    {
        //fprintf(stderr,"deletion pos 1 %d \t and pos %zu\n",i,pos-seq_len);
        //std::cout << str[pos] << " " << str[pos+1] << std::endl;
        //std::cout << seq_len << " " << pos << " " << pos - seq_len << std::endl;
        memmove(&str[pos], &str[pos+1], pos-seq_len);
        //memmove(&str[pos-1], &str[pos], pos-seq_len-1);
        //fprintf(stderr,"deletion pos 2 %d\n",seq_len - pos);
        seq_len--;
    }
}

void delete_seq_ins(char *str, int seq_len, int del_len, size_t pos){
    // after insertion it deleted the reference allele, e.g. REF A -> ALT AGGGGGG, which creates a 6 bp insertion
    for (int i = 0; i < del_len; i++)
    {
        memmove(&str[pos], &str[pos+1], seq_len - pos);
        seq_len--;
    }
}

void instert_seq(char *str, int len, char insert_seq[],int ins_len, size_t pos){
    for (int i = 0; i < ins_len; i++)
    {
        memmove(&str[pos+1], &str[pos], len - pos + 1);
        str[pos] = insert_seq[i];
        len++;
    }
    //the insertion doesn't overlap so i remove the nucleotide which the insertion replace
    delete_seq_ins(str, (int) strlen(str),1,pos+ins_len);
}

void DNA_CAPITAL(char seq[]){
  while (*seq) {
    switch(*seq) {
      case 'A':
      case 'a':
        *seq = 'A';
        break;
      case 'G':
      case 'g':
        *seq = 'G';
        break;
      case 'C':
      case 'c':
        *seq = 'C';
        break;
      case 'T':
      case 't':
        *seq = 'T';
        break;
      case 'N':
      case 'n':
        *seq = 'N';
        break;  
    }
    ++seq;
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
      case 'N':
      case 'n':
        *seq = 'N';
        break;  
    }
    ++seq;
  }
}

void reverseChar(char* str,int length) {
    std::reverse(str, str + length);
}

void ReversComplement(char seq[]){
  // generates the reverse complementary sequence from an input sequence
  unsigned char nuc2int[255];
  nuc2int['a'] = nuc2int['A'] = nuc2int[0] = 0;
  nuc2int['t'] = nuc2int['T'] = nuc2int[1] = 1;
  nuc2int['g'] = nuc2int['G'] = nuc2int[2] = 2;
  nuc2int['c'] = nuc2int['C'] = nuc2int[3] = 3;
  nuc2int['n'] = nuc2int['N'] = nuc2int[4] = 4;

  char ntdeam[4] = {'T', 'A', 'C', 'G'};
  char seq_intermediate[1024] = {0};
  strcpy(seq_intermediate,seq);

  //Complementing sequence
  for(int i=0;i<strlen(seq);i++){
    seq_intermediate[i] = ntdeam[nuc2int[seq_intermediate[i]]];
  }

  //reverse complement
  for(int i=strlen(seq)-1;i>-1;i--){
    seq[strlen(seq)-i-1] = seq_intermediate[i];
  }

  //just to ensure no issues arise in case of not clearing out the intermediate sequence
  memset(seq_intermediate, 0, sizeof seq_intermediate);
}

void deletechar(char* array,int seq_len, size_t index_to_remove, int del_len){
  /*memmove(&array[index_to_remove],
        &array[index_to_remove + 5],
        seq_len - index_to_remove - 5);*/
  std::string s(array);
  s.erase(s.begin()+index_to_remove, s.end()-(seq_len-index_to_remove-del_len));
  strcpy(array, s.c_str());
}

void InsertChar(char* array,std::string ins,int index){
  std::string s(array);
  s.insert(index,ins);  
  strcpy(array, s.c_str());    
}

void Header_func(htsFormat *fmt_hts,const char *outfile_nam,samFile *outfile,sam_hdr_t *header,faidx_t *seq_ref,int chr_total,int chr_idx_arr[],size_t genome_len,char CommandArray[1024],const char* version){
  // Creates a header for the bamfile. The header is initialized before the function is called //

  if (header == NULL) { fprintf(stderr, "sam_hdr_init");}
    
  // Creating header information
  
  char genome_len_buf[1024];
  for(int i=0;i<chr_total;i++){
    const char *name = faidx_iseq(seq_ref,chr_idx_arr[i]);
    //fprintf(stderr,"chromosomes added to bam header %d\n",chr_idx_arr[i]);

    int name_len =  faidx_seq_len(seq_ref,name);
    snprintf(genome_len_buf,1024,"%d", name_len);
    
    // reference part of the header, int r variable ensures the header is added
    int r = sam_hdr_add_line(header, "SQ", "SN", name, "LN", genome_len_buf, NULL);
    if (r < 0) { fprintf(stderr,"sam_hdr_add_line");}
    memset(genome_len_buf,0, sizeof(genome_len_buf));
    //snprintf(genome_len_buf,1024,"COMMAND ./ngsngs");
    //int r = sam_hdr_add_line(header, "TEST COMMAND", genome_len_buf, NULL);
    memset(genome_len_buf,0, sizeof(genome_len_buf));
  }
  // Adding PG tag
  sam_hdr_add_pg(header,"NGSNGS","VN",version,"CL",CommandArray,NULL);
  // saving the header to the file
  if (sam_hdr_write(outfile, header) < 0) fprintf(stderr,"writing headers to %s", outfile_nam); //outfile
}

char* HaploGenome(char* genome,char genome_data1[],char genome_data2[],int chr_sizes,const char* bcf_file,const char *chr_names[],const char* VarType,FILE *VarInfoFile,const char* HeaderIndiv){

  DNA_CAPITAL(genome_data1);DNA_CAPITAL(genome_data2);

  htsFile *bcf_obj = NULL;
  bcf_hdr_t *bcf_head = NULL;
  // "/home/wql443/WP1/data/Chr14_subset.vcf  "/home/wql443/WP1/data/test2.vcf"   "/home/wql443/WP1/data/chr14_full.vcf"
  bcf_obj = bcf_open(bcf_file, "r");
  hts_idx_t *idx = bcf_index_load(bcf_file);
  bcf_head = bcf_hdr_read(bcf_obj);
  fprintf(stderr, "File contains before %i samples\n",bcf_hdr_nsamples(bcf_head)); //bcf.nsamples(fp)
  int nseq = 0;

  const char **seqnames = NULL;
  seqnames = bcf_hdr_seqnames(bcf_head, &nseq);

  // check if the headers match
  int chr_exit = 0;
  for (int vcf_i = 0; vcf_i < nseq; vcf_i++){
    if(strcasecmp(seqnames[vcf_i],chr_names[0])==0){continue;}//fprintf(stderr,"number of chromosomes %d and names %s\n",chr_i,chr_names[chr_i]);
    else{chr_exit++;}
  }
  if(chr_exit == nseq){
    fprintf(stderr,"Reference file chromosome %s is not identified in vcf/bcf header\n",chr_names[0]);
    exit(0);
  }
  
  // Set sample number 
  //const char* HeaderIndiv = "HG00096";
  
  //bcf_hdr_t *bcf_head2 =bcf_hdr_subset(const bcf_hdr_t *h0, int n, char *const* samples, int *imap);
  
  //int indiv = bcf_hdr_set_samples(bcf_head, "HG00099", 0);
  //fprintf(stderr, "File contains after %i samples %d\n",bcf_hdr_nsamples(bcf_head),indiv); //bcf.nsamples(fp)
  //bcf_hdr_parse_sample_line(bcf_head,htxt)
  //initializing variation sampling 
  bcf1_t *bcf_records = bcf_init();
  if (bcf_records == NULL){fprintf(stderr,"WARNING NO VARIANTS IN VCF FILE"); exit(0);}
  else if(bcf_records != NULL && bcf_read(bcf_obj, bcf_head, bcf_records)==0){
    // else if{bcf_records != NULL && bcf_read(bcf_obj, bcf_head, bcf_records)==0; 
    //bcf_read(bcf_obj, bcf_head, bcf_records)==0; //ignoring return value (int)
    
    /*kstring_t shared,indiv = bcf_records->indiv;
    std::cout << shared << std::endl;
    fprintf(stderr,"test %s \n",indiv.s);
    exit(0);*/

    size_t chr_pos_before;
    size_t chr_pos_after;
    size_t insert_total = 0;
    size_t del_total = 0;
      
    hts_itr_t *itr = bcf_itr_querys(idx, bcf_head,chr_names[0]);
    int record_indiv;
    size_t old_pos = 0;
    while ((record_indiv = bcf_itr_next(bcf_obj, itr, bcf_records)) == 0){
      //fprintf(stderr,"WHILE LOOP\n");
      if(strcasecmp(bcf_hdr_id2name(bcf_head, bcf_records->rid),chr_names[0])==0){
        bcf_unpack((bcf1_t*)bcf_records, BCF_UN_ALL);

        //https://github.com/samtools/htslib/blob/226c1a813bc5d0582f7e0b0bdb4b3ea9e3ee4ce4/htslib/vcf.h#L1019
        int nsamples = bcf_hdr_nsamples(bcf_head);

        //fprintf(stderr,"number of samples %d\n",nsamples);

        // gt data for each call
        int32_t ngt_arr = 0;     
        int32_t *gt_arr = NULL;
        int ngt = bcf_get_genotypes(bcf_head, bcf_records, &gt_arr, &ngt_arr);
        int max_ploidy = ngt/nsamples;

        for (int i =0; i<nsamples; i++){
          char* haplotype1;char* haplotype2;
          //iterates through all samples
          // fprintf(stderr,"--------\nTHE SAMPLE INDEX IS %d and sample name is %s\n",i,bcf_head->samples[i]);
          // match the names in the header with the input indivduals
          if(strcasecmp(bcf_head->samples[i],HeaderIndiv)==0){
            // fprintf(stderr,"--------\nTHE SAMPLE INDEX IS %d and sample name is %s and input %s\n",i,bcf_head->samples[i],HeaderIndiv);
            int32_t *ptr = gt_arr + i*max_ploidy;
            int haplotype1_int; int haplotype2_int; 
            //for (int j=0; j<max_ploidy; j++){fprintf(stderr,"the ploidy index is %d \t alelle values %d\n",j,bcf_gt_allele(ptr[j]));}
            haplotype1_int = bcf_gt_allele(ptr[0]); haplotype2_int = bcf_gt_allele(ptr[1]);

            //Extract actual allele chars;
            if(haplotype1_int == 0 && haplotype2_int == 0){
              //homozygous for reference
              //fprintf(stderr,"homozygous for reference\n");
              haplotype1 = bcf_records->d.allele[bcf_gt_allele(gt_arr[0])];
              haplotype2 = bcf_records->d.allele[bcf_gt_allele(gt_arr[0])];
            }
            else if(haplotype1_int == 1 && haplotype2_int == 1){
              //homozygous for alternative
              //fprintf(stderr,"homozygous for alternative\n");
              haplotype1 = bcf_records->d.allele[bcf_gt_allele(gt_arr[1])];
              haplotype2 = bcf_records->d.allele[bcf_gt_allele(gt_arr[1])];
            }
            else if(haplotype1_int != haplotype2_int){
              //heterozygous
              //fprintf(stderr,"heterozygous\n");
              haplotype1 = bcf_records->d.allele[bcf_gt_allele(gt_arr[0])];
              haplotype2 = bcf_records->d.allele[bcf_gt_allele(gt_arr[1])];
            }

            size_t pos = (int) bcf_records->pos;
            //fprintf(stderr,"Position %zu \t The haplotypes values are %d \t %d and characters are %c \t %c \n",pos,haplotype1_int,haplotype2_int,*haplotype1,*haplotype2);
            //Current SNP location
            //fprintf(stderr,"before alterations %c%c%c\t%c%c%c\n",genome_data1[pos-1],genome_data1[pos],genome_data1[pos+1],genome_data2[pos-1],genome_data2[pos],genome_data2[pos+1]);
            genome_data1[pos] = *haplotype1;genome_data2[pos] = *haplotype2;
            //fprintf(stderr,"after alterations %c%c%c\t%c%c%c\n",genome_data1[pos-1],genome_data1[pos],genome_data1[pos+1],genome_data2[pos-1],genome_data2[pos],genome_data2[pos+1]);
            
          }
        }
      }
      else{
        fprintf(stderr,"Reference chromosome %s has no variations in vcf file at chromosomes %s\n",chr_names[0],bcf_hdr_id2name(bcf_head, bcf_records->rid));
      }
    }
  }
  //valgrind ./ngsngs -i ../vcfdata/chr14.fa -r 10000000 -t1 4 -l 100 -seq SE -f fq -ne -q1 Test_Examples/Qual_profiles/AccFreqL150R1.txt -bcf ../vcfdata/chr14_indiv_3.bcf -v snp -indiv HG00099 -o Chr14_HG00099
  bcf_hdr_destroy(bcf_head);
  bcf_destroy(bcf_records); 
  bcf_close(bcf_obj);
  //fprintf(stderr,"Genome 1    example %c%c%c%c%c\n",genome[19000014],genome[19000015],genome[19000016], genome[19000017], genome[19000018]);
  sprintf(genome+strlen(genome_data1),"%s",genome_data1);
  strcat(genome,genome_data2);
  /*fprintf(stderr,"Haplotype 1 example %c%c%c%c%c\n",genome_data1[19000014],genome_data1[19000015],genome_data1[19000016], genome_data1[19000017], genome_data1[19000018]);
  fprintf(stderr,"Haplotype 2 example %c%c%c%c%c\n",genome_data2[19000014],genome_data2[19000015],genome_data2[19000016], genome_data2[19000017], genome_data2[19000018]);
  fprintf(stderr,"Genome II   example %c%c%c%c%c\n",genome[19000014],genome[19000015],genome[19000016], genome[19000017], genome[19000018]);
  fprintf(stderr,"Genome III  example %c%c%c%c%c\n",genome[19000014+107349540],genome[19000015+107349540],genome[19000016+107349540], genome[19000017+107349540], genome[19000018+107349540]);*/

  return genome;
}

char* full_genome_create(faidx_t *seq_ref,int chr_total,int chr_sizes[],const char *chr_names[],size_t chr_size_cumm[]){
  size_t genome_size = 0;
  chr_size_cumm[0] = 0;
  /*std::string Nstr = std::string(300, 'N');
  std::cout << Nstr << std::endl;
  const char* Ndata =  Nstr.c_str();
  std::cout << Ndata << std::endl;
  std::cout << strlen(Ndata) << std::endl;*/
  for (int i = 0; i < chr_total; i++){
    const char *chr_name = faidx_iseq(seq_ref,i);
    int chr_len = faidx_seq_len(seq_ref,chr_name);
    chr_sizes[i] = chr_len;
    chr_names[i] = chr_name;
    genome_size += chr_len;// + strlen(Ndata);
    chr_size_cumm[i+1] = genome_size;
  }
  
  char* genome = (char*) malloc(sizeof(char) * (genome_size+chr_total+1));//(strlen(Ndata)*chr_total)
  genome[0] = 0; //Init to create proper C string before strcat
  //chr_total
  for (int i = 0; i < chr_total; i++){

    const char *data = fai_fetch(seq_ref,chr_names[i],&chr_sizes[i]);
    //sprintf(&genome[strlen(genome)],data);
    //strcat(genome,data);  //Both gives conditional jump or move error
    if (data != NULL){
      //std::cout << strlen(genome) << std::endl;
      sprintf(genome+strlen(genome),"%s",data);
      //strcat(genome,Ndata);
    }
    // several of the build in functions allocates memory without freeing it again.
    free((char*)data); //Free works on const pointers, so we have to cast into a const char pointer
  }
  return genome;
}

char* partial_genome_create(faidx_t *seq_ref,int chr_total,int chr_sizes[],const char *chr_names[],size_t chr_size_cumm[]){
  
  size_t genome_size = 0;
  chr_size_cumm[0] = 0;
  for (int i = 0; i < chr_total; i++){
    int chr_len = faidx_seq_len(seq_ref,chr_names[i]);
    fprintf(stderr,"chr len %d\n",chr_len);
    chr_sizes[i] = chr_len;
    genome_size += chr_len;
    chr_size_cumm[i+1] = genome_size;
  }
  fprintf(stderr,"OUT OF FOR\n");
  char* genome = (char*) malloc(sizeof(char) * (genome_size+chr_total+1));
  genome[0] = 0; //Init to create proper C string before strcat
  //chr_total
  for (int i = 0; i < chr_total; i++){
    const char *data = fai_fetch(seq_ref,chr_names[i],&chr_sizes[i]);
    if (data != NULL){
      sprintf(genome+strlen(genome),"%s",data);  
    }
    // several of the build in functions allocates memory without freeing it again.
    free((char*)data); //Free works on const pointers, so we have to cast into a const char pointer
  }
  return genome;
}

char* full_vcf_genome_create(faidx_t *seq_ref,int chr_total,int chr_sizes[],const char *chr_names[],size_t chr_size_cumm[],const char* bcf_file,const char* VarType,const char* HeaderIndiv){
  fprintf(stderr,"BCF FILE %s\n",bcf_file);
  size_t genome_size = 0;
  chr_size_cumm[0] = 0;
  
  const char *chr_name = faidx_iseq(seq_ref,0);
  int chr_len = faidx_seq_len(seq_ref,chr_name);
  chr_sizes[0] = chr_len;chr_sizes[1] = chr_len;
  chr_names[0] = chr_name;chr_names[1] = chr_name;
  genome_size += chr_len;
  chr_size_cumm[1] = genome_size;
  genome_size += chr_len;
  chr_size_cumm[2] = genome_size;

  char* genome = (char*) malloc(sizeof(char) * (genome_size+(chr_total*2))*2);
  genome[0] = 0; //Init to create proper C string before strcat
  //chr_total
  FILE *VarInfoFile = fopen("Variant_pos_log.txt", "w");
  fprintf(VarInfoFile,"Type\tChr\tRef_Pos\tRef\tRead_Pos\tAlt\n");
  for (int i = 0; i < chr_total; i++){

    char *Chr_hapl1 = fai_fetch(seq_ref,chr_names[i],&chr_sizes[i]);
    char *Chr_hapl2 = fai_fetch(seq_ref,chr_names[i],&chr_sizes[i]);

    if (Chr_hapl1 != NULL){
      fprintf(stderr,"INSIDE NULL DATA \t WITH VARTYPE %s and chr %s\n",VarType,chr_names[i]);
      //char* HaploGenome(char genome_data[],int chr_sizes,const char* bcf_file,const char *chr_names[]){
      
      HaploGenome(genome,Chr_hapl1, Chr_hapl2,chr_sizes[i],bcf_file,&chr_names[i],VarType,VarInfoFile,HeaderIndiv); //Chr14_subset4_bcf
    }
    //std::cout<< genome[19000016] << genome[19000017] << genome[19000018] << genome[19000016+107349540] << genome[19000017+107349540] << genome[19000018+107349540] << std::endl;
    //std::cout << strlen(genome)<<std::endl;
    // several of the build in functions allocates memory without freeing it again.
    free((char*)Chr_hapl1); //Free works on const pointers, so we have to cast into a const char pointer
    free((char*)Chr_hapl2); //Free works on const pointers, so we have to cast into a const char pointer
  }
  /*fprintf(stderr,"Genome IIII example %c%c%c%c%c\n",genome[19000014],genome[19000015],genome[19000016], genome[19000017], genome[19000018]);
  fprintf(stderr,"Genome IV   example %c%c%c%c%c\n",genome[19000014+107349540],genome[19000015+107349540],genome[19000016+107349540], genome[19000017+107349540], genome[19000018+107349540]);
  fprintf(stderr,"DONE WITH DATA FOR LOOP \n");*/
  fclose(VarInfoFile);
  return genome;
}

ransampl_ws ***ReadQuality(char *ntqual, double *ErrProb, int ntcharoffset,const char *freqfile,unsigned long &readcycle){
  ransampl_ws ***dists = new ransampl_ws**[5];

  std::vector<char *> all_lines;
  gzFile gz = Z_NULL;
  assert(((gz = gzopen(freqfile,"rb")))!=Z_NULL);
  char buf[LENS];
  while(gzgets(gz,buf,LENS))
    all_lines.push_back(strdup(buf));
  gzclose(gz);

  //fprintf(stderr,"All lines: %lu\n",all_lines.size());

  unsigned long readcyclelength = (all_lines.size()-2)/5; //all_lines.size()-1

  fprintf(stderr,"\t-> Inferred read cycle lengths: %lu\n",readcyclelength);
  readcycle = readcyclelength; 
  //loop over inputdata
  int nbins = -1;
  double probs[MAXBINS];
  for(int b=0;b<5;b++){
    dists[b] = new ransampl_ws *[readcyclelength];
    for(unsigned long pos = 0 ; pos<readcyclelength;pos++){
      int at = 0;
      probs[at++] = atof(strtok(all_lines[2+b*readcyclelength+pos],"\n\t ")); //1+b*readcyclelength+pos
      char *tok = NULL;
      while(((tok=strtok(NULL,"\n\t ")))){
	      probs[at++] = atof(tok);
	      assert(at<MAXBINS);
      }
      if(nbins==-1){
	      nbins = at;
	      //fprintf(stderr,"Number of qualities/bins in inputfile: %d\n",nbins);
      }
      if(nbins!=at){
	      fprintf(stderr,"Problems, number of columns is different nbins: %d at: %d\n",nbins,at);
	      exit(0);
      }
      dists[b][pos] =  ransampl_alloc( nbins );
      ransampl_set(dists[b][pos],probs);
    }
  }
  //fprintf(stderr,"\t-> ransampl_ws done\n");

  //printf(all_lines[0]);
  int qualidx = 1;
  //extract the first token
  ntqual[0] = (char) (atoi(strtok(all_lines[0],"\n\t "))+ntcharoffset);
  char *qualtok = NULL;
  //extract the next
  while(((qualtok=strtok(NULL,"\n\t ")))){ntqual[qualidx++] = (char) (atoi(qualtok)+ntcharoffset);}
  
  int Err_idx = 1;
  ErrProb[0] = atof(strtok(all_lines[1],"\n\t"));
  char *Errtok = NULL;
  while ((Errtok=strtok (NULL,"\n\t"))){ErrProb[Err_idx++] = (double) atof(Errtok);}
  
  //strdup function allocate necessary memory to store the sourcing string implicitly, i need to free the returned string
  for (unsigned long i = 0; i < all_lines.size(); i++){free(all_lines[i]);}
  //fprintf(stderr,"\t-> Before return in READQUAL FUNC\n");
  return dists;
}