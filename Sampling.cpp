#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>

#include <pthread.h>

#include "NGSNGS_func.h"

pthread_mutex_t Fq_write_mutex;

// ---------------------- SINGLE-END ---------------------- //
struct Parsarg_for_Fafq_se_thread{
  kstring_t *fqresult_r1;
  kstring_t *fqresult_r2;
  char *genome; // The actual concatenated genome
  int chr_no;
  int threadno;
  int *size_cumm;
  const char **names;

  int* FragLen;
  double* FragFreq;
  int No_Len_Val;

  double* Qualfreq;
  int threadseed;
  size_t reads;
  
  BGZF *bgzf_fp1;
  BGZF *bgzf_fp2;

  samFile *SAMout;
  sam_hdr_t *SAMHeader;
  bam1_t **list_of_reads;
  int l;
  int m;

  const char* Adapter_flag;
  const char* Adapter_1;
  const char* OutputFormat;
  const char* SeqType;
};
      
void* Fafq_thread_se_run(void *arg){
  //casting my struct as arguments for the thread creation
  Parsarg_for_Fafq_se_thread *struct_obj = (Parsarg_for_Fafq_se_thread*) arg;

  fprintf(stderr,"%s",struct_obj->SeqType);
  // creating random objects for all distributions.
  unsigned int loc_seed = struct_obj->threadseed+struct_obj->threadno;
  
  size_t genome_len = strlen(struct_obj->genome);

  //coverage method2
  char seq_r1[1024] = {0};
  char seq_r1_mod[1024] = {0};
  char read[1024] = {0};
  char readadapt[1024] = {0};

  char seq_r2[1024] = {0};
  char seq_r2_mod[1024] = {0};
  char read2[1024] = {0};
  char readadapt2[1024] = {0};
  
  // for the coverage examples
  int reads = struct_obj -> reads;

  //float cov_current = 0;
  size_t rand_start;
  //int nread = 0;

  char qual[1024] = "\0"; // {0};
  //char *qual = (char*) malloc(sizeof(char) * (151));
  //int D_i = 0;
  int localread = 0;
  int iter = 0;
  int current_reads_atom = 0;

  while (current_reads_atom < reads){
    double rand_val = ((double) rand_r(&loc_seed)/ RAND_MAX);
    double rand_val2 = rand_val*RAND_MAX;
    double rand_val3 = myrand((unsigned int) rand_val2); //((double) rand_r(&test)/ RAND_MAX);// 

    rand_start = rand_val3 * (genome_len-300); //genome_len-100000;

    int lengthbin = BinarySearch_fraglength(struct_obj->FragFreq,0, struct_obj->No_Len_Val - 1, rand_val);

    int fraglength = struct_obj->FragLen[lengthbin];//75; //struct_obj->FragLen[lengthbin];
    
    int chr_idx = 0;
    while (rand_start > struct_obj->size_cumm[chr_idx+1]){chr_idx++;}

    if (fraglength > 150){strncpy(seq_r1,struct_obj->genome+rand_start-1,150);}   // case 1
    else {strncpy(seq_r1,struct_obj->genome+rand_start-1,fraglength);}  // case 2
    
    if(strcasecmp("PE",struct_obj->SeqType)==0){
      if (fraglength > 150){strncpy(seq_r2,struct_obj->genome+rand_start+fraglength-1-150,150);}
      else {strncpy(seq_r2,struct_obj->genome+rand_start-1,fraglength);}  // case 2
    }

    int rand_id = rand_val * fraglength-1; //100

    //removes reads with NNN
    char * pch;
    char * pch2;
    pch = strchr(seq_r1,'N');
    pch2 = strrchr(seq_r1,'N');
    //if (pch != NULL){continue;}
    int seqlen = strlen(seq_r1);

    // FOR BAM FILES WE NEED THE FOLLOWING
    size_t n_cigar; uint32_t cigar_bitstring; const uint32_t *cigar;
    
    if ((int )(pch-seq_r1+1) == 1 && (int)(pch2-seq_r1+1)  == seqlen){memset(seq_r1, 0, sizeof seq_r1);}
    else{
      //memset(qual, '\0', seqlen);

      //for (int j = 0; j < fraglength; j++){D_total += 1;}
      //std::time(nullptr)
      SimBriggsModel(seq_r1, seq_r1_mod, fraglength, 0.024, 0.36, 0.68, 0.0097,loc_seed);
      
      int strand = (int) rand_r(&loc_seed)%2;//1;//rand() % 2;
      
      if (strcasecmp("SE",struct_obj->SeqType)==0 && strand == 0){
        DNA_complement(seq_r1);
        reverseChar(seq_r1);
      }

      if (strcasecmp("PE",struct_obj->SeqType)==0 && strand == 0){
        DNA_complement(seq_r1);
      }
      else if (strcasecmp("PE",struct_obj->SeqType)==0 && strand == 1){
        DNA_complement(seq_r2);
        reverseChar(seq_r2);
      }

      char READ_ID[1024]; int read_id_length;
      read_id_length = sprintf(READ_ID,"T%d_RID%d_S%d_%s:%d-%d_length:%d", struct_obj->threadno, rand_id,strand,
        struct_obj->names[chr_idx],rand_start-struct_obj->size_cumm[chr_idx],rand_start+fraglength-1-struct_obj->size_cumm[chr_idx],
        fraglength);
      
      if(strcasecmp(struct_obj->Adapter_flag,"true")==0){
        strcpy(read, seq_r1);
        strcat(read,struct_obj->Adapter_1);
        strncpy(readadapt, read, 150);

        if (strcasecmp("PE",struct_obj->SeqType)==0){
          strcpy(read2, seq_r2);
          strcat(read2,struct_obj->Adapter_1);
          strncpy(readadapt2, read2, 150);
        }

        if (strcasecmp(struct_obj -> OutputFormat,"fa")==0|| strcasecmp(struct_obj -> OutputFormat,"fa.gz")==0){
          ksprintf(struct_obj->fqresult_r1,">%s\n%s\n",READ_ID,readadapt);
          if (strcasecmp("PE",struct_obj->SeqType)==0){ksprintf(struct_obj->fqresult_r2,">%s\n%s\n",READ_ID,readadapt2);}
        }
        if (strcasecmp(struct_obj -> OutputFormat,"fq")==0|| strcasecmp(struct_obj -> OutputFormat,"fq.gz")==0){
          Read_Qual_new(readadapt,qual,loc_seed,struct_obj->Qualfreq,33);
          ksprintf(struct_obj->fqresult_r1,"@%s\n%s\n+\n%s\n",READ_ID,readadapt,qual);
          if (strcasecmp("PE",struct_obj->SeqType)==0){ksprintf(struct_obj->fqresult_r2,"@%s\n%s\n+\n%s\n",READ_ID,readadapt2,qual);}
        }
        if (struct_obj->SAMout){
          Read_Qual_new(readadapt,qual,loc_seed,struct_obj->Qualfreq,0);        
          ksprintf(struct_obj->fqresult_r1,"%s",readadapt);
          if (strcasecmp("PE",struct_obj->SeqType)==0){ksprintf(struct_obj->fqresult_r2,"%s",readadapt2);}
        }
      }
      else{
        if (strcasecmp(struct_obj -> OutputFormat,"fa")==0|| strcasecmp(struct_obj -> OutputFormat,"fa.gz")==0){
          ksprintf(struct_obj->fqresult_r1,">%s\n%s\n",READ_ID,seq_r1);
          if (strcasecmp("PE",struct_obj->SeqType)==0){ksprintf(struct_obj->fqresult_r2,">%s\n%s\n",READ_ID,seq_r2);}
        }
        if (strcasecmp(struct_obj -> OutputFormat,"fq")==0|| strcasecmp(struct_obj -> OutputFormat,"fq.gz")==0){
          Read_Qual_new(seq_r1,qual,loc_seed,struct_obj->Qualfreq,33);
          ksprintf(struct_obj->fqresult_r1,"@%s\n%s\n+\n%s\n",READ_ID,seq_r1,qual);
          if (strcasecmp("PE",struct_obj->SeqType)==0){ksprintf(struct_obj->fqresult_r2,"@%s\n%s\n+\n%s\n",READ_ID,seq_r2,qual);}
        }
        if (struct_obj->SAMout){
          Read_Qual_new(seq_r1,qual,loc_seed,struct_obj->Qualfreq,0);
          ksprintf(struct_obj->fqresult_r1,"%s",seq_r1);
          if (strcasecmp("PE",struct_obj->SeqType)==0){ksprintf(struct_obj->fqresult_r2,"%s",seq_r2);}
        }
      }
      if (struct_obj->bgzf_fp1){
        if (struct_obj->fqresult_r1->l > 30000000){
          //fprintf(stderr,"\t Buffer mutex with thread no %d\n", struct_obj->threadno);fflush(stderr);
          pthread_mutex_lock(&Fq_write_mutex);
          // bgzf_write(struct_obj->bgzf,struct_obj->fqresult_r1->s,struct_obj->fqresult_r1->l);
          assert(bgzf_write(struct_obj->bgzf_fp1,struct_obj->fqresult_r1->s,struct_obj->fqresult_r1->l)!=0);
          if (strcasecmp("PE",struct_obj->SeqType)==0){assert(bgzf_write(struct_obj->bgzf_fp2,struct_obj->fqresult_r2->s,struct_obj->fqresult_r2->l)!=0);}
          pthread_mutex_unlock(&Fq_write_mutex);
          struct_obj->fqresult_r1->l =0;
          struct_obj->fqresult_r2->l =0;
        }
      }
      if (struct_obj->SAMout){
        size_t l_aux = 0; uint8_t mapq = 60;
        hts_pos_t min_beg; //max_end, insert;
        cigar_bitstring = bam_cigar_gen(strlen(struct_obj->fqresult_r1->s), BAM_CMATCH); 
        n_cigar = 1; // Number of cigar operations, 1 since we only have matches
        uint32_t cigar_arr[] = {cigar_bitstring}; //converting uint32_t {aka unsigned int} to const uint32_t* 
        cigar = cigar_arr;
        min_beg = rand_start-struct_obj->size_cumm[chr_idx] - 1;
        uint16_t flag;
        if (strand == 0){flag = 16;}
        else{flag = 0;}

        bam_set1(struct_obj->list_of_reads[struct_obj->l++],read_id_length,READ_ID,flag,chr_idx,min_beg,mapq,n_cigar,cigar,-1,-1,0,seqlen,struct_obj->fqresult_r1->s,qual,l_aux);
        
        if (struct_obj->l < struct_obj->m){   
          pthread_mutex_lock(&Fq_write_mutex);
          for (int k = 0; k < struct_obj->l; k++){
            assert(sam_write1(struct_obj->SAMout,struct_obj->SAMHeader,struct_obj->list_of_reads[k]) != 1);
            //sam_write1(struct_obj->SAMout,struct_obj->SAMHeader,struct_obj->list_of_reads[k]);
          }
          //fprintf(stderr,"\n sam_write works\n");
          pthread_mutex_unlock(&Fq_write_mutex);
          struct_obj->l = 0;

          //fprintf(stderr,"%d\n",struct_obj->l);
          //sam_write1(struct_obj->SAMout,struct_obj->SAMHeader,struct_obj->list_of_reads[struct_obj->l]);
        }
        struct_obj->fqresult_r1->l =0;
        struct_obj->fqresult_r2->l =0;
      }

      memset(qual, 0, sizeof qual);  
      memset(seq_r1, 0, sizeof seq_r1);
      memset(seq_r1_mod, 0, sizeof seq_r1_mod);

      memset(seq_r2, 0, sizeof seq_r2);
      memset(seq_r2_mod, 0, sizeof seq_r2_mod);
      
      chr_idx = 0;
      //fprintf(stderr,"start %d, fraglength %d\n",rand_start,fraglength);
      iter++;
      localread++;
      current_reads_atom++;
    }
  }
  if (struct_obj->bgzf_fp1){
    if (struct_obj->fqresult_r1->l > 0){
      pthread_mutex_lock(&Fq_write_mutex);
      assert(bgzf_write(struct_obj->bgzf_fp1,struct_obj->fqresult_r1->s,struct_obj->fqresult_r1->l)!=0);
      if (strcasecmp("PE",struct_obj->SeqType)==0){assert(bgzf_write(struct_obj->bgzf_fp2,struct_obj->fqresult_r2->s,struct_obj->fqresult_r2->l)!=0);}
      pthread_mutex_unlock(&Fq_write_mutex);
      struct_obj->fqresult_r1->l =0;
      struct_obj->fqresult_r2->l =0;
    } 
  }
  //Freeing allocated memory
  free(struct_obj->size_cumm);
  free(struct_obj->names);
  //bam_destroy1(struct_obj->list_of_reads[0]);
  for(int j=0; j<struct_obj->m;j++){bam_destroy1(struct_obj->list_of_reads[j]);}

  fprintf(stderr,"\t number of reads generated by this thread %d \n",localread);

  return NULL;
}

void* Create_se_threads(faidx_t *seq_ref,int thread_no, int seed, int reads,const char* OutputName,const char* Adapt_flag,const char* Adapter_1,const char* OutputFormat,const char* SeqType){
  //creating an array with the arguments to create multiple threads;
  int nthreads=thread_no;
  pthread_t mythreads[nthreads];
  
  int chr_total = faidx_nseq(seq_ref);
  const char *chr_names[chr_total];
  int chr_sizes[chr_total];
  int chr_size_cumm[chr_total+1];
  char *genome_data = full_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm);
  size_t genome_size = strlen(genome_data);

  if (genome_data != NULL){
    fprintf(stderr,"\t-> Full genome function run!\n");
    fprintf(stderr,"\t-> Full genome size %lu \n",genome_size);
  
    //std::cout << " genome length " << genome_len << std::endl;
    Parsarg_for_Fafq_se_thread struct_for_threads[nthreads];

    // declare files and headers
    BGZF *bgzf_fp1 = NULL;
    BGZF *bgzf_fp2 = NULL;

    samFile *SAMout = NULL;
    sam_hdr_t *SAMHeader;
    htsFormat *fmt_hts =(htsFormat*) calloc(1,sizeof(htsFormat));

    char file1[80];
    char file2[80];
    const char* fileprefix = OutputName; //"chr22_out";
    strcpy(file1,fileprefix);
    strcpy(file2,fileprefix);

    const char* suffix1;
    const char* suffix2;

    const char *mode;
    if(strcasecmp("fa",OutputFormat)==0){
      mode = "wu";
      if(strcasecmp("SE",SeqType)==0){suffix1 = ".fa";}
      else{suffix1 = "_r1.fa";suffix2 = "_r2.fa";}
    }
    else if(strcasecmp("fa.gz",OutputFormat)==0){
      mode = "wb";
      if(strcasecmp("SE",SeqType)==0){suffix1 = ".fa.gz";}
      else{suffix1 = "_r1.fa.gz";suffix2 = "_r2.fa.gz";}
    }
    else if(strcasecmp("fq",OutputFormat)==0){
      mode = "wu";
      if(strcasecmp("SE",SeqType)==0){suffix1 = ".fq";}
      else{suffix1 = "_r1.fq";suffix2 = "_r2.fq";}
    }
    else if(strcasecmp("fq.gz",OutputFormat)==0){
      mode = "wb";
      if(strcasecmp("SE",SeqType)==0){suffix1 = ".fq.gz";}
      else{suffix1 = "_r1.fq.gz";suffix2 = "_r2.fq.gz";}
    }
    else if(strcasecmp("bam",OutputFormat)==0){
      mode = "wb";
      suffix1 = ".bam";
    }
    else{fprintf(stderr,"\t-> Fileformat is currently not supported \n");}
    strcat(file1,suffix1);

    fprintf(stderr,"\t-> File output name is %s\n",file1);
    const char* filename1 = file1;
    const char* filename2 = NULL;

    if(strcasecmp("bam",OutputFormat)!=0){
      int mt_cores = 1;
      int bgzf_buf = 256;
      
      bgzf_fp1 = bgzf_open(filename1,mode);
      bgzf_mt(bgzf_fp1,mt_cores,bgzf_buf);

      if(strcasecmp("PE",SeqType)==0){
        strcat(file2,suffix2);
        filename2 = file2;
        bgzf_fp2 = bgzf_open(filename2,mode);
        bgzf_mt(bgzf_fp2,mt_cores,bgzf_buf);
      }

      //fprintf(stderr,"\t-> Number of cores for bgzf_mt: %d\n",mt_cores); 
    }
    else{
      SAMout = sam_open_format(filename1, mode, fmt_hts);
      SAMHeader = sam_hdr_init();
      Header_func(fmt_hts,filename1,SAMout,SAMHeader,seq_ref,chr_total,genome_size);
    }
  
    // READ QUAL ARRAY
    double* Qual_freq_array = new double[6000];
    Qual_freq_array = Qual_array(Qual_freq_array,"Qual_profiles/Acc_freq1.txt");
  
    // FRAGMENT LENGTH CREATING ARRAY
    int* Frag_len = new int[4096];
    double* Frag_freq = new double[4096];
    int number;

    FragArray(number,Frag_len,Frag_freq,"Size_dist/Size_dist_sampling.txt");

    int maxsize = 5;
    //initialzie values that should be used for each thread

    for (int i = 0; i < nthreads; i++){
      struct_for_threads[i].fqresult_r1 =new kstring_t;
      struct_for_threads[i].fqresult_r1 -> l = 0;
      struct_for_threads[i].fqresult_r1 -> m = 0;
      struct_for_threads[i].fqresult_r1 -> s = NULL;

      struct_for_threads[i].fqresult_r2 =new kstring_t;
      struct_for_threads[i].fqresult_r2 -> l = 0;
      struct_for_threads[i].fqresult_r2 -> m = 0;
      struct_for_threads[i].fqresult_r2 -> s = NULL;

      struct_for_threads[i].threadno = i;
      struct_for_threads[i].genome = genome_data;
      struct_for_threads[i].chr_no = chr_total;
      struct_for_threads[i].threadseed = seed;

      struct_for_threads[i].FragLen =Frag_len;
      struct_for_threads[i].FragFreq = Frag_freq;
      struct_for_threads[i].No_Len_Val = number; 

      struct_for_threads[i].Qualfreq = Qual_freq_array;
      struct_for_threads[i].reads = reads;
      
      struct_for_threads[i].bgzf_fp1 = bgzf_fp1;
      struct_for_threads[i].bgzf_fp2 = bgzf_fp2;
      struct_for_threads[i].SAMout = SAMout;
      struct_for_threads[i].SAMHeader = SAMHeader;
      struct_for_threads[i].l = 0;
      struct_for_threads[i].m = maxsize;
      struct_for_threads[i].list_of_reads = (bam1_t**) malloc(sizeof(bam1_t)*maxsize); // need to free this space

      for(int j=0; j<maxsize;j++){struct_for_threads[i].list_of_reads[j]=bam_init1();} // but also destroy the bam_init1 objects    

      struct_for_threads[i].Adapter_flag = Adapt_flag;
      struct_for_threads[i].Adapter_1 = Adapter_1;
      struct_for_threads[i].OutputFormat = OutputFormat;
      struct_for_threads[i].SeqType = SeqType;

      //declaring the size of the different arrays
      struct_for_threads[i].size_cumm = (int*)malloc(sizeof(int) * (struct_for_threads[i].chr_no+1));
      //fprintf(stderr,"struct_for_threads[i]: %p\n",struct_for_threads[i]);
      struct_for_threads[i].size_cumm[0] = 0;
      memcpy(struct_for_threads[i].size_cumm, chr_size_cumm, sizeof(chr_size_cumm));
      
      struct_for_threads[i].names = (const char**)malloc(sizeof(const char*) * struct_for_threads[i].chr_no+1);
      struct_for_threads[i].names[0] = 0;
      memcpy(struct_for_threads[i].names, chr_names, sizeof(chr_names));
    }
    
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    
    //fprintf(stderr,"Creating a bunch of threads\n"); 
    for (int i = 0; i < nthreads; i++){
      pthread_create(&mythreads[i],&attr,Fafq_thread_se_run,&struct_for_threads[i]);
    }
    // fprintf(stderr,"Done Creating a bunch of threads\n");
    

    for (int i = 0; i < nthreads; i++)
    {  
      //fprintf(stderr,"joing threads\n");fflush(stderr);
      // PRINT 
      pthread_join(mythreads[i],NULL);
      //fprintf(stderr, "\t[ANDET STED] walltime used for join =  %.2f sec\n", (float)(time(NULL) - t3));  
    }

    if(strcasecmp("bam",OutputFormat)!=0){
      bgzf_close(bgzf_fp1);
      if(strcasecmp("PE",SeqType)==0){bgzf_close(bgzf_fp2);}
    }
    else{
      sam_hdr_destroy(SAMHeader);
      sam_close(SAMout);
    } 
    
    for(int i=0;i<nthreads;i++){
      free(struct_for_threads[i].fqresult_r1 -> s);//4 errors from 4 contexts (suppressed: 0 eventhough that delete goes with new and not free
      free(struct_for_threads[i].list_of_reads);
      delete struct_for_threads[i].fqresult_r1;

      free(struct_for_threads[i].fqresult_r2 -> s);//4 errors from 4 contexts (suppressed: 0 eventhough that delete goes with new and not free
      delete struct_for_threads[i].fqresult_r2;      
    }
    
    //struct_for_threads[i].list_of_reads = (bam1_t**) malloc(sizeof(bam1_t)*maxsize); // need to free this space
    //for(int j=0; j<maxsize;j++){struct_for_threads[i].list_of_reads[j]=bam_init1();} // but also destroy the bam_init1 objects   

    free(fmt_hts);
    delete[] Frag_freq;
    delete[] Frag_len;
    delete[] Qual_freq_array;
    free(genome_data);
    fflush(stderr);
  }
  return NULL;
}