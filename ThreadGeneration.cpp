#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <math.h>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>
#include <htslib/thread_pool.h>

#include <pthread.h>

#include "mrand.h"
#include "Briggs.h"
#include "NtSubModels.h"
#include "NGSNGS_func.h"
#include "RandSampling.h"
#include "getFragmentLength.h"
#include "Sampling.h"

#define LENS 4096
#define MAXBINS 100

unsigned char nuc2intThread[255];

void* Create_se_threads(faidx_t *seq_ref,int thread_no, int seed, size_t reads,const char* OutputName,const char* Adapt_flag,const char* Adapter_1,
                        const char* Adapter_2,const char* OutputFormat,const char* SeqType,float BriggsParam[4],const char* Briggs_flag,
                        const char* Sizefile,int FixedSize,int SizeDistType, double val1, double val2,
                        int qualstringoffset,const char* QualProfile1,const char* QualProfile2, int threadwriteno,
                        const char* QualStringFlag,const char* Polynt,const char* ErrorFlag,const char* Specific_Chr[1024],const char* FastaFileName,
                        const char* MisMatchFlag,const char* SubProfile,int MisLength,int RandMacro,const char *VCFformat,char* Variant_flag,const char *VarType,
                        char CommandArray[1024],const char* version,const char* HeaderIndiv,const char* NoAlign,size_t BufferLength){
  //creating an array with the arguments to create multiple threads;
  //fprintf(stderr,"Random MacIntType %d\n",MacroRandType);
  //fprintf(stderr,"\t-> Command 3 : %s \n",CommandArray);
  int nthreads=thread_no;
  pthread_t mythreads[nthreads];

  nuc2intThread['a'] = nuc2intThread['A'] = nuc2intThread[0] = 0;
  nuc2intThread['t'] = nuc2intThread['T'] = nuc2intThread[1] = 1;
  nuc2intThread['g'] = nuc2intThread['G'] = nuc2intThread[2] = 2;
  nuc2intThread['c'] = nuc2intThread['C'] = nuc2intThread[3] = 3;
  nuc2intThread['n'] = nuc2intThread['N'] = nuc2intThread[4] = 4; 

  int chr_total = 0;
  char *genome_data;
  if (Specific_Chr[0] != NULL){
    while (Specific_Chr[chr_total]){
      std::cout << " lol " << chr_total << std::endl;
      chr_total++;
      }
  }
  else{
    chr_total = faidx_nseq(seq_ref);
  }

  const char *chr_names[chr_total];
  int chr_sizes[chr_total];
  int chr_idx_arr[chr_total];
  size_t chr_size_cumm[chr_total+1];
  /*fprintf(stderr,"Chromosome count %d\n",chr_total);
  fprintf(stderr,"DONE WITH LOOP\n");*/
  
  if (chr_total < faidx_nseq(seq_ref)){
    for (int j = 0; j < faidx_nseq(seq_ref); j++){
      for (int i = 0; i < chr_total; i++){
        if(strcasecmp(faidx_iseq(seq_ref, j),Specific_Chr[i])==0){
          chr_idx_arr[i] = j;
        }
      } 
    }
  }
  else
  {
    for (int j = 0; j < faidx_nseq(seq_ref); j++){;chr_idx_arr[j] = j;}
  }
  
  if(VCFformat != NULL && strcasecmp(Variant_flag,"bcf")==0){
    //const char* HeaderIndiv = "HG00097";
    genome_data = full_vcf_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm,VCFformat,VarType,HeaderIndiv);
  }
  else{
    if (chr_total == faidx_nseq(seq_ref)){
      genome_data = full_genome_create(seq_ref,chr_total,chr_sizes,chr_names,chr_size_cumm);
    }  
    else{
      fprintf(stderr,"Generating partial\n");
      std::cout << chr_total << " " << faidx_nseq(seq_ref) << std::endl;
      fprintf(stderr,"chr total %d \t specific chr %s\n",chr_total,Specific_Chr[0]);
      genome_data = partial_genome_create(seq_ref,chr_total-1,chr_sizes,Specific_Chr,chr_size_cumm);
      fprintf(stderr,"Generating done\n");
      for (int i = 0; i < chr_total; i++){
        chr_names[i] = Specific_Chr[i];
      }
    }
  }

  size_t genome_size = strlen(genome_data);;
  if (genome_data != NULL){
    fprintf(stderr,"\t-> Creating the large concatenated contig, with size of %lu bp\n",genome_size);
  
    Parsarg_for_Sampling_thread struct_for_threads[nthreads];

    // declare files and headers
    BGZF *bgzf_fp1 = NULL;
    BGZF *bgzf_fp2 = NULL;

    samFile *SAMout = NULL;
    sam_hdr_t *SAMHeader;
    htsFormat *fmt_hts =(htsFormat*) calloc(1,sizeof(htsFormat));
    htsThreadPool p = {NULL, 0};

    char file1[80];
    char file2[80];
    const char* fileprefix = OutputName; //"chr22_out";
    strcpy(file1,fileprefix);
    strcpy(file2,fileprefix);

    const char* suffix1;
    const char* suffix2;
    const char *mode;
    int alnformatflag = 0;
    if(strcasecmp("fa",OutputFormat)==0){
      //fprintf(stderr,"\t-> FA SE FILE\n");
      mode = "wu";
      if(strcasecmp("SE",SeqType)==0){suffix1 = ".fa";}
      else{suffix1 = "_R1.fa";suffix2 = "_R2.fa";}
    }
    else if(strcasecmp("fa.gz",OutputFormat)==0){
      mode = "wb";
      if(strcasecmp("SE",SeqType)==0){suffix1 = ".fa.gz";}
      else{suffix1 = "_R1.fa.gz";suffix2 = "_R2.fa.gz";}
    }
    else if(strcasecmp("fq",OutputFormat)==0){
      //fprintf(stderr,"\t-> FQ SE FILE\n");
      mode = "wu";
      if(strcasecmp("SE",SeqType)==0){suffix1 = ".fq";}
      else{suffix1 = "_R1.fq";suffix2 = "_R2.fq";}
    }
    else if(strcasecmp("fq.gz",OutputFormat)==0){
      mode = "w";
      if(strcasecmp("SE",SeqType)==0){suffix1 = ".fq.gz";}
      else{suffix1 = "_R1.fq.gz";suffix2 = "_R2.fq.gz";}
    }
    else if(strcasecmp("sam",OutputFormat)==0){
      //fprintf(stderr,"\t-> SAM SE FILE\n");
      mode = "ws";
      suffix1 = ".sam";
      alnformatflag++;
    }
    else if(strcasecmp("bam",OutputFormat)==0){
      //fprintf(stderr,"\t-> BAM SE FILE\n");
      mode = "wb";//"wc";
      suffix1 = ".bam"; //".cram";
      alnformatflag++;
    }
    else if(strcasecmp("cram",OutputFormat)==0){
      mode = "wc";
      suffix1 = ".cram";
      alnformatflag++;
    }
    else{fprintf(stderr,"\t-> Fileformat is currently not supported \n");}
    strcat(file1,suffix1);
    //fprintf(stderr,"\t-> Compression level for file %s and writing mode %s \n",OutputFormat,mode);

    fprintf(stderr,"\t-> File output name is %s\n",file1);
    const char* filename1 = file1;
    const char* filename2 = NULL;

    if(alnformatflag == 0){
      //fprintf(stderr,"not bam loop \n");
      int mt_cores = threadwriteno;
      int bgzf_buf = 256;
      
      bgzf_fp1 = bgzf_open(filename1,mode); //w
      //fprintf(stderr,"Number of writing cores %d\n",mt_cores);
      bgzf_mt(bgzf_fp1,mt_cores,bgzf_buf); //
      
      //fprintf(stderr,"\t-> BGZF FILE\n");
      if(strcasecmp("PE",SeqType)==0){
        strcat(file2,suffix2);
        filename2 = file2;
        bgzf_fp2 = bgzf_open(filename2,mode);
        bgzf_mt(bgzf_fp2,mt_cores,bgzf_buf);
      }
    }
    else{
      //fprintf(stderr,"Fasta input file name is %s\n",FastaFileName);
      char *ref =(char*) malloc(10 + strlen(FastaFileName) + 1);
      sprintf(ref, "reference=%s", FastaFileName);
      //char *ref =(char*) malloc(10 + strlen("Test_Examples/Mycobacterium_leprae.fa.gz") + 1);
      //sprintf(ref, "reference=%s", "Test_Examples/Mycobacterium_leprae.fa.gz");
      hts_opt_add((hts_opt **)&fmt_hts->specific,ref);
      //fprintf(stderr,"Writing mode is %s\n",mode);
      SAMout = sam_open_format(filename1, mode, fmt_hts);
      SAMHeader = sam_hdr_init();

      if(threadwriteno>0){
        if (!(p.pool = hts_tpool_init(threadwriteno))) {
          fprintf(stderr, "Error creating thread pool\n");
          exit(0);
        }
        hts_set_opt(SAMout, HTS_OPT_THREAD_POOL, &p);
      }
      // generate header
      //hts_set_threads(SAMout, 4);
      Header_func(fmt_hts,filename1,SAMout,SAMHeader,seq_ref,chr_total,chr_idx_arr,genome_size,CommandArray,version);
      free(ref);
      hts_opt_free((hts_opt *)fmt_hts->specific);
    }
    //fprintf(stderr,"\t-> AFTER OUTPUT FORMAT\n");

    //generate file array before creating the threads
    int no_elem;double* Frag_freq;int* Frag_len;
    if(SizeDistType==1){
        Frag_len = new int[LENS];Frag_freq = new double[LENS];
        ReadLengthFile(no_elem,Frag_len,Frag_freq,Sizefile);
    }
    else{no_elem = -1;}
    
    const char *freqfile_r1; //"Qual_profiles/AccFreqL150R1.txt";
    const char *freqfile_r2;
    int outputoffset = qualstringoffset;
    unsigned long readcyclelength;
    //fprintf(stderr,"\t-> FRAG ARRAY LE\n");
    ransampl_ws ***QualDist;
    char nt_qual_r1[1024];
    ransampl_ws ***QualDist2;
    char nt_qual_r2[1024];
    //fprintf(stderr,"\t-> QUAL POINTER POINTER POINTER\n");
    double ErrArray_r1[1024];
    double ErrArray_r2[1024];
    //fprintf(stderr,"\t-> BEFORE QUAL STRING IF\n");
    if(strcasecmp("true",QualStringFlag)==0){ //|| strcasecmp("bam",OutputFormat)==0
      freqfile_r1 = QualProfile1;
      QualDist = ReadQuality(nt_qual_r1,ErrArray_r1,outputoffset,freqfile_r1,readcyclelength);
      //fprintf(stderr,"\t-> CREATING QUALDIST\n");
      if(strcasecmp("PE",SeqType)==0){
        //fprintf(stderr,"\t-> PE LOOP\n");
        freqfile_r2 = QualProfile2;
        QualDist2 = ReadQuality(nt_qual_r2,ErrArray_r2,outputoffset,freqfile_r2,readcyclelength);
      }
    }
    //fprintf(stderr,"\t-> AFTER QUAL STRING IF\n");
    
    int maxsize = 20;
    char polynucleotide;
    //fprintf(stderr,"\t-> BEFORE POLY\n");
    if (Polynt != NULL && strlen(Polynt) == 1){polynucleotide = (char) Polynt[0];}
    else{polynucleotide = 'F';}
    //fprintf(stderr,"\t-> AFTER POLY\n");

    double* MisMatchFreqArray;
    int mismatchcyclelength = 0;
    if (SubProfile != NULL){
      MisMatchFreqArray = new double[LENS];
      MisMatchFreqArray = MisMatchFileArray(MisMatchFreqArray,SubProfile,mismatchcyclelength);
    }

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
      struct_for_threads[i].chr_idx_array = chr_idx_arr;
      struct_for_threads[i].chr_no = chr_total;
      struct_for_threads[i].threadseed = seed;
      struct_for_threads[i].RandMacro = RandMacro;

      struct_for_threads[i].FragLen = Frag_len;
      struct_for_threads[i].FragFreq = Frag_freq;
      struct_for_threads[i].No_Len_Val = no_elem;
      struct_for_threads[i].FixedSize = FixedSize;
      struct_for_threads[i].distparam1 = val1;
      struct_for_threads[i].distparam2 = val2;
      struct_for_threads[i].LengthType = SizeDistType;

      struct_for_threads[i].NtQual_r1 = nt_qual_r1;
      struct_for_threads[i].NtQual_r2 = nt_qual_r2;
      struct_for_threads[i].QualDist_r1 = QualDist;
      struct_for_threads[i].QualDist_r2 = QualDist2;
      struct_for_threads[i].NtErr_r1 = ErrArray_r1;
      struct_for_threads[i].NtErr_r2 = ErrArray_r2;

      struct_for_threads[i].MisMatch = MisMatchFreqArray;
      struct_for_threads[i].SubFlag = MisMatchFlag;
      struct_for_threads[i].MisLength = (int) mismatchcyclelength;
      struct_for_threads[i].readcycle = (int) readcyclelength;
      struct_for_threads[i].reads = reads;
      struct_for_threads[i].BufferLength = BufferLength;

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
      struct_for_threads[i].Adapter_2 = Adapter_2;
      struct_for_threads[i].Briggs_flag = Briggs_flag;
      struct_for_threads[i].BriggsParam = BriggsParam;
      struct_for_threads[i].OutputFormat = OutputFormat;
      struct_for_threads[i].SeqType = SeqType;
      struct_for_threads[i].QualFlag = QualStringFlag;
      struct_for_threads[i].PolyNt = polynucleotide;
      struct_for_threads[i].ErrorFlag = (char) ErrorFlag[0];
      struct_for_threads[i].NoAlign = (char) NoAlign[0];
      struct_for_threads[i].Variant_flag = Variant_flag;


      //declaring the size of the different arrays
      struct_for_threads[i].size_cumm = (size_t*)malloc(sizeof(size_t) * (struct_for_threads[i].chr_no+1));
      struct_for_threads[i].size_cumm[0] = 0;
      memcpy(struct_for_threads[i].size_cumm, chr_size_cumm, sizeof(chr_size_cumm));
      
      struct_for_threads[i].names = (char**)malloc(sizeof(char*) * struct_for_threads[i].chr_no+1);
      struct_for_threads[i].names[0] = 0;
      memcpy(struct_for_threads[i].names, chr_names, sizeof(chr_names));
    }
    
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    
    for (int i = 0; i < nthreads; i++){
      pthread_create(&mythreads[i],&attr,Sampling_threads,&struct_for_threads[i]);
    }
    
    for (int i = 0; i < nthreads; i++){  
      pthread_join(mythreads[i],NULL);
    }
        
    if(alnformatflag == 0){
      bgzf_close(bgzf_fp1);
      if(strcasecmp("PE",SeqType)==0){bgzf_close(bgzf_fp2);}
    }
    else{
      sam_hdr_destroy(SAMHeader);
      sam_close(SAMout);
    } 
    
    for(int i=0;i<nthreads;i++){
      free(struct_for_threads[i].fqresult_r1 -> s);
      free(struct_for_threads[i].list_of_reads);
      delete struct_for_threads[i].fqresult_r1;

      free(struct_for_threads[i].fqresult_r2 -> s);
      delete struct_for_threads[i].fqresult_r2;      
    }
    
    if(strcasecmp("true",QualStringFlag)==0){
      for(int b=0;b<5;b++){
        for(int pos = 0 ; pos< (int) readcyclelength;pos++){
          ransampl_free(QualDist[b][pos]);
        }
        delete[] QualDist[b];
      }
      delete[] QualDist;

      if(strcasecmp("PE",SeqType)==0){
        for(int b=0;b<5;b++){
          for(int pos = 0 ; pos< (int) readcyclelength;pos++){
            ransampl_free(QualDist2[b][pos]);
          }
          delete[] QualDist2[b];
        }
        delete[] QualDist2;
      }
    }
    
    free(fmt_hts);
    if(SizeDistType==1){
      delete[] Frag_freq;
      delete[] Frag_len;
    }
    
    if(SubProfile != NULL){delete[] MisMatchFreqArray;}
    
    free(genome_data);
    fflush(stderr);
  }
  return NULL;
}