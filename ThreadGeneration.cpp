#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>
#include <htslib/thread_pool.h>
#include <iostream>
#include <cmath>
#include <math.h>

#include <pthread.h>

#include "mrand.h"
#include "Briggs.h"
#include "NtSubModels.h"
#include "RandSampling.h"
#include "getFragmentLength.h"
#include "Sampling.h"
#include "sample_qscores.h"
#include "NGSNGS_cli.h"
#include "fasta_sampler.h"
#include "add_variants.h"

#define LENS 4096

void Header_func(htsFormat *fmt_hts,const char *outfile_nam,samFile *outfile,sam_hdr_t *header,fasta_sampler *fs,char CommandArray[1024],const char* version){
  // Creates a header for the bamfile. The header is initialized before the function is called //
  
  char genome_len_buf[1024];
  sam_hdr_add_line(header, "HD", "VN",version, "SO", "unsorted", NULL);
  for(int i=0;i<fs->nref;i++){
    snprintf(genome_len_buf,1024,"%d", fs->seqs_l[i]);
    
    // reference part of the header, int r variable ensures the header is added
    int r = sam_hdr_add_line(header, "SQ", "SN", fs->seqs_names[i], "LN", genome_len_buf, NULL);
    if (r < 0) { fprintf(stderr,"sam_hdr_add_line");}
   
    memset(genome_len_buf,0, sizeof(genome_len_buf));
  
  }
  // Adding PG tag specifying the command used for simulations
  sam_hdr_add_pg(header,"NGSNGS","VN",version,"CL",CommandArray,NULL);
  // saving the header to the file
  if (sam_hdr_write(outfile, header) < 0) fprintf(stderr,"writing headers to %s", outfile_nam); //outfile
}

void* ThreadInitialization(const char* refSseq,int thread_no, int seed, size_t reads,const char* OutputName,int AddAdapt,const char* Adapter_1,
                        const char* Adapter_2,outputformat_e OutputFormat,seqtype_e SeqType,float BriggsParam[4],int DoBriggs,int DoBriggsBiotin,
                        const char* Sizefile,int FixedSize,int SizeDistType, double val1, double val2,int readcycle,
                        int qualstringoffset,const char* QualProfile1,const char* QualProfile2, int threadwriteno,
                        const char* QualStringFlag,const char* Polynt,int DoSeqErr,const char* Specific_Chr,
                        int doMisMatchErr,const char* SubProfile,int MisLength,int RandMacro,const char *VariantFile,float IndelFuncParam[4],int DoIndel,
                        char CommandArray[1024],const char* version,int HeaderIndiv,int NoAlign,size_t BufferLength,const char* FileDump){
  //creating an array with the arguments to create multiple threads;

  int nthreads=thread_no;
  pthread_t *mythreads = new pthread_t[nthreads]; //pthread_t mythreads[nthreads];
  
  //allocate for reference file
  fasta_sampler *reffasta = fasta_sampler_alloc(refSseq,Specific_Chr);
  
  if(VariantFile){
    fprintf(stderr,"specific chr %s \n",Specific_Chr);
    fprintf(stderr,"Variant file %s and individual in Header %d \n",VariantFile,HeaderIndiv);
    add_variants(reffasta,VariantFile,HeaderIndiv);
  }

  //const char* file = "FileDump.txt";
  if (VariantFile && FileDump != NULL){
    fprintf(stderr,"BEFORE \n");
    char dumpfile1[80];
    const char* dumpfile1prefix = FileDump;
    const char* dumpfile1suffix = ".fa";

    strcpy(dumpfile1,dumpfile1prefix);
    strcat(dumpfile1,dumpfile1suffix);
    const char* dumpfilefull = dumpfile1;
    dump_internal(reffasta,dumpfilefull);
    fprintf(stderr,"AFTER \n");
  }
 
  //fasta_sampler_print(stderr,reffasta);
  
  fprintf(stderr,"\t-> Allocated memory for %d chromosomes/contigs/scaffolds from input reference genome\n",reffasta->nref);
  fprintf(stderr,"\t-> Chromosome name first %s and length %d and full length %zu\n",reffasta->seqs_names[0],reffasta->seqs_l[0],reffasta->seq_l_total);


  if (reffasta->seqs != NULL){
  
    Parsarg_for_Sampling_thread *struct_for_threads = new Parsarg_for_Sampling_thread[nthreads]; //Parsarg_for_Sampling_thread struct_for_threads[nthreads];

    // declare files and headers
    BGZF **bgzf_fp = (BGZF **) calloc(2,sizeof(BGZF *));

    samFile *SAMout = NULL;
    sam_hdr_t *SAMHeader = NULL;
    htsFormat *fmt_hts =(htsFormat*) calloc(1,sizeof(htsFormat));
    htsThreadPool p = {NULL, 0};

    char file1[80];
    char file2[80];
    const char* fileprefix = OutputName; //"chr22_out";
    strcpy(file1,fileprefix);
    strcpy(file2,fileprefix);

    const char* suffix1 = NULL;
    const char* suffix2 = NULL;
    const char* mode = NULL;
    int alnformatflag = 0;
    switch(OutputFormat){
    case faT:
      mode = "wu";
      if(SE==SeqType)
	      suffix1 = ".fa";
      else{
	      suffix1 = "_R1.fa";
	      suffix2 = "_R2.fa";
      }
      break;
    case fagzT:
      mode = "wb";
      if(SE== SeqType)
	      suffix1 = ".fa.gz";
      else{
        suffix1 = "_R1.fa.gz";
        suffix2 = "_R2.fa.gz";
      }
      break;
    case fqT:
      mode = "wu";
      if(SE ==SeqType)
	      suffix1 = ".fq";
      else{
        suffix1 = "_R1.fq";
        suffix2 = "_R2.fq";
      }
      break;
    case fqgzT:
      mode = "w";
      if(SE==SeqType)
	      suffix1 = ".fq.gz";
      else{
        suffix1 = "_R1.fq.gz";
        suffix2 = "_R2.fq.gz";
      }
      break;

    case samT:
      mode = "ws";
      suffix1 = ".sam";
      alnformatflag++;
      break;

    case bamT:
      mode = "wb";//"wc";
      suffix1 = ".bam"; //".cram";
      alnformatflag++;
      break;
    case cramT:
      mode = "wc";
      suffix1 = ".cram";
      alnformatflag++;
      break;
    default:
      fprintf(stderr,"\t-> Fileformat is currently not supported \n");
      break;
    }
     
    strcat(file1,suffix1);

    fprintf(stderr,"\t-> File output name is %s\n",file1);
    const char* filename1 = file1;
    const char* filename2 = NULL;

    if(alnformatflag == 0){
      //fprintf(stderr,"not bam loop \n");
      int mt_cores = threadwriteno;
      int bgzf_buf = 256;
      
      bgzf_fp[0] = bgzf_open(filename1,mode); //w
      bgzf_mt(bgzf_fp[0],mt_cores,bgzf_buf); //
      
      if(PE==SeqType){
        strcat(file2,suffix2);
        filename2 = file2;
        bgzf_fp[1] = bgzf_open(filename2,mode);
        bgzf_mt(bgzf_fp[1],mt_cores,bgzf_buf);
      }
    }
    else{
      //strlen(".fastq.gz") = longest suffix for fasta file + 2 in case of Null string terminators
      char *ref =(char*) malloc(strlen(".fasta.gz") + strlen(refSseq) + 2);
      sprintf(ref, "reference=%s", refSseq);
      
      // Save reference file name for header creation of the sam output
      //  hts_opt_add((hts_opt **)&fmt_hts->specific,ref);
      SAMout = sam_open_format(filename1, mode, fmt_hts);
      SAMHeader = sam_hdr_init();

      if(threadwriteno>0){
        if (!(p.pool = hts_tpool_init(threadwriteno))) {
          fprintf(stderr, "Error creating thread pool\n");
          exit(0);
        }
        hts_set_opt(SAMout, HTS_OPT_THREAD_POOL, &p);
      }
      hts_set_opt(SAMout, CRAM_OPT_REFERENCE, refSseq);
      // generate header
      Header_func(fmt_hts,filename1,SAMout,SAMHeader,reffasta,CommandArray,version);
      //void Header_func( int chr_total,char* chr_names[],int arra[],char CommandArray[1024],const char* version){

      free(ref);
      // hts_opt_free((hts_opt *)fmt_hts->specific);
    }

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
    //fprintf(stderr,"\t-> FRAG ARRAY LE\n");
    ransampl_ws ***QualDist = NULL;
    char nt_qual_r1[1024];
    ransampl_ws ***QualDist2 = NULL;
    char nt_qual_r2[1024];
    //fprintf(stderr,"\t-> QUAL POINTER POINTER POINTER\n");
    double ErrArray_r1[1024];
    double ErrArray_r2[1024];
    //fprintf(stderr,"\t-> BEFORE QUAL STRING IF\n");
    if(strcasecmp("true",QualStringFlag)==0){ //|| strcasecmp("bam",OutputFormat)==0
      freqfile_r1 = QualProfile1;
      QualDist = ReadQuality(nt_qual_r1,ErrArray_r1,outputoffset,freqfile_r1);
      //fprintf(stderr,"\t-> CREATING QUALDIST\n");
      if(PE==SeqType){
        //fprintf(stderr,"\t-> PE LOOP\n");
        freqfile_r2 = QualProfile2;
        QualDist2 = ReadQuality(nt_qual_r2,ErrArray_r2,outputoffset,freqfile_r2);
      }
    }
    //fprintf(stderr,"\t-> AFTER QUAL STRING IF\n");
    
    int maxsize = 20;
    char polynucleotide;

    //std::string Polystr = std::string(1000, (char) Polynt[0]);
    //fprintf(stderr,"AGAGAGA%s\n",(const char*) Polystr.substr(0, 8).c_str());
    if (Polynt != NULL && strlen(Polynt) == 1){polynucleotide = (char) Polynt[0];}
    else{polynucleotide = 'F';}

    //generating mismatch matrix to parse for each string
    double* MisMatchFreqArray = new double[LENS]; //moved new up here
    int mismatchcyclelength = 0;
    if (SubProfile != NULL){
      MisMatchFreqArray = MisMatchFileArray(MisMatchFreqArray,SubProfile,mismatchcyclelength);
    }

    for (int i = 0; i < nthreads; i++){
      struct_for_threads[i].reffasta = reffasta;

      // The output format, output files, and structural elements for SAM outputs
      struct_for_threads[i].OutputFormat = OutputFormat;
      struct_for_threads[i].SeqType = SeqType;
      struct_for_threads[i].bgzf_fp = bgzf_fp;
      struct_for_threads[i].SAMout = SAMout;
      struct_for_threads[i].SAMHeader = SAMHeader;
      struct_for_threads[i].LengthData = 0;
      struct_for_threads[i].MaximumLength = maxsize;
      struct_for_threads[i].list_of_reads = (bam1_t**) malloc(sizeof(bam1_t)*maxsize); // need to free this space
      for(int j=0; j<maxsize;j++){struct_for_threads[i].list_of_reads[j]=bam_init1();} // but also destroy the bam_init1 objects    

      // Thread generation and sampling specific information
      struct_for_threads[i].threadno = i;
      struct_for_threads[i].totalThreads = nthreads;
      struct_for_threads[i].threadseed = seed;
      struct_for_threads[i].rng_type = RandMacro;

      // Sequence alteration models
      // 1) nucleotide quality score and sequencing errors,  
      struct_for_threads[i].QualFlag = QualStringFlag;
      struct_for_threads[i].DoSeqErr = DoSeqErr;
      struct_for_threads[i].NtQual_r1 = nt_qual_r1;
      struct_for_threads[i].NtQual_r2 = nt_qual_r2;
      struct_for_threads[i].QualDist_r1 = QualDist;
      struct_for_threads[i].QualDist_r2 = QualDist2;
      struct_for_threads[i].NtErr_r1 = ErrArray_r1;
      struct_for_threads[i].NtErr_r2 = ErrArray_r2;
      struct_for_threads[i].maxreadlength = (int) readcycle;
      struct_for_threads[i].IndelFuncParam = IndelFuncParam;
      struct_for_threads[i].DoIndel = DoIndel;
      
      // 2) briggs model
      struct_for_threads[i].MisMatch = MisMatchFreqArray;
      struct_for_threads[i].doMisMatchErr = doMisMatchErr;
      struct_for_threads[i].MisLength = (int) mismatchcyclelength;

      // 3) misincorporation matrix
      struct_for_threads[i].DoBriggs = DoBriggs;
      struct_for_threads[i].DoBriggsBiotin = DoBriggsBiotin;
      struct_for_threads[i].BriggsParam = BriggsParam;
      
      // Fragment lengths 
      struct_for_threads[i].FragLen = Frag_len;
      struct_for_threads[i].FragFreq = Frag_freq;
      struct_for_threads[i].No_Len_Val = no_elem;
      struct_for_threads[i].FixedSize = FixedSize;
      struct_for_threads[i].distparam1 = val1;
      struct_for_threads[i].distparam2 = val2;
      struct_for_threads[i].LengthType = SizeDistType;

      // Sequence output specific
      struct_for_threads[i].BufferLength = BufferLength;

      // Additional information for sequence reads
      struct_for_threads[i].AddAdapt = AddAdapt;
      struct_for_threads[i].Adapter_1 = Adapter_1;
      struct_for_threads[i].Adapter_2 = Adapter_2;
      struct_for_threads[i].PolyNt = polynucleotide;
      struct_for_threads[i].NoAlign = NoAlign;
    }

    size_t ThreadReads = (size_t) floor( reads / (double) thread_no);
    for (int i = 0; i < nthreads-1; i++){
      struct_for_threads[i].reads = ThreadReads;
    }
    struct_for_threads[nthreads-1].reads = reads - (ThreadReads*(nthreads-1));
    
    //for (int i = 0; i < nthreads; i++){fprintf(stderr,"The number of threads reads for threads %zu \n",struct_for_threads[i].reads);}    
    //fprintf(stderr,"The number of threads %d and the number of reads %zu\n",thread_no,reads);

    pthread_attr_t attr;
    pthread_attr_init(&attr);
    if(nthreads==1){
      Sampling_threads(struct_for_threads);
    }
    else{
      for (int i = 0; i < nthreads; i++){
	    pthread_create(&mythreads[i],&attr,Sampling_threads,&struct_for_threads[i]);
      }
      
      for (int i = 0; i < nthreads; i++){  
	    pthread_join(mythreads[i],NULL);
      }
    }
    if(bgzf_fp[0]!=NULL)
      bgzf_close(bgzf_fp[0]);
     if(bgzf_fp[1]!=NULL)
      bgzf_close(bgzf_fp[1]);
     free(bgzf_fp); //free the calloc
     
     if(SAMHeader)
      sam_hdr_destroy(SAMHeader);
     if(SAMout)
       sam_close(SAMout);
     if (p.pool)
       hts_tpool_destroy(p.pool);
    
    fasta_sampler_destroy(reffasta);

    for(int i=0;i<nthreads;i++)
      free(struct_for_threads[i].list_of_reads);

    
    delete[] mythreads; //pthread_t *mythreads = new pthread_t[nthreads]; 

    if(strcasecmp("true",QualStringFlag)==0){
      for(int base=0;base<5;base++){
        for(int pos = 0 ; pos< (int) readcycle;pos++){
          ransampl_free(QualDist[base][pos]);
        }
        delete[] QualDist[base];
      }
      delete[] QualDist;

      if(PE==SeqType){
        for(int base=0;base<5;base++){
          for(int pos = 0 ; pos< (int) readcycle;pos++){
            ransampl_free(QualDist2[base][pos]);
          }
          delete[] QualDist2[base];
        }
        delete[] QualDist2;
      }
    }
    
    free(fmt_hts);
    if(SizeDistType==1){
      delete[] Frag_freq;
      delete[] Frag_len;
    }

    delete[] struct_for_threads; //Parsarg_for_Sampling_thread *struct_for_threads = new Parsarg_for_Sampling_thread[nthreads];

    delete[] MisMatchFreqArray; //if(SubProfile != NULL){delete[] MisMatchFreqArray;}
    
    fflush(stderr);
  }
  return NULL;
}

