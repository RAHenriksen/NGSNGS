
#include "fasta_sampler.h"
#include "RandSampling.h"
#include "mrand.h"
#include "NGSNGS_func.h"
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#define MAXBINS 100

// fasta_sampler_alloc(refSseq,Specific_Chr,VariantFile,VarType,HeaderIndiv);
fasta_sampler *fasta_sampler_alloc(const char *fa,const char *SpecificChr,const char *VariantFile,const char *VarType,const char *HeaderIndiv){
  
  char *chrtok = NULL;
  std::vector<char *> SubsetChr;
  if (SpecificChr != NULL){
    //fprintf(stderr,"INSIDE SPECIFIC SUBSET CHR\n");
    SubsetChr.push_back(strtok(strdup(SpecificChr),"\", \t"));
    char *chrtok = NULL;
    while(((chrtok=strtok(NULL,"\", \t"))))
      SubsetChr.push_back(strdup(chrtok));
  }

  fasta_sampler *fs = new fasta_sampler;
  fs->fai = NULL;
  assert(((fs->fai =fai_load(fa)))!=NULL);
  
  if(SubsetChr.size()==0)
    fs->nref = faidx_nseq(fs->fai);
  else if(VariantFile != NULL){
    fs->nref = SubsetChr.size()*2;
  }
  else
    fs->nref = SubsetChr.size();

  fs->seqs = new char* [fs->nref];
  fs->seqs_l = new int[fs->nref];
  fs->seqs_names = new char *[fs->nref];
  
  if(SubsetChr.size()==0)
    for(int i=0;i<fs->nref;i++){
      fs->seqs_names[i] = strdup(faidx_iseq(fs->fai,i));
      // fs->seqs_l[i] = faidx_seq_len(fs->fai,fs->seqs_names[i]);
      fs->seqs[i] = fai_fetch(fs->fai,fs->seqs_names[i],fs->seqs_l+i);
      //  fprintf(stderr,"%d)\tchr:%s\tlen:%d\n",i,fs->seqs_names[i],fs->seqs_l[i]);
    }
  else{
    if(VariantFile == NULL){
      //subset of chromosomes
      for(int i=0;i<SubsetChr.size();i++){
        fs->seqs_names[i] = strdup(SubsetChr[i]);
        fs->seqs[i] = fai_fetch(fs->fai,fs->seqs_names[i],fs->seqs_l+i);
        fprintf(stderr,"%d)\tchr:%s\tlen:%d\n",i,fs->seqs_names[i],fs->seqs_l[i]);
      }
    }
    else{
      if(SubsetChr.size()>1){
        fprintf(stderr,"Specify a single chromsome matching the chromosome from the input bcf");
      }
      else{
        htsFile *bcf_obj = bcf_open(VariantFile, "r");
        bcf_hdr_t *bcf_head = bcf_hdr_read(bcf_obj);
        hts_idx_t *bcf_idx = bcf_index_load(VariantFile);
        bcf1_t *bcf_records = bcf_init();
        
        char SeqName[1024] = {0};
        char SeqName2[1024] = {0};
        fs->seqs_names[0] = SeqName;
        fs->seqs_names[1] = SeqName2;
        
        //size_t pos = 19000017;

        for(int i=0;i<fs->nref;i++){
          //fs->seqs_names[i] = strdup(SubsetChr[0]);
          snprintf(fs->seqs_names[i],1024,"%sH%d",strdup(SubsetChr[0]),i+1);
          fs->seqs[i] = fai_fetch(fs->fai,strdup(SubsetChr[0]),fs->seqs_l+i);
          fprintf(stderr,"%d)\tchr:%s\tlen:%d\n",i,fs->seqs_names[i],fs->seqs_l[i]);
          //fprintf(stderr,"%zu POSITION \t before alterations %c%c%c\n",pos,fs->seqs[i][pos-1],fs->seqs[i][pos],fs->seqs[i][pos+1]);
          HaploType(fs->seqs[i],bcf_obj,bcf_head,bcf_idx,bcf_records,"14",i,HeaderIndiv);
          //fprintf(stderr,"%zu POSITION \t after alterations %c%c%c\n",pos,fs->seqs[i][pos-1],fs->seqs[i][pos],fs->seqs[i][pos+1]);
        }
        //fprintf(stderr,"--------------------\n");
        fprintf(stderr,"AFTER FOR \n%d)\tchr: %s\tlen:%d\n",0,fs->seqs_names[0],fs->seqs_l[0]);
        fprintf(stderr,"AFTER FOR \n%d)\tchr: %s\tlen:%d\n",1,fs->seqs_names[1],fs->seqs_l[1]);
        
        //exit(0);
        bcf_hdr_destroy(bcf_head);
        bcf_destroy(bcf_records); 
        bcf_close(bcf_obj);
      }
    }
  }
  fprintf(stderr,"number of nref %d\n",fs->nref);
  fs->ws = ransampl_alloc(fs->nref);
  double *p = new double[fs->nref];
  fs->seq_l_total =0;
  for(int i=0;i<fs->nref;i++){
    fs->seq_l_total += fs->seqs_l[i];
    fprintf(stderr,"Total length %d\t%zu\n",i,fs->seq_l_total);
  }
  for(int i=0;i<fs->nref;i++)
    p[i] = ((double) fs->seqs_l[i])/fs->seq_l_total;
  ransampl_set(fs->ws,p);
  return fs;
}

//functions returns head of chromosome, with posB, posE, chr_idx and fraglength set accordingly
char *sample(fasta_sampler *fs,mrand_t *mr,char **chromoname,int &chr_idx,int &posB,int &posE,int &fraglength){
  chr_idx = ransampl_draw2(fs->ws,mrand_pop(mr),mrand_pop(mr));
  *chromoname = fs->seqs_names[chr_idx];
  //fprintf(stderr,"CHROMOSOME NAME XX  %s \t %s and index %d\n",fs->seqs_names[0],fs->seqs_names[1],chr_idx);
  posB = abs(mrand_pop_long(mr)) % fs->seqs_l[chr_idx]+1;
  posE = posB +fraglength;
  if(posE>fs->seqs_l[chr_idx]){
    posE=fs->seqs_l[chr_idx];
    fraglength = posE-posB;
  }
  return fs->seqs[chr_idx];
}

void fasta_sampler_destroy(fasta_sampler *fs){
  fai_destroy(fs->fai);
  for(int i=0;i<fs->nref;i++){
    free(fs->seqs[i]);
    free(fs->seqs_names[i]);
  }
  delete [] fs->seqs_names;
  delete [] fs->seqs_l;
  delete [] fs->seqs;

  ransampl_free(fs->ws);
  delete fs;
}

#ifdef __WITH_MAIN__
int main(int argc,char**argv){

  const char* SubsetChr = NULL;// = "chr12,chr14,chr15";
      

  int seed = 101;
  mrand_t *mr = mrand_alloc(3,seed);
  char ref[] = "/projects/lundbeck/scratch/wql443/reference_files/Hg19/chr10_15.fa";
  fasta_sampler *fs = fasta_sampler_alloc(ref,SubsetChr);
  //fasta_sampler *fs = fasta_sampler_alloc("Test_Examples/fire.fa");
  
  char *chr; //this is an unallocated pointer to a chromosome name, eg chr1, chrMT etc
  int posB,posE;//this is the first and last position of our fragment
  char *seq;//actual sequence, this is unallocated
  int fraglength = 100;

  size_t nit =0;
  size_t ngen = 30e6;
  char seq_r1[1024] = {0};
  while(nit++<100){
    //fraglength = abs(mrand_pop_long(mr)) % 1000;
    seq = sample(fs,mr,&chr,posB,posE,fraglength);
    strncpy(seq_r1,seq,fraglength);
    fprintf(stdout,"nit:%lu\tchromo:%s\tposB:%d\tposE:%d\tfraglength:%d\texample:%s\n",nit,chr,posB,posE,fraglength,seq_r1);
  }
  fprintf(stderr,"after\n");
  return 0;
}
#endif
//g++ fasta_sampler.cpp mrand.o RandSampling.o /projects/lundbeck/scratch/wql443/WP1/htslib/libhts.a -std=c++11 -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto -D__WITH_MAIN__ -o Fatest
