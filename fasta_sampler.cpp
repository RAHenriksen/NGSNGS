#include "fasta_sampler.h"
#include "RandSampling.h"
#include "mrand.h"
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <htslib/bgzf.h>

#define MAXBINS 100


void fasta_sampler_setprobs(fasta_sampler *fs){
  //add sampling operation
  if(fs->ws)
    ransampl_free(fs->ws);
  fs->ws = ransampl_alloc(fs->nref);
  double *p = new double[fs->nref];
  fs->seq_l_total =0;
  for(int i=0;i<fs->nref;i++){
    fs->seq_l_total += fs->seqs_l[i];
  }
  for(int i=0;i<fs->nref;i++)
    p[i] = ((double) fs->seqs_l[i])/fs->seq_l_total;
  
  ransampl_set(fs->ws,p);
  delete [] p;
}
void fasta_sampler_print(FILE *fp,fasta_sampler *fs){
  fprintf(fp,"--------------\n[%s]\tfai: %p nref: %d\n",__FUNCTION__,fs->fai,fs->nref);
  for(int i=0;i<fs->nref;i++)
    fprintf(fp,"[%s]\tidx:%d) name:%s example:%.15s length: %d realidx: %d\n",__FUNCTION__,i,fs->seqs_names[i],fs->seqs[i],fs->seqs_l[i],fs->realnameidx[i]);

  for(char2int::iterator it=fs->char2idx.begin();it!=fs->char2idx.end();it++)
    fprintf(fp,"[%s]\tkey: %s val:%d\n",__FUNCTION__,it->first,it->second);
  fprintf(fp,"------------\n");
}

fasta_sampler *fasta_sampler_alloc(const char *fa,const char *SpecificChr){
  
  std::vector<char *> SubsetChr;
  if (SpecificChr != NULL){
    SubsetChr.push_back(strtok(strdup(SpecificChr),"\", \t"));
    char *chrtok = NULL;
    while(((chrtok=strtok(NULL,"\", \t"))))
      SubsetChr.push_back(strdup(chrtok));
    fprintf(stderr,"\t-> Number of entries in chromosome list: %lu\n",SubsetChr.size());
  }

  fasta_sampler *fs = new fasta_sampler;
  fs->ws = NULL;
  fs->fai = NULL;
  assert(((fs->fai =fai_load(fa)))!=NULL);
  
  if(SubsetChr.size()==0)
    fs->nref = faidx_nseq(fs->fai);
  else
    fs->nref = SubsetChr.size();

  fs->seqs = new char* [fs->nref];
  fs->seqs_l = new int[fs->nref];
  fs->seqs_names = new char *[fs->nref];
  fs->realnameidx = new int[fs->nref];
  if(SubsetChr.size()==0)
    for(int i=0;i<fs->nref;i++){
      fs->seqs_names[i] = strdup(faidx_iseq(fs->fai,i));
      fs->char2idx[fs->seqs_names[i]] = i;
      fs->seqs[i] = fai_fetch(fs->fai,fs->seqs_names[i],fs->seqs_l+i);
    }
  else{
    int at = 0;
    for(int i=0;i<(int)SubsetChr.size();i++){
      if( faidx_has_seq(fs->fai, SubsetChr[i])){
        fs->seqs_names[at] = strdup(SubsetChr[i]);
        fs->seqs[at] = fai_fetch(fs->fai,fs->seqs_names[i],fs->seqs_l+at);
        fs->char2idx[fs->seqs_names[i]] = at;
        at++;
      }
    }
    fs->nref = at;
  }
  for(int i=0;i<fs->nref;i++)
    fs->realnameidx[i] = i;

  fprintf(stderr,"\t-> Number of nref %d in file: \'%s\'\n",fs->nref,fa);
  if(fs->nref==0){
    fprintf(stderr,"\t-> Possible error, no sequences loaded\n");
    exit(1);
  }
  fasta_sampler_setprobs(fs);
  for(int i=0;i<(int)SubsetChr.size();i++)
    free(SubsetChr[i]);
  return fs;
}

//functions returns head of chromosome, with posB, posE, chr_idx and fraglength set accordingly
char *sample(fasta_sampler *fs,mrand_t *mr,char **chromoname,int &chr_idx,int &posB,int &posE,int &fraglength){
  chr_idx = 0;
  if(fs->nref>1)
    chr_idx = ransampl_draw2(fs->ws,mrand_pop(mr),mrand_pop(mr));
  *chromoname = fs->seqs_names[chr_idx];
  posB = mrand_pop(mr)*fs->seqs_l[chr_idx]; //abs(mrand_pop_long(mr)) % fs->seqs_l[chr_idx];
  posE = posB +fraglength;
  
  if(posE>=fs->seqs_l[chr_idx]){
    posE=fs->seqs_l[chr_idx];
    posB=posE-fraglength;
  }
  char *ret =fs->seqs[chr_idx]; 
  chr_idx = fs->realnameidx[chr_idx];
  return ret;
}

void dump_internal(fasta_sampler *fs,const char* filename){
  kstring_t *fa_s[1];
  fa_s[0] =(kstring_t*) calloc(1,sizeof(kstring_t));
  fa_s[0]->s = NULL;
  fa_s[0]->l = fa_s[0]->m = 0;

  BGZF **bgzf_fp = (BGZF **) calloc(1,sizeof(BGZF *));

  int mt_cores = 1;
  int bgzf_buf = 256;
  const char* mode = "wu";

  bgzf_fp[0] = bgzf_open(filename,mode);
  bgzf_mt(bgzf_fp[0],mt_cores,bgzf_buf);
  
  for(int i=0;i<fs->nref;i++){
    ksprintf(fa_s[0],">%s\n%s\n",fs->seqs_names[i],fs->seqs[i]);
  }
  
  assert(bgzf_write(bgzf_fp[0],fa_s[0]->s,fa_s[0]->l)!=0);

  free(fa_s[0]->s);
  free(fa_s[0]);

  bgzf_close(bgzf_fp[0]);
  free(bgzf_fp);
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
  delete [] fs->realnameidx;
  for(ploidymap::iterator it=fs->pldmap.begin();it!=fs->pldmap.end();it++){
    delete [] it->second;
  }
  ransampl_free(fs->ws);
  delete fs;
}

#ifdef __WITH_MAIN__
int main(int argc,char**argv){

  const char* SubsetChr = NULL;
      

  int seed = 101;
  mrand_t *mr = mrand_alloc(3,seed);
  char ref[] = "Test_Examples/Mycobacterium_leprae.fa.gz";
  fasta_sampler *fs = fasta_sampler_alloc(ref,SubsetChr);
  
  char *chr; //this is an unallocated pointer to a chromosome name, eg chr1, chrMT etc
  int chr_idx;
  int posB,posE;//this is the first and last position of our fragment
  char *seq;//actual sequence, this is unallocated
  int fraglength = 100;

  size_t nit =0;
  size_t ngen = 30e6;
  char seq_r1[1024] = {0};
  while(nit++<100){
    //fraglength = abs(mrand_pop_long(mr)) % 1000;
    seq = sample(fs,mr,&chr,chr_idx,posB,posE,fraglength);
    strncpy(seq_r1,seq,fraglength);
    fprintf(stdout,"nit:%lu\tchromo:%s\tposB:%d\tposE:%d\tfraglength:%d\texample:%s\n",nit,chr,posB,posE,fraglength,seq_r1);
  }
  return 0;
}
#endif
//g++ fasta_sampler.cpp mrand.o RandSampling.o htslib/libhts.a -std=c++11 -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto -D__WITH_MAIN__ -o fasta
