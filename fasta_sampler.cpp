
#include "fasta_sampler.h"

fasta_sampler *fasta_sampler_alloc(char *fa){
  fasta_sampler *fs = new fasta_sampler;
  fs->fai = NULL;
  assert(((fs->fai =fai_load(fa)))!=NULL);
  fs->nref = faidx_nseq(fs->fai);
  fs->seqs = new char* [fs->nref];
  fs->seqs_l = new int[fs->nref];
  fs->seqs_names = new char *[fs->nref];
  for(int i=0;i<fs->nref;i++){
    fs->seqs_names[i] = strdup(faidx_iseq(fs->fai,i));
    // fs->seqs_l[i] = faidx_seq_len(fs->fai,fs->seqs_names[i]);
    fs->seqs[i] = fai_fetch(fs->fai,fs->seqs_names[i],fs->seqs_l+i);
    //  fprintf(stderr,"%d)\tchr:%s\tlen:%d\n",i,fs->seqs_names[i],fs->seqs_l[i]);
  }
  fs->ws = ransampl_alloc(fs->nref);
  double *p = new double[fs->nref];
  double sum =0;
  for(int i=0;i<fs->nref;i++)
    sum += fs->seqs_l[i];
  for(int i=0;i<fs->nref;i++)
    p[i] = ((double) fs->seqs_l[i])/sum;
  ransampl_set(fs->ws,p);
  return fs;
}

char *sample(fasta_sampler *fs,mrand_t *mr,char **chromoname,int &posB,int &posE,int &fraglength){
  int idx = ransampl_draw2(fs->ws,mrand_pop(mr),mrand_pop(mr));
  *chromoname = fs->seqs_names[idx];
  posB =  abs(mrand_pop_long(mr)) % fs->seqs_l[idx];
  posE = posB +fraglength;
  if(posE>fs->seqs_l[idx]){
    posE=fs->seqs_l[idx];
    fraglength = posE-posB;
  }
  return fs->seqs[idx] + posB;
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
  int seed = 101;
  mrand_t *mr = mrand_alloc(3,seed);
  fasta_sampler *fs = fasta_sampler_alloc("../hs37d5.fa.gz");
  //fasta_sampler *fs = fasta_sampler_alloc("Test_Examples/fire.fa");
  
  char *chr; //this is an unallocated pointer to a chromosome name, eg chr1, chrMT etc
  int posB,posE;//this is the first and last position of our fragment
  char *seq;//actual sequence, this is unallocated
  int fraglength;

  size_t nit =0;
  size_t ngen = 30e6;
  while(nit++<ngen){
    fraglength = abs(mrand_pop_long(mr)) % 1000;
    seq = sample(fs,mr,&chr,posB,posE,fraglength);
    fprintf(stdout,"nit:%lu\tchromo:%s\tposB:%d\tposE:%d\tfraglength:%d\texample:%.10s\n",nit,chr,posB,posE,fraglength,seq);
  }
  fprintf(stderr,"after\n");
  return 0;
}
#endif
