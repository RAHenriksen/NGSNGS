
#include "fasta_sampler.h"
#include "RandSampling.h"
#include "mrand.h"
#define MAXBINS 100

fasta_sampler *fasta_sampler_alloc(const char *fa,const char *SpecificChr){
  
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
  else
    for(int i=0;i<SubsetChr.size();i++){
      fs->seqs_names[i] = strdup(SubsetChr[i]);
      fs->seqs[i] = fai_fetch(fs->fai,fs->seqs_names[i],fs->seqs_l+i);
      fprintf(stderr,"%d)\tchr:%s\tlen:%d\n",i,fs->seqs_names[i],fs->seqs_l[i]);
    }
  fs->ws = ransampl_alloc(fs->nref);
  double *p = new double[fs->nref];
  fs->seq_l_total =0;
  for(int i=0;i<fs->nref;i++)
    fs->seq_l_total += fs->seqs_l[i];
  for(int i=0;i<fs->nref;i++)
    p[i] = ((double) fs->seqs_l[i])/fs->seq_l_total;
  ransampl_set(fs->ws,p);
  return fs;
}

char *sample(fasta_sampler *fs,mrand_t *mr,char **chromoname,int &posB,int &posE,int fraglength){
  int idx = ransampl_draw2(fs->ws,mrand_pop(mr),mrand_pop(mr));
  *chromoname = fs->seqs_names[idx];
  //std::cout << fs->seqs_names[0] << fs->seqs_names[1] << fs->seqs_names[2] << std::endl;
  posB = abs(mrand_pop_long(mr)) % fs->seqs_l[idx]+1;
  posE = posB +fraglength;
  if(posE>fs->seqs_l[idx]){
    posE=fs->seqs_l[idx]-posB;
  }
  return fs->seqs[idx]+posB;
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



/*


fasta_sampler *fasta_sampler_alloc(const char *fa,std::vector<char *> chromoname){
  fasta_sampler *fs = new fasta_sampler;
  fs->fai = NULL;
  assert(((fs->fai =fai_load(fa)))!=NULL);
  if(chromoname.size()==0)
  fs->nref = faidx_nseq(fs->fai);
  else
  fs-> chromoname.size();

  fs->seqs = new char* [fs->nref];
  fs->seqs_l = new int[fs->nref];
  fs->seqs_names = new char *[fs->nref];
  if(chromoname.size()==0)
  for(int i=0;i<fs->nref;i++){
    fs->seqs_names[i] = strdup(faidx_iseq(fs->fai,i));
    // fs->seqs_l[i] = faidx_seq_len(fs->fai,fs->seqs_names[i]);
    fs->seqs[i] = fai_fetch(fs->fai,fs->seqs_names[i],fs->seqs_l+i);
    //  fprintf(stderr,"%d)\tchr:%s\tlen:%d\n",i,fs->seqs_names[i],fs->seqs_l[i]);
  }
  else
  for(int i=0;i<crhomoname.size();i++){
    fs->seqs_names[i] = strdup(chromoname(i));
    // fs->seqs_l[i] = faidx_seq_len(fs->fai,fs->seqs_names[i]);
    fs->seqs[i] = fai_fetch(fs->fai,fs->seqs_names[i],fs->seqs_l+i);
    //  fprintf(stderr,"%d)\tchr:%s\tlen:%d\n",i,fs->seqs_names[i],fs->seqs_l[i]);
  }
  else
  fs->ws = ransampl_alloc(fs->nref);
  double *p = new double[fs->nref];
  fs->seq_l_total =0;
  for(int i=0;i<fs->nref;i++)
    fs->seq_l_total += fs->seqs_l[i];
  for(int i=0;i<fs->nref;i++)
    p[i] = ((double) fs->seqs_l[i])/fs->seq_l_total;
  ransampl_set(fs->ws,p);
  return fs;
}

*/