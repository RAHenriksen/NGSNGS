
#include "fasta_sampler.h"
#include "NGSNGS_misc.h"
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

void fasta_sampler_print2(fasta_sampler *fs){
  fprintf(stderr,"--------------\n[%s]\tfai: %p nref: %d\n",__FUNCTION__,fs->fai,fs->nref);
  for(int i=0;i<fs->nref;i++)
    fprintf(stderr,"[%s]\tidx:%d) name:%s example:%.15s length: %d realidx: %d\n",__FUNCTION__,i,fs->seqs_names[i],fs->seqs[i],fs->seqs_l[i],fs->realnameidx[i]);

  for(char2int::iterator it=fs->char2idx.begin();it!=fs->char2idx.end();it++)
    fprintf(stderr,"[%s]\tkey: %s val:%d\n",__FUNCTION__,it->first,it->second);
  fprintf(stderr,"------------\n");
}

fasta_sampler *fasta_sampler_alloc_full(const char *fa){

  fasta_sampler *fs = new fasta_sampler;
  fs->ws = NULL;
  fs->fai = NULL;
  assert(((fs->fai =fai_load(fa)))!=NULL);
  fs->BedReferenceEntries = NULL;
  fs->BedReferenceCount = 0;
  
  fs->nref = faidx_nseq(fs->fai);
  fs->seqs = new char* [fs->nref];
  fs->seqs_l = new int[fs->nref];
  fs->seqs_names = new char *[fs->nref];
  fs->realnameidx = new int[fs->nref];

  for(int i=0;i<fs->nref;i++){
    fs->seqs_names[i] = strdup(faidx_iseq(fs->fai,i));
    fs->char2idx[fs->seqs_names[i]] = i;
    fs->seqs[i] = fai_fetch(fs->fai,fs->seqs_names[i],fs->seqs_l+i);
    /*
    abstract    Fetch the sequence in a region.
    param  fai  Pointer to the faidx_t struct
    param  reg  Region in the format "chr2:20,000-30,000"
    param  len  Length of the region; -2 if seq not present, -1 general error
    return      Pointer to the sequence; null on failure

    @discussion The returned sequence is allocated by malloc family
    and should be destroyed by end users by calling free() on it.
    char *fai_fetch(const faidx_t *fai, const char *reg, int *len);
     */
  }

  for(int i=0;i<fs->nref;i++)
    fs->realnameidx[i] = i;

  fprintf(stderr,"\t-> Number of nref %d in file: \'%s\'\n",fs->nref,fa);
  if(fs->nref==0){
    fprintf(stderr,"\t-> Possible error, no sequences loaded\n");
    exit(1);
  }
  fasta_sampler_setprobs(fs);

  return fs;
}

fasta_sampler *fasta_sampler_alloc_subset(const char *fa,const char *SpecificChr){
  
  std::vector<char *> SubsetChr;
  
  SubsetChr.push_back(strtok(strdup(SpecificChr),"\", \t"));
  char *chrtok = NULL;
  while(((chrtok=strtok(NULL,"\", \t"))))
    SubsetChr.push_back(strdup(chrtok));
  fprintf(stderr,"\t-> Number of entries in chromosome list: %lu\n",SubsetChr.size());

  fasta_sampler *fs = new fasta_sampler;
  fs->ws = NULL;
  fs->fai = NULL;
  assert(((fs->fai =fai_load(fa)))!=NULL);
  fs->BedReferenceEntries = NULL;
  fs->BedReferenceCount = 0;
  
  fs->nref = SubsetChr.size();

  fs->seqs = new char* [fs->nref];
  fs->seqs_l = new int[fs->nref];
  fs->seqs_names = new char *[fs->nref];
  fs->realnameidx = new int[fs->nref];
  
  int at = 0;
  for(int i=0;i<(int)SubsetChr.size();i++){
    if(faidx_has_seq(fs->fai, SubsetChr[i])){
      fs->seqs_names[at] = strdup(SubsetChr[i]);
      fs->seqs[at] = fai_fetch(fs->fai,fs->seqs_names[i],fs->seqs_l+at);
      fs->char2idx[fs->seqs_names[i]] = at;
      at++;
    }
  }
  fs->nref = at;
  
  for(int i=0;i<fs->nref;i++)
    fs->realnameidx[i] = i;

  fprintf(stderr,"\t-> Number of nref %d in file: \'%s\' and selected %d chromosomes/contigs/scaffolds \n",faidx_nseq(fs->fai),fa,fs->nref);
  if(fs->nref==0){
    fprintf(stderr,"\t-> Possible error, no sequences loaded\n");
    exit(1);
  }
  fasta_sampler_setprobs(fs);

  //couldn't i in theory close down the reference here and now?
  for(int i=0;i<(int)SubsetChr.size();i++)
    free(SubsetChr[i]);
  return fs;
}

fasta_sampler *fasta_sampler_alloc_bedentry(const char *fa,const char *bedfilename,size_t flanking){

  fasta_sampler *fs = new fasta_sampler;
  fs->ws = NULL;
  fs->fai = NULL;
  assert(((fs->fai =fai_load(fa)))!=NULL);
  fs->BedReferenceEntries = NULL;
  fs->BedReferenceCount = 0;

  int BedEntryCount = 0;
  BedEntry* BedEntries = NULL;

  //read in the bed entries
  BedEntries = readBedFile(bedfilename,&BedEntryCount);
  //sort the coordinates
  sortBedEntries(BedEntries, BedEntryCount);

  //merge the coordinates
  int mergedCount = 0;
  BedEntry* mergedEntries = mergeOverlappingRegions(BedEntries, BedEntryCount, &mergedCount);
  free(BedEntries);

  fs->BedReferenceEntries = checkbedentriesfasta(fs,mergedEntries,mergedCount,&fs->BedReferenceCount);

  //for (int i = 0; i < fs->BedReferenceCount; i++){fprintf(stderr,"merge bed file information bed entry \t%s\t%d\t%d\n", fs->BedReferenceEntries[i].chromosome, fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);}

  free(mergedEntries);

  //fprintf(stderr,"\t-> Input bed file had %d regions, after merging the overlapping regions there is %d and post-filtering there is %d\n",BedEntryCount,mergedCount,fs->BedReferenceCount);
 
  fs->nref = fs->BedReferenceCount; 
  fs->seqs = new char* [fs->nref];
  fs->seqs_l = new int[fs->nref];
  fs->seqs_names = new char *[fs->nref];
  fs->realnameidx = new int[fs->nref];
     
  //create the names to extract the region defined in the bed entry when using fai_fetch
  for(int i=0;i<fs->nref;i++){
    //ReferenceEntries[i].chromosome should still be in accordance with the reference genome
    size_t refchr_length = faidx_seq_len(fs->fai, fs->BedReferenceEntries[i].chromosome);
    size_t bedflankstart; 
    size_t bedflankend;   
    
    if((int)(fs->BedReferenceEntries[i].start-flanking) < 0){
      //fprintf(stderr,"case 1 \t %d\n",(ReferenceEntries[i].start-flanking));
      bedflankstart = 1;
    }
    else if((int)(fs->BedReferenceEntries[i].start-flanking) > 0){
      //fprintf(stderr,"case 2 \t %d\n",(ReferenceEntries[i].start-flanking));
      bedflankstart = fs->BedReferenceEntries[i].start-flanking; 
    }
    
    if((fs->BedReferenceEntries[i].end+flanking) < refchr_length){
      //fprintf(stderr,"case 3 \t %d\n",(ReferenceEntries[i].end+flanking));
      bedflankend = fs->BedReferenceEntries[i].end + flanking;
    }
    if((fs->BedReferenceEntries[i].end+flanking) > refchr_length){
      //fprintf(stderr,"case 4 \t %d\n",(ReferenceEntries[i].end+flanking));
      bedflankend = refchr_length-1;
    }

    fs->BedReferenceEntries[i].start = bedflankstart;
    fs->BedReferenceEntries[i].end = bedflankend;

    int length = snprintf(NULL, 0, "%s:%d-%d", fs->BedReferenceEntries[i].chromosome, bedflankstart, bedflankend);
    fs->seqs_names[i] = (char*) malloc((length + 1) * sizeof(char));
    sprintf(fs->seqs_names[i], "%s:%d-%d", fs->BedReferenceEntries[i].chromosome, bedflankstart, bedflankend);
  }

  for(int i=0;i<fs->nref;i++){
    fs->char2idx[fs->seqs_names[i]] = i;
    //fprintf(stderr,"names are %s\n",fs->seqs_names[i]);
    fs->seqs[i] = fai_fetch(fs->fai,fs->seqs_names[i],fs->seqs_l+i);
  }

  for(int i=0;i<fs->nref;i++)
    fs->realnameidx[i] = i;

  fprintf(stderr,"\t-> Number of nref %d in file: \'%s\'\n",fs->nref,fa);
  if(fs->nref==0){
    fprintf(stderr,"\t-> Possible error, no sequences loaded\n");
    exit(1);
  }
  fasta_sampler_setprobs(fs);

  return fs;
}

fasta_sampler *fasta_sampler_alloc_maskbedentry(const char *fa,const char *bedfilename,size_t flanking){

  fasta_sampler *fs = new fasta_sampler;
  fs->ws = NULL;
  fs->fai = NULL;
  assert(((fs->fai =fai_load(fa)))!=NULL);
  fs->BedReferenceEntries = NULL;
  fs->BedReferenceCount = 0;

  int BedEntryCount = 0;
  BedEntry* BedEntries = NULL;

  //read in the bed entries
  BedEntries = readBedFile(bedfilename,&BedEntryCount);
  //for (int i = 0; i < BedEntryCount; i++){fprintf(stderr,"read bed file information bed entry \t%s\t%d\t%d\n", BedEntries[i].chromosome, BedEntries[i].start, BedEntries[i].end);}
  //fprintf(stderr,"DONE readBedFile\n");

  //sort the coordinates
  sortBedEntries(BedEntries, BedEntryCount);
  //for (int i = 0; i < BedEntryCount; i++){fprintf(stderr,"sort bed file information bed entry \t%s\t%d\t%d\n", BedEntries[i].chromosome, BedEntries[i].start, BedEntries[i].end);}
  //fprintf(stderr,"DONE sortBedEntries\n");
  
  //merge the coordinates
  int mergedCount = 0;
  BedEntry* mergedEntries = mergeOverlappingRegions(BedEntries, BedEntryCount, &mergedCount);
  //for (int i = 0; i < mergedCount; i++){fprintf(stderr,"merge bed file information bed entry \t%s\t%d\t%d\n", mergedEntries[i].chromosome, mergedEntries[i].start, mergedEntries[i].end);}
  free(BedEntries);

  int ReferenceCount = 0;
  //fprintf(stderr,"DONE mergeOverlappingRegions\n");

  fs->BedReferenceEntries = maskbedentriesfasta(fs,mergedEntries,mergedCount,&fs->BedReferenceCount);
  
  //fprintf(stderr,"\t-> Input bed file had %d regions, after merging the overlapping regions there is %d and post-filtering there is %d\n",BedEntryCount,mergedCount,fs->BedReferenceCount);
  //for (int i = 0; i < fs->BedReferenceCount; i++){fprintf(stderr,"merge bed file information mask bed entry \t%s\t%d\t%d\n", fs->BedReferenceEntries[i].chromosome, fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);}

  free(mergedEntries);
  //fprintf(stderr,"DONE maskbedentriesfasta\n");

  fs->nref = fs->BedReferenceCount; 
  fs->seqs = new char* [fs->nref];
  fs->seqs_l = new int[fs->nref];
  fs->seqs_names = new char *[fs->nref];
  fs->realnameidx = new int[fs->nref];


  //create the names to extract the region defined in the bed entry
  for(int i=0;i<fs->nref;i++){
    int length = snprintf(NULL, 0, "%s:%d-%d", fs->BedReferenceEntries[i].chromosome,  fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);
    fs->seqs_names[i] = (char*) malloc((length + 1) * sizeof(char));
    sprintf(fs->seqs_names[i], "%s:%d-%d", fs->BedReferenceEntries[i].chromosome, fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);
  }

  for(int i=0;i<fs->nref;i++){
    fs->char2idx[fs->seqs_names[i]] = i;
    //fprintf(stderr,"names are %s\n",fs->seqs_names[i]);
    fs->seqs[i] = fai_fetch(fs->fai,fs->seqs_names[i],fs->seqs_l+i);
  }

  for(int i=0;i<fs->nref;i++)
    fs->realnameidx[i] = i;

  fprintf(stderr,"\t-> Number of nref %d in file: \'%s\'\n",fs->nref,fa);
  if(fs->nref==0){
    fprintf(stderr,"\t-> Possible error, no sequences loaded\n");
    exit(1);
  }
  fasta_sampler_setprobs(fs);

  return fs;
}

//functions returns head of chromosome, with posB, posE, chr_idx and fraglength set accordingly
char *sample(fasta_sampler *fs,mrand_t *mr,char **chromoname,int &chr_idx,int &posB,int &posE,int &fraglength, size_t &chr_end, int simmode){
  chr_idx = 0;
  if(fs->nref>1)
    chr_idx = ransampl_draw2(fs->ws,mrand_pop(mr),mrand_pop(mr));
  *chromoname = fs->seqs_names[chr_idx];
  chr_end = fs->seqs_l[chr_idx];
  posB = mrand_pop(mr)*fs->seqs_l[chr_idx]; //abs(mrand_pop_long(mr)) % fs->seqs_l[chr_idx];
  posE = posB +fraglength;
  
  if(posE>=fs->seqs_l[chr_idx] && simmode == 0){
    // linear simulations, move the end point to ensure no crossing of the breakpoint.
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

  free(fs->BedReferenceEntries);
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
