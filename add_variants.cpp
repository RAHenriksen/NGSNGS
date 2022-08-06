#include <htslib/vcf.h>
#include "fasta_sampler.h"

int verbose2000 = 0;
int *mapper(bcf_hdr_t *bcf_hdr,char2int &fai2idx,int &bongo){
  int seqnames_l;
  const char **bcf_chrom_names = bcf_hdr_seqnames(bcf_hdr,&seqnames_l);
  fprintf(stderr,"Raw number of contigs from bcf: %d\n",seqnames_l);

  std::map<int,int> res;
  int max_l = -1;
  for (int i=0; i<seqnames_l;i++){
    char2int::iterator it = fai2idx.find(bcf_chrom_names[i]);
    if(it==fai2idx.end()){
      if(verbose2000++<5)
	fprintf(stderr,"\t-> Chromesome id from bcf: %s does not exists in fastafile this message is printed out %d more\n", bcf_chrom_names[i],5-verbose2000);
      continue;
    }else{
      res[i] = it->second;
      if(i>max_l)
	max_l = i;
    }
  }
  
  fprintf(stderr,"\t-> From header in bcf we observe %lu different scaffolds/chromosomes we will allocate lookup table with dim: %d\n",res.size(),max_l);
  if(res.size()==0)
    return NULL;
      
  int *lookup = new int[max_l];
  for(int i=0;i<max_l;i++)
    lookup[i] = -1;

  for(std::map<int,int>::iterator it=res.begin();it!=res.end();it++)
    lookup[it->first] = it->second;
  
  return lookup;
}

//fasta sampler struct, index for chromosomenaem, position, the alleles, the genotypes and the ploidy. 
void add_variant(fasta_sampler *fs, int chr_idx,int pos,char **alleles, int32_t *gts, int ploidy){
  //  fprintf(stderr,"Adding genotype informaiton for chromosome: %s pos: %d\n",fs->seqs_names[chr_idx],pos);
  char ref = fs->seqs[chr_idx][pos];
  for(int i=0;0&&i<ploidy;i++){
    fprintf(stderr,"%d) gt: %d and GT alleles are  '%s' (first one should be reference?)\n",i,gts[i],alleles[i]);
  }
  //first check if all parental chromosomes exists
  int isThere  = 1;
  char buf[1024];
  for(int i=1;i<ploidy;i++){
    snprintf(buf,1024,"%s_ngsngs%d",fs->seqs_names[chr_idx],i);
    char2int::iterator it= fs->char2idx.find(buf);
    if(it==fs->char2idx.end())
      isThere  = 0;
  }
  //its not there lets add it
  if(isThere == 0){
    int nref = fs->nref+ploidy-1;
    char **seqs = (char**) new char*[nref];
    char **seqs_names = (char**) new char*[nref];
    int *seqs_l = new int[nref];
    for(int i=0;i<fs->nref;i++){
      seqs[i] = fs->seqs[i];
      seqs_names[i] = fs->seqs_names[i];
      seqs_l[i] =fs->seqs_l[i];
    }
    delete [] fs->seqs;
    delete [] fs->seqs_names;
    delete [] fs->seqs_l;

    fs->seqs = seqs;
    fs->seqs_names = seqs_names;
    fs->seqs_l = seqs_l;
    for(int i=1;i<ploidy;i++) {
      snprintf(buf,1024,"%s_ngsngs%d",fs->seqs_names[chr_idx],i);
      fprintf(stderr,"Allocating extra space space space\n");
      fs->seqs[fs->nref+i-1] = strdup(seqs[chr_idx]);
      fs->seqs_names[fs->nref+i-1] = strdup(buf);
      fs->seqs_l[fs->nref+i-1] = fs->seqs_l[chr_idx];
      fs->char2idx[fs->seqs_names[fs->nref+i-1]] =fs->nref+i-1;
    }
    fs->nref = nref;
  }

  //now lets add snp
  for(int i=0;i<ploidy;i++){
    if(i==0)
      fs->seqs[chr_idx][pos] = alleles[gts[i]][0];
    else{
      snprintf(buf,1024,"%s_ngsngs%d",fs->seqs_names[chr_idx],i);
      char2int::iterator it = fs->char2idx.find(buf);
      assert(it!=fs->char2idx.end());
      fs->seqs[it->second][pos] = alleles[gts[i]][0];
    }
  }
  
}


int add_variants(fasta_sampler *fs,const char *bcffilename){
  htsFile *bcf = bcf_open(bcffilename, "r");
  bcf_hdr_t *bcf_head = bcf_hdr_read(bcf);
  bcf1_t *brec = bcf_init();
  int whichsample=-1;//if minus one then ref and alt fields are used, if nonnegative then it is used as offset to which genotype to use from the GT fields
  //map chromosomenames in bcf to index in fastafile
  int max_l;
  int *bcf_idx_2_fasta_idx = mapper(bcf_head,fs->char2idx,max_l);
  int ret = -1;
  int nsamples = bcf_hdr_nsamples(bcf_head);
  int32_t ngt_arr = 0;     
  int32_t *gt_arr = NULL;
  int ngt;
  int inferred_ploidy =2;//assume ploidy is two.
  int32_t nodatagt[2] = {0,1};
  
  while(((ret=bcf_read(bcf,bcf_head,brec)))==0){
    bcf_unpack((bcf1_t*)brec, BCF_UN_ALL);
    int fai_chr = bcf_idx_2_fasta_idx[brec->rid];
    
    if(fai_chr==-1){
      fprintf(stderr,"chrname: %s does not exists in fasta reference\n",bcf_hdr_id2name(bcf_head,brec->rid));
      exit(0);
    }
    
    if(whichsample!=-1){
      ngt = bcf_get_genotypes(bcf_head, brec, &gt_arr, &ngt_arr);
      assert((ngt %nsamples)==0);
      inferred_ploidy = ngt/nsamples;
    }
    
    if(whichsample!=-1)
      add_variant(fs,brec->rid,brec->pos,brec->d.allele,gt_arr+inferred_ploidy*whichsample,inferred_ploidy);
    else
      add_variant(fs,brec->rid,brec->pos,brec->d.allele,nodatagt,inferred_ploidy);
  }
  bcf_hdr_destroy(bcf_head);
  bcf_destroy(brec); 
  bcf_close(bcf);
  return 0;
}


#ifdef __WITH_MAIN__

int main(int argc,char **argv){
  const char* subsetchr = "1,2,3,4,5,MT";
  
  int seed = 101;
  mrand_t *mr = mrand_alloc(3,seed);
  char *ref = "../hs37d5.fa.gz";
  char *vcf = "sub2.bcf";
  fasta_sampler *fs = fasta_sampler_alloc(ref,subsetchr);
  fprintf(stderr,"Done adding fasta, will now add variants\n");
  fasta_sampler_print(stderr,fs);
  add_variants(fs,vcf);
  fasta_sampler_setprobs(fs);
  fasta_sampler_print(stderr,fs);
  char *chr; //this is an unallocated pointer to a chromosome name, eg chr1, chrMT etc
  int chr_idx;
  int posB,posE;//this is the first and last position of our fragment
  char *seq;//actual sequence, this is unallocated
  int fraglength = 100;

  size_t nit =0;
  size_t ngen = 30e6;
  char seq_r1[1024] = {0};
  while(nit++<4){
    //fraglength = abs(mrand_pop_long(mr)) % 1000;
    seq = sample(fs,mr,&chr,chr_idx,posB,posE,fraglength);
    strncpy(seq_r1,seq,fraglength);
    fprintf(stdout,"nit:%lu\tchromo:%s\tposB:%d\tposE:%d\tfraglength:%d\texample:%s\n",nit,chr,posB,posE,fraglength,seq_r1);
  }
  fprintf(stderr,"after\n");

  return 0;
  // g++ add_variants.cpp -D__WITH_MAIN__ -lhts mrand.o  fasta_sampler.o RandSampling.o   -ggdb
}

#endif
