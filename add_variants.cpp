#include <map>
#include <htslib/vcf.h>
#include "fasta_sampler.h"
#include "add_variants.h"

int verbose2000 = 0;
int *mapper(bcf_hdr_t *bcf_hdr,char2int &fai2idx,int &bongo){
  int seqnames_l;
  const char **bcf_chrom_names = bcf_hdr_seqnames(bcf_hdr,&seqnames_l);
  fprintf(stderr,"\t-> Raw number of contigs from bcf: %d\n",seqnames_l);

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
  free(bcf_chrom_names);
  max_l = max_l +1;
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

int extend_fasta_sampler(fasta_sampler *fs,int fs_chr_idx,int ploidy){
  if(ploidy ==1){
    //ploidy is haploid, will perform updownstream sorting of the reverse positions
    ploidymap::iterator it = fs->pldmap.find(fs_chr_idx);
    if(it!=fs->pldmap.end()){
      return 0;
    }
    int *pldmap = new int[5];
    pldmap[0]=fs_chr_idx;
    pldmap[1]=pldmap[2]=pldmap[3]=pldmap[4]=-1;
    fs->pldmap[fs_chr_idx] = pldmap;
    return 0;
  }
  int isThere  = 1;
  char buf[1024];
  for(int i=1;i<ploidy;i++){
    snprintf(buf,1024,"%s_ngsngs%d",fs->seqs_names[fs_chr_idx],i);
    char2int::iterator it= fs->char2idx.find(buf);
    if(it==fs->char2idx.end())
      isThere  = 0;
  }
  if(isThere == 0){
    // Will duplicate chromosomes to emulate parental chromosomes
    int nref = fs->nref+ploidy-1;
    char **seqs = (char**) new char*[nref];
    char **seqs_names = (char**) new char*[nref];
    int *seqs_l = new int[nref];
    int *realnameidx = new int[nref];
    int *pldmap = new int[5];
    for(int i=0;i<fs->nref;i++){
      seqs[i] = fs->seqs[i];
      seqs_names[i] = fs->seqs_names[i];
      seqs_l[i] =fs->seqs_l[i];
      realnameidx[i] = fs->realnameidx[i];
    }
    delete [] fs->seqs;
    delete [] fs->seqs_names;
    delete [] fs->seqs_l;
    delete [] fs->realnameidx;

    fs->seqs = seqs;
    fs->seqs_names = seqs_names;
    fs->seqs_l = seqs_l;
    fs->realnameidx = realnameidx;
    ploidymap::iterator it = fs->pldmap.find(fs_chr_idx);
    assert(it==fs->pldmap.end());
    pldmap[0] = fs_chr_idx;
    for(int i=1;i<ploidy;i++) {
      snprintf(buf,1024,"%s_ngsngs%d",fs->seqs_names[fs_chr_idx],i);
      fs->seqs[fs->nref+i-1] = strdup(seqs[fs_chr_idx]);
      fs->seqs_names[fs->nref+i-1] = strdup(buf);
      fs->seqs_l[fs->nref+i-1] = fs->seqs_l[fs_chr_idx];
      fs->char2idx[fs->seqs_names[fs->nref+i-1]] =fs->nref+i-1;
      fs->realnameidx[fs->nref+i-1] = fs_chr_idx;
      pldmap[i] = fs->nref+i-1;
    }
    fs->pldmap[fs_chr_idx]  = pldmap;
    fs->nref = nref;
    
  }
  return 0;
}

/*
  Assumption, input map is sorted according to chr and position.
  reset is set to one when change of chromosome is detected. This will flush remainder of chrosomes.
  Everytime a new indel is found, program will write all bp from previous event until (including) the current position.
  event is then updated by incrementing reference bp for indel. 
*/

void add_indels(fasta_sampler *fs,bcfmap &mybcfmap,bcf_hdr_t *hdr,int ploidy){
  int last[ploidy];
  int last_fid = -1;
  int reset = 1;
  int maxsize = -1;
  int *fsoffsets = NULL;
  for(int i=0;i<fs->nref;i++)
    if(fs->seqs_l[i]>maxsize)
      maxsize= fs->seqs_l[i];

  maxsize += 1000;

  char **elephant =new char*[ploidy];
  for(int i=0;i<ploidy;i++)
    elephant[i] =(char*) calloc(maxsize,sizeof(char));

  #if 0  
  for(bcfmap::iterator it=mybcfmap.begin();0&&it!=mybcfmap.end();it++){
    fprintf(stderr,"%d %d\n",it->first.rid,it->first.pos);
  }
  #endif
  for(bcfmap::iterator it=mybcfmap.begin();it!=mybcfmap.end();it++){
    #if 0
    for(int i=0;0&&i<ploidy;i++)
      fprintf(stderr,"[%s] gt[%d]: %d\n",__FUNCTION__,i,it->first.gt[i]);
    #endif
    if(it->first.rid!=last_fid){
      last_fid = it->first.rid;
      reset = 1;
    }
    
    if(reset){
      if(fsoffsets!=NULL){
        for(int i=0;i<ploidy;i++){

          strcat(elephant[i],fs->seqs[fsoffsets[i]]+last[i]);
          fs->seqs_l[fsoffsets[i]] = strlen(elephant[i]);
          free(fs->seqs[i]);
          fs->seqs[i] = elephant[i];
        }
      }

      ploidymap::iterator it = fs->pldmap.find(last_fid);
      assert(it!=fs->pldmap.end());
      fsoffsets = it->second;
      for(int i=0;i<ploidy;i++)
	      last[i] = 0;
        reset = 0;
      }
      #if 0
      for(int i=0;i<ploidy;i++)
        fprintf(stderr,"last[%d]: %d\n",i,last[i]);
      #endif
      bcf1_t *brec = it->second;
      bcf_unpack(brec, BCF_UN_ALL);
      //we now loop over mother, father and other parental chromosomes
      for(int i=0;i<ploidy;i++) {
        //first copy contigous block from last operator
        if(last[i]>it->first.pos){
          last[i] = it->first.pos;
      }
      int nitems2copy = brec->pos-last[i]+1;
      assert(strlen(elephant[i])+brec->pos-last[i]< (size_t)maxsize );//funky assert
      strncat(elephant[i],fs->seqs[fsoffsets[i]]+last[i],nitems2copy);
      char *allele = NULL;
      allele = it->second->d.allele[it->first.gt[i]];
      assert(allele!=NULL);
            
      int turbocharge =0;
      #if 0
      int isdel = 0;
      int isins = 0;
      // allele -> alternatie
      fprintf(stderr,"Position: %lld ref %s alt %s and length of ref %zu and alt: %zu\n",it->first.pos,it->second->d.allele[0],allele,strlen(it->second->d.allele[0]),strlen(allele));
      if(it->first.gt[i]>0 && (strlen(it->second->d.allele[0]) > strlen(allele))){
        isdel = strlen(it->second->d.allele[0])-strlen(allele);
        fprintf(stderr,"Deletion of length: %d \n",isdel);
      }
      else if(it->first.gt[i]>0 && (strlen(it->second->d.allele[0]) < strlen(allele))){
        isins = strlen(allele)-strlen(it->second->d.allele[0]);
        fprintf(stderr,"Insertion of length: %d \n",isdel);
      }
      #endif
      if(it->first.gt[i]>0)
	      turbocharge =  strlen(allele)-strlen(it->second->d.allele[0]);
	      
      if(turbocharge<0){
        //is deletion then we just skip the number of bases
        last[i] = brec->pos+abs(turbocharge)+1;
      }
      else if(turbocharge>0){
        assert(strlen(elephant[i])+strlen(allele)<strlen(elephant[i])+maxsize);
        strncat(elephant[i],allele+1,strlen(allele));
        last[i] = brec->pos+1;
      }
      else{
        last[i] = brec->pos+1;
      }
    }
    bcf_destroy(brec);
    delete [] it->first.gt;
  }
  if(fsoffsets!=NULL){
    for(int i=0;i<ploidy;i++){
      
      strcat(elephant[i],fs->seqs[fsoffsets[i]]+last[i]);
      fs->seqs_l[fsoffsets[i]] = strlen(elephant[i]);
      free(fs->seqs[i]);
      fs->seqs[i] = elephant[i];
    }
  }
  delete[] elephant;
}

//fasta sampler struct, index for chromosomenaem, position, the alleles, the genotypes and the ploidy. 
void add_variant(fasta_sampler *fs, int fs_chr_idx,int pos,char **alleles, int32_t *gts, int ploidy){
  // Adding genotype information for chromosome,pos,ploidy
  char buf[1024];
  #if 0
  for(int i=0;0&&i<ploidy;i++){
    fprintf(stderr,"%d) gt: %d and GT alleles are  '%s' (first one should be reference?)\n",i,bcf_gt_allele(gts[i]),alleles[i]);
  }
  #endif

  //now lets add snp
  for(int i=0;i<ploidy;i++){
    if(i==0)
      fs->seqs[fs_chr_idx][pos] = alleles[bcf_gt_allele(gts[i])][0];
    else{
      snprintf(buf,1024,"%s_ngsngs%d",fs->seqs_names[fs_chr_idx],i);
      char2int::iterator it = fs->char2idx.find(buf);
      assert(it!=fs->char2idx.end());
      fs->seqs[it->second][pos] = alleles[bcf_gt_allele(gts[i])][0];
      
    }
  }
}

//if minus one then ref and alt fields are used, if nonnegative then it is used as offset to which genotype to use from the GT fields
int add_variants(fasta_sampler *fs,const char *bcffilename,int id,const char* Name){
  if(bcffilename==NULL)
    return 0;
  htsFile *bcf = bcf_open(bcffilename, "r");
  bcf_hdr_t *bcf_head = bcf_hdr_read(bcf);
  bcf1_t *brec = bcf_init();
  // So this could be the parameter to change for which individual
  int max_l;
  int *bcf_idx_2_fasta_idx = mapper(bcf_head,fs->char2idx,max_l);
  int ret = -1;
  int nsamples = bcf_hdr_nsamples(bcf_head);
  int32_t ngt_arr = 0;     
  int32_t *gt_arr = NULL;
  int ngt;
  int inferred_ploidy =5;//assume max ploidy is five, this will be adjusted when parsing data
  int32_t nodatagt[2]={bcf_gt_unphased(0),bcf_gt_unphased(1)};

  //int32_t nodatagt[2]={bcf_gt_unphased(0),bcf_gt_unphased(1)};
  int32_t *mygt=NULL;
  
  int whichsample = -1;
  if((int)id < 0){
    for (int i = 0; i < nsamples; i++){
      const char *sample_name = bcf_hdr_int2id(bcf_head, BCF_DT_SAMPLE, i);
      //fprintf(stderr,"SAMPLE NAME IS %s \n",sample_name);
      //fprintf(stderr,"PROVIDE INPUT IS %s \n",Name);
      if(strcasecmp(sample_name,Name)==0){
        //fprintf(stderr,"provided index is %d with name %s and option name %s\n",i,sample_name,Name);
        whichsample = (int)i;
        break;
      }
    }
  }
  else{
    whichsample = (int)id; //index of individual
  }
  
  //fprintf(stderr,"-------\nwhichsample is %d-------\n",whichsample);
  
  if ((int)whichsample < 0 ) {
    fprintf(stderr, "Error: Sample index out of range\n");
    bcf_hdr_destroy(bcf_head);
    bcf_destroy(brec);
    bcf_close(bcf);
    exit(1);
  }
  
  bcfmap mybcfmap;
  
  while(((ret=bcf_read(bcf,bcf_head,brec)))==0){
    bcf_unpack((bcf1_t*)brec, BCF_UN_ALL);
    int fai_chr = bcf_idx_2_fasta_idx[brec->rid];
    
    if (fai_chr == -1) {
      fprintf(stderr, "Error: Chromosome name not found in the reference\n");
      bcf_hdr_destroy(bcf_head);
      bcf_destroy(brec);
      bcf_close(bcf);
      exit(1);
    }

    ngt = bcf_get_genotypes(bcf_head, brec, &gt_arr, &ngt_arr);
    assert(ngt>0);
    assert((ngt %nsamples)==0);
    inferred_ploidy = ngt/nsamples;

    mygt = gt_arr+inferred_ploidy*whichsample;
    if(bcf_gt_is_missing(mygt[0])){
      fprintf(stderr,"\t-> Genotype is missing for pos: %ld will skip\n",brec->pos+1);
      continue;
    }

    //check if parental chromosomes exists otherwise add them
    extend_fasta_sampler(fs,fai_chr,inferred_ploidy);

    //figure out if its an indel
    int isindel  = 0;
    bcfkey key;
    key.gt = new int[inferred_ploidy];
    for(int i=0;i<inferred_ploidy;i++){
      key.gt[i] = bcf_gt_allele(mygt[i]);
      if(strlen(brec->d.allele[bcf_gt_allele(mygt[i])])>1) 
	      isindel =1;
      if(strlen(brec->d.allele[0])>1)//this is confusion
	      isindel =1;
    }
    if(isindel==0){
      add_variant(fs,fai_chr,brec->pos,brec->d.allele,mygt,inferred_ploidy);
      //fprintf(stderr,"\t-> Found an variant as rid: %d pos:%lld\n",brec->rid,brec->pos+1);
    }
    else{
      //fprintf(stderr,"\t-> Found an indel as rid: %d pos:%lld\n",brec->rid,brec->pos+1);
      key.rid=fai_chr;
      key.pos= (int) brec->pos;
      bcf1_t *duped= bcf_dup(brec);
      mybcfmap[key] =duped;
    }
  }
  fprintf(stderr,"\t-> Done adding snp variants\n");
  if(mybcfmap.size()>0){
    //fprintf(stderr,"\t-> Found some indels, these will now be added to internal datastructures\n");
    add_indels(fs,mybcfmap,bcf_head,inferred_ploidy);

  }
  
  free(gt_arr);
  delete[] bcf_idx_2_fasta_idx;
  fasta_sampler_setprobs(fs);
  bcf_hdr_destroy(bcf_head);
  bcf_destroy(brec); 
  bcf_close(bcf);
  return 0;
}



#ifdef __WITH_MAIN__

int main(int argc,char **argv){
  const char* subsetchr = "14"; //"NZ_CP029543.1";
  
  int seed = 101;
  mrand_t *mr = mrand_alloc(3,seed);
  char *ref = "chr14.fa";
  char *vcf = "Chr14_indiv4.vcf";//"DiploidTest.vcf"; //"DiploidTest.vcf";//"sub3.vcf";
  fasta_sampler *fs = fasta_sampler_alloc(ref,subsetchr);
  fprintf(stderr,"Done adding fasta, will now add variants\n");
  fasta_sampler_print(stderr,fs);

  add_variants(fs,vcf,1);
  fasta_sampler_print(stderr,fs);
  char *chr; //this is an unallocated pointer to a chromosome name, eg chr1, chrMT etc
  int chr_idx;
  char *seq;//actual sequence, this is unallocated
  int fraglength = 100;

  size_t nit =0;
  size_t ngen = 30e6;
  char seq_r1[1024] = {0};
  /*seq = sample(fs,mr,&chr,chr_idx,posB,posE,fraglength);
  strncpy(seq_r1,seq+34999998,30);
  fprintf(stderr,"sequence %s\n",seq_r1);*/
  while(nit++<1){
    int posB,posE;
    //this is the first and last position of our fragment
    //fraglength = abs(mrand_pop_long(mr)) % 1000;
    seq = sample(fs,mr,&chr,chr_idx,posB,posE,fraglength);
    strncpy(seq_r1,seq+(35000000-2),30);
    fprintf(stdout,"nit:%lu\tchromo:%s\tposB:%d\tposE:%d\tfraglength:%d\texample:%s\n",nit,chr,posB,posE,fraglength,seq_r1);
  }
  fprintf(stderr,"after\n");

  return 0;
  // g++ add_variants.cpp -D__WITH_MAIN__ mrand.o  fasta_sampler.o RandSampling.o htslib/libhts.a -std=c++11 -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto -ggdb
}

#endif
