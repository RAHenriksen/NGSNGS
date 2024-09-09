
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
  /*
  fasta_sampler_setprobs - calculates the total length of all sequences, computes the probability of selecting each sequence based on its length, 
  and updates the sampling weights.
  
  @param fs: A pointer to a `fasta_sampler` struct. The function modifies the `ws` field to allocate new sampling weights and update the probabilities.
  */

  //add sampling operation
  if(fs->ws)
    ransampl_free(fs->ws);  // Free previously allocated sampling weights

  fs->ws = ransampl_alloc(fs->nref); // Allocate new sampling weights

  double *p = new double[fs->nref]; // Create an array to store probabilities

  // Calculate the total length of all sequences
  fs->seq_l_total =0;
  for(int i=0;i<fs->nref;i++){
    fs->seq_l_total += fs->seqs_l[i];
  }

  // Compute the probability for each sequence based on its length
  for(int i=0;i<fs->nref;i++)
    p[i] = ((double) fs->seqs_l[i])/fs->seq_l_total;
  
  // Update sampling weights with the computed probabilities
  ransampl_set(fs->ws,p);
  delete [] p;
}

void fasta_sampler_print(FILE *fp,fasta_sampler *fs){
  /*
  fasta_sampler_print - Prints information about the `fasta_sampler` object
  
  @param fp: A file pointer where the information will be printed.
  @param fs: A pointer to a `fasta_sampler` struct. The function retrieves and prints its internal details.
  */
  
  fprintf(fp,"--------------\n[%s]\tfai: %p nref: %d\n",__FUNCTION__,fs->fai,fs->nref);
  for(int i=0;i<fs->nref;i++)
    fprintf(fp,"[%s]\tidx:%d) name:%s example:%.15s length: %d realidx: %d\n",__FUNCTION__,i,fs->seqs_names[i],fs->seqs[i],fs->seqs_l[i],fs->realnameidx[i]);

  for(char2int::iterator it=fs->char2idx.begin();it!=fs->char2idx.end();it++)
    fprintf(fp,"[%s]\tkey: %s val:%d\n",__FUNCTION__,it->first,it->second);
  fprintf(fp,"------------\n");
}

void fasta_sampler_print2(fasta_sampler *fs){
  /*
  fasta_sampler_print2 - Prints detailed information about the `fasta_sampler` object to the standard error stream.
  
  This function outputs information including the file index, number of references, sequence details, and character-to-integer mappings, 
  similar to `fasta_sampler_print`, but prints to the standard error stream instead.
  
  @param fs: A pointer to a `fasta_sampler` struct. The function retrieves and prints its internal details.
  */

  fprintf(stderr,"--------------\n[%s]\tfai: %p nref: %d\n",__FUNCTION__,fs->fai,fs->nref);
  for(int i=0;i<fs->nref;i++)
    fprintf(stderr,"[%s]\tidx:%d) name:%s example:%.15s length: %d realidx: %d\n",__FUNCTION__,i,fs->seqs_names[i],fs->seqs[i],fs->seqs_l[i],fs->realnameidx[i]);

  for(char2int::iterator it=fs->char2idx.begin();it!=fs->char2idx.end();it++)
    fprintf(stderr,"[%s]\tkey: %s val:%d\n",__FUNCTION__,it->first,it->second);
  fprintf(stderr,"------------\n");
}

fasta_sampler *fasta_sampler_alloc_full(const char *fa){
  /*
  fasta_sampler_alloc_full - Allocates and initializes a `fasta_sampler` object based on the entire input reference FASTA file (-i).
    This function creates a new `fasta_sampler` instance, loads sequence data from a FASTA file, and initializes internal fields such as sequence names,
    lengths, and sampling weights. It also verifies that sequences are successfully loaded and sets up initial probabilities for sampling.
  
  @param fa: A string representing the path to the FASTA file from which sequences will be loaded.
  @return: A pointer to a newly allocated `fasta_sampler` struct with initialized fields and loaded sequences. Returns NULL on failure.
  */

  // Allocate memory for the fasta_sampler object and load index
  fasta_sampler *fs = new fasta_sampler;
  fs->ws = NULL;
  fs->fai = NULL;
  assert(((fs->fai =fai_load(fa)))!=NULL);

  // Initialize the number of references and allocate memory for sequence data,
  fs->BedReferenceEntries = NULL; // the bedreference is part of the struct but remains unused in this function
  fs->BedReferenceCount = 0;
  
  fs->nref = faidx_nseq(fs->fai);
  fs->seqs = new char* [fs->nref];
  fs->seqs_l = new int[fs->nref];
  fs->seqs_names = new char *[fs->nref];
  fs->realnameidx = new int[fs->nref];

  // Fetch and store sequence names, lengths, and actual sequences
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

  // initialize names with indices
  for(int i=0;i<fs->nref;i++)
    fs->realnameidx[i] = i;

  fprintf(stderr,"\t-> Number of nref %d in file: \'%s\'\n",fs->nref,fa);
  if(fs->nref==0){
    fprintf(stderr,"\t-> Possible error, no sequences loaded\n");
    exit(1);
  }

  // Set initial probabilities for sampling
  fasta_sampler_setprobs(fs);

  // the initialized fasta_sampler object from which our sequence read is sampled
  return fs;
}

fasta_sampler *fasta_sampler_alloc_subset(const char *fa,const char *SpecificChr){
  /*
  fasta_sampler_alloc_subset - Allocates and initializes a `fasta_sampler` object based on a desired subset of the contigs present within the input reference FASTA file (-i).
  
  @param fa: A string representing the path to the FASTA file from which sequences will be loaded.
  @param SpecificChr: A comma-separated string of chromosome or contig names to be included in the subset.
  
  similar to fasta_sampler_alloc_full function
  */

  // initialize vector to store subset which is extracted from comma seperated list by tokenizing the string
  std::vector<char *> SubsetChr;
  SubsetChr.push_back(strtok(strdup(SpecificChr),"\", \t"));
  char *chrtok = NULL;
  while(((chrtok=strtok(NULL,"\", \t"))))
    SubsetChr.push_back(strdup(chrtok));
  fprintf(stderr,"\t-> Number of entries in chromosome list: %lu\n",SubsetChr.size());

  // Allocate memory for the fasta_sampler object and load index similar to fasta_sampler_alloc_full function
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
  // Fetch and store sequence names, lengths, and actual sequences for the desired subset defined in the vector
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

  // Set initial probabilities for sampling
  fasta_sampler_setprobs(fs);

  // the fasta_sampler object is succesfully created so free memory for the names of subset contigs in the vector
  for(int i=0;i<(int)SubsetChr.size();i++)
    free(SubsetChr[i]);
  
  // the initialized fasta_sampler object with desired subset from which our sequence read is sampled
  return fs;
}

fasta_sampler *fasta_sampler_alloc_bedentry(const char *fa,const char *bedfilename,size_t flanking){
  /*
  fasta_sampler_alloc_bedentry - Allocates and initializes a `fasta_sampler` object based on capture simulations extracting regions of interest from provided bed files from an input reference FASTA file (-i).
  
  @param fa: A string representing the path to the FASTA file from which sequences will be loaded.
  @param bedfilename: A string representing the path to the BED file containing the regions of interest
  @param flanking: A size_t value specifying the number of bases to include as flanking regions around each BED entry.

  similar to fasta_sampler_alloc_full function
  */

  // Allocate memory for the fasta_sampler object and load index
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

  //sort the coordinates of the bed entries
  sortBedEntries(BedEntries, BedEntryCount);

  //merge the bed entries with overlapping regions of interest
  int mergedCount = 0;
  BedEntry* mergedEntries = mergeOverlappingRegions(BedEntries, BedEntryCount, &mergedCount);
  free(BedEntries);

  //Check if merged regions of interest is present within provided reference file
  fs->BedReferenceEntries = checkbedentriesfasta(fs,mergedEntries,mergedCount,&fs->BedReferenceCount);

  //for (int i = 0; i < fs->BedReferenceCount; i++){fprintf(stderr,"merge bed file information bed entry \t%s\t%d\t%d\n", fs->BedReferenceEntries[i].chromosome, fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);}

  free(mergedEntries);

  // Initialize the number of references and allocate memory for sequence data,
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
  
    // Adjust start and end positions with flanking regions, ensuring they are within sequence bounds

    if((int)(fs->BedReferenceEntries[i].start-flanking) < 0){
      //Case 1: region of interest with flanking region to start coordinate extends outside chromosome start
      bedflankstart = 1;
    }
    else if((int)(fs->BedReferenceEntries[i].start-flanking) > 0){
      //Case 2: region of interest has flanking nucleotides substracted from start coordinate
      bedflankstart = fs->BedReferenceEntries[i].start-flanking; 
    }
    
    if((fs->BedReferenceEntries[i].end+flanking) < refchr_length){
      //Case 3: region of interest with flanking nucleotides added to end coordinate
      bedflankend = fs->BedReferenceEntries[i].end + flanking;
    }
    if((fs->BedReferenceEntries[i].end+flanking) > refchr_length){
      //Case 4: region of interest with flanking nucleotides added to end coordinate extends outside chromosome end
      bedflankend = refchr_length-1;
    }

    // addjust merge bed entries coordinates according to the new flanking coordinates
    fs->BedReferenceEntries[i].start = bedflankstart;
    fs->BedReferenceEntries[i].end = bedflankend;

    // Create a string representation of the sequence region for fetching to match string when using samtools faidx chr1:10000-100100
    int length = snprintf(NULL, 0, "%s:%d-%d", fs->BedReferenceEntries[i].chromosome, bedflankstart, bedflankend);
    fs->seqs_names[i] = (char*) malloc((length + 1) * sizeof(char));
    sprintf(fs->seqs_names[i], "%s:%d-%d", fs->BedReferenceEntries[i].chromosome, bedflankstart, bedflankend);
  }

  // Fetch and store sequence names, lengths, and actual sequences for regions of interest
  for(int i=0;i<fs->nref;i++){
    fs->char2idx[fs->seqs_names[i]] = i;
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

  // the initialized fasta_sampler object with regions of interest from which our sequence read is sampled
  return fs;
}

fasta_sampler *fasta_sampler_alloc_maskbedentry(const char *fa,const char *bedfilename,size_t flanking){
  /*
  fasta_sampler_alloc_maskbedentry - Allocates and initializes a `fasta_sampler` object based on capture simulations excluding genomic regions from a provided bed files from an input reference FASTA file (-i).
  
  @param fa: A string representing the path to the FASTA file from which sequences will be loaded.
  @param bedfilename: A string representing the path to the BED file containing the regions to be excluded
  @param flanking: A size_t value specifying the number of bases to include as flanking regions around each BED entry.

  similar to fasta_sampler_alloc_bedentry function
  */

  // Allocate memory for the fasta_sampler object and load index
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

  //sort the coordinates
  sortBedEntries(BedEntries, BedEntryCount);
  //for (int i = 0; i < BedEntryCount; i++){fprintf(stderr,"sort bed file information bed entry \t%s\t%d\t%d\n", BedEntries[i].chromosome, BedEntries[i].start, BedEntries[i].end);}
  
  //merge the bed entries with overlapping regions of interest
  int mergedCount = 0;
  BedEntry* mergedEntries = mergeOverlappingRegions(BedEntries, BedEntryCount, &mergedCount);
  //for (int i = 0; i < mergedCount; i++){fprintf(stderr,"merge bed file information bed entry \t%s\t%d\t%d\n", mergedEntries[i].chromosome, mergedEntries[i].start, mergedEntries[i].end);}
  free(BedEntries);

  int ReferenceCount = 0;

  //Exclude (mask) merged genomic regions from the provided reference file
  fs->BedReferenceEntries = maskbedentriesfasta(fs,mergedEntries,mergedCount,&fs->BedReferenceCount);
  
  //fprintf(stderr,"\t-> Input bed file had %d regions, after merging the overlapping regions there is %d and post-filtering there is %d\n",BedEntryCount,mergedCount,fs->BedReferenceCount);
  //for (int i = 0; i < fs->BedReferenceCount; i++){fprintf(stderr,"merge bed file information mask bed entry \t%s\t%d\t%d\n", fs->BedReferenceEntries[i].chromosome, fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);}

  free(mergedEntries);

  // Initialize the number of references and allocate memory for sequence data,
  fs->nref = fs->BedReferenceCount; 
  fs->seqs = new char* [fs->nref];
  fs->seqs_l = new int[fs->nref];
  fs->seqs_names = new char *[fs->nref];
  fs->realnameidx = new int[fs->nref];

  //create the names to extract the region outside of the coordinates defined in the bed entry
  for(int i=0;i<fs->nref;i++){
    int length = snprintf(NULL, 0, "%s:%d-%d", fs->BedReferenceEntries[i].chromosome,  fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);
    fs->seqs_names[i] = (char*) malloc((length + 1) * sizeof(char));
    sprintf(fs->seqs_names[i], "%s:%d-%d", fs->BedReferenceEntries[i].chromosome, fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);
  }

  // Fetch and store sequence names, lengths, and actual sequences for regions of interest
  for(int i=0;i<fs->nref;i++){
    fs->char2idx[fs->seqs_names[i]] = i;
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

  // the initialized fasta_sampler object without the provided genomic regions from which our sequence read is sampled
  return fs;
}

fasta_sampler *fasta_sampler_alloc_vcf(const char *fa, const char *bcffilename,int id,const char* Name,size_t flanking){
   /*
  fasta_sampler_alloc_vcf - Allocates and initializes a `fasta_sampler` object from an input reference FASTA file (-i) based on capture simulations 
    extracting regions of interest in close proximity to the variant positions from a provided vcf
  
  @param fa: A string representing the path to the FASTA file from which sequences will be loaded.
  @param bcffilename: A string representing the path to the VCF file containing the variant information.
  @param id: An integer representing the identifier of an individual defined in the VCF header.
  @param Name: A string representing the name of an individual defined in the VCF header.
  @param flanking: A size_t value specifying the number of bases to include as flanking regions around the position in each VCF entry.
 
  similar to fasta_sampler_alloc_bedentry function
  */
  
  // Allocate memory for the fasta_sampler object and load index
  fasta_sampler *fs = new fasta_sampler;
  fs->ws = NULL;
  fs->fai = NULL;
  assert(((fs->fai =fai_load(fa)))!=NULL);
  fs->BedReferenceEntries = NULL;
  fs->BedReferenceCount = 0;
  int ploidy = 5;
  
  //Convert positions of vcf entries into a bed format with the provided flanking regions
  fs->BedReferenceEntries = vcftobedentries(bcffilename,id,Name,flanking, &fs->BedReferenceCount,&ploidy);

  // Initialize the number of references and allocate memory for sequence data,
  fs->nref = fs->BedReferenceCount*ploidy; 
  fs->seqs = new char* [fs->nref];
  fs->seqs_l = new int[fs->nref];
  fs->seqs_names = new char *[fs->nref];
  fs->realnameidx = new int[fs->nref];

  int nref_entry=0;

  //create the names to extract the region obtained from variant positions using fai_fetch by also creating internal copies of the contigs depending on the ploidy
  for(int i = 0; i < fs->BedReferenceCount; i++) {
    for(int j = 0; j < ploidy; j++) {
      int length = snprintf(NULL, 0, "%s:%d-%d", fs->BedReferenceEntries[i].chromosome, fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);
      fs->seqs_names[nref_entry] = (char*) malloc((length + 1) * sizeof(char));
      if (fs->seqs_names[nref_entry] == NULL) {
        fprintf(stderr, "Error allocating memory for seqs_names\n");
        exit(1);
      }
      memset(fs->seqs_names[nref_entry], 0, (length + 1) * sizeof(char));
      
      //sprintf(fs->seqs_names[nref_entry], "%s:%d-%d", fs->BedReferenceEntries[i].chromosome, tmp_start, fs->BedReferenceEntries[i].end);

      char chr_reg_tmp[128];
      snprintf(chr_reg_tmp,sizeof(chr_reg_tmp),"%s:%d-%d", fs->BedReferenceEntries[i].chromosome, fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);      

      fs->seqs[nref_entry] = fai_fetch(fs->fai,chr_reg_tmp,fs->seqs_l+nref_entry);

      if (fs->seqs[nref_entry] == NULL) {
          fprintf(stderr, "Error fetching sequence from fai\n");
          exit(1);
      }

      // Modify the fasta_sampler contig sequence to include variant
      fs->seqs[nref_entry][fs->BedReferenceEntries[i].overlappositions[0]-fs->BedReferenceEntries[i].start] = *fs->BedReferenceEntries[i].variants[j];

      // Update sequence name to reflect the allele and region
      int new_length = snprintf(NULL, 0, "%sallele%d:%d-%d", fs->BedReferenceEntries[i].chromosome, j+1, fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);
      fs->seqs_names[nref_entry] = (char*) realloc(fs->seqs_names[nref_entry], (new_length + 1) * sizeof(char));
      fs->char2idx[fs->seqs_names[nref_entry]] = nref_entry;

      snprintf(fs->seqs_names[nref_entry], new_length + 1, "%sallele%d:%d-%d", fs->BedReferenceEntries[i].chromosome, j+1, fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);

      nref_entry++;
    }
  }

  for(int i=0;i<fs->nref;i++){
    fs->realnameidx[i] = i;
  }
  if(fs->nref==0){
    exit(1);
  }
  fasta_sampler_setprobs(fs);

  // the initialized fasta_sampler object with regions solely overlapping the desired variant entries within a vcf file, each contig contains only one vcf entry

  return fs;
}

fasta_sampler *fasta_sampler_alloc_vcf_LD(const char *fa, const char *bcffilename,int id,const char* Name,size_t flanking){
  /*
  fasta_sampler_alloc_vcf_LD - Allocates and initializes a `fasta_sampler` object from an input reference FASTA file (-i) based on capture simulations 
    extracting regions of interest in close proximity to the variant positions from a provided vcf, assuming those regions with variants which overlap is in Linkage Disequilibrium
    such that the copies of the contigs can contain numerous variants
  
  @param fa: A string representing the path to the FASTA file from which sequences will be loaded.
  @param bcffilename: A string representing the path to the VCF file containing the variant information.
  @param id: An integer representing the identifier of an individual defined in the VCF header.
  @param Name: A string representing the name of an individual defined in the VCF header.
  @param flanking: A size_t value specifying the number of bases to include as flanking regions around the position in each VCF entry.
 
  similar to fasta_sampler_alloc_vcf function
  */

  // Allocate memory for the fasta_sampler object and load index
  fasta_sampler *fs = new fasta_sampler;
  fs->ws = NULL;
  fs->fai = NULL;
  assert(((fs->fai =fai_load(fa)))!=NULL);
  fs->BedReferenceEntries = NULL;
  fs->BedReferenceCount = 0;
  int ploidy = 5;

  int BedEntryCount = 0;
  BedEntry* BedEntries = NULL;

  //Convert positions of vcf entries into a bed format with the provided flanking regions
  BedEntries = vcftobedentries(bcffilename,id,Name,flanking,&BedEntryCount,&ploidy);

  //merge those regions after flanking nucleotides has been added to variant position - assuming variants in LD
  int mergedCount = 0;
  fs->BedReferenceEntries = VCFLinkageDisequilibrium(BedEntries,BedEntryCount, &fs->BedReferenceCount,ploidy,flanking);
  
  // Free memory created from the vcf file
  for(int i = 0; i<BedEntryCount;i++){
    for(int j = 0; j < BedEntries[i].ploidy; j++){
      free(BedEntries[i].variants[j]);
    }
    free(BedEntries[i].overlappositions);
    free(BedEntries[i].variants);
  }
  free(BedEntries);

  // Initialize the number of references and allocate memory for sequence data,
  fs->nref = fs->BedReferenceCount*ploidy; 
  //fprintf(stderr,"total entries %d\n",fs->nref);
  fs->seqs = new char* [fs->nref];
  fs->seqs_l = new int[fs->nref];
  fs->seqs_names = new char *[fs->nref];
  fs->realnameidx = new int[fs->nref];

  int nref_entry=0;

  //create the names to extract the region obtained from variant positions using fai_fetch by also creating internal copies of the contigs depending on the ploidy
  for(int i = 0; i < fs->BedReferenceCount; i++) {
    for(int j = 0; j < ploidy; j++) {

      int length = snprintf(NULL, 0, "%s:%d-%d", fs->BedReferenceEntries[i].chromosome,fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);
      fs->seqs_names[nref_entry] = (char*) malloc((length + 1) * sizeof(char));
      if (fs->seqs_names[nref_entry] == NULL) {
        fprintf(stderr, "Error allocating memory for seqs_names\n");
        exit(1);
      }
      memset(fs->seqs_names[nref_entry], 0, (length + 1) * sizeof(char));
      //sprintf(fs->seqs_names[nref_entry], "%s:%d-%d", fs->BedReferenceEntries[i].chromosome, tmp_start, fs->BedReferenceEntries[i].end);
      char chr_reg_tmp[128];
      snprintf(chr_reg_tmp,sizeof(chr_reg_tmp),"%s:%d-%d", fs->BedReferenceEntries[i].chromosome,fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);

      fs->seqs[nref_entry] = fai_fetch(fs->fai,chr_reg_tmp,fs->seqs_l+nref_entry);
      
      if (fs->seqs[nref_entry] == NULL) {
        fprintf(stderr, "Error fetching sequence from fai\n");
        exit(1);
      }
      //fprintf(stderr,"fs->seqs_names[i] %d \t %s\n",i,chr_reg_tmp);
      
      for (int v=0;v<(fs->BedReferenceEntries[i].variantCount/ploidy);v++){
        // adding the variations stored within each merged bed entry
        int var_tmp_idx = v*ploidy+j;
        int var_pos_tmp = fs->BedReferenceEntries[i].overlappositions[var_tmp_idx] - fs->BedReferenceEntries[i].start+1;
        
        // Modify the fasta_sampler contig sequence to include variant
        if (var_pos_tmp >= 0 && var_pos_tmp < fs->seqs_l[nref_entry]) {
          fs->seqs[nref_entry][var_pos_tmp] = *fs->BedReferenceEntries[i].variants[var_tmp_idx];
        }
        else {
          fprintf(stderr, "Error: var_pos_tmp out of bounds.\n");
        }

        fs->seqs[nref_entry][var_pos_tmp] = *fs->BedReferenceEntries[i].variants[var_tmp_idx];
      }
      fs->seqs[nref_entry][fs->BedReferenceEntries[i].overlappositions[0]-1] = *fs->BedReferenceEntries[i].variants[j];
      
      // Update sequence name to reflect the allele and region
      int new_length = snprintf(NULL, 0, "%sallele%d:%d-%d", fs->BedReferenceEntries[i].chromosome, j+1, fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);
      fs->seqs_names[nref_entry] = (char*) realloc(fs->seqs_names[nref_entry], (new_length + 1) * sizeof(char));
      fs->char2idx[fs->seqs_names[nref_entry]] = nref_entry;

      snprintf(fs->seqs_names[nref_entry], new_length + 1, "%sallele%d:%d-%d", fs->BedReferenceEntries[i].chromosome, j+1, fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);   
      nref_entry++;
    }
  }
  
  // Initialize the realnameidx with indices
  for(int i=0;i<fs->nref;i++){
    fs->realnameidx[i] = i;
  }
  if(fs->nref==0){
    fprintf(stderr,"\t-> Possible error, no sequences loaded\n");
    exit(1);
  }
  fasta_sampler_setprobs(fs);
  
  // the initialized fasta_sampler object with regions solely overlapping the desired variant entries within a vcf file, assuming variants in LD are present within the same contig 

  return fs;
}

//functions returns head of chromosome, with posB, posE, chr_idx and fraglength set accordingly
char *sample(fasta_sampler *fs,mrand_t *mr,char **chromoname,int &chr_idx,int &posB,int &posE,int &fraglength, size_t &chr_end, int simmode){
  /*
  sample - Samples a chromosome randomly within a FASTA sampler and extract a fragment from a random region wiht a specified length.
 
  @param fs: A pointer to the `fasta_sampler` struct containing the sequence data.
  @param mr: A pointer to the random number generator.
  @param chromoname: A pointer to a string where the name of the selected chromosome will be stored.
  @param chr_idx: An integer reference where the index of the selected chromosome will be stored.
  @param posB: An integer reference where the starting position of the fragment will be stored.
  @param posE: An integer reference where the ending position of the fragment will be stored.
  @param fraglength: An integer representing the length of the fragment to sample.
  @param chr_end: A size_t reference where the length of the selected chromosome will be stored.
  @param simmode: An integer indicating the simulation mode (0 for linear simulations and 1 for circular simulations).
   
  @return: A pointer to the sequence of the selected chromosome.
  */

  chr_idx = 0;

  // Randomly select a chromosome index if more than one chromosome is available
  if(fs->nref>1)
    chr_idx = ransampl_draw2(fs->ws,mrand_pop(mr),mrand_pop(mr));

  // Set the chromosome name, its length and the regional positions
  *chromoname = fs->seqs_names[chr_idx];
  chr_end = fs->seqs_l[chr_idx];
  posB = mrand_pop(mr)*fs->seqs_l[chr_idx]; //abs(mrand_pop_long(mr)) % fs->seqs_l[chr_idx];
  posE = posB +fraglength;
  
  // linear simulations, move the end point to ensure no crossing of the breakpoint.
  if(posE>=fs->seqs_l[chr_idx] && simmode == 0){
    posE=fs->seqs_l[chr_idx];
    posB=posE-fraglength;
  }
  
  // extract and return the sequence of the selected chromosome
  char *ret =fs->seqs[chr_idx]; 
  chr_idx = fs->realnameidx[chr_idx];
  return ret;
}

void dump_internal(fasta_sampler *fs,const char* filename){
  /*
  dump_internal - Saves the sequences from a `fasta_sampler` to a BGZF-compressed file.
  
  @param fs: A pointer to the `fasta_sampler` struct containing the sequences.
  @param filename: A string representing the path to the output file where the sequences will be saved.
 */

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
  /*
  fasta_sampler_destroy - cleans up allocated memory for the fasta_sampler object
  
  @param fs: A pointer to the `fasta_sampler` struct containing the sequences.
  */

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

void fasta_sampler_destroy_captureLD(fasta_sampler *fs){
    /*
  fasta_sampler_destroy_captureLD - cleans up allocated memory for the fasta_sampler object created when assuming variants are in LD
  
  @param fs: A pointer to the `fasta_sampler` struct containing the sequences.
  */

  if (fs == NULL) return;

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

  //handle the capture memory leak
  if(fs->BedReferenceEntries[0].variantCount == 0){
    for(int i = 0; i < (fs->nref/fs->BedReferenceEntries[0].ploidy); i++){
      for(int j = 0; j < fs->BedReferenceEntries[i].ploidy; j++){
        free(fs->BedReferenceEntries[i].variants[j]);
      }
      free(fs->BedReferenceEntries[i].overlappositions);
      free(fs->BedReferenceEntries[i].variants);
    }
  }
  else{
    //fprintf(stderr,"entry count %d \t nref %d \t ploidy %d\n",fs->BedReferenceCount,fs->nref,fs->BedReferenceEntries[0].ploidy);
    for(int i = 0; i < fs->BedReferenceCount; i++){
      if (fs->BedReferenceEntries[i].variants){
        for(int j = 0; j < fs->BedReferenceEntries[i].variantCount; j++){ //variantCount
          free(fs->BedReferenceEntries[i].variants[j]);
        }
        free(fs->BedReferenceEntries[i].variants);
      }
      if (fs->BedReferenceEntries[i].overlappositions) {
        free(fs->BedReferenceEntries[i].overlappositions); // Free the positions array
      }
      //fprintf(stderr,"TEST %d \n",fs->BedReferenceEntries[i].variantCount);
    }
  }

  free(fs->BedReferenceEntries);

  //bedentries[variant_number].overlappositions = (int*) malloc(sizeof(int));

  //fs->BedReferenceCount;

  delete fs;
}

int extend_fasta_sampler(fasta_sampler *fs,int fs_chr_idx,int ploidy,const char* sample_name){
  /*
  extend_fasta_sampler - Extends the `fasta_sampler` instance to accommodate additional ploidy levels for a given chromosome. 
    By modifying the object in place, creating internal copies of the chromosome containing the alleles from a provided vcf file 

  @param fs: A pointer to the `fasta_sampler` struct that will be extended.
  @param fs_chr_idx: The index of the chromosome to be extended in the `fasta_sampler`.
  @param ploidy: The desired ploidy level for the chromosome.

 */

  if(ploidy ==1){
    // Ploidy is haploid; check if chromosome is already in the map
    ploidymap::iterator it = fs->pldmap.find(fs_chr_idx);
    if(it!=fs->pldmap.end()){
      // Chromosome is already present, no extension needed
      return 0;
    }

    // Create a new ploidy map entry for haploid chromosomes
    int *pldmap = new int[5];
    pldmap[0]=fs_chr_idx;
    pldmap[1]=pldmap[2]=pldmap[3]=pldmap[4]=-1;
    fs->pldmap[fs_chr_idx] = pldmap;
    return 0;
  }

  // Check if all required ploidy levels are already present
  int isThere  = 1;
  char buf[1024];
  for(int i=1;i<ploidy;i++){
    snprintf(buf,1024,"%s_%s_allele_%d",fs->seqs_names[fs_chr_idx],sample_name,i);
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

    // Copy existing sequences and metadata
    for(int i=0;i<fs->nref;i++){
      seqs[i] = fs->seqs[i];
      seqs_names[i] = fs->seqs_names[i];
      seqs_l[i] =fs->seqs_l[i];
      realnameidx[i] = fs->realnameidx[i];
    }
  
    // Free original fasta_sampler object
    delete [] fs->seqs;
    delete [] fs->seqs_names;
    delete [] fs->seqs_l;
    delete [] fs->realnameidx;

    // Update `fasta_sampler` with new arrays
    fs->seqs = seqs;
    fs->seqs_names = seqs_names;
    fs->seqs_l = seqs_l;
    fs->realnameidx = realnameidx;

    // Create and extend the fasta_sampler object with the new ploidy levels
    ploidymap::iterator it = fs->pldmap.find(fs_chr_idx);
    assert(it==fs->pldmap.end());
    pldmap[0] = fs_chr_idx;
    for(int i=1;i<ploidy;i++) {
      snprintf(buf,1024,"%s_%s_allele_%d",fs->seqs_names[fs_chr_idx],sample_name,i);
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
