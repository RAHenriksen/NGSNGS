#include <map>
#include <htslib/vcf.h>
#include "fasta_sampler.h"
#include "add_variants.h"

int ChromCounterErrMsg = 0;
int *fabcflookup(bcf_hdr_t *bcf_hdr,char2int &fai2idx,int &maxIndex){
  /*
    fabcflookup - Creates a lookup table for chromosome indices between a BCF header and a FASTA index, and compares if chromosome names in the BCF header with those in the FASTA index
  
    @param bcf_hdr: A pointer to the BCF header containing chromosome information.
    @param fai2idx: A reference to a map that associates chromosome names (from the FASTA index) with their corresponding indices.
    @param maxIndex: An integer reference that may be used to track additional information outside the scope of this function.
  
    @return: A pointer to an integer array that serves as a lookup. If no matching chromosomes are found, returns NULL.
  */

  // extract chromosome names from bcf header
  int seqnames_l;
  const char **bcf_chrom_names = bcf_hdr_seqnames(bcf_hdr,&seqnames_l);
  fprintf(stderr,"\t-> Raw number of contigs from bcf: %d\n",seqnames_l);

  std::map<int,int> res; //map storing matching indeces between bcf header and fasta index
  int max_l = -1;
  for (int i=0; i<seqnames_l;i++){
    // Check if the current chromosome name exists in the FASTA index map
    char2int::iterator it = fai2idx.find(bcf_chrom_names[i]);

    if(it==fai2idx.end()){
      // if not found move onto next
      if(ChromCounterErrMsg++<5)
	      fprintf(stderr,"\t-> Chromesome id from bcf: %s does not exists in fastafile this message is printed out %d more\n", bcf_chrom_names[i],5-ChromCounterErrMsg);
      continue;
    }else{
      // else map the bcf index to fasta index
      res[i] = it->second;
      if(i>max_l)
	    max_l = i;
    }
  }
  free(bcf_chrom_names); // free memory
  max_l = max_l +1;
  fprintf(stderr,"\t-> From header in bcf we observe %lu different scaffolds/chromosomes we will allocate lookup table with dim: %d\n",res.size(),max_l);
  if(res.size()==0)
    return NULL; // for no matches

  // max_l is maximum index so initialize lookup table according to highest number, and set all values to -1  
  int *lookup = new int[max_l];
  for(int i=0;i<max_l;i++)
    lookup[i] = -1;

  // insert mapped indeces
  for(std::map<int,int>::iterator it=res.begin();it!=res.end();it++)
    lookup[it->first] = it->second;
  
  maxIndex = max_l;
  return lookup;
}

/*
  Assumption, input map is sorted according to chr and position.
  reset is set to one when change of chromosome is detected. This will flush remainder of chrosomes.
  Everytime a new indel is found, program will write all bp from previous event until (including) the current position.
  event is then updated by incrementing reference bp for indel. 
*/

void add_indels_simple(fasta_sampler *fs,bcfmap &mybcfmap,bcf_hdr_t *hdr,int ploidy){
  /*
    add_indels_simple - Inserts simple indels (insertions and deletions) into a set of sequences based on BCF data.

    Simple Insertions
    #CHROM	POS	ID	REF	ALT
    MT	10	.	T	TGTCC

    Simple deletions
    #CHROM	POS	ID	REF	ALT
    MT	10	.	TCTATC	T


    @param fs: A pointer to the `fasta_sampler` structure, which contains the sequences of each chromosome present within the mapped indeces.
    @param mybcfmap: A map that associates BCF data with fasta chromosomes.
    @param hdr: A pointer to the BCF header, which provides context for the BCF data.
    @param ploidy: The ploidy level of the organism.
  */

  int last[ploidy];     //the last processed position
  int last_chrid = -1;  // last processed chromosome ID
  int reset = 1;        //flag to reset when to process next chromosome
  int maxsize = -1;     // 
  int *fsoffsets = NULL;//ofsett into fasta_sampler structure

  // Determine the maximum sequence length in the `fasta_sampler` structure
  for(int i=0;i<fs->nref;i++)
    if(fs->seqs_l[i]>maxsize)
      maxsize= fs->seqs_l[i];

  maxsize += 1000; // Increase the maximum size to accommodate indels

  // Allocate memory for storing indels for each chromosome set (ploidy)
  char **indels =new char*[ploidy];
  for(int i=0;i<ploidy;i++)
    indels[i] =(char*) calloc(maxsize,sizeof(char));

  #if 0  
  // Debugging code to print BCF map entries (currently disabled)
  for(bcfmap::iterator it=mybcfmap.begin();0&&it!=mybcfmap.end();it++){
    fprintf(stderr,"%d %d\n",it->first.rid,it->first.pos);
  }
  #endif

  // iterate over each mapped bcf chromosome, position and genotype to process indels
  for(bcfmap::iterator it=mybcfmap.begin();it!=mybcfmap.end();it++){
    #if 0
    for(int i=0;0&&i<ploidy;i++)
      fprintf(stderr,"[%s] gt[%d]: %d\n",__FUNCTION__,i,it->first.gt[i]);
    #endif
    
    // check if processing new chromosome
    if(it->first.rid!=last_chrid){
      last_chrid = it->first.rid;
      reset = 1;
    }
    
    if(reset){
      // process previous chromosome from final position
      if(fsoffsets!=NULL){
        for(int i=0;i<ploidy;i++){
          // Append sequence data after the last processed position
          strcat(indels[i],fs->seqs[fsoffsets[i]]+last[i]);
          fs->seqs_l[fsoffsets[i]] = strlen(indels[i]);
          free(fs->seqs[i]);
          fs->seqs[i] = indels[i]; // replace fasta reference sequence with updated indel sequence
        }
      }

      // Find the sequence offsets for the new chromosome
      ploidymap::iterator it = fs->pldmap.find(last_chrid);
      assert(it!=fs->pldmap.end());
      fsoffsets = it->second;
      for(int i=0;i<ploidy;i++)
        last[i] = 0; // Reset the last processed position
        
      reset = 0;
    }
    #if 0
    //Debugging code to print the last processed positions 
    for(int i=0;i<ploidy;i++)
      fprintf(stderr,"last[%d]: %d\n",i,last[i]);
    #endif
    
    bcf1_t *brec = it->second; // Retrieve the BCF record for the current position
    bcf_unpack(brec, BCF_UN_ALL); // extract all information stored within bcf structure

    //loop over each chromosome and extract allele information (determined by ploidy)
    for(int i=0;i<ploidy;i++) {
      // Update reference sequences with the modified indels
      if(last[i]>it->first.pos){
        last[i] = it->first.pos;
      }

      int nitems2copy = brec->pos-last[i]+1; //number of bases to copy
      assert(strlen(indels[i])+brec->pos-last[i]< (size_t)maxsize ); // ensure enough buffer for indels
      strncat(indels[i],fs->seqs[fsoffsets[i]]+last[i],nitems2copy); // append contigous block sequence
      
      // retrieve genotype information for both reference and alternative allele
      char *allele = NULL;
      allele = it->second->d.allele[it->first.gt[i]];
      assert(allele!=NULL);
      
      int ref_length = strlen(it->second->d.allele[0]);
      int alt_length = strlen(it->second->d.allele[it->first.gt[i]]);
      int allele_length_diff = strlen(allele)-ref_length;

      #if 0
      // Debugging code to print indel information (currently disabled)
      int isdel = 0;
      int isins = 0;
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

      if(allele_length_diff<0){
        // If deletion, skip the number of bases
        last[i] = brec->pos+abs(allele_length_diff)+1;
      }
      else if(allele_length_diff>0){
        // If insertion, append number of bases 
        assert(strlen(indels[i])+strlen(allele)<strlen(indels[i])+maxsize);
        strncat(indels[i],allele+1,strlen(allele));
        last[i] = brec->pos+1;
      }
      else{
        // no indels move on to next position
        last[i] = brec->pos+1;
      }
    }
    bcf_destroy(brec);
    delete [] it->first.gt;
  }
  
  if(fsoffsets!=NULL){
    // Finalize the sequences for the last processed chromosome
    for(int i=0;i<ploidy;i++){
      
      strcat(indels[i],fs->seqs[fsoffsets[i]]+last[i]);
      fs->seqs_l[fsoffsets[i]] = strlen(indels[i]);
      free(fs->seqs[i]);
      fs->seqs[i] = indels[i];
    }
  }
  delete[] indels;
}

void add_ins_complex(fasta_sampler *fs,bcfmap &mybcfmap,bcf_hdr_t *hdr,int ploidy){
  /*
    add_ins_complex - Inserts complex indels (insertions and deletions) into a set of sequences based on BCF data.

    Multiple insertions of different lengths - useful for STR analysis (artificial examples below not empirically obtained)
    #CHROM	POS	ID	REF	ALT
    MT	1	.	GAT	GATGAT,GATGATGAT,GATGATGATGAT

    or both insertions and deletions of different lengths
    MT	1	.	GATGAT	G,GAT,GATGATGAT

    @param fs: A pointer to the `fasta_sampler` structure, which contains the sequences of each chromosome present within the mapped indeces.
    @param mybcfmap: A map that associates BCF data with fasta chromosomes.
    @param hdr: A pointer to the BCF header, which provides context for the BCF data.
    @param ploidy: The ploidy level of the organism.
  */

  int last[ploidy]; // the last processed position
  int last_chrid = -1; // last processed chromosome ID
  int reset = 1;
  int maxsize = -1; //flag to reset when to process next chromosome
  int *fsoffsets = NULL; //ofsett into fasta_sampler structure
  
  // Determine the maximum sequence length in the `fasta_sampler` structure
  for(int i=0;i<fs->nref;i++)
    if(fs->seqs_l[i]>maxsize)
      maxsize= fs->seqs_l[i];

  maxsize += 1000; // Increase the maximum size to accommodate indels
  
  // Allocate memory for storing indels for each chromosome set (ploidy)
  char **indels =new char*[ploidy];
  for(int i=0;i<ploidy;i++)
    indels[i] =(char*) calloc(maxsize,sizeof(char));

  // iterate over each mapped bcf chromosome, position and genotype to process indels
  for(bcfmap::iterator it=mybcfmap.begin();it!=mybcfmap.end();it++){

    // check if processing new chromosome
    if(it->first.rid!=last_chrid){
      last_chrid = it->first.rid;
      reset = 1;
    }

    if(reset){
      // process previous chromosome from final position
      if(fsoffsets!=NULL){
        for(int i=0;i<ploidy;i++){
          // Update reference sequences with the modified indels
          strcat(indels[i],fs->seqs[fsoffsets[i]]+last[i]);
          fs->seqs_l[fsoffsets[i]] = strlen(indels[i]);
          free(fs->seqs[i]);
          fs->seqs[i] = indels[i];
        }
      }

      ploidymap::iterator it = fs->pldmap.find(last_chrid);
      assert(it!=fs->pldmap.end());
      fsoffsets = it->second;
      for(int i=0;i<ploidy;i++){
	      last[i] = 0;
      }
      
      reset = 0;
    }

    bcf1_t *brec = it->second; // Retrieve the BCF record for the current position
    bcf_unpack(brec, BCF_UN_ALL); // extract all information stored within bcf structure

    //loop over each chromosome and extract allele information (determined by ploidy)
    for(int i=0;i<ploidy;i++) {
      //first copy contigous block from last operator
      if(last[i]>it->first.pos){
        last[i] = it->first.pos;
      }

      // retrieve genotype information for both reference and alternative allele   
      char *allele = NULL;
      allele = it->second->d.allele[it->first.gt[i]];
      assert(allele!=NULL);

      int ref_length = strlen(it->second->d.allele[0]);
      int alt_length = strlen(it->second->d.allele[it->first.gt[i]]);    
      int allele_length_diff = strlen(allele)-ref_length;

      int pos = brec->pos;
      
      int GT_indiv = it->first.gt[i];
      int nitems2copy;
      if(GT_indiv == 0){  
        // If the genotype matches the reference, include the position itself
        nitems2copy = pos - last[i] + 1;
      }
      else{
        // For alternate alleles, do not include the position itself as we alter the sequence directly and overwrite this position
        nitems2copy =  pos - last[i];
      }

      // Copy the contiguous sequence block before the indel
      if (strlen(indels[i]) + nitems2copy < (size_t)maxsize) {
        strncat(indels[i], fs->seqs[fsoffsets[i]] + last[i], nitems2copy);
      }

      if(allele_length_diff < 0){
        // If deletion, skip the number of bases
        last[i] = pos+abs(allele_length_diff);
      }
      else if(allele_length_diff > 0){
        // If insertion, append number of bases 
        if (it->first.gt[i] > 0){
          if (strlen(indels[i]) + alt_length < maxsize) {
            strncat(indels[i], allele, alt_length);
          }
          last[i] = pos + ref_length; // Update the last position to skip over the ref allele
        }
        else{
          // no indels move on to next position
          last[i] = pos + 1;
        }
      }
      else {
        // no indels move on to next position (SNP's)
        last[i] = pos + 1;
      }
    }
    bcf_destroy(brec);
    delete [] it->first.gt;
  }
  if(fsoffsets!=NULL){
    // Finalize the sequences for the last processed chromosome
    for(int i=0;i<ploidy;i++){
      
      strcat(indels[i],fs->seqs[fsoffsets[i]]+last[i]);
      fs->seqs_l[fsoffsets[i]] = strlen(indels[i]);
      free(fs->seqs[i]);
      fs->seqs[i] = indels[i];
    }
  }
  delete[] indels;
}

//fasta sampler struct, index for chromosomenaem, position, the alleles, the genotypes and the ploidy. 
void add_snp(fasta_sampler *fs, int fs_chr_idx,int pos,char **alleles, int32_t *gts, int ploidy,const char* sample_name){
  /*
    add_snp - Adds single nucleotide polymorphisms (SNPs) to sequences based on genotype information.
      
    @param fs: A pointer to the `fasta_sampler` structure containing the chromosome sequences.
    @param fs_chr_idx: The index of the chromosome in the `fasta_sampler` structure.
    @param pos: The position in the sequence where the SNP will be incorporated.
    @param alleles: An array of allele strings where alleles[0] is the reference allele and other indices are alternate alleles.
    @param gts: An array of genotypes where each value represents an allele index.
    @param ploidy: The number of chromosome sets, indicating how many genotypes there are.
  */

  char buf[1024]; // Buffer for creating sequence names with suffixes

  // Iterate over each genotype to add the SNP
  for(int i=0;i<ploidy;i++){
    if(i==0){
      // For haploid organism, update the SNP directly
      fs->seqs[fs_chr_idx][pos] = alleles[bcf_gt_allele(gts[i])][0];
    }
    else{
      // create internal copies of the allele and alter that creating another entry in the fasta_sampler structure representing the reference contig with alternative genotypes
      snprintf(buf,1024,"%s_%s_allele_%d",fs->seqs_names[fs_chr_idx],sample_name,i);
      char2int::iterator it = fs->char2idx.find(buf);
      assert(it!=fs->char2idx.end());

      // Update the nucleotide in the sequence at the position
      fs->seqs[it->second][pos] = alleles[bcf_gt_allele(gts[i])][0];
      
    }
  }
}

//if minus one then ref and alt fields are used, if nonnegative then it is used as offset to which genotype to use from the GT fields
int add_vcf_variants(fasta_sampler *fs,const char *bcffilename,int id,const char* Name){
  /*
    add_vcf_variants - Processes VCF (Variant Call Format) file to add SNPs, insertions, and deletions to the input reference contig sequence creating internal copies to accomodate the alternatives.
  
    @param fs: A pointer to the `fasta_sampler` structure containing the chromosome sequences.
    @param bcffilename: The path to the VCF file containing variant information.
    @param id: The index of the sample in the VCF file. If negative, the sample will be identified by name.
    @param Name: The name of the sample to be used if `id` is negative.
    
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00096	HG00097	HG00099

    -id  | --indiv 0 = HG00096
    -name  | --headername HG00099 = id 2
  */

  if(bcffilename==NULL)
    return 0;

  // store the original fasta file entries, as the number is altered during sequence alteration and duplication
  int OrigFastaEntry = fs->nref;

  // Open VCF file and read header and initialize all information stored within

  htsFile *bcf = bcf_open(bcffilename, "r");
  bcf_hdr_t *bcf_head = bcf_hdr_read(bcf);
  bcf1_t *brec = bcf_init();

  int max_l; //maximum index of matched chromosomes - as of 27-08-2024 it remains unused but has during testing been useful.
  int *bcf_idx_2_fasta_idx = fabcflookup(bcf_head,fs->char2idx,max_l);
  int ret = -1;
  int nsamples = bcf_hdr_nsamples(bcf_head);
  int32_t ngt_arr = 0;     
  int32_t *gt_arr = NULL;
  int ngt;
  int inferred_ploidy =5;//assume max ploidy is five, this will be adjusted when parsing data
  int32_t nodatagt[2]={bcf_gt_unphased(0),bcf_gt_unphased(1)};

  //int32_t nodatagt[2]={bcf_gt_unphased(0),bcf_gt_unphased(1)};
  int32_t *mygt=NULL;
  
  // Determine the sample index based on provided ID or name
  int whichsample = -1;
  const char *sample_name;
  if((int)id < 0){
    for (int i = 0; i < nsamples; i++){
      sample_name = bcf_hdr_int2id(bcf_head, BCF_DT_SAMPLE, i);
      //fprintf(stderr,"sample name %s and provide input %s \n",sample_name,Name);
      if(strcasecmp(sample_name,Name)==0){
        //fprintf(stderr,"provided index is %d with name %s and option name %s\n",i,sample_name,Name);
        whichsample = (int)i;
        break;
      }
    }
  }
  else{
    whichsample = (int)id; //index of individual
    sample_name = bcf_hdr_int2id(bcf_head, BCF_DT_SAMPLE, whichsample);
  }
  
  // Check if the sample index is valid
  if ((int)whichsample < 0 ) {
    fprintf(stderr, "Error: Sample index out of range\n");
    bcf_hdr_destroy(bcf_head);
    bcf_destroy(brec);
    bcf_close(bcf);
    exit(1);
  }
  
  bcfmap simple_indels; // To store simple indels extracted with function void add_indels_simple(fasta_sampler *fs,bcfmap &mybcfmap,bcf_hdr_t *hdr,int ploidy);
  bcfmap complex_insertions; // To store complex indels extracted with function void add_ins_complex(fasta_sampler *fs,bcfmap &mybcfmap,bcf_hdr_t *hdr,int ploidy);

  // read the variants and information present within vcf file
  while(((ret=bcf_read(bcf,bcf_head,brec)))==0){
    bcf_unpack((bcf1_t*)brec, BCF_UN_ALL); // extract info
    int fai_chr = bcf_idx_2_fasta_idx[brec->rid];
    
    if (fai_chr == -1) {
      fprintf(stderr, "Error: Chromosome name not found in the reference\n");
      bcf_hdr_destroy(bcf_head);
      bcf_destroy(brec);
      bcf_close(bcf);
      exit(1);
    }

    // extract genotype information
    ngt = bcf_get_genotypes(bcf_head, brec, &gt_arr, &ngt_arr);
    assert(ngt>0);
    assert((ngt %nsamples)==0);
    inferred_ploidy = ngt/nsamples;

    mygt = gt_arr+inferred_ploidy*whichsample;
    if(bcf_gt_is_missing(mygt[0])){
      fprintf(stderr,"\t-> Genotype is missing for pos: %ld will skip\n",brec->pos+1);
      continue;
    }

    //check if internal copy of chromosomes with alternative alleles exists otherwise extend the fasta_sampler structure them
    extend_fasta_sampler(fs,fai_chr,inferred_ploidy,sample_name);

    bcfkey key;
    key.gt = new int[inferred_ploidy];
    for(int i=0;i<inferred_ploidy;i++){
      int issnp = 0;
      int isindel  = 0;
      key.gt[i] = bcf_gt_allele(mygt[i]);

      int ref_length = strlen(brec->d.allele[0]);
      int alt_length = strlen(brec->d.allele[bcf_gt_allele(mygt[i])]);
      //fprintf(stderr,"GT %d | %d \t Ref %s \t Len %d \t ALT %s \t Len %d \n",key.gt[0],key.gt[1],brec->d.allele[0],ref_length,brec->d.allele[bcf_gt_allele(mygt[i])],alt_length);
      
      if(ref_length == 1 && alt_length == 1){
        //snp
        issnp = 1;
        add_snp(fs,fai_chr,brec->pos,brec->d.allele,mygt,inferred_ploidy,sample_name);
      }
      else if(ref_length > alt_length){
        //deletion
        if(alt_length == 1){
          // simple deletion -> MT	10	.	TCTATC	T
          //fprintf(stderr,"Simple deletion with ref length %d \t alt length %d\n",ref_length,alt_length);
          isindel = 1;
          key.rid=fai_chr;
          key.pos= (int) brec->pos;
          bcf1_t *simpl_del_dup= bcf_dup(brec);
          simple_indels[key] =simpl_del_dup;
        }
      } 
      else if(alt_length > ref_length){
        if(ref_length == 1){
          // simple insertions ->  MT	10	.	T	TGTCC
          //fprintf(stderr,"Simple insertion with ref length %d \t alt length %d\n",ref_length,alt_length);
          isindel = 1;
          key.rid=fai_chr;
          key.pos= (int) brec->pos;
          bcf1_t *simpl_ins_dup= bcf_dup(brec);
          simple_indels[key] =simpl_ins_dup;
        }
        else if(ref_length > 1){
          // complex insertions -> MT	1	.	GAT	GATGAT,GATGATGAT,GATGATGATGAT

          //fprintf(stderr,"Complex insertion with ref length %d \t alt length %d\n",ref_length,alt_length);
          isindel = 1;
          key.rid=fai_chr;
          key.pos= (int) brec->pos;
          bcf1_t *complex_ins_dup= bcf_dup(brec);
          complex_insertions[key] = complex_ins_dup;
        }
      }
    }
  }

  fprintf(stderr,"\t-> Done adding the provided variants from the -vcf\n");
  if(simple_indels.size()>0){
    // add the simple indels
    //fprintf(stderr,"\t-> Found some indels, these will now be added to internal datastructures\n");
    add_indels_simple(fs,simple_indels,bcf_head,inferred_ploidy);
  }
  if(complex_insertions.size()>0){
    // add the complex insertions
    //fprintf(stderr,"\t-> Found some indels, these will now be added to internal datastructures\n");
    add_ins_complex(fs,complex_insertions,bcf_head,inferred_ploidy);
  }

  
  /*
  rename chromosomes
  int new_length = snprintf(NULL, 0, "%sallele%d:%d-%d", fs->BedReferenceEntries[i].chromosome, j+1, fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);
  fs->seqs_names[nref_entry] = (char*) realloc(fs->seqs_names[nref_entry], (new_length + 1) * sizeof(char));
  fs->char2idx[fs->seqs_names[nref_entry]] = nref_entry;

  snprintf(fs->seqs_names[nref_entry], new_length + 1, "%sallele%d:%d-%d", fs->BedReferenceEntries[i].chromosome, j+1, fs->BedReferenceEntries[i].start, fs->BedReferenceEntries[i].end);   
  nref_entry++;
  */
  
  // Update sequence names with sample information
  for(int i = 0; i < OrigFastaEntry;i++){
    char chr_reg_tmp[128];
    snprintf(chr_reg_tmp,sizeof(chr_reg_tmp),"%s",fs->seqs_names[i]);      

    int new_length = snprintf(NULL, 0, "%s_%s_allele_0",fs->seqs_names[i],sample_name);
    //fprintf(stderr,"number of ref %d \t ref name %s_%s_allele_%d \t pos %s\n",i,fs->seqs_names[j],sample_name,j,position);
    fs->seqs_names[i] = (char*) realloc(fs->seqs_names[i], (new_length + 1) * sizeof(char));
    fs->char2idx[fs->seqs_names[i]] = i;
    snprintf(fs->seqs_names[i], new_length + 1, "%s_%s_allele_0",chr_reg_tmp,sample_name);
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

  add_vcf_variants(fs,vcf,1);
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