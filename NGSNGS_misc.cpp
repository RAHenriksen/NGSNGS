#include <algorithm>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <iostream>

#include "NGSNGS_misc.h"
#include <htslib/kstring.h>
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include "fasta_sampler.h"

#define INITIAL_BED_CAPACITY 100

/*
refToInt - A lookup table mapping ASCII values to integer codes for nucleotides.
  The table maps A, C, G, T to 0, 1, 2, 3, respectively, and other characters to 4.
*/
int refToInt[256] = {
  0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//15
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//31
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//47
  0,1,2,3,4,4,4,4,4,4,4,4,4,4,4,4,//63
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//79
  4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//95
  4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,//111
  4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,//127
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//143
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//159
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//175
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//191
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//207
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//223
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,//239
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4//255
};

/*
  intToRef - A lookup table mapping integer codes to nucleotide characters.
             Maps 0, 1, 2, 3 to 'A', 'C', 'G', 'T'.
*/
char intToRef[4] = {'A','C','G','T'};

/*
  intToRef - A lookup table mapping integer codes to nucleotide characters.
             Maps 0, 1, 2, 3 to 'A', 'C', 'G', 'T', respectively.
*/

char NtComp[5] = {'T', 'G', 'C', 'A','N'};

const char *bass = "ACGTN";

void Complement(char* seq) {
  // Complement input sequence in place
  for (; *seq; ++seq) {
    *seq = NtComp[refToInt[(unsigned char)*seq]];
  }
}

void reverseChar(char* str,int length) {
  // Reverse the input sequence in place
  std::reverse(str, str + length);
}

void ReversComplement(char* seq) {
  // Reverse and complement input sequence in place
  reverseChar(seq, strlen(seq));
  Complement(seq);
}

void Complement_k(kstring_t* seq) {
  // Complement input sequence in place with kstring_t structure
  char *s = seq->s;
  for (; *s; ++s) {
    *s = NtComp[refToInt[(unsigned char)*s]];
  }
}

void ReversComplement_k(kstring_t* seq) {
  // Reverse and complement input sequence in place with kstring_t structure
  reverseChar(seq->s, seq->l);
  Complement_k(seq);
}

char* PrintCigarBamSet1(size_t n_cigar,const uint32_t *cigar){
  /*
  PrintCigarBamSet1 - Converts a BAM CIGAR array into a readable CIGAR string.
  
  @param n_cigar: The number of CIGAR operations.
  @param cigar: A pointer to an array of CIGAR operations.
  */
  
  // Allocate a buffer to store the CIGAR string
  char *cigar_str = (char *)malloc(n_cigar * 2); // a safe overestimate
  char *ptr = cigar_str;

  /*
  Bam structure in general:
    size_t l_qname, const char *qname,
    uint16_t flag, int32_t tid, hts_pos_t pos, uint8_t mapq,
    size_t n_cigar, const uint32_t *cigar,
    int32_t mtid, hts_pos_t mpos, hts_pos_t isize,
    size_t l_seq, const char *seq, const char *qual,
    size_t l_aux
  */

  // Convert CIGAR to string
  for (int i = 0; i < n_cigar; ++i) {
    int len = sprintf(ptr, "%d%c", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i]));
    ptr += len;
  }

  return cigar_str;
}

void CreateSeqQualKString(bam1_t *aln, kstring_t *Sequence, kstring_t *Quality, int offset){
  /*
  CreateSeqQualKString - Extracts sequence and quality strings from a BAM alignment and stores them in kstring_t structures.

  @param aln: A pointer to a BAM alignment structure (`bam1_t`) containing sequence and quality data.
  @param Sequence: A pointer to a kstring_t structure where the extracted sequence string will be stored.
  @param Quality: A pointer to a kstring_t structure where the extracted quality string will be stored.
  @param offset: An integer value to adjust the quality scores by adding the offset to each quality score,  33 only for fastq it should be 0 for bam.
  */

  //sequence length
  int len = aln->core.l_qseq;

  // Allocate memory for sequence and quality string +1 for null terminator
  Sequence->s = (char *) realloc(Sequence->s, len + 1); 
  Sequence->l = len;
  Sequence->m = len + 1;

  Quality->s = (char *) realloc(Quality->s, len + 1); 
  Quality->l = len;
  Quality->m = len + 1;

  // Extract raw sequence and quality data from the BAM alignment
  uint8_t *seq_data = bam_get_seq(aln);
  uint8_t *qual_data = bam_get_qual(aln);

  // Convert sequence and quality data to their respective string representations
  for (int i = 0; i < len; i++) {
    Sequence->s[i] = seq_nt16_str[bam_seqi(seq_data, i)]; // convert nucleotide id to IUPAC id
    Quality->s[i] = qual_data[i]+offset;
  }
  
  Sequence->s[len] = '\0'; // null-terminate the string
  Quality->s[len] = '\0';
}

BedEntry* readBedFile(const char* filename, int* entryCount){
  /*
  readBedFile - Reads a BED file and stores each entry in an dynamically allocated array of BedEntry structures containing the BED entries with chromosome, start, and end fields populated. 
  The other fields (variants, overlappositions, variantCount, ploidy) are uninitialized when reading in the BED file and populated later one

  @param filename: The name of the BED file to be read.
  @param entryCount: A pointer to an integer where the function will store the number of entries read from the file.

  */

  // open bed file
  FILE* file = fopen(filename, "r");
  if (!file) {
    // Check if the file opened successfully
    fprintf(stderr,"Error opening bed file");
    exit(1);
  }

  // Initilize the capacity for the BED entries array
  int region_capacity = INITIAL_BED_CAPACITY;
  BedEntry* bedentries = (BedEntry*) malloc(region_capacity * sizeof(BedEntry));
  if (!bedentries) {
    // Check if memory allocation was successful
    fprintf(stderr,"Error allocating memory");
    fclose(file);
    return NULL;
  }

  char line[256]; // Buffer to store each line from the BED file
  int tmpcount = 0; // Counter to track the number of BED entries read
  
  // Read Bed entries from the file
  while (fgets(line, sizeof(line), file)) {
    // dynamically reallocate the entry capacity
    if (tmpcount >= region_capacity) {
      region_capacity *= 2; // Double the capacity
      BedEntry* temp = (BedEntry*) realloc(bedentries, region_capacity * sizeof(BedEntry));
      if (!temp) {
        fprintf(stderr, "Error reallocating memory for bedentries");
        free(bedentries);
        fclose(file);
        return NULL;
      }
      bedentries = temp; //update pointer to the new allocated memory for the entry
    }

    // Parse the line and store in the structure
    sscanf(line, "%s\t%d\t%d", bedentries[tmpcount].chromosome, &bedentries[tmpcount].start, &bedentries[tmpcount].end);
    //fprintf(stderr,"%s\t%d\t%d\t%d\n", bedentries[*entryCount].chromosome, bedentries[*entryCount].start, bedentries[*entryCount].end,*entryCount);
    
    // Initialize unparsed fields to default values (consider if this is necessary)
    bedentries[tmpcount].variants = NULL;
    bedentries[tmpcount].overlappositions = NULL;
    bedentries[tmpcount].variantCount = 0;
    bedentries[tmpcount].ploidy = 0;
    
    tmpcount++;
  }

  *entryCount = tmpcount;
  fclose(file);
  return bedentries;
}

int compareBedEntries(const void* a, const void* b) {
  /*
  compareBedEntries - Function used by qsort to compare two BedEntry structures.

  @param a: Pointer to the first BedEntry structure.
  @param b: Pointer to the second BedEntry structure.
  */

  BedEntry* entryA = (BedEntry*)a;
  BedEntry* entryB = (BedEntry*)b;

  int chromCompare = strcmp(entryA->chromosome, entryB->chromosome);
  if (chromCompare == 0) {
    return (entryA->start - entryB->start);
  }
  else {
    return chromCompare;
  }
}

void sortBedEntries(BedEntry* entries, int entryCount) {
  /*
  sortBedEntries - Sorts (in place) an array of BedEntry structures based on chromosome names and start positions.

  @param entries: A pointer to an array of BedEntry structures to be sorted.
  @param entryCount: The number of entries in the array.
  */
  qsort(entries, entryCount, sizeof(BedEntry), compareBedEntries);
}

BedEntry* mergeOverlappingRegions(BedEntry* entries, int entryCount, int* mergedCount){
  /*
  mergeOverlappingRegions - Merges overlapping BED regions in a sorted array of BedEntry structures.

  @param entries: A pointer to an array of BedEntry structures, assumed to be sorted by chromosome and start position.
  @param entryCount: The number of entries in the array.
  @param mergedCount: A pointer to an integer where the function will store the number of merged entries.
  */
  
  //check empty input
  if (entryCount == 0) {
    *mergedCount = 0;
    return NULL;
  }

  // Allocate memory for the merged entries array
  BedEntry* mergedEntries = (BedEntry*) malloc(entryCount * sizeof(BedEntry));
  if (!mergedEntries) {
    fprintf(stderr,"Error allocating memory for merged entries");
    return NULL;
  }

  int j = 0;
  mergedEntries[j] = entries[0];

  // Iterate through the entries and merge overlapping regions
  for (int i = 1; i < entryCount; i++) {
    // Check for overlap: same chromosome and overlapping regions
    if (strcmp(mergedEntries[j].chromosome, entries[i].chromosome) == 0 && mergedEntries[j].end >= entries[i].start) {
      // Merge overlapping regions by extending the end position if necessary
      if (entries[i].end > mergedEntries[j].end) {
        mergedEntries[j].end = entries[i].end;
      }
    }
    else {
      // No overlap, move to the next entry in the merged array
      j++;
      mergedEntries[j] = entries[i];
    }
  }

  *mergedCount = j + 1;
  return mergedEntries;
}

BedEntry* checkbedentriesfasta(fasta_sampler *fs,BedEntry* entries,int entryCount,int* ReferenceCount) {
  /*
  checkbedentriesfasta - dynamically allocated array of BedEntry structures that match sequences in the reference FASTA.

  @param fs: A pointer to the fasta_sampler structure, which contains the reference genome information (-i)
  @param entries: A pointer to an array of BedEntry structures representing the BED file entries.
  @param entryCount: The number of BED entries to be checked.
  @param ReferenceCount: A pointer to integer that will be updated with the count of BED entries found in the reference genome.
  */
  
  //If no BED entries
  if (entryCount == 0) {
    *ReferenceCount = 0;
    return NULL;
  }
    
  
  /*for(int i=0;i<n_seqs;i++){
    const char *chrname = faidx_iseq(fs->fai,i);
    fprintf(stderr,"nref val %d \t and name %s\n",i,chrname);
  }*/


  // Allocate memory to store the BED entries that are present in the reference genome
  BedEntry* ReferenceEntries = (BedEntry*) malloc(entryCount * sizeof(BedEntry));

  if (!ReferenceEntries) {
    fprintf(stderr, "Error allocating memory for ReferenceEntries");
    *ReferenceCount = 0;
    return NULL;
  }

  int warning = 0;
  int n_seqs = faidx_nseq(fs->fai);
  int reference_int = 0;

  // Iterate over each sequence in the reference FASTA index
  for(int i=0;i<n_seqs;i++){
    const char *chrname = faidx_iseq(fs->fai,i); //chromosome name

    //fprintf(stderr,"nref val %d \t and name %s\n",i,chrname);

    // Iterate through all BED entries to identify potential matches with reference genome
    for (int j = 0; j < entryCount; j++) {
      char region[1024];
      // Format the current BED entry as a genomic region string
      sprintf(region, "%s:%d-%d", entries[j].chromosome,entries[j].start,entries[j].end);
      if (strcmp(chrname, entries[j].chromosome) == 0) {
        //fprintf(stderr,"ref seq %d \t entry val %d \n",i,j);
        //fprintf(stderr,"Bed entry %s \t and name %s \t %d \t %d\n",entries[j].chromosome,chrname,entries[j].start,entries[j].end);

        //Copy matching bed entries to new BedEntry struct - ReferenceEntries
        strcpy(ReferenceEntries[reference_int].chromosome,entries[j].chromosome);
        ReferenceEntries[reference_int].start = entries[j].start;
        ReferenceEntries[reference_int].end = entries[j].end;
        reference_int++; 
      }
      else{
        // Set the warning flag if there is a mismatch
        warning =1 ;
        continue;
      }
    }
  }
  
  if(warning==1)
    fprintf(stderr,"\t-> Discrepancy between input referenDce genome (-i) and input bed file, some regions in bed is not present in reference or some chromosomes in reference is not present in bed file, and these remains unused for further simulation \n");

   // Set the count of valid BED entries found in the reference genome
  *ReferenceCount = reference_int;

  return ReferenceEntries;
}

BedEntry* maskbedentriesfasta(fasta_sampler *fs,BedEntry* entries,int entryCount,int* ReferenceCount) {
  /*
  maskbedentriesfasta - A dynamically allocated array of BedEntry structures representing regions in the reference genome not covered by BED entries.

  @param fs: A pointer to the fasta_sampler structure, which contains the reference genome information (-i)
  @param entries: An array of BedEntry structures representing regions in the BED file.
  @param entryCount: The number of BED entries.
  @param ReferenceCount: A pointer to an integer updated with the count of masked regions found.
  */

  // no bed entries
  if (entryCount == 0) {
    *ReferenceCount = 0;
    return NULL;
  }
  
  // Initial capacity for the array of masked regions
  int region_capacity = INITIAL_BED_CAPACITY;

  int warning = 0;
  int n_seqs = faidx_nseq(fs->fai);
  int reference_int = 0;

  // Allocate memory for masked regions
  BedEntry* ReferenceEntries = (BedEntry*) malloc(region_capacity * sizeof(BedEntry)); //think about how many i actually need for allocation
  if (!ReferenceEntries) {
    fprintf(stderr, "Error allocating memory for ReferenceEntries\n");
    *ReferenceCount = 0;
    return NULL;
  }

  // Iterate over each sequence in the reference FASTA index
  for(int i=0;i<n_seqs;i++){
    const char *chrname = faidx_iseq(fs->fai,i); // Get chromosome name
    size_t seq_len = faidx_seq_len(fs->fai, chrname); // Get length of the sequence

    // initialize potentially masekd region
    int tmp_start = 1;
    int chr_not_present = 0;
    int capture_present_chr_end = 0;

    // Iterate through each BED entry to find regions within the current reference sequence
    for (int j = 0; j < entryCount; j++){
      //fprintf(stderr,"reference %d \t entrycount %d\n",i,j);

      // Check for matches
      if((strcmp(chrname, entries[j].chromosome) == 0)){
        capture_present_chr_end = 1;
        chr_not_present = 1;

        // Check if there is an uncovered region before the current BED entry
        if(tmp_start < seq_len){
          //fprintf(stderr,"Bed entry %s \t and name %s \t %d \t %d\n",entries[j].chromosome,chrname,tmp_start,entries[j].start);

          //ensure enough capacity to create regions from which not to sample
          if (reference_int >= region_capacity) {
            region_capacity *= 2;
            BedEntry* temp =  (BedEntry*) realloc(ReferenceEntries, region_capacity * sizeof(BedEntry));
            if (!temp) {
              fprintf(stderr, "Error reallocating memory for ReferenceEntries\n");
              free(ReferenceEntries);
              *ReferenceCount = 0;
              return NULL;
            }
            ReferenceEntries = temp;
          }

          // Store the coordinates of uncovered region
          strcpy(ReferenceEntries[reference_int].chromosome,entries[j].chromosome);
          ReferenceEntries[reference_int].start = tmp_start;
          ReferenceEntries[reference_int].end = entries[j].start;
          
          tmp_start = entries[j].end;
          reference_int++;
        }
      }
    }

    // Capture the end of the chromosome if no further entries cover it
    if (capture_present_chr_end || !chr_not_present) {
      //fprintf(stderr,"%s\t%d\t%d\n",chrname,tmp_start, seq_len);
      if (reference_int >= region_capacity) {
        region_capacity *= 2;
        BedEntry* temp =  (BedEntry*)realloc(ReferenceEntries, region_capacity * sizeof(BedEntry));
        if (!temp) {
          fprintf(stderr, "Error reallocating memory for ReferenceEntries\n");
          free(ReferenceEntries);
          *ReferenceCount = 0;
          return NULL;
        }
        ReferenceEntries = temp;
      }

      // Store the uncovered region at the end of the chromosome
      strcpy(ReferenceEntries[reference_int].chromosome,chrname);
      ReferenceEntries[reference_int].start = tmp_start;
      ReferenceEntries[reference_int].end = seq_len;
      
      reference_int++;
    }

  }
  
  *ReferenceCount = reference_int;

  return ReferenceEntries;
}

BedEntry* vcftobedentries(const char* bcffilename, int id,const char* Name,size_t flanking, int* entryCount,int* ploidy){
  /*
  vcftobedentries - Reads a BCF/VCF file and extracts variant information to create BED entries with flanking regions in each end, to create pseudo-contigs from which we can sample.

  @param bcffilename: The name of the BCF/VCF file to read.
  @param id: The index of the sample; if negative, the sample name is used instead. (-id 0)
  @param Name: The name of the sample to extract; used if id is negative. (-name HG00096)
  @param flanking: The size of the flanking region to add around each variant position.
  @param entryCount: A pointer to an integer updated with the count of BED entries created.
  @param ploidy: A pointer to an integer with the inferred ploidy of the sample.
  */

  // initialize bed entries
  int region_capacity = INITIAL_BED_CAPACITY;
  BedEntry* bedentries = (BedEntry*) malloc(region_capacity * sizeof(BedEntry));
  if (!bedentries) {
    fprintf(stderr,"Error allocating memory");
    exit(1);
  }

  //open the vcf file and store the header information
  htsFile *bcf = bcf_open(bcffilename, "r");

  if (!bcf) {
    fprintf(stderr, "Error: Unable to open BCF file: %s\n", bcffilename);
    free(bedentries);
    exit(1);
  }

  // initialize and read header
  bcf_hdr_t *bcf_head = bcf_hdr_read(bcf);
  if (!bcf_head) {
    fprintf(stderr, "Error: Unable to read BCF header\n");
    bcf_close(bcf);
    free(bedentries);
    exit(1);
  }

  // Initialize a BCF record structure
  bcf1_t *brec = bcf_init();
  /*
  https://github.com/samtools/htslib/blob/develop/htslib/vcf.h

  typedef struct bcf1_t {
    hts_pos_t pos;  // POS
    hts_pos_t rlen; // length of REF
    int32_t rid;  // CHROM
    float qual;   // QUAL
    uint32_t n_info:16, n_allele:16;
    uint32_t n_fmt:8, n_sample:24;
    kstring_t shared, indiv;
    bcf_dec_t d;            // lazy evaluation: $d is not generated by bcf_read(), but by explicitly calling bcf_unpack()
    int max_unpack;         // Set to BCF_UN_STR, BCF_UN_FLT, or BCF_UN_INFO to boost performance of vcf_parse when some of the fields won't be needed
    int unpacked;           // remember what has been unpacked to allow calling bcf_unpack() repeatedly without redoing the work
    int unpack_size[3];     // the original block size of ID, REF+ALT and FILTER
    int errcode;    // one of BCF_ERR_* codes
  } bcf1_t;
  */

  if (!brec) {
    fprintf(stderr, "Error: Unable to initialize BCF record\n");
    bcf_hdr_destroy(bcf_head);
    bcf_close(bcf);
    free(bedentries);
    exit(1);
  }

  // extracting sample information
  int nsamples = bcf_hdr_nsamples(bcf_head);
  int32_t ngt_arr = 0;     
  int32_t *gt_arr = NULL;
  int ngt;
  int inferred_ploidy =5;//assume max ploidy is five, this will be adjusted when parsing data
  //int32_t nodatagt[2]={bcf_gt_unphased(0),bcf_gt_unphased(1)};
  int32_t *mygt=NULL;
  int seqnames_l;
  int whichsample = -1;
  int variant_number = 0;
  
  if(id < 0){
    //if sample name is provided - identify sample index
    for (int i = 0; i < nsamples; i++){
      const char *sample_name = bcf_hdr_int2id(bcf_head, BCF_DT_SAMPLE, i);
      //fprintf(stderr,"Current sample name is %s \t and index %d \t while provide -name is \n",sample_name,i,Name);
      if(strcasecmp(sample_name,Name)==0){
        whichsample = (int)i;
        break;
      }
    }
  }
  else{
    whichsample = (int)id; //index of individual
  }
    
  if ((int)whichsample < 0 ) {
    fprintf(stderr, "Error: Sample index out of range\n");
    bcf_hdr_destroy(bcf_head);
    bcf_destroy(brec);
    bcf_close(bcf);
    free(bedentries);
    exit(1);
  }
  //fprintf(stderr,"NGSNGS_misc.cpp \t The number of individuals present in the bcf header is %d with the chosen individual being %s\n",nsamples,sample_name);

  //const char **bcf_chrom_names = bcf_hdr_seqnames(bcf_head,&seqnames_l);
  bcf_idpair_t *chr_info = bcf_head->id[BCF_DT_CTG];

  // Iterate through each record in the BCF file to create BedEntry struct

  while (bcf_read(bcf, bcf_head, brec) == 0) {
    bcf_unpack((bcf1_t *)brec, BCF_UN_ALL);

    //extract genotype information
    ngt = bcf_get_genotypes(bcf_head, brec, &gt_arr, &ngt_arr);
    assert(ngt>0);
    assert((ngt %nsamples)==0);
    inferred_ploidy = ngt/nsamples;
      
    mygt = gt_arr+inferred_ploidy*whichsample;
    if(bcf_gt_is_missing(mygt[0])){
      fprintf(stderr,"\t-> Genotype is missing for pos: %ld will skip\n",brec->pos+1);
      continue;
    }

    
    // create the coordinate information in our internal bed struct from which we sample data
    int fai_chr = brec->rid;
    if (fai_chr == -1) {
      fprintf(stderr, "Error: Chromosome name not found in the reference\n");
      bcf_hdr_destroy(bcf_head);
      bcf_destroy(brec);
      bcf_close(bcf);
      free(bedentries);
      exit(1);
    }
    int chr_end = chr_info[fai_chr].val->info[0];

    // Ensure enough memory for BedEntry depending on variants
    if (variant_number >= region_capacity) {
      region_capacity *= 2;
      BedEntry* temp = (BedEntry*) realloc(bedentries, region_capacity * sizeof(BedEntry));

      if (!temp) {
        fprintf(stderr,"Error reallocating memory for bed entries\n");
        free(bedentries);
        bcf_hdr_destroy(bcf_head);
        bcf_destroy(brec);
        bcf_close(bcf);
        exit(1);
      }
      bedentries = temp;
      //free(temp);
    }
    
    //initialize the variables for the struct despite remaining unused for the vcf capture
    bedentries[variant_number].variantCount = 0; 
    bedentries[variant_number].variants = (char**) malloc(inferred_ploidy * sizeof(char*));
    
    if (!bedentries[variant_number].variants){
      fprintf(stderr, "Error allocating memory for variants\n");
      free(bedentries[variant_number].overlappositions);
      free(bedentries);
      bcf_hdr_destroy(bcf_head);
      bcf_destroy(brec);
      bcf_close(bcf);
      exit(1);
    }

    int n = 1;
    bedentries[variant_number].overlappositions = (int*) malloc(n * sizeof(int));
      
    if (!bedentries[variant_number].overlappositions) {
      fprintf(stderr, "Error allocating memory for overlappositions\n");
      free(bedentries);
      bcf_hdr_destroy(bcf_head);
      bcf_destroy(brec);
      bcf_close(bcf);
      exit(1);
    }

    bedentries[variant_number].chromosome[0] = '\0'; // Ensure it's empty initially

    //Assign initial values to struct
    bedentries[variant_number].overlappositions[0] = (int)(brec->pos +1);
    strcpy(bedentries[variant_number].chromosome,bcf_hdr_id2name(bcf_head, brec->rid));

    //sscanf(line, "%s %d %d", bedentries[*entryCount].chromosome, &bedentries[*entryCount].start, &bedentries[*entryCount].end);  
    int tmp_start = (size_t) brec->pos - (size_t) flanking; //size_t is unsigned integer thus non-negative values
    size_t tmp_end = (size_t) brec->pos + (size_t) flanking; //size_t is unsigned integer thus non-negative values

    if(tmp_start < 1){
      // if flanking extend beyond chromosome start coordinate 
      bedentries[variant_number].start = 1;
      bedentries[variant_number].end = brec->pos + flanking;
    }
    else if(tmp_end > chr_end){
      // if flanking extend beyond chromosome end coordinate  
      bedentries[variant_number].start = brec->pos - flanking;
      bedentries[variant_number].end = chr_end;
    }
    else{
      bedentries[variant_number].start = brec->pos - flanking;
      bedentries[variant_number].end = brec->pos + flanking;
    }

    /*
    fprintf(stderr,"flanking start %d \t flanking end %d \t variation %d \n",bedentries[variant_number].start,bedentries[variant_number].end,bedentries[variant_number].overlappositions[0]);
    fprintf(stderr,"postition %d \t allele 1 is %s \t allele 2 is %s \n",brec->pos+1,brec->d.allele[bcf_gt_allele(mygt[0])],brec->d.allele[bcf_gt_allele(mygt[1])]);
    */

    //bedentries[variant_number].variants = new char *[inferred_ploidy];

    // Iterate through alleles from the inferred ploidy
    for (int v = 0; v < inferred_ploidy; v++){
      bedentries[variant_number].variants[v] = NULL;

      const char* allele_str = brec->d.allele[bcf_gt_allele(mygt[v])];

      if (allele_str != NULL) {
        bedentries[variant_number].variants[v] = strdup(allele_str); // Duplicate the variant
        if (bedentries[variant_number].variants[v] == NULL) {
          fprintf(stderr, "Error duplicating allele string\n");
          free(bedentries[variant_number].overlappositions);
          free(bedentries[variant_number].variants);
          free(bedentries);
          bcf_hdr_destroy(bcf_head);
          bcf_destroy(brec);
          bcf_close(bcf);
          exit(1);
        }
      }

      //fprintf(stderr,"NGSNGS_misc.cpp \t bed entry allele is %d \t %d \t %d \t %d \t %s \n",variant_number,v,bedentries[variant_number].start,bedentries[variant_number].end,bedentries[variant_number].variants[v]);
      //brec->d.allele[bcf_gt_allele(mygt[v])]
    }
    //fprintf(stderr,"postition %d \t allele 1 is %s \t allele 2 is %s \n",brec->pos+1,brec->d.allele[bcf_gt_allele(mygt[0])],brec->d.allele[bcf_gt_allele(mygt[1])]);
    bedentries[variant_number].ploidy = inferred_ploidy;
    variant_number++;
  }

  *ploidy = inferred_ploidy;
  *entryCount = variant_number;

  bcf_hdr_destroy(bcf_head);
  bcf_destroy(brec);
  bcf_close(bcf);

  return bedentries;
}

void addVariant(BedEntry* entry, const char* variant,int position) {
  /*
  addVariant - Adds a new variant to its position in the BedEntry struct.

  @param entry: Pointer to the BedEntry to which the variant will be added.
  @param variant: The allele to be added.
  @param position: The genomic position extracted from the bcf file for the variant.
  */
  
  // Reallocate memory for the list of variants
  char** newVariants = (char**) realloc(entry->variants, (entry->variantCount + 1) * sizeof(char*));
  int* newPositions = (int*) realloc(entry->overlappositions, (entry->variantCount + 1) * sizeof(int));

  if (!newVariants || !newPositions) {
    fprintf(stderr, "Error reallocating memory for variants or positions\n");
    
    // If either allocation failed, avoid memory leaks by ensuring consistent state
    if (newVariants) {
      entry->variants = newVariants; // Ensure old memory is not lost
    }
    if (newPositions) {
      entry->overlappositions = newPositions; // Ensure old memory is not lost
    }
    return;
  }

  // Update entries after reallocation
  entry->variants = newVariants;
  entry->overlappositions = newPositions;

  // Add the new variant to the array
  entry->variants[entry->variantCount] = strdup(variant);
  if (!entry->variants[entry->variantCount]) {
    fprintf(stderr, "Error allocating memory for variant string\n");
    return;
  }
  
  // Store the position in the overlap positions array
  entry->overlappositions[entry->variantCount] = position; // Store the position

  // Update the variant count
  entry->variantCount++;
}

BedEntry* VCFLinkageDisequilibrium(BedEntry* entries, int entryCount, int* mergedCount,int ploidy,size_t flanking) {
  /*
  VCFLinkageDisequilibrium - Merges overlapping BedEntry regions with inserted variants to mimic variants in LD inherited together.

  @param entries: Pointer to an array of BedEntry structures representing VCF regions.
  @param entryCount: The number of entries in the input array.
  @param mergedCount: Pointer to an integer to store the count of merged entries.
  @param ploidy: The ploidy level of the variants.
  @param flanking: The flanking size added to positions when adding variants.
  */
  
  // Check for empty input
  if (entryCount == 0) {
    *mergedCount = 0;
    return NULL;
  }

  // Allocate memory for the merged entries array
  BedEntry* mergedEntries = (BedEntry*) malloc(entryCount * sizeof(BedEntry));
  if (!mergedEntries) {
    fprintf(stderr,"Error allocating memory for merged entries");
    return NULL;
  }
  
  // Initialize the first entry in the merged Entries
  int j = 0;
  mergedEntries[j] = entries[0];
  mergedEntries[j].variantCount = 0;
  mergedEntries[j].variants = NULL;
  mergedEntries[j].overlappositions = NULL;

  // Add initial variants from the first entry
  for (int k = 0; k < ploidy; k++) {
    addVariant(&mergedEntries[j], entries[0].variants[k],entries[0].start+flanking);
  }

  // Iterate through the entries and merge where applicable
  for (int i = 1; i < entryCount; i++) {

    // Check if the current merged entry overlaps with previous
    if (strcmp(mergedEntries[j].chromosome, entries[i].chromosome) == 0 && mergedEntries[j].end >= entries[i].start) {
      // Merge regions which overlaps and extend end coordinate to accomodate both regions
      if (entries[i].end > mergedEntries[j].end) {
        mergedEntries[j].end = entries[i].end;
      }
      // Add variants
      for (int k = 0; k < ploidy; k++) {
        addVariant(&mergedEntries[j], entries[i].variants[k],entries[i].start+flanking);
      }
    }
    else {
      // No overlap, move to the next entry
      j++;
      mergedEntries[j] = entries[i];
      mergedEntries[j].variantCount = 0;
      mergedEntries[j].variants = NULL;
      mergedEntries[j].overlappositions = NULL;

      // Add variants for the new non-overlapping entry
      for (int k = 0; k < ploidy; k++) {
        addVariant(&mergedEntries[j], entries[i].variants[k], entries[i].start+flanking);
      }
    }
  }

  *mergedCount = j + 1;
  return mergedEntries;
}

#ifdef __WITH_MAIN__
//g++ NGSNGS_misc.cpp -std=c++11 -lz -lm -lbz2 -llzma -lpthread -lcurl -lcrypto ../htslib/libhts.a -D __WITH_MAIN__ -o ngsngsmisc

int main(){
  // Example usage of kstring_t and Complement_k
  kstring_t seq;
  char sequence[] = "ACGTN";
  seq.s = sequence;
  seq.l = strlen(sequence);
  seq.m = seq.l + 1; // Just for safety, usually seq.m is set to the allocated size
  std::cout << "Original sequence: " << seq.s << std::endl;

  Complement_k(&seq);

  fprintf(stderr,"----------------------\n");
  std::cout << "Complemented sequence: " << seq.s << std::endl;
  
  ReversComplement_k(&seq);
  fprintf(stderr,"----------------------\n");

  std::cout << "Reverse complemented sequence: " << seq.s << std::endl;

  reverseChar(seq.s, seq.l);
  fprintf(stderr,"----------------------\n");

  std::cout << "Reverse sequence: " << seq.s << std::endl;

  const char* bcffilename = "Test_Examples/ChrMtSubSNPDiploid.vcf";
  int range = 30;
  int sample_index = 0;
  
  BedEntry* bedentries = NULL;
  int entryCount = 0;
  fprintf(stderr,"entrycount before vcftobed %d \n",entryCount);
  int result = VCFtoBED(bcffilename, sample_index, range, &bedentries, &entryCount);
  fprintf(stderr,"entrycount after vcftobed %d \n",entryCount);

  for (int i = 0; i < entryCount; i++) {
    fprintf(stderr,"vcf information as bed entry \t%s\t%d\t%d\n", bedentries[i].chromosome, bedentries[i].start, bedentries[i].end);
  }
  free(bedentries);
  fprintf(stderr,"----------------------\n");

  int entryCount2 = 0;
  const char* bedfilename = "Coord.bed";
  BedEntry* bedentries2 = NULL;

  fprintf(stderr,"entrycount before read bed file %d \n",entryCount2);
  bedentries2 = readBedFile(bedfilename,&entryCount2);
  fprintf(stderr,"entrycount after read bed file %d \n",entryCount2);

  for (int i = 0; i < entryCount2; i++) {
    fprintf(stderr,"bed file information bed entry \t%s\t%d\t%d\n", bedentries2[i].chromosome, bedentries2[i].start, bedentries2[i].end);
  }
  free(bedentries2);
  fprintf(stderr,"---------------------- sort ----------------------\n");
  
  int entryCount3 = 0;
  const char* bedfilename2 = "Coord2.bed";
  BedEntry* bedentries3 = NULL;

  fprintf(stderr,"entrycount before read bed file %d \n",entryCount3);
  bedentries3 = readBedFile(bedfilename2,&entryCount3);
  fprintf(stderr,"entrycount after read bed file %d \n",entryCount3);
  fprintf(stderr,"----\n");

  for (int i = 0; i < entryCount3; i++) {
    fprintf(stderr,"bed file information bed entry \t%s\t%d\t%d\n", bedentries3[i].chromosome, bedentries3[i].start, bedentries3[i].end);
  }
  fprintf(stderr,"----\n");
  sortBedEntries(bedentries3, entryCount3);

  for (int i = 0; i < entryCount3; i++) {
    fprintf(stderr,"bed file information bed entry \t%s\t%d\t%d\n", bedentries3[i].chromosome, bedentries3[i].start, bedentries3[i].end);
  }
  fprintf(stderr,"---------------------- merge ----------------------\n");

  int mergedCount = 0;
  BedEntry* mergedEntries = mergeOverlappingRegions(bedentries3, entryCount3, &mergedCount);

  for (int i = 0; i < mergedCount; i++) {
    fprintf(stderr,"bed file information bed entry \t%s\t%d\t%d\n", mergedEntries[i].chromosome, mergedEntries[i].start, mergedEntries[i].end);
  }

  return 0;
}

#endif