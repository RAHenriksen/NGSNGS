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

char intToRef[4] = {'A','C','G','T'};

char refToChar[256] = {
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

char NtComp[5] = {'T', 'G', 'C', 'A','N'};

const char *bass = "ACGTN";

void ReversComplement(char* seq) {
  // Reverse and complement input sequence in place
  reverseChar(seq, strlen(seq));
  Complement(seq);
}

void Complement(char* seq) {
  // Complemen input sequence in place
  for (; *seq; ++seq) {
    *seq = NtComp[refToInt[(unsigned char)*seq]];
  }
}

void reverseChar(char* str,int length) {
    std::reverse(str, str + length);
}

void Complement_k(kstring_t* seq) {
  // Complement input sequence in place
  char *s = seq->s;
  for (; *s; ++s) {
    *s = NtComp[refToInt[(unsigned char)*s]];
  }
}

void ReversComplement_k(kstring_t* seq) {
  // Reverse and complement input sequence in place
  reverseChar(seq->s, seq->l);
  Complement_k(seq);
}

char* PrintCigarBamSet1(size_t n_cigar,const uint32_t *cigar){
   // Allocate a buffer to store the CIGAR string
  char *cigar_str = (char *)malloc(n_cigar * 2); // a safe overestimate
  char *ptr = cigar_str;

  /*size_t l_qname, const char *qname,
  uint16_t flag, int32_t tid, hts_pos_t pos, uint8_t mapq,
  size_t n_cigar, const uint32_t *cigar,
  int32_t mtid, hts_pos_t mpos, hts_pos_t isize,
  size_t l_seq, const char *seq, const char *qual,
  size_t l_aux)*/

  // Convert CIGAR to string
  for (int i = 0; i < n_cigar; ++i) {
    int len = sprintf(ptr, "%d%c", bam_cigar_oplen(cigar[i]), bam_cigar_opchr(cigar[i]));
    ptr += len;
  }

  return cigar_str;
}

void CreateSeqQualKString(bam1_t *aln, kstring_t *Sequence, kstring_t *Quality, int offset){
  int len = aln->core.l_qseq;

  Sequence->s = (char *) realloc(Sequence->s, len + 1); // allocate memory for the sequence string (+1 for null terminator)
  Sequence->l = len;
  Sequence->m = len + 1;

  Quality->s = (char *) realloc(Quality->s, len + 1); // allocate memory for the sequence string (+1 for null terminator)
  Quality->l = len;
  Quality->m = len + 1;

  uint8_t *seq_data = bam_get_seq(aln);
  uint8_t *qual_data = bam_get_qual(aln);

  for (int i = 0; i < len; i++) {
    Sequence->s[i] = seq_nt16_str[bam_seqi(seq_data, i)]; // convert nucleotide id to IUPAC id
    Quality->s[i] = qual_data[i]+offset;
  }
  
  Sequence->s[len] = '\0'; // null-terminate the string
  Quality->s[len] = '\0';
}

BedEntry* readBedFile(const char* filename, int* entryCount){
  // Function to read BED file and store entries in an array of BedEntry
  FILE* file = fopen(filename, "r");
  if (!file) {
    fprintf(stderr,"Error opening bed file");
    exit(1);
  }

  int region_capacity = INITIAL_BED_CAPACITY;
  BedEntry* bedentries = (BedEntry*) malloc(region_capacity * sizeof(BedEntry));
  if (!bedentries) {
    fprintf(stderr,"Error allocating memory");
    fclose(file);
    return NULL;
  }

  char line[256];
  int tmpcount = 0;
  while (fgets(line, sizeof(line), file)) {
    // dynamically reallocate the entry capacity
    if (tmpcount >= region_capacity) {
      region_capacity *= 2;
      BedEntry* temp = (BedEntry*) realloc(bedentries, region_capacity * sizeof(BedEntry));
      bedentries = temp;
    }

    // Parse the line and store in the structure
    sscanf(line, "%s\t%d\t%d", bedentries[tmpcount].chromosome, &bedentries[tmpcount].start, &bedentries[tmpcount].end);
    //fprintf(stderr,"%s\t%d\t%d\t%d\n", bedentries[*entryCount].chromosome, bedentries[*entryCount].start, bedentries[*entryCount].end,*entryCount);
    tmpcount++;
  }

  *entryCount = tmpcount;
  fclose(file);
  return bedentries;
}

// Comparison function for qsort
int compareBedEntries(const void* a, const void* b) {
    BedEntry* entryA = (BedEntry*)a;
    BedEntry* entryB = (BedEntry*)b;

    int chromCompare = strcmp(entryA->chromosome, entryB->chromosome);
    if (chromCompare == 0) {
        return (entryA->start - entryB->start);
    } else {
        return chromCompare;
    }
}

void sortBedEntries(BedEntry* entries, int entryCount) {
    qsort(entries, entryCount, sizeof(BedEntry), compareBedEntries);
}

BedEntry* mergeOverlappingRegions(BedEntry* entries, int entryCount, int* mergedCount) {
    if (entryCount == 0) {
        *mergedCount = 0;
        return NULL;
    }

    BedEntry* mergedEntries = (BedEntry*) malloc(entryCount * sizeof(BedEntry));
    if (!mergedEntries) {
        fprintf(stderr,"Error allocating memory for merged entries");
        return NULL;
    }

    int j = 0;
    mergedEntries[j] = entries[0];

    for (int i = 1; i < entryCount; i++) {
        if (strcmp(mergedEntries[j].chromosome, entries[i].chromosome) == 0 && mergedEntries[j].end >= entries[i].start) {
            // There is an overlap, so merge the regions
            if (entries[i].end > mergedEntries[j].end) {
                mergedEntries[j].end = entries[i].end;
            }
        } else {
            // No overlap, move to the next entry
            j++;
            mergedEntries[j] = entries[i];
        }
    }

    *mergedCount = j + 1;
    return mergedEntries;
}

BedEntry* checkbedentriesfasta(fasta_sampler *fs,BedEntry* entries,int entryCount,int* ReferenceCount) {
  if (entryCount == 0) {
    *ReferenceCount = 0;
    return NULL;
  }
    
  int n_seqs = faidx_nseq(fs->fai);
  /*for(int i=0;i<n_seqs;i++){
    const char *chrname = faidx_iseq(fs->fai,i);
    fprintf(stderr,"nref val %d \t and name %s\n",i,chrname);
  }*/

  int warning = 0;

  BedEntry* ReferenceEntries = (BedEntry*) malloc(entryCount * sizeof(BedEntry));
  int reference_int = 0;
  for(int i=0;i<n_seqs;i++){
    const char *chrname = faidx_iseq(fs->fai,i);
    //fprintf(stderr,"nref val %d \t and name %s\n",i,chrname);
    for (int j = 0; j < entryCount; j++) {
      char region[1024];
      sprintf(region, "%s:%d-%d", entries[j].chromosome,entries[j].start,entries[j].end);
      if (strcmp(chrname, entries[j].chromosome) == 0) {
        //fprintf(stderr,"ref seq %d \t entry val %d \n",i,j);
        //fprintf(stderr,"Bed entry %s \t and name %s \t %d \t %d\n",entries[j].chromosome,chrname,entries[j].start,entries[j].end);
        strcpy(ReferenceEntries[reference_int].chromosome,entries[j].chromosome);
        ReferenceEntries[reference_int].start = entries[j].start;
        ReferenceEntries[reference_int].end = entries[j].end;
        reference_int++; 
      }
      else{
        warning =1 ;
        continue;
      }
    }
  }
  
  if(warning==1)
    fprintf(stderr,"\t-> iscrepancy between input referenDce genome (-i) and input bed file, some regions in bed is not present in reference or some chromosomes in reference is not present in bed file, and these remains unused for further simulation \n");

  *ReferenceCount = reference_int;

  return ReferenceEntries;
}

BedEntry* maskbedentriesfasta(fasta_sampler *fs,BedEntry* entries,int entryCount,int* ReferenceCount) {
  if (entryCount == 0) {
    *ReferenceCount = 0;
    return NULL;
  }

  int region_capacity = INITIAL_BED_CAPACITY;

  int n_seqs = faidx_nseq(fs->fai);

  int warning = 0;

  BedEntry* ReferenceEntries = (BedEntry*) malloc(region_capacity * sizeof(BedEntry)); //think about how many i actually need for allocation
  int reference_int = 0;
  
  for(int i=0;i<n_seqs;i++){
    const char *chrname = faidx_iseq(fs->fai,i);
    size_t seq_len = faidx_seq_len(fs->fai, chrname);

    int tmp_start = 1;
    int chr_not_present = 0;
    int capture_present_chr_end = 0;

    for (int j = 0; j < entryCount; j++){
      //fprintf(stderr,"reference %d \t entrycount %d\n",i,j);
      if((strcmp(chrname, entries[j].chromosome) == 0)){
        capture_present_chr_end = 1;
        chr_not_present = 1;
        if(tmp_start < seq_len){
          //fprintf(stderr,"Bed entry %s \t and name %s \t %d \t %d\n",entries[j].chromosome,chrname,tmp_start,entries[j].start);
          strcpy(ReferenceEntries[reference_int].chromosome,entries[j].chromosome);
          ReferenceEntries[reference_int].start = tmp_start;
          ReferenceEntries[reference_int].end = entries[j].start;
          
          tmp_start = entries[j].end;
          reference_int++;
        }
      }
    }

    if (capture_present_chr_end || !chr_not_present) {
      //fprintf(stderr,"%s\t%d\t%d\n",chrname,tmp_start, seq_len);
      strcpy(ReferenceEntries[reference_int].chromosome,chrname);
      ReferenceEntries[reference_int].start = tmp_start;
      ReferenceEntries[reference_int].end = seq_len;
      
      reference_int++;
    }

  }
  
  *ReferenceCount = reference_int;

  return ReferenceEntries;
}

BedEntry* vcftobedentries(const char* bcffilename, int id,size_t flanking, int* entryCount,int* ploidy){
  //fprintf(stderr,"NGSNGS_misc.cpp \t INSIDE VCF BED ENTRY\n");
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

  // initialize header
  bcf_hdr_t *bcf_head = bcf_hdr_read(bcf);
  if (!bcf_head) {
    fprintf(stderr, "Error: Unable to read BCF header\n");
    bcf_close(bcf);
    free(bedentries);
    exit(1);
  }

  // initialize vcf information
  bcf1_t *brec = bcf_init();
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
  int32_t nodatagt[2]={bcf_gt_unphased(0),bcf_gt_unphased(1)};
  int32_t *mygt=NULL;
  int whichsample = id; //index of individual
  
  if (whichsample < 0 || whichsample >= nsamples) {
    fprintf(stderr, "Error: Sample index out of range\n");
    bcf_hdr_destroy(bcf_head);
    bcf_destroy(brec);
    bcf_close(bcf);
    free(bedentries);
    exit(1);
  }

  const char *sample_name = bcf_hdr_int2id(bcf_head, BCF_DT_SAMPLE, whichsample);
  //fprintf(stderr,"NGSNGS_misc.cpp \t The number of individuals present in the bcf header is %d with the chosen individual being %s\n",nsamples,sample_name);

  int seqnames_l;
  const char **bcf_chrom_names = bcf_hdr_seqnames(bcf_head,&seqnames_l);
  bcf_idpair_t *chr_info = bcf_head->id[BCF_DT_CTG];

  int variant_number = 0;
  // extract the bcf file information to create the bedstruct
  while (bcf_read(bcf, bcf_head, brec) == 0) {
    bcf_unpack((bcf1_t *)brec, BCF_UN_ALL);

    ngt = bcf_get_genotypes(bcf_head, brec, &gt_arr, &ngt_arr);
    assert(ngt>0);
    assert((ngt %nsamples)==0);
    inferred_ploidy = ngt/nsamples;
      
    mygt = gt_arr+inferred_ploidy*whichsample;
    if(bcf_gt_is_missing(mygt[0])){
      fprintf(stderr,"\t-> Genotype is missing for pos: %ld will skip\n",brec->pos+1);
      continue;
    }
    int fai_chr = brec->rid;
    if (fai_chr == -1) {
      fprintf(stderr, "Error: Chromosome name not found in the reference\n");
      exit(1);
    }
    // create the coordinate information in our internal bed struct from which we sample data
    int chr_end = chr_info[fai_chr].val->info[0];
      
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
    }
    
    //sscanf(line, "%s %d %d", bedentries[*entryCount].chromosome, &bedentries[*entryCount].start, &bedentries[*entryCount].end);
    
    strcpy(bedentries[variant_number].chromosome,bcf_hdr_id2name(bcf_head, brec->rid));

    int tmp_start = (size_t) brec->pos - (size_t) flanking; //size_t is unsigned integer thus non-negative values
    size_t tmp_end = (size_t) brec->pos + (size_t) flanking; //size_t is unsigned integer thus non-negative values

    // i dont add the flanking to the start yet, i need to store the position of the variant first to alter the sequence
    int n = 1;
    bedentries[variant_number].overlappositions = (int*) malloc(n * sizeof(int));;
    if(tmp_start < 1){
      bedentries[variant_number].start = 1;
      bedentries[variant_number].end = brec->pos + flanking;
    }
    else if(tmp_end > chr_end){
      bedentries[variant_number].start = brec->pos - flanking;
      bedentries[variant_number].end = chr_end;
    }
    else{
      bedentries[variant_number].start = brec->pos - flanking;
      bedentries[variant_number].end = brec->pos + flanking;
    }
    //bedentries[variant_number].overlappositions[0] = brec->pos;
    bedentries[variant_number].overlappositions[0] = (int)(brec->pos +1);
    /*fprintf(stderr,"flanking start %d \t flanking end %d \t variation %d \n",bedentries[variant_number].start,bedentries[variant_number].end,bedentries[variant_number].overlappositions[0]);
    exit(1);*/
    //fprintf(stderr,"postition %d \t allele 1 is %s \t allele 2 is %s \n",brec->pos+1,brec->d.allele[bcf_gt_allele(mygt[0])],brec->d.allele[bcf_gt_allele(mygt[1])]);

    bedentries[variant_number].variants = new char *[inferred_ploidy];
    for (int v = 0; v < inferred_ploidy; v++){
      const char* allele_str = brec->d.allele[bcf_gt_allele(mygt[v])];
      bedentries[variant_number].variants[v] = strdup(allele_str); // Duplicate the string
      //fprintf(stderr,"NGSNGS_misc.cpp \t bed entry allele is %d \t %d \t %d \t %d \t %s \n",variant_number,v,bedentries[variant_number].start,bedentries[variant_number].end,bedentries[variant_number].variants[v]);
      //brec->d.allele[bcf_gt_allele(mygt[v])]
    }
    //fprintf(stderr,"postition %d \t allele 1 is %s \t allele 2 is %s \n",brec->pos+1,brec->d.allele[bcf_gt_allele(mygt[0])],brec->d.allele[bcf_gt_allele(mygt[1])]);
    variant_number++;
  }

  /*
  fprintf(stderr,"NGSNGS_misc.cpp \t BEDENTRIES with %d entries and ploidy %d\n",variant_number,inferred_ploidy);
  fprintf(stderr,"NGSNGS_misc.cpp \t bed entry allele is \t %d \t %d \t %s \t %s  \n",bedentries[0].start,bedentries[0].end,bedentries[0].variants[0],bedentries[0].variants[0]);
  fprintf(stderr,"NGSNGS_misc.cpp \t bed entry allele is \t %d \t %d \t %s \t %s  \n",bedentries[1].start,bedentries[1].end,bedentries[1].variants[0],bedentries[1].variants[1]);
  fprintf(stderr,"NGSNGS_misc.cpp \t The inferred ploidy being %d\n",inferred_ploidy);
  */
  
  *ploidy = inferred_ploidy;
  *entryCount = variant_number;

  bcf_hdr_destroy(bcf_head);
  bcf_destroy(brec);
  bcf_close(bcf);

  //fprintf(stderr,"NGSNGS_misc.cpp \t DONE INSIDE VCF BED ENTRY\n");
  return bedentries;
}

void addVariant(BedEntry* entry, const char* variant,int position) {
    // Reallocate memory for the list of variants
    char** newVariants = (char**) realloc(entry->variants, (entry->variantCount + 1) * sizeof(char*));
    int* newPositions = (int*) realloc(entry->overlappositions, (entry->variantCount + 1) * sizeof(int));
    if (!newVariants || !newPositions) {
        fprintf(stderr, "Error reallocating memory for variants or positions\n");
        return;
    }
    entry->variants = newVariants;
    entry->overlappositions = newPositions;

    // Add the new variant to the array
    entry->variants[entry->variantCount] = strdup(variant);
    if (!entry->variants[entry->variantCount]) {
        fprintf(stderr, "Error allocating memory for variant string\n");
        return;
    }
    entry->overlappositions[entry->variantCount] = position; // Store the position

    // Update the variant count
    entry->variantCount++;
}

BedEntry* VCFLinkageDisequilibrium(BedEntry* entries, int entryCount, int* mergedCount,int ploidy,size_t flanking) {
    fprintf(stderr,"INSIDE MERGE VCF\n");
    if (entryCount == 0) {
        *mergedCount = 0;
        return NULL;
    }

    BedEntry* mergedEntries = (BedEntry*) malloc(entryCount * sizeof(BedEntry));
    if (!mergedEntries) {
        fprintf(stderr,"Error allocating memory for merged entries");
        return NULL;
    }

    /*for(int i = 0; i < entryCount; i++) {
      for(int j = 0; j < ploidy; j++) {
        fprintf(stderr,"entry %d\t ploidy %d\t%s\t%d\t%d\t%s\n",i,j,entries[i].chromosome,entries[i].start,entries[i].end,entries[i].variants[j]);
      }
    }
    fprintf(stderr,"BEFORE MERGE\n");
    */

    // Initialize the first entry
    int j = 0;
    mergedEntries[j] = entries[0];
    mergedEntries[j].variantCount = 0;
    mergedEntries[j].variants = NULL;
    mergedEntries[j].overlappositions = NULL;

    // Add initial variants
    for (int k = 0; k < ploidy; k++) {
        addVariant(&mergedEntries[j], entries[0].variants[k],entries[0].start+flanking);
    }

    for (int i = 1; i < entryCount; i++) {
      if (strcmp(mergedEntries[j].chromosome, entries[i].chromosome) == 0 && mergedEntries[j].end >= entries[i].start) {
        // There is an overlap, so merge the regions
        if (entries[i].end > mergedEntries[j].end) {
          mergedEntries[j].end = entries[i].end;
        }
        // Add new variants
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

        // Add variants for new entry
        for (int k = 0; k < ploidy; k++) {
          addVariant(&mergedEntries[j], entries[i].variants[k], entries[i].start+flanking);
        }
      }
    }
    
    fprintf(stderr,"AFTER MERGE\n");

    /*for (int i = 0; i < j + 1; i++) {
      for (int k = 0; k < mergedEntries[i].variantCount; k++) {
        fprintf(stderr,"entry %d\t%s\t%d\t%d\t%s\t%d\n", i, mergedEntries[i].chromosome, mergedEntries[i].start, mergedEntries[i].end,mergedEntries[i].variants[k],mergedEntries[i].overlappositions[k]);
      }
      fprintf(stderr,"\n---------\n");
    }*/
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