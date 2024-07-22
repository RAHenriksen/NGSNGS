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

int VCFtoBED(const char* bcffilename, int id,int range,BedEntry** bedentries, int* entryCount){

  //create the bed entry struct
  int region_capacity = INITIAL_BED_CAPACITY;
  *bedentries = (BedEntry*) malloc(region_capacity * sizeof(BedEntry));  // Cast to BedEntry*

  //open the vcf file and store the header information
  htsFile *bcf = bcf_open(bcffilename, "r");

  if (!bcf) {
    fprintf(stderr, "Error: Unable to open BCF file: %s\n", bcffilename);
    free(*bedentries);
    exit(1);
  }

  bcf_hdr_t *bcf_head = bcf_hdr_read(bcf);
  if (!bcf_head) {
    fprintf(stderr, "Error: Unable to read BCF header\n");
    bcf_close(bcf);
    free(*bedentries);
    exit(1);
  }

  // initialize vcf information
  bcf1_t *brec = bcf_init();
  if (!brec) {
    fprintf(stderr, "Error: Unable to initialize BCF record\n");
    bcf_hdr_destroy(bcf_head);
    bcf_close(bcf);
    free(*bedentries);
    exit(1);
  }

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
  fprintf(stderr,"The number of individuals present in the bcf header is %d with the chosen individual being %s\n",nsamples,sample_name);

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
      return 1;
    }

    // create the coordinate information in our internal bed struct from which we sample data
    variant_number++;
    int chr_end = chr_info[fai_chr].val->info[0];
      
    if (variant_number >= region_capacity) {
      region_capacity *= 2;
      BedEntry* temp = (BedEntry*) realloc(*bedentries, region_capacity * sizeof(BedEntry));

      if (!temp) {
        fprintf(stderr,"Error reallocating memory for bed entries\n");
        free(*bedentries);
        bcf_hdr_destroy(bcf_head);
        bcf_destroy(brec);
        bcf_close(bcf);
        exit(1);
      }
      *bedentries = temp;
    }

    //sscanf(line, "%s %d %d", bedentries[*entryCount].chromosome, &bedentries[*entryCount].start, &bedentries[*entryCount].end);
      
    strncpy((*bedentries)[variant_number-1].chromosome,bcf_hdr_id2name(bcf_head, brec->rid),MAX_CHROM_NAME_LEN - 1);
    (*bedentries)[variant_number-1].chromosome[MAX_CHROM_NAME_LEN - 1] = '\0'; //null termination for chromosome name

    if((brec->pos - range) < 1){
      (*bedentries)[variant_number-1].start = 1;
      (*bedentries)[variant_number-1].end = brec->pos + range;
    }
    else if((brec->pos + range) > chr_end){
      (*bedentries)[variant_number-1].start = brec->pos - range;
      (*bedentries)[variant_number-1].end = chr_end;
    }
    else{
      (*bedentries)[variant_number-1].start = brec->pos - range;
      (*bedentries)[variant_number-1].end = brec->pos + range;
    }
  }
  fprintf(stderr,"The inferred ploidy being %d\n",inferred_ploidy);

  *entryCount = variant_number;

  bcf_hdr_destroy(bcf_head);
  bcf_destroy(brec);
  bcf_close(bcf);

  return 0;
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
  fprintf(stderr,"----------------------\n");
  free(bedentries2);

  return 0;
}

#endif