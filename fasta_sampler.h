#ifndef FASTA_SAMPLER_H
#define FASTA_SAMPLER_H
#include <cstring>//for strcmp
#include <cstdio>//for stderr
#include <cassert>//for assert
#include <map> //for map
#include "RandSampling.h"
#include "mrand.h"
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#define MAX_CHROM_NAME_LEN 100

struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};

typedef std::map<const char *,int,ltstr> char2int;

typedef std::map<int,int*> ploidymap;

// Declaration of BedEntry
typedef struct {
    char chromosome[MAX_CHROM_NAME_LEN];
    int start;
    int end;
    char** variants;
    int* overlappositions;     // Array of original positions for each variant - relevant for merged regions
    int variantCount; // Number of variants
    int ploidy;
} BedEntry;

typedef struct{
  faidx_t *fai;
  int nref;
  char **seqs;
  char **seqs_names;
  int *seqs_l;
  size_t seq_l_total;
  ransampl_ws *ws;
  char2int char2idx;
  int *realnameidx;
  ploidymap  pldmap;
  BedEntry* BedReferenceEntries;
  int BedReferenceCount;
}fasta_sampler;

fasta_sampler *fasta_sampler_alloc_full(const char *fa);
fasta_sampler *fasta_sampler_alloc_subset(const char *fa,const char *SpecificChr);
fasta_sampler *fasta_sampler_alloc_bedentry(const char *fa,const char *bedfilename,size_t flanking);
fasta_sampler *fasta_sampler_alloc_maskbedentry(const char *fa,const char *bedfilename,size_t flanking);
fasta_sampler *fasta_sampler_alloc_vcf(const char *fa, const char *bcffilename,int id,const char* Name,size_t flanking);
fasta_sampler *fasta_sampler_alloc_vcf_LD(const char *fa, const char *bcffilename,int id,const char* Name,size_t flanking);

int extend_fasta_sampler(fasta_sampler *fs,int fs_chr_idx,int ploidy,const char* Name);

void fasta_sampler_destroy(fasta_sampler *fs);
void fasta_sampler_destroy_captureLD(fasta_sampler *fs);

char* sample(fasta_sampler *fs,mrand_t *mr,char **chromoname,int &chr_idx,int &posB,int &posE,int &fraglength,size_t& chr_end,int simmode);
void fasta_sampler_setprobs(fasta_sampler *fs);
void fasta_sampler_print(FILE *fp,fasta_sampler *fs);
void fasta_sampler_print2(fasta_sampler *fs);

void dump_internal(fasta_sampler *fs,const char* filename);
#endif
