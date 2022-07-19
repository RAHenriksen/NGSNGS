#ifndef NGSNGSFUNC_H
#define NGSNGSFUNC_H
#include "mrand.h"
#include "RandSampling.h"

void FragDistArray(int& number,int*& Length, double*& Frequency,int SizeDistType,int seed,int val1, int val2);

void FragArray(int& number,int*& Length, double*& Frequency,const char* filename);

void delete_seq(char *str, int seq_len, int del_len, size_t pos,int alt_len);

void delete_seq_ins(char *str, int seq_len, int del_len, size_t pos);

void instert_seq(char *str, int len, char insert_seq[],int ins_len, size_t pos);

void DNA_CAPITAL(char seq[]);

void DNA_complement(char seq[]);

void reverseChar(char* str,int length);

void deletechar(char* str,int seq_len, size_t index_to_remove,int del_len);

void InsertChar(char* array,std::string ins,int index);

void Header_func(htsFormat *fmt_hts,const char *outfile_nam,samFile *outfile,sam_hdr_t *header,faidx_t *seq_ref,int chr_total, int chr_idx[], size_t genome_len,char CommandArray[1024],const char* version);

char* HaploGenome(char* genome,char genome_data[],char genome_data2[],int chr_sizes,const char* bcf_file,const char *chr_names[],const char* VarType,FILE *VarInfoFile,const char* HeaderIndiv);

char* full_genome_create(faidx_t *seq_ref,int chr_total,int chr_sizes[],const char *chr_names[],size_t chr_size_cumm[]);

char* partial_genome_create(faidx_t *seq_ref,int chr_total,int chr_sizes[],const char *chr_names[],size_t chr_size_cumm[]);

char* full_vcf_genome_create(faidx_t *seq_ref,int chr_total,int chr_sizes[],const char *chr_names[],size_t chr_size_cumm[],const char* bcf_file,const char* VarType,const char* HeaderIndiv);

ransampl_ws ***ReadQuality(char *ntqual, double *ErrProb,int ntcharoffset,const char *freqfile,unsigned long &readcycle);




#endif /* NGSNGSFUNC_H */