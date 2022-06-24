#ifndef NGSNGSFUNC_H
#define NGSNGSFUNC_H
#include "mrand.h"

static void delete_seq(char *str, size_t seq_len, size_t del_len, size_t pos,int alt_len);

static void delete_seq_ins(char *str, size_t seq_len, size_t del_len, size_t pos);

static void instert_seq(char *str, size_t len, char insert_seq[],size_t ins_len, size_t pos);

void DNA_CAPITAL(char seq[]);

void DNA_complement(char seq[]);

void reverseChar(char* str,int length);

int Random_geometric_k(unsigned int  seed,const double p);

double myrand(unsigned int persistent);

int myrandgenmodulo(unsigned int seed, int modulo);

double uniform();

void SimBriggsModel(char* reffrag, char* frag, int L, double nv, double lambda, double delta_s, double delta, unsigned int seed,mrand_t *mr);

double* DeamFileArray(double* freqval,const char* filename,int &deamcyclelength);

void Deam_File(char seq[],mrand_t *mr,double* freqval,int LEN);

const char* Error_lookup(double a,double err[6000],int nt_offset, int read_pos,int outputoffset);

double* Qual_array(double* freqval,const char* filename);

int BinarySearch_fraglength(double* SearchArray,int low, int high, double key);

void FragArray(int& number,int*& Length, double*& Frequency,const char* filename);

void printTime(FILE *fp);

void Header_func(htsFormat *fmt_hts,const char *outfile_nam,samFile *outfile,sam_hdr_t *header,faidx_t *seq_ref,int chr_total, int chr_idx[], size_t genome_len);

char* full_genome_create(faidx_t *seq_ref,int chr_total,int chr_sizes[],const char *chr_names[],size_t chr_size_cumm[]);

char* partial_genome_create(faidx_t *seq_ref,int chr_total,int chr_sizes[],const char *chr_names[],size_t chr_size_cumm[]);

char* full_vcf_genome_create(faidx_t *seq_ref,int chr_total,int chr_sizes[],const char *chr_names[],size_t chr_size_cumm[],const char* bcf_file,const char* VarType);

char* vcf_info(char genome_data[],int chr_sizes,const char* bcf_file,const char *chr_names[],const char* VarType,FILE *VarInfoFile);

typedef struct{
    int n;
    int* alias;
    double* prob;
} ransampl_ws;

ransampl_ws* ransampl_alloc(int n);

void ransampl_set( ransampl_ws *ws, const double *p );

int ransampl_draw2( ransampl_ws *ws,double r1, double r2); //added below function to make it threadsafe tsk 23dec 2021

void ransampl_free( ransampl_ws *ws );

ransampl_ws ***ReadQuality(char *ntqual, double *ErrProb,int ntcharoffset,const char *freqfile,unsigned long &readcycle);

void deletechar(char* str,int seq_len, size_t index_to_remove,int del_len);

void InsertChar(char* array,std::string ins,int index);

void ErrorSub(double randval,char seqchar[], int pos);


#endif /* NGSNGSFUNC_H */