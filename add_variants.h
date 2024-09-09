#ifndef ADD_VARIANTS_H
#define ADD_VARIANTS_H


typedef struct{
  int rid;
  int pos;
  int *gt;
}bcfkey;

struct bcfkey_cmp
{
  bool operator()(bcfkey s1,bcfkey s2) const
  {
    if(s1.rid==s2.rid)
      return s1.pos<s2.pos;
    return s1.rid<s2.rid;
  }
};

//NB notice that proper rid is NOT bcf1_t->rid but bcfkey->rid which has been shifted relative to fs
typedef std::map<bcfkey,bcf1_t*,bcfkey_cmp> bcfmap;

int *fabcflookup(bcf_hdr_t *bcf_hdr,char2int &fai2idx,int &maxIndex);

void add_snp(fasta_sampler *fs, int fs_chr_idx,int pos,char **alleles, int32_t *gts, int ploidy,const char* Name);

int add_vcf_variants(fasta_sampler *fs,const char *bcffilename,int HeaderIndiv,const char* Name);

void add_indels_simple(fasta_sampler *fs,bcfmap &mybcfmap,bcf_hdr_t *hdr,int ploidy);

void add_ins_complex(fasta_sampler *fs,bcfmap &mybcfmap,bcf_hdr_t *hdr,int ploidy);

#endif
