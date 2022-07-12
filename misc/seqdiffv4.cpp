#include <map>//for map
#include <cstring>//for strcmp
#include <cstdio>//for stderr
#include <cassert>//for assert
#include <htslib/bgzf.h> //for bgzf
#include <htslib/faidx.h>//for faidx
#include <htslib/kstring.h>

struct ltstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return strcmp(s1, s2) < 0;
  }
};

typedef struct{
  char *seq;
  int seq_l;
}mydata;

#define RLEN 1024 //very long readlengths
size_t maxl = 0;//contains our maximum readlength
size_t MM[2][RLEN][5][5] = {0};


typedef std::map<char *,mydata,ltstr> aMap;

aMap mymap;
faidx_t *fai;

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

char intToRef[5] = {'A','C','G','T','N'};

int main(int argc,char**argv){
  char *fq = NULL;
  char *fa = NULL;
  if(argc!=3){
    fprintf(stderr,"\t-> Supply ./a.out fq fa\n");
    return 0;
  }
  fq = argv[1];
  fa = argv[2];
  fprintf(stderr,"fq: %s fa: %s\n",fq,fa);

  BGZF *fq_fp = NULL;
  assert(((fq_fp = bgzf_open(fq,"r"))!=NULL));
  fai = NULL;
  assert(((fai =fai_load(fa)))!=NULL);

  kstring_t *kstr =(kstring_t*) calloc(1,sizeof(kstring_t));

  int verbose[2] = {0,0};
  char *contig = (char*) calloc(2048,1);//maxlength 1024, should be enough
  char *readname = (char*) calloc(2048,1);//maxlength 1024, should be enough
  while(bgzf_getline(fq_fp,'\n',kstr)>0){
    //    fprintf(stderr,"%s\n",kstr->s);
    strcpy(readname,kstr->s);
    char *tok = strtok(kstr->s,"_");

    int threadid;
    int readid;
    int strand;//assumes zero or one 0/1    
    int posB;
    int posE;
    int length;
    assert(sscanf(tok,"@T%d",&threadid)==1);
    tok = strtok(NULL,"_");
    assert(sscanf(tok,"RID%d",&readid)==1);
    tok = strtok(NULL,"_");
    assert(sscanf(tok,"S%d",&strand)==1);
    tok = strtok(NULL,":");
    strcpy(contig,tok);
    tok = strtok(NULL,"\n");
    assert(sscanf(tok,"%d-%d_length:%d",&posB,&posE,&length)==3);
    //    fprintf(stderr,"T=%d R=%d S=%d C=%s P1=%d P2=%d L=%d\n",threadid,readid,strand,contig,posB,posE,length);
    aMap::iterator it = mymap.find(contig);
    if(it==mymap.end()){
      fprintf(stderr,"contig: \'%s\' has not been loaded: will load\n",contig);
      mydata md;
      md.seq = fai_fetch(fai,contig,&md.seq_l);
      fprintf(stderr,"contig: %s length: %d\n",contig,md.seq_l);
      mymap[strdup(contig)] = md;
    }
    it = mymap.find(contig);
    assert(it!=mymap.end());
    mydata md = it->second;

    //read in the read, kstr->s contains read
    assert(bgzf_getline(fq_fp,'\n',kstr)>0);
    if(kstr->l>maxl)
      maxl = kstr->l;
    
    //strand = -1 is used for signalling not to do anything
    if(strand!=-1) {
      if(kstr->l+posB!=posE+1) {
	if(verbose[0]++<5)
	  fprintf(stderr,"[%s] Problem with offsets of reads vs end of read?: %lu %d, kstr->l: %lu, this message is printed %d time more\n",readname,kstr->l+posB,posE,kstr->l,5-verbose[0]);
	strand = -1;
      }
    }
    
    if(strand!=-1){
      if(kstr->l+posB>=md.seq_l||posB>md.seq_l){
	if(verbose[1]++<5)
	  fprintf(stderr,"[%s] Problem with end of read and length of contig: %d_%d and  %lu fai_l: %d, this message is printed %d time more\n",readname,posB,posE,kstr->l+posB,md.seq_l,5-verbose[1]);
	strand = -1;
      }

    }

    if(strand!=-1){
      for(unsigned i=0;i<kstr->l;i++){
	MM[strand][i][refToInt[md.seq[posB+i]]][refToInt[kstr->s[i]]]++;
      }
    }

    //forget about + and qscores
    assert(bgzf_getline(fq_fp,'\n',kstr)>0);
    assert(bgzf_getline(fq_fp,'\n',kstr)>0);
    kstr->l = 0;
  }
  //printout results
  fprintf(stdout,"Strand\tPosition");
  for(int r=0;r<5;r++)
    for(int o=0;o<5;o++)
      fprintf(stdout,"\t%c->%c",intToRef[r],intToRef[o]);
  fprintf(stdout,"\n");

  for(int s=0;s<2;s++){
    for(size_t i=0;i<maxl;i++){
      int tmpsum = 0;
      for(int r=0;r<5;r++)
	for(int o=0;o<5;o++)
	  tmpsum += MM[s][i][r][o];
      if(tmpsum==0)
	continue;
      fprintf(stdout,"%d\t%lu",s,i);
      for(int r=0;r<5;r++)
	for(int o=0;o<5;o++)
	  fprintf(stdout,"\t%lu",MM[s][i][r][o]);
      fprintf(stdout,"\n");
    }
  }

  //cleanup
  free(contig);
  free(readname);
  for(aMap::iterator it=mymap.begin();it!=mymap.end();it++){
    mydata md = it->second;
    free(md.seq);
    free(it->first);
  }
  fai_destroy(fai);
  bgzf_close(fq_fp);
  free(kstr->s);
  free(kstr);
  return 0;
}
