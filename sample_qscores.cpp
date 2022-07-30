#include <zlib.h>
#include <cstring>
#include <cassert>
#include <zlib.h>
#include "sample_qscores.h"

#define LENS 10000
#define MAXBINS 100

char int2nuc[5] = {'A','C','G','T','N'};
double phred2Prob(int qscore){
  double d = qscore;
  return pow(10,-d/10);
}
ransampl_ws ***ReadQuality(char *ntqual, double *ErrProb, int ntcharoffset,const char *freqfile,int &readcycle){
  ransampl_ws ***dists = new ransampl_ws**[5];

  std::vector<char *> all_lines;
  gzFile gz = Z_NULL;
  assert(((gz = gzopen(freqfile,"rb")))!=Z_NULL);
  char buf[LENS];
  while(gzgets(gz,buf,LENS))
    all_lines.push_back(strdup(buf));
  gzclose(gz);

  //fprintf(stderr,"All lines: %lu\n",all_lines.size());

  readcycle = (all_lines.size()-2)/5; //all_lines.size()-1

  fprintf(stderr,"\t-> Inferred read cycle lengths: %d\n",readcycle);

  //loop over inputdata
  int nbins = -1;
  double probs[MAXBINS];
  for(int b=0;b<5;b++){
    dists[b] = new ransampl_ws *[readcycle];
    for(unsigned long pos = 0 ; pos<readcycle;pos++){
      int at = 0;
      probs[at++] = atof(strtok(all_lines[2+b*readcycle+pos],"\n\t ")); //1+b*readcyclelength+pos
      char *tok = NULL;
      while(((tok=strtok(NULL,"\n\t ")))){
	      probs[at++] = atof(tok);
	      assert(at<MAXBINS);
      }
      if(nbins==-1){
	      nbins = at;
	      //fprintf(stderr,"Number of qualities/bins in inputfile: %d\n",nbins);
      }
      if(nbins!=at){
	      fprintf(stderr,"Problems, number of columns is different nbins: %d at: %d\n",nbins,at);
	      exit(0);
      }
      dists[b][pos] =  ransampl_alloc( nbins );
      ransampl_set(dists[b][pos],probs);
    }
  }
  //fprintf(stderr,"\t-> ransampl_ws done\n");

  //printf(all_lines[0]);
  int qualidx = 1;
  //extract the first token
  ntqual[0] = (char) (atoi(strtok(all_lines[0],"\n\t "))+ntcharoffset);
  char *qualtok = NULL;
  //extract the next
  while(((qualtok=strtok(NULL,"\n\t ")))){ntqual[qualidx++] = (char) (atoi(qualtok)+ntcharoffset);}
  
  int Err_idx = 1;
  ErrProb[0] = atof(strtok(all_lines[1],"\n\t"));
  char *Errtok = NULL;
  while ((Errtok=strtok (NULL,"\n\t"))){ErrProb[Err_idx++] = (double) atof(Errtok);}
  
  //strdup function allocate necessary memory to store the sourcing string implicitly, i need to free the returned string
  for (unsigned long i = 0; i < all_lines.size(); i++){free(all_lines[i]);}
  //fprintf(stderr,"\t-> Before return in READQUAL FUNC\n");
  return dists;
}
void sample_qscores(char *bases, char *qscores,int len,ransampl_ws ***ws,char *NtQuals,mrand_t *mr,int simError){
  unsigned char nuc2int[256];
  memset(nuc2int,4,256);
  nuc2int['a'] = nuc2int['A'] = nuc2int[0] = 0;
  nuc2int['t'] = nuc2int['T'] = nuc2int[1] = 1;
  nuc2int['g'] = nuc2int['G'] = nuc2int[2] = 2;
  nuc2int['c'] = nuc2int['C'] = nuc2int[3] = 3;
  nuc2int['n'] = nuc2int['N'] = nuc2int[4] = 4;
  
  for(int i = 0;i<len;i++){
    double dtemp1 = mrand_pop(mr);
    double dtemp2 = mrand_pop(mr);
    
    char inbase = bases[i];// <- is this ACGT,or 01234
    int qscore = ransampl_draw2(ws[nuc2int[inbase]][i],dtemp1,dtemp2);
    qscores[i] = NtQuals[qscore]+33;
      if (simError){
      if ( mrand_pop(mr) < phred2Prob(qscore)){
	int outbase;
	int inbase = nuc2int[bases[i]];
	while ((outbase=((int)floor(4*mrand_pop(mr)))) == inbase);
	bases[i] = int2nuc[outbase];
      }
    }
  }
}


#ifdef __WITH_MAIN__
int main(int argc, char **argv){
  const char *profile_fname = "Test_Examples/Qual_profiles/AccFreqL150R1.txt";
  char ntquals[1024];
  double errorArray[1024];
  int maxreadcycles;
  ransampl_ws ***ws = ReadQuality(ntquals,errorArray,33,profile_fname,maxreadcycles);
  mrand_t *mr = mrand_alloc(3,88);
  char bases[30];
  char qscores[30];
  memset(bases,'\0',30);
  memset(qscores,'\0',30);
  for(int i=0;i<30;i++)
    bases[i] = int2nuc[(int)floor(drand48()*4)];
  sample_qscores(bases,qscores,30,ws,ntquals,mr,1);
  fprintf(stderr,"@readname\n%s\n+\n%s\n",bases,qscores);
  return 0;
}
#endif
