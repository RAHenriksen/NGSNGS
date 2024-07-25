#include <zlib.h>
#include <cstring>
#include <cassert>
#include <htslib/kstring.h>

#include "RandSampling.h"
#include "sample_qscores.h"
#include "mrand.h"
#include "NGSNGS_misc.h"

#define LENS 10000
#define MAXBINS 100

extern int refToInt[256];
extern char intToRef[4];
extern char refToChar[256];

double phred2Prob[256];

ransampl_ws ***ReadQuality(char *ntqual, double *ErrProb, int ntcharoffset,const char *freqfile,int &inferredreadcycle){
  for(int qscore =0 ;qscore<256;qscore++){
    double d = qscore;
    phred2Prob[qscore] = pow(10,((-d)/10));
  }
  ransampl_ws ***dists = new ransampl_ws**[5];

  std::vector<char *> all_lines;
  gzFile gz = Z_NULL;
  assert(((gz = gzopen(freqfile,"rb")))!=Z_NULL);
  char buf[LENS];
  while(gzgets(gz,buf,LENS))
    all_lines.push_back(strdup(buf));
  gzclose(gz);

  inferredreadcycle = (all_lines.size()-2)/5;
  //loop over inputdata
  int nbins = -1;
  double probs[MAXBINS];
  for(int b=0;b<5;b++){
    dists[b] = new ransampl_ws *[inferredreadcycle];
    for(int pos = 0 ; pos<inferredreadcycle;pos++){
      int at = 0;
      probs[at++] = atof(strtok(all_lines[2+b*inferredreadcycle+pos],"\n\t "));
      char *tok = NULL;
      while(((tok=strtok(NULL,"\n\t ")))){
	      probs[at++] = atof(tok);
	      assert(at<MAXBINS);
      }
      if(nbins==-1){
	      nbins = at;
      }
      if(nbins!=at){
	      fprintf(stderr,"Problems, number of columns is different nbins: %d at: %d\n",nbins,at);
	      exit(1);
      }
      dists[b][pos] =  ransampl_alloc( nbins );
      ransampl_set(dists[b][pos],probs);
    }
  }

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
  return dists;
}

int sample_qscores(char *bases, char *qscores,int len,ransampl_ws ***ws,char *NtQuals,mrand_t *mr,int simError, int ntcharoffset){
  int seq_err = 0;
  for(int i = 0;i<len;i++){
    double dtemp1 = mrand_pop(mr);
    double dtemp2 = mrand_pop(mr);
    
    char inbase = refToInt[bases[i]];
    int qscore_idx =  ransampl_draw2(ws[inbase][i],dtemp1,dtemp2);//ransampl_draw2(ws[inbase][i],1,1); //ransampl_draw2(ws[inbase][i],dtemp1,dtemp2);

    qscores[i] = NtQuals[qscore_idx];
    if (simError){

      double tmprand = mrand_pop(mr);
      if ( tmprand < phred2Prob[qscores[i]-ntcharoffset]){

        int outbase=(int)floor(4.0*phred2Prob[qscores[i]-ntcharoffset]*tmprand);//DRAGON
        while (((outbase=((int)floor(4*mrand_pop(mr))))) == inbase);
	      bases[i] = intToRef[outbase];
        seq_err = 1;
      }
    }
  }

  return seq_err;
}

int sample_qscores_amplicon(kstring_t* seq,kstring_t* qual,ransampl_ws ***ws,char *NtQuals,mrand_t *mr,int simError, int ntcharoffset){
  int seq_err = 0;

  for(int i = 0;i<seq->l;i++){
    double dtemp1 = mrand_pop(mr);
    double dtemp2 = mrand_pop(mr);

    int qscore_idx =  ransampl_draw2(ws[refToInt[seq->s[i]]][i],dtemp1,dtemp2);//ransampl_draw2(ws[inbase][i],1,1); //ransampl_draw2(ws[inbase][i],dtemp1,dtemp2);
    qual->s[i] = NtQuals[qscore_idx];
  
    char inbase = refToInt[seq->s[i]];

    if (simError){
      double tmprand = mrand_pop(mr);
      if ( tmprand < phred2Prob[qual->s[i]-ntcharoffset]){

        int outbase=(int)floor(4.0*phred2Prob[qual->s[i]-ntcharoffset]*tmprand);//DRAGON
        while (((outbase=((int)floor(4*mrand_pop(mr))))) == inbase);
        seq->s[i] = intToRef[outbase];
        seq_err = 1;
      }
    }  
  }


  return seq_err;
}

int sample_qscores_fix(char *bases, char *qscores, int qscoresval,int len,mrand_t *mr,int simError, int ntcharoffset){
  //fprintf(stderr,"INSIDE THE SAMPLE QSCORE FIX VALUE\n");
  int seq_err = 0;
      
  if (simError){
    for(int qscore =0 ;qscore<256;qscore++){
      double d = qscore;
      phred2Prob[qscore] = pow(10,((-d)/10));
    }
  }

  char fixedscore = (char)(qscoresval + ntcharoffset);
  //fprintf(stderr,"Fixedscore %c\tprob %f\t length %d\n",fixedscore,phred2Prob[fixedscore-ntcharoffset],len);
  
  for(int i = 0;i<len;i++){
    qscores[i] = fixedscore;
    char inbase = refToInt[bases[i]];

    if (simError){
      double tmprand = mrand_pop(mr);
      if ( tmprand < phred2Prob[fixedscore-ntcharoffset]){
        int outbase=(int)floor(4.0*phred2Prob[fixedscore-ntcharoffset]*tmprand);
        while (((outbase=((int)floor(4*mrand_pop(mr))))) == inbase);
        //fprintf(stderr,"bases before %c \n",bases[i]);
	      bases[i] = intToRef[outbase];
        //fprintf(stderr,"bases after %c \n",bases[i]);
        seq_err = 1;
      }
    }
  }
  return seq_err;
}

//void add_indel_amplicon_fa(mrand_t *mr,kstring_t* seq,double *pars,int* ops){
//add_indel_amplicon_fqbam(mr,&FQseq->seq,&FQseq->qual,pars,ops,ErrProbTypeOffset);

int sample_qscores_fix_amplicon(mrand_t *mr,kstring_t* seq,int qscoresval,int ntcharoffset){
  //std::cout << seq->s << std::endl;
  int seq_err = 0;
  
  kstring_t seq_intermediate;
  seq_intermediate.l = seq->l;
  seq_intermediate.m = seq->l;
  seq_intermediate.s = (char *)malloc((seq->l + 1) * sizeof(char)); // Allocate memory
  strcpy(seq_intermediate.s, seq->s); // Copy seq->s to seq_intermediate.s
  

  for(int qscore =0 ;qscore<256;qscore++){
    double d = qscore;
    phred2Prob[qscore] = pow(10,((-d)/10));
  }

  for(int i = 0;i<seq->l;i++){
    char inbase = refToInt[seq->s[i]];

    double tmprand = mrand_pop(mr);
    //fprintf(stderr,"inside simerror %f \n",phred2Prob[qscoresval-33]);
    //exit(1);
    if ( tmprand < phred2Prob[qscoresval-ntcharoffset]){
      int outbase=(int)floor(4.0*phred2Prob[qscoresval]*tmprand);
      while (((outbase=((int)floor(4*mrand_pop(mr))))) == inbase);
	    seq_intermediate.s[i] = intToRef[outbase];
      seq->s[i] = seq_intermediate.s[i];
      seq_err = 1;
    }
  }

  seq->s[seq->l] = '\0';

  free(seq_intermediate.s);
  seq_intermediate.s = NULL;
  seq_intermediate.l = seq_intermediate.m = 0;

  return seq_err;
}

#ifdef __WITH_MAIN__
int main(int argc, char **argv){
  const char *profile_fname = "Test_Examples/AccFreqL150R1.txt";
  char ntquals[1024];
  double errorArray[1024];
  int maxreadcycles = 150;
  //ransampl_ws ***ws = ReadQuality(ntquals,errorArray,33,profile_fname);
  mrand_t *mr = mrand_alloc(3,88);
  char bases[30];
  char qscores[30];
  memset(bases,'\0',30);
  memset(qscores,'\0',30);

  char specific_char = 'I';
  char* char_array = (char*)malloc((30 + 1) * sizeof(char));
  for(int i=0;i<30;i++){
    bases[i] = intToRef[(int)floor(drand48()*4)];
    char_array[i] = specific_char;
    std::cout << " i " << i << " bases " << bases[i] << std::endl;
  }
  sample_qscores_fix(bases,qscores,10,30,mr,1,33);
  //sample_qscores(bases,qscores,30,ws,ntquals,mr,0,0);
  fprintf(stderr,"@readname\n%s\n+\n%s\n",bases,qscores);
  //gtggTAGAGATAAAGCACATTCTTTAGGAGTGAATATGGNNTNNCTGCNCGCANANTGNNATTGNNTTGCNNNTNNANCGNNNCNNTNNNGNTTNGCNACAGCNANGNNA
  return 0;
}
#endif

//g++ sample_qscores.cpp RandSampling.o mrand.o NGSNGS_misc.o -std=c++11 -lm -lz -D__WITH_MAIN__ -o Scores