#include <zlib.h>
#include <cstring>
#include <cassert>
#include "sample_qscores.h"
#include "RandSampling.h"
#include "mrand.h"
#include "NGSNGS_misc.h"

#define LENS 10000
#define MAXBINS 100

extern int refToInt[256];
extern char intToRef[4];
extern char refToChar[256];

double phred2Prob[256];

ransampl_ws ***ReadQuality(char *ntqual, double *ErrProb, int ntcharoffset,const char *freqfile){
  for(int qscore =0 ;qscore<256;qscore++){
    double d = qscore;
    phred2Prob[qscore] = pow(10,((-d)/10));
    //fprintf(stderr,"qscore prob cal %lf \t %lf\n",qscore,phred2Prob[qscore]);
    //std::cout << qscore << " " << phred2Prob[qscore] << std::endl;
  }
  ransampl_ws ***dists = new ransampl_ws**[5];

  std::vector<char *> all_lines;
  gzFile gz = Z_NULL;
  assert(((gz = gzopen(freqfile,"rb")))!=Z_NULL);
  char buf[LENS];
  while(gzgets(gz,buf,LENS))
    all_lines.push_back(strdup(buf));
  gzclose(gz);

  //fprintf(stderr,"All lines: %lu\n",all_lines.size());

  int readcycle = (all_lines.size()-2)/5; //all_lines.size()-1

  //loop over inputdata
  int nbins = -1;
  double probs[MAXBINS];
  for(int b=0;b<5;b++){
    dists[b] = new ransampl_ws *[readcycle];
    for(int pos = 0 ; pos<readcycle;pos++){
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
  //std::cout << "nt qual" << ntqual[1] << std::endl;
  int Err_idx = 1;
  ErrProb[0] = atof(strtok(all_lines[1],"\n\t"));
  char *Errtok = NULL;
  while ((Errtok=strtok (NULL,"\n\t"))){ErrProb[Err_idx++] = (double) atof(Errtok);}
  //std::cout << "error prob " << ErrProb[1] << std::endl;

  //strdup function allocate necessary memory to store the sourcing string implicitly, i need to free the returned string
  for (unsigned long i = 0; i < all_lines.size(); i++){free(all_lines[i]);}
  //fprintf(stderr,"\t-> Before return in READQUAL FUNC\n");
  return dists;
}
int sample_qscores(char *bases, char *qscores,int len,ransampl_ws ***ws,char *NtQuals,mrand_t *mr,int simError, int ntcharoffset){
  int seq_err = 0;
  for(int i = 0;i<len;i++){
    double dtemp1 = mrand_pop(mr);
    double dtemp2 = mrand_pop(mr);
    
    char inbase = refToInt[bases[i]];
    //fprintf(stdout,"dtemp %lf\t%lf\n",dtemp1,dtemp2);
    //fflush(stdout);
    int qscore_idx =  ransampl_draw2(ws[inbase][i],dtemp1,dtemp2);//ransampl_draw2(ws[inbase][i],1,1); //ransampl_draw2(ws[inbase][i],dtemp1,dtemp2);
    //fprintf(stderr,"qscore_idx %d \n",qscore_idx);
    //fprintf(stderr,"DID NOT FAIL HURRA TK HAPPY \n");
    //std::cout << qscore << std::endl;
    qscores[i] = NtQuals[qscore_idx];
    //std::cout << " qscore" << qscore_idx << " " << qscores[i]<<std::endl;
    //std::cout << " qscore " << qscores << " qscore i " << qscores[i] << " phred2prob " << phred2Prob[qscore] << std::endl; 
    if (simError){
      //fprintf(stderr,"PERFORMING SEQUENCING SUBSITUTIONS\n");
      //fprintf(stderr,"ntcharoffset %d\n",ntcharoffset);
      double tmprand = mrand_pop(mr);
      //std::cout << phred2Prob[qscore] << std::endl;
      //std::cout << phred2Prob[qscore] << " X " << " qscore i " << qscores[i] << " " << phred2Prob[qscores[i]-33] << std::endl;
      //phred2Prob[qscores[i]-33]
      if ( tmprand < phred2Prob[qscores[i]-ntcharoffset]){
        /*fprintf(stderr,"SIMULATION ERROR sequence length %d\n",len);
        std::cout << phred2Prob[qscores[i]]<< std::endl;
        std::cout << bases[i] << std::endl;*/
        int outbase=(int)floor(4.0*phred2Prob[qscores[i]-ntcharoffset]*tmprand);//DRAGON
        while (((outbase=((int)floor(4*mrand_pop(mr))))) == inbase);
	      bases[i] = intToRef[outbase];
        seq_err = 1;
        //std::cout << bases[i] << std::endl;
      }
    }
  }
  return seq_err;
}


#ifdef __WITH_MAIN__
int main(int argc, char **argv){
  const char *profile_fname = "Test_Examples/AccFreqL150R1.txt";
  char ntquals[1024];
  double errorArray[1024];
  int maxreadcycles = 150;
  ransampl_ws ***ws = ReadQuality(ntquals,errorArray,33,profile_fname);
  mrand_t *mr = mrand_alloc(3,88);
  char bases[30];
  char qscores[30];
  memset(bases,'\0',30);
  memset(qscores,'\0',30);

  for(int i=0;i<30;i++){
    bases[i] = intToRef[(int)floor(drand48()*4)];
    std::cout << " i " << i << " bases " << bases[i] << std::endl;
  }
  sample_qscores(bases,qscores,30,ws,ntquals,mr,0,0);
  //sample_qscores(seq_r1,qual_r1,strlen(seq_r1),struct_obj->QualDist_r1,struct_obj->NtQual_r1,rand_alloc,struct_obj->DoSeqErr,ErrProbTypeOffset);
  fprintf(stderr,"@readname\n%s\n+\n%s\n",bases,qscores);
  //gtggTAGAGATAAAGCACATTCTTTAGGAGTGAATATGGNNTNNCTGCNCGCANANTGNNATTGNNTTGCNNNTNNANCGNNNCNNTNNNGNTTNGCNACAGCNANGNNA
  return 0;
}
#endif

//g++ sample_qscores.cpp RandSampling.o mrand.o -std=c++11 -lm -lz -D__WITH_MAIN__ -o Scores