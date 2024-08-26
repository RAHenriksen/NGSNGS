#include <zlib.h>
#include <cassert>
#include <cstring>
#include "NtSubModels.h"
#include "mrand.h"
#include <map>
#include "htslib/bgzf.h"

#define LENS 4096
#define MAXBINS 100

void ErrorSub(double randval,char seqchar[], int pos){
  // Generates nucleotide substitutions
  if (seqchar[pos] == 'A' || seqchar[pos] == 'a'){
    if (0 < randval && randval <= 1.0/3.0){seqchar[pos] = 'C';}
    else if (1.0/3.0 < randval && randval <= 2.0/3.0){seqchar[pos] = 'G';}
    else if (2.0/3.0 < randval && randval <= 1){seqchar[pos]  = 'T';}
  }
  else if (seqchar[pos] == 'C'|| seqchar[pos] == 'c'){
    if (0 < randval && randval <= 1.0/3.0){seqchar[pos] = 'G';}
    else if (1.0/3.0 < randval && randval <= 2.0/3.0){seqchar[pos]  = 'T';}
    else if (2.0/3.0 < randval && randval <= 1){seqchar[pos]  = 'A';}
  }
  else if (seqchar[pos] == 'G'|| seqchar[pos] == 'g'){
    if (0 < randval && randval <= 1.0/3.0){seqchar[pos] = 'T';}
    else if (1.0/3.0 < randval && randval <= 2.0/3.0){seqchar[pos]  = 'A';}
    else if (2.0/3.0 < randval && randval <= 1){seqchar[pos]  = 'C';}
  }
  else if (seqchar[pos] == 'T'|| seqchar[pos] == 't'){
    if (0 < randval && randval <= 1.0/3.0){seqchar[pos] = 'A';}
    else if (1.0/3.0 < randval && randval <= 2.0/3.0){seqchar[pos]  = 'C';}
    else if (2.0/3.0 < randval && randval <= 1){seqchar[pos]  = 'G';}
  }
}

void MisMatchFileArray(double* freqval,const char* filename,int &mismatchcyclelength,int &elements){
    char buf[LENS];
    int i = 0;
    gzFile gz = Z_NULL;
    gz = gzopen(filename,"r");
    assert(gz!=Z_NULL);
    while(gzgets(gz,buf,LENS)){
        double val1;double val2;double val3;double val4;
        sscanf(buf,"%lf\t%lf\t%lf\t%lf\n",&val1,&val2,&val3,&val4);
        freqval[i*4] = val1; freqval[i*4+1] = val2; freqval[i*4+2] = val3; freqval[i*4+3] = val4;
        i++;
        elements += 4;
    }
    gzclose(gz);
    mismatchcyclelength = i/8;
}

double* MisMatchMatrixMetagenomic(double* freqval,const char* filename,int &mismatchcyclelength){
  fprintf(stderr,"INSIDE MisMatchMatrixMetagenomic \n");
  return freqval;
}


int MisMatchFile(char seq[],mrand_t *mr,double* freqval,int LEN){
  int IsMis5 = 0;
  int IsMis3 = 0;
  int IsMis = 0; //if 1 then only 5' change, if 2 only 3' change, if 3 change in both ends

  char ntdeam[4] = {'A', 'T', 'G', 'C'};
  double dtemp1;
  // 5' moving downwards from the 1 postion in the sequence 
  int Astart = 0;
  int Tstart = LEN*4;
  int Gstart = LEN*8;
  int Cstart = LEN*12;

  // moving upstream from the last position
  int Aend3 = LEN*20;
  int Tend3 = LEN*24;
  int Gend3 = LEN*28;
  int Cend3 = LEN*32;
  
  /*
  fprintf(stderr,"LEN IS %d\n",LEN);
  fprintf(stderr,"ASTART IS %d \t TSTART %d \t GSTART %d \t CSTART %d\n",Astart,Tstart,Gstart,Cstart);
  fprintf(stderr,"AEND IS %d \t TEND %d \t GEND %d \t CEND %d\n",Aend3,Tend3,Gend3,Cend3);
  */
  
  int seqlen = strlen(seq);
  //fprintf(stderr,"SEQUENCIN BEFORE %d \n %s \n",IsMis,seq);
  //5'
  for (int row_idx = 0; row_idx < LEN;row_idx++){
    dtemp1 = mrand_pop(mr);
    if (seq[row_idx] == 'A' || seq[row_idx] == 'a'){
      if (dtemp1 <= freqval[Astart+(row_idx*4)]){
        seq[row_idx] = ntdeam[0];
      }
      else if (freqval[Astart+(row_idx*4)] < dtemp1 && dtemp1 <= freqval[Astart+(row_idx*4)+1]){
        seq[row_idx] = ntdeam[1];
        IsMis5 = 1;
      }
      else if (freqval[Astart+(row_idx*4)+1] < dtemp1 && dtemp1 <= freqval[Astart+(row_idx*4)+2]){
        seq[row_idx] = ntdeam[2];
        IsMis5 = 1;
      }
      else if (freqval[Astart+(row_idx*4)+2] < dtemp1 && dtemp1 <= freqval[Astart+(row_idx*4)+3]){
        seq[row_idx] = ntdeam[3];
        IsMis5 = 1;
      }
    }
    else if (seq[row_idx] == 'T' || seq[row_idx] == 't'){
      if (dtemp1 <= freqval[Tstart+(row_idx*4)]){
        seq[row_idx] = ntdeam[0];
        IsMis5 = 1;
      }
      else if (freqval[Tstart+(row_idx*4)] < dtemp1 && dtemp1 <= freqval[Tstart+(row_idx*4)+1]){
        seq[row_idx] = ntdeam[1];
      }
      else if (freqval[Tstart+(row_idx*4)+1] < dtemp1 && dtemp1 <= freqval[Tstart+(row_idx*4)+2]){
        seq[row_idx] = ntdeam[2];
        IsMis5 = 1;
      }
      else if (freqval[Tstart+(row_idx*4)+2] < dtemp1 && dtemp1 <= freqval[Tstart+(row_idx*4)+3]){
        seq[row_idx] = ntdeam[3];
        IsMis5 = 1;
      }
    }
    else if (seq[row_idx] == 'G' || seq[row_idx] == 'g'){
      if (dtemp1 <= freqval[Gstart+(row_idx*4)]){
        seq[row_idx] = ntdeam[0];
        IsMis5 = 1;
      }
      else if (freqval[Gstart+(row_idx*4)] < dtemp1 && dtemp1 <= freqval[Gstart+(row_idx*4)+1]){
        seq[row_idx] = ntdeam[1];
        IsMis5 = 1;
      }
      else if (freqval[Gstart+(row_idx*4)+1] < dtemp1 && dtemp1 <= freqval[Gstart+(row_idx*4)+2]){
        seq[row_idx] = ntdeam[2];
      }
      else if (freqval[Gstart+(row_idx*4)+2] < dtemp1 && dtemp1 <= freqval[Gstart+(row_idx*4)+3]){
        seq[row_idx] = ntdeam[3];
        IsMis5 = 1;
      }
    }
    else if (seq[row_idx] == 'C' || seq[row_idx] == 'c'){
      if (dtemp1 <= freqval[Cstart+(row_idx*4)]){
        seq[row_idx] = ntdeam[0];
        IsMis5 = 1;
      }
      else if (freqval[Cstart+(row_idx*4)] < dtemp1 && dtemp1 <= freqval[Cstart+(row_idx*4)+1]){
        seq[row_idx] = ntdeam[1];
        IsMis5 = 1;
      }
      else if (freqval[Cstart+(row_idx*4)+1] < dtemp1 && dtemp1 <= freqval[Cstart+(row_idx*4)+2]){
        seq[row_idx] = ntdeam[2];
        IsMis5 = 1;
      }
      else if (freqval[Cstart+(row_idx*4)+2] < dtemp1 && dtemp1 <= freqval[Cstart+(row_idx*4)+3]){
        seq[row_idx] = ntdeam[3];
      }
    }
  }
  //3'
  for (int row_idx = 0; row_idx < LEN;row_idx++){
    int row_idx_3p = seqlen-(LEN-row_idx);
    if (seq[row_idx_3p] == 'A' || seq[row_idx_3p] == 'a'){
      if (dtemp1 <= freqval[Aend3-((row_idx)*4)-4]){
        seq[row_idx_3p] = ntdeam[0];
      }
      else if (freqval[Aend3-((row_idx)*4)-4] < dtemp1 && dtemp1 <= freqval[Aend3-((row_idx)*4)-3]){
        seq[row_idx_3p] = ntdeam[1];
        IsMis3 = 2;
      }
      else if (freqval[Aend3-((row_idx)*4)-3] < dtemp1 && dtemp1 <= freqval[Aend3-((row_idx)*4)-2]){
        seq[row_idx_3p] = ntdeam[2];
        IsMis3 = 2;
      }
      else if (freqval[Aend3-((row_idx)*4)-2] < dtemp1 && dtemp1 <= freqval[Aend3-((row_idx)*4)-1]){
        seq[row_idx_3p] = ntdeam[3];
        IsMis3 = 2;
      }
    }
    else if (seq[row_idx_3p] == 'T' || seq[row_idx_3p] == 't'){
      if (dtemp1 <= freqval[Tend3-((row_idx)*4)-4]){
        seq[row_idx_3p] = ntdeam[0];
        IsMis3 = 2;
      }
      else if (freqval[Tend3-((row_idx)*4)-4] < dtemp1 && dtemp1 <= freqval[Tend3-((row_idx)*4)-3]){
        seq[row_idx_3p] = ntdeam[1];
      }
      else if (freqval[Tend3-((row_idx)*4)-3] < dtemp1 && dtemp1 <= freqval[Tend3-((row_idx)*4)-2]){
        seq[row_idx_3p] = ntdeam[2];
        IsMis3 = 2;
      }
      else if (freqval[Tend3-((row_idx)*4)-2] < dtemp1 && dtemp1 <= freqval[Tend3-((row_idx)*4)-1]){
        seq[row_idx_3p] = ntdeam[3];
        IsMis3 = 2;
      }
    }
    else if (seq[row_idx_3p] == 'G' || seq[row_idx_3p] == 'g'){
      if (dtemp1 <= freqval[Gend3-((row_idx)*4)-4]){
        seq[row_idx_3p] = ntdeam[0];
      }
      else if (freqval[Gend3-((row_idx)*4)-4] < dtemp1 && dtemp1 <= freqval[Gend3-((row_idx)*4)-3]){
        seq[row_idx_3p] = ntdeam[1];
        IsMis3 = 2;
      }
      else if (freqval[Gend3-((row_idx)*4)-3] < dtemp1 && dtemp1 <= freqval[Gend3-((row_idx)*4)-2]){
        seq[row_idx_3p] = ntdeam[2];
      }
      else if (freqval[Gend3-((row_idx)*4)-2] < dtemp1 && dtemp1 <= freqval[Gend3-((row_idx)*4)-1]){
        seq[row_idx_3p] = ntdeam[3];
        IsMis3 = 2;
      }
    }
    else if (seq[row_idx_3p] == 'C' || seq[row_idx_3p] == 'c'){
      if (dtemp1 <= freqval[Cend3-((row_idx)*4)-4]){
        seq[row_idx_3p] = ntdeam[0];
        IsMis3 = 2;
      }
      else if (freqval[Cend3-((row_idx)*4)-4] < dtemp1 && dtemp1 <= freqval[Cend3-((row_idx)*4)-3]){
        seq[row_idx_3p] = ntdeam[1];
        IsMis3 = 2;
      }
      else if (freqval[Cend3-((row_idx)*4)-3] < dtemp1 && dtemp1 <= freqval[Cend3-((row_idx)*4)-2]){
        seq[row_idx_3p] = ntdeam[2];
        IsMis3 = 2;
      }
      else if (freqval[Cend3-((row_idx)*4)-2] < dtemp1 && dtemp1 <= freqval[Cend3-((row_idx)*4)-1]){
        seq[row_idx_3p] = ntdeam[3];
      }
    }
  }

  IsMis = IsMis5 + IsMis3;
  //fprintf(stderr," %s \n SEQUENCIN AFTER %d \n-----",seq,IsMis);

  return IsMis;
}


int MisMatchFile_kstring(kstring_t* seq,mrand_t *mr,double* freqval,int LEN){
  int IsMis5 = 0;
  int IsMis3 = 0;
  int IsMis = 0; //if 1 then only 5' change, if 2 only 3' change, if 3 change in both ends

  char ntdeam[4] = {'A', 'T', 'G', 'C'};
  double dtemp1;
  // 5' moving downwards from the 1 postion in the sequence 
  int Astart = 0;
  int Tstart = LEN*4;
  int Gstart = LEN*8;
  int Cstart = LEN*12;

  // moving upstream from the last position
  int Aend3 = LEN*20;
  int Tend3 = LEN*24;
  int Gend3 = LEN*28;
  int Cend3 = LEN*32;

  int seqlen = seq->l;


  kstring_t seq_intermediate;
  seq_intermediate.l = seq->l;
  seq_intermediate.m = seq->l;
  seq_intermediate.s = (char *)malloc((seq->l + 1) * sizeof(char)); // Allocate memory

  strcpy(seq_intermediate.s, seq->s); // Copy seq->s to seq_intermediate.s

  //fprintf(stderr,"SEQUENCIN BEFORE %d \n %s \n",IsMis,seq);
  //5'
  for (int row_idx = 0; row_idx < LEN;row_idx++){
    dtemp1 = mrand_pop(mr);
    if (seq_intermediate.s[row_idx] == 'A' || seq_intermediate.s[row_idx] == 'a'){
      if (dtemp1 <= freqval[Astart+(row_idx*4)]){
        seq->s[row_idx] = ntdeam[0];
      }
      else if (freqval[Astart+(row_idx*4)] < dtemp1 && dtemp1 <= freqval[Astart+(row_idx*4)+1]){
        seq->s[row_idx] = ntdeam[1];
        IsMis5 = 1;
      }
      else if (freqval[Astart+(row_idx*4)+1] < dtemp1 && dtemp1 <= freqval[Astart+(row_idx*4)+2]){
        seq->s[row_idx] = ntdeam[2];
        IsMis5 = 1;
      }
      else if (freqval[Astart+(row_idx*4)+2] < dtemp1 && dtemp1 <= freqval[Astart+(row_idx*4)+3]){
        seq->s[row_idx] = ntdeam[3];
        IsMis5 = 1;
      }
    }
    else if (seq_intermediate.s[row_idx] == 'T' || seq_intermediate.s[row_idx] == 't'){
      if (dtemp1 <= freqval[Tstart+(row_idx*4)]){
        seq->s[row_idx] = ntdeam[0];
        IsMis5 = 1;
      }
      else if (freqval[Tstart+(row_idx*4)] < dtemp1 && dtemp1 <= freqval[Tstart+(row_idx*4)+1]){
        seq->s[row_idx] = ntdeam[1];
      }
      else if (freqval[Tstart+(row_idx*4)+1] < dtemp1 && dtemp1 <= freqval[Tstart+(row_idx*4)+2]){
        seq->s[row_idx] = ntdeam[2];
        IsMis5 = 1;
      }
      else if (freqval[Tstart+(row_idx*4)+2] < dtemp1 && dtemp1 <= freqval[Tstart+(row_idx*4)+3]){
        seq->s[row_idx] = ntdeam[3];
        IsMis5 = 1;
      }
    }
    else if (seq_intermediate.s[row_idx] == 'G' || seq_intermediate.s[row_idx] == 'g'){
      if (dtemp1 <= freqval[Gstart+(row_idx*4)]){
        seq->s[row_idx] = ntdeam[0];
        IsMis5 = 1;
      }
      else if (freqval[Gstart+(row_idx*4)] < dtemp1 && dtemp1 <= freqval[Gstart+(row_idx*4)+1]){
        seq->s[row_idx] = ntdeam[1];
        IsMis5 = 1;
      }
      else if (freqval[Gstart+(row_idx*4)+1] < dtemp1 && dtemp1 <= freqval[Gstart+(row_idx*4)+2]){
        seq->s[row_idx] = ntdeam[2];
      }
      else if (freqval[Gstart+(row_idx*4)+2] < dtemp1 && dtemp1 <= freqval[Gstart+(row_idx*4)+3]){
        seq->s[row_idx] = ntdeam[3];
        IsMis5 = 1;
      }
    }
    else if (seq_intermediate.s[row_idx] == 'C' || seq_intermediate.s[row_idx] == 'c'){
      if (dtemp1 <= freqval[Cstart+(row_idx*4)]){
        seq->s[row_idx] = ntdeam[0];
        IsMis5 = 1;
      }
      else if (freqval[Cstart+(row_idx*4)] < dtemp1 && dtemp1 <= freqval[Cstart+(row_idx*4)+1]){
        seq->s[row_idx] = ntdeam[1];
        IsMis5 = 1;
      }
      else if (freqval[Cstart+(row_idx*4)+1] < dtemp1 && dtemp1 <= freqval[Cstart+(row_idx*4)+2]){
        seq->s[row_idx] = ntdeam[2];
        IsMis5 = 1;
      }
      else if (freqval[Cstart+(row_idx*4)+2] < dtemp1 && dtemp1 <= freqval[Cstart+(row_idx*4)+3]){
        seq->s[row_idx] = ntdeam[3];
      }
    }
  }
  //3'
  for (int row_idx = 0; row_idx < LEN;row_idx++){
    int row_idx_3p = seqlen-(LEN-row_idx);
    if (seq_intermediate.s[row_idx_3p] == 'A' || seq_intermediate.s[row_idx_3p] == 'a'){
      if (dtemp1 <= freqval[Aend3-((row_idx)*4)-4]){
        seq->s[row_idx_3p] = ntdeam[0];
      }
      else if (freqval[Aend3-((row_idx)*4)-4] < dtemp1 && dtemp1 <= freqval[Aend3-((row_idx)*4)-3]){
        seq->s[row_idx_3p] = ntdeam[1];
        IsMis3 = 2;
      }
      else if (freqval[Aend3-((row_idx)*4)-3] < dtemp1 && dtemp1 <= freqval[Aend3-((row_idx)*4)-2]){
        seq->s[row_idx_3p] = ntdeam[2];
        IsMis3 = 2;
      }
      else if (freqval[Aend3-((row_idx)*4)-2] < dtemp1 && dtemp1 <= freqval[Aend3-((row_idx)*4)-1]){
        seq->s[row_idx_3p] = ntdeam[3];
        IsMis3 = 2;
      }
    }
    else if (seq_intermediate.s[row_idx_3p] == 'T' || seq_intermediate.s[row_idx_3p] == 't'){
      if (dtemp1 <= freqval[Tend3-((row_idx)*4)-4]){
        seq->s[row_idx_3p] = ntdeam[0];
        IsMis3 = 2;
      }
      else if (freqval[Tend3-((row_idx)*4)-4] < dtemp1 && dtemp1 <= freqval[Tend3-((row_idx)*4)-3]){
        seq->s[row_idx_3p] = ntdeam[1];
      }
      else if (freqval[Tend3-((row_idx)*4)-3] < dtemp1 && dtemp1 <= freqval[Tend3-((row_idx)*4)-2]){
        seq->s[row_idx_3p] = ntdeam[2];
        IsMis3 = 2;
      }
      else if (freqval[Tend3-((row_idx)*4)-2] < dtemp1 && dtemp1 <= freqval[Tend3-((row_idx)*4)-1]){
        seq->s[row_idx_3p] = ntdeam[3];
        IsMis3 = 2;
      }
    }
    else if (seq_intermediate.s[row_idx_3p] == 'G' || seq_intermediate.s[row_idx_3p] == 'g'){
      if (dtemp1 <= freqval[Gend3-((row_idx)*4)-4]){
        seq->s[row_idx_3p] = ntdeam[0];
      }
      else if (freqval[Gend3-((row_idx)*4)-4] < dtemp1 && dtemp1 <= freqval[Gend3-((row_idx)*4)-3]){
        seq->s[row_idx_3p] = ntdeam[1];
        IsMis3 = 2;
      }
      else if (freqval[Gend3-((row_idx)*4)-3] < dtemp1 && dtemp1 <= freqval[Gend3-((row_idx)*4)-2]){
        seq->s[row_idx_3p] = ntdeam[2];
      }
      else if (freqval[Gend3-((row_idx)*4)-2] < dtemp1 && dtemp1 <= freqval[Gend3-((row_idx)*4)-1]){
        seq->s[row_idx_3p] = ntdeam[3];
        IsMis3 = 2;
      }
    }
    else if (seq_intermediate.s[row_idx_3p] == 'C' || seq_intermediate.s[row_idx_3p] == 'c'){
      if (dtemp1 <= freqval[Cend3-((row_idx)*4)-4]){
        seq->s[row_idx_3p] = ntdeam[0];
        IsMis3 = 2;
      }
      else if (freqval[Cend3-((row_idx)*4)-4] < dtemp1 && dtemp1 <= freqval[Cend3-((row_idx)*4)-3]){
        seq->s[row_idx_3p] = ntdeam[1];
        IsMis3 = 2;
      }
      else if (freqval[Cend3-((row_idx)*4)-3] < dtemp1 && dtemp1 <= freqval[Cend3-((row_idx)*4)-2]){
        seq->s[row_idx_3p] = ntdeam[2];
        IsMis3 = 2;
      }
      else if (freqval[Cend3-((row_idx)*4)-2] < dtemp1 && dtemp1 <= freqval[Cend3-((row_idx)*4)-1]){
        seq->s[row_idx_3p] = ntdeam[3];
      }
    }
  }

  IsMis = IsMis5 + IsMis3;
  //fprintf(stderr," %s \n SEQUENCIN AFTER %d \n-----",seq,IsMis);

  // Cleanup
  free(seq_intermediate.s);
  seq_intermediate.s = NULL;
  seq_intermediate.l = seq_intermediate.m = 0;

  return IsMis;
}

std::map<int, mydataD> load_mismatch(const char* fname,int &printlength){
  //  fprintf(stderr,"./metadamage print file.bdamage.gz [-names file.gz -bam file.bam]\n");
  const char *infile = fname;
  //  fprintf(stderr,"infile: %s howmany: %d \n",infile,howmany);
  
  BGZF *bgfp = NULL;

  if(((bgfp = bgzf_open(infile, "r")))== NULL){
    fprintf(stderr,"Could not open input BAM file: %s\n",infile);
    exit(0);
  }

  std::map<int,mydataD> retmap;
  printlength =0;
  assert(sizeof(int)==bgzf_read(bgfp,&printlength,sizeof(int)));

  int ref_nreads[2];
 
  while(1){
    int nread=bgzf_read(bgfp,ref_nreads,2*sizeof(int));
    if(nread==0)
      break;
    assert(nread==2*sizeof(int));
    mydataD md;
    
    md.fwD = new double[16*printlength];
    md.bwD = new double[16*printlength];
    md.nreads = ref_nreads[1];

    float tmp[16];
    for(int i=0;i<printlength;i++){
      assert(16*sizeof(float)==bgzf_read(bgfp,tmp,sizeof(float)*16));
      for(int ii=0;ii<16;ii++)
	md.fwD[i*16+ii] = tmp[ii];
    }
  
    for(int i=0;i<printlength;i++){
      assert(16*sizeof(float)==bgzf_read(bgfp,tmp,sizeof(float)*16));
      for(int ii=0;ii<16;ii++)
	md.bwD[i*16+ii] = tmp[ii];
    }
    retmap[ref_nreads[0]] = md;
  }

  if(bgfp)
    bgzf_close(bgfp);
  fprintf(stderr,"\t-> Done loading binary bdamage.gz file. It contains: %lu\n",retmap.size());
    #if 0
    for(std::map<int,mydata>::iterator it = retmap.begin();it!=retmap.end();it++)
        fprintf(stderr,"it->second:%p\n",it->second);
    #endif
  return retmap;
}

void parse_mismatch_Data(mydataD &md,double **dat,int howmany) {
    
    for (int i = 0; i < 15; i++) {
        // Populate dat with counts for forward direction
        for (int j = 0; j < 16; j++) {
            dat[i][j] = md.fwD[i *  16 + j]; // Assuming dat contains counts for each nucleotide
            //std::cout << dat[0][j] << std::endl;
        }
        
        // Populate dat with counts for backward direction
        for (int j = 0; j < 16; j++) {
            dat[i+15][j] = md.bwD[i * 16 + j]; // Assuming dat contains counts for each nucleotide
        }
    }
}

void MisMatchMetaFileArray(double* freqval,const char* filename,int &mismatchcyclelength,int &num_elem,const char* fileoutname){    
    double** mm5p, **mm3p;

    // Load data using load_bdamage_full
    std::map<int, mydataD> retmap = load_mismatch(filename,mismatchcyclelength);
    fprintf(stderr, "\t-> %lu mismatch matrices read for %d base pairs\n", retmap.size(), mismatchcyclelength);

    for (std::map<int, mydataD>::iterator it = retmap.begin(); it != retmap.end(); it++) {
        int taxid = it->first;
        mydataD md = it->second;
        if (it->second.nreads == 0)
        continue;
        
        // allocating memory for my mm5p and mm3p such that i can incorporate the tables 
        mm5p = (double**) malloc(mismatchcyclelength * sizeof(double*));
        mm3p = (double**) malloc(mismatchcyclelength * sizeof(double*));
        for (int i = 0; i < mismatchcyclelength-1; i++){
            mm5p[i] =(double *) malloc(16 * sizeof(double));
            mm3p[i] =(double *) malloc(16 * sizeof(double));
        }
        mm5p[mismatchcyclelength-1] =(double *) malloc(16 * sizeof(double));
        mm3p[mismatchcyclelength-1] =(double *) malloc(16 * sizeof(double));
        for (int i=0; i<mismatchcyclelength-1;i++){
            for (int j=0; j<16;j++){
                mm5p[i][j]=0;
                mm3p[i][j]=0;
            }
        }
        for (int j=0; j<16;j++){
            mm5p[mismatchcyclelength-1][j]=0;
            mm3p[mismatchcyclelength-1][j]=0;
        }
        
        int numpos = mismatchcyclelength*2+1;
        int numcolumn = 16;
        double ** Table = (double **) malloc(numpos*(sizeof(double *))); /*I allocate memory here.  If this function is called many times it may be better to move the memmory allocation out of this function*/
        for (int i=0; i<numpos; i++){
            Table[i]=(double *) malloc(numcolumn*(sizeof(double)));
        }  

        parse_mismatch_Data(md, Table, mismatchcyclelength);

        for (int pos=0; pos<mismatchcyclelength; pos++){
            for (int nt_count=0; nt_count<numcolumn;nt_count++){
                mm5p[pos][nt_count] = Table[pos][nt_count];
                mm3p[pos][nt_count] = Table[pos+mismatchcyclelength][nt_count];
            }
        }
        
        double* A_A = (double *)malloc(2*mismatchcyclelength * sizeof(double)+1);
        double* A_C = (double *)malloc(2*mismatchcyclelength * sizeof(double)+1);
        double* A_G = (double *)malloc(2*mismatchcyclelength * sizeof(double)+1);
        double* A_T = (double *)malloc(2*mismatchcyclelength * sizeof(double)+1);
        
        double* C_A = (double *)malloc(2*mismatchcyclelength * sizeof(double)+1);
        double* C_C = (double *)malloc(2*mismatchcyclelength * sizeof(double)+1);
        double* C_G = (double *)malloc(2*mismatchcyclelength * sizeof(double)+1);
        double* C_T = (double *)malloc(2*mismatchcyclelength * sizeof(double)+1);
        
        double* G_A = (double *)malloc(2*mismatchcyclelength * sizeof(double)+1);
        double* G_C = (double *)malloc(2*mismatchcyclelength * sizeof(double)+1);
        double* G_G = (double *)malloc(2*mismatchcyclelength * sizeof(double)+1);
        double* G_T = (double *)malloc(2*mismatchcyclelength * sizeof(double)+1);
        
        double* T_A = (double *)malloc(2*mismatchcyclelength * sizeof(double)+1);
        double* T_C = (double *)malloc(2*mismatchcyclelength * sizeof(double)+1);
        double* T_G = (double *)malloc(2*mismatchcyclelength * sizeof(double)+1);
        double* T_T = (double *)malloc(2*mismatchcyclelength * sizeof(double)+1);
        
        for (int i=0; i<mismatchcyclelength;i++){
            // 0 -> 3 A>A A>C A>G A>T
            A_A[i] = mm5p[i][0]/(mm5p[i][0]+mm5p[i][1]+mm5p[i][2]+mm5p[i][3]);
            A_A[i+mismatchcyclelength] = mm3p[i][0]/(mm3p[i][0]+mm3p[i][1]+mm3p[i][2]+mm3p[i][3]);
            //So in A_A the first index i (0 to 14) is the counts of position 5' 1 to 15, index i+MAXLENGTH (15 to 29) is position 3' 1 to 15
            
            A_T[i] = A_A[i] + mm5p[i][3]/(mm5p[i][0]+mm5p[i][1]+mm5p[i][2]+mm5p[i][3]);
            A_T[i+mismatchcyclelength] = A_A[i+mismatchcyclelength] + mm3p[i][3]/(mm3p[i][0]+mm3p[i][1]+mm3p[i][2]+mm3p[i][3]);
            
            A_G[i] = A_T[i] + mm5p[i][2]/(mm5p[i][0]+mm5p[i][1]+mm5p[i][2]+mm5p[i][3]);
            A_G[i+mismatchcyclelength] = A_T[i+mismatchcyclelength] + mm3p[i][2]/(mm3p[i][0]+mm3p[i][1]+mm3p[i][2]+mm3p[i][3]);
            
            A_C[i] = A_G[i] + mm5p[i][1]/(mm5p[i][0]+mm5p[i][1]+mm5p[i][2]+mm5p[i][3]);
            A_C[i+mismatchcyclelength] = A_G[i+mismatchcyclelength] + mm3p[i][1]/(mm3p[i][0]+mm3p[i][1]+mm3p[i][2]+mm3p[i][3]);

            /*
                AA	AC	AG	AT
            5' 861933	1713	6988	2831
            3' 851392  1724    4909    2217
            */
            /*
            fprintf(stderr,"position 5p %d \t AA %f\t AT %f\t AG %f\t AC %f\n",i,A_A[i],A_T[i],A_G[i],A_C[i]);
            fprintf(stderr,"position 5p %d \t AA %d\t AT %d\t AG %d\t AC %d\n",i,(int)mm5p[i][0],(int)mm5p[i][3],(int)mm5p[i][2],(int)mm5p[i][1]);
            fprintf(stderr,"position 3p %d \t AA %f\t AT %f\t AG %f\t AC %f\n",i,A_A[i+mismatchcyclelength],A_T[i+mismatchcyclelength],A_G[i+mismatchcyclelength],A_C[i+mismatchcyclelength]);
            fprintf(stderr,"position 3p %d \t AA %d\t AT %d\t AG %d\t AC %d\n",i,(int)mm3p[i][0],(int)mm3p[i][3],(int)mm3p[i][2],(int)mm3p[i][1]);
            */
            
            T_A[i] = mm5p[i][12]/(mm5p[i][12]+mm5p[i][13]+mm5p[i][14]+mm5p[i][15]);
            T_A[i+mismatchcyclelength] = mm3p[i][12]/(mm3p[i][12]+mm3p[i][13]+mm3p[i][14]+mm3p[i][15]);

            T_T[i] = T_A[i] + mm5p[i][15]/(mm5p[i][12]+mm5p[i][13]+mm5p[i][14]+mm5p[i][15]);
            T_T[i+mismatchcyclelength] = T_A[i+mismatchcyclelength] + mm3p[i][15]/(mm3p[i][12]+mm3p[i][13]+mm3p[i][14]+mm3p[i][15]);

            T_G[i] = T_T[i] + mm5p[i][14]/(mm5p[i][12]+mm5p[i][13]+mm5p[i][14]+mm5p[i][15]);
            T_G[i+mismatchcyclelength] = T_T[i+mismatchcyclelength] + mm3p[i][14]/(mm3p[i][12]+mm3p[i][13]+mm3p[i][14]+mm3p[i][15]);

            T_C[i] = T_G[i] + mm5p[i][13]/(mm5p[i][12]+mm5p[i][13]+mm5p[i][14]+mm5p[i][15]);
            T_C[i+mismatchcyclelength] = T_G[i+mismatchcyclelength] + mm3p[i][13]/(mm3p[i][12]+mm3p[i][13]+mm3p[i][14]+mm3p[i][15]);

            /*
                TA	  TC	  TG  	TT
            5' 2496	  5012	1885	863057
            3' 2920   6217  2047  865570
            */

            /*
            fprintf(stderr,"position 5p %d \t TA %f\t TT %f\t TG %f\t TC %f\n",i,T_A[i],T_T[i],T_G[i],T_C[i]);
            fprintf(stderr,"position 5p %d \t TA %d\t TT %d\t TG %d\t TC %d\n",i,(int)mm5p[i][12],(int)mm5p[i][15],(int)mm5p[i][14],(int)mm5p[i][13]);
            fprintf(stderr,"position 3p %d \t TA %f\t TT %f\t TG %f\t TC %f\n",i,T_A[i+mismatchcyclelength],T_T[i+mismatchcyclelength],T_G[i+mismatchcyclelength],T_C[i+mismatchcyclelength]);
            fprintf(stderr,"position 3p %d \t TA %d\t TT %d\t TG %d\t TC %d\n",i,(int)mm3p[i][12],(int)mm3p[i][15],(int)mm3p[i][14],(int)mm3p[i][13]);
            */

            G_A[i] = mm5p[i][8]/(mm5p[i][8]+mm5p[i][9]+mm5p[i][10]+mm5p[i][11]);
            G_A[i+mismatchcyclelength] = mm3p[i][8]/(mm3p[i][8]+mm3p[i][9]+mm3p[i][10]+mm3p[i][11]);

            G_T[i] = G_A[i] + mm5p[i][11]/(mm5p[i][8]+mm5p[i][9]+mm5p[i][10]+mm5p[i][11]);
            G_T[i+mismatchcyclelength] = G_A[i+mismatchcyclelength] + mm3p[i][11]/(mm3p[i][8]+mm3p[i][9]+mm3p[i][10]+mm3p[i][11]);

            G_G[i] = G_T[i] + mm5p[i][10]/(mm5p[i][8]+mm5p[i][9]+mm5p[i][10]+mm5p[i][11]);
            G_G[i+mismatchcyclelength] = G_T[i+mismatchcyclelength] + mm3p[i][10]/(mm3p[i][8]+mm3p[i][9]+mm3p[i][10]+mm3p[i][11]);

            G_C[i] = G_G[i] + mm5p[i][9]/(mm5p[i][8]+mm5p[i][9]+mm5p[i][10]+mm5p[i][11]);
            G_C[i+mismatchcyclelength] = G_G[i+mismatchcyclelength] + mm3p[i][9]/(mm3p[i][8]+mm3p[i][9]+mm3p[i][10]+mm3p[i][11]);

            /*
              GA	GC	GG	GT
            5' 6780	629	589369	2895
            3' 87327	1242	492873	2377
            */
            /*
            fprintf(stderr,"position 5p %d \t GA %f\t GT %f\t GG %f\t GC %f\n",i,G_A[i],G_T[i],G_G[i],G_C[i]);
            fprintf(stderr,"position 5p %d \t GA %d\t GT %d\t GG %d\t GC %d\n",i,(int)mm5p[i][8],(int)mm5p[i][11],(int)mm5p[i][10],(int)mm5p[i][9]);
            fprintf(stderr,"position 3p %d \t GA %f\t GT %f\t GG %f\t GC %f\n",i,G_A[i+mismatchcyclelength],G_T[i+mismatchcyclelength],G_G[i+mismatchcyclelength],G_C[i+mismatchcyclelength]);
            fprintf(stderr,"position 3p %d \t GA %d\t GT %d\t GG %d\t GC %d\n",i,(int)mm3p[i][8],(int)mm3p[i][11],(int)mm3p[i][10],(int)mm3p[i][9]);
            */

            C_A[i] = mm5p[i][4]/(mm5p[i][4]+mm5p[i][5]+mm5p[i][6]+mm5p[i][7]);
            C_A[i+mismatchcyclelength] = mm3p[i][4]/(mm3p[i][4]+mm3p[i][5]+mm3p[i][6]+mm3p[i][7]);

            C_T[i] = C_A[i] + mm5p[i][7]/(mm5p[i][4]+mm5p[i][5]+mm5p[i][6]+mm5p[i][7]);
            C_T[i+mismatchcyclelength] = C_A[i+mismatchcyclelength] + mm3p[i][7]/(mm3p[i][4]+mm3p[i][5]+mm3p[i][6]+mm3p[i][7]);
            
            C_G[i] = C_T[i] + mm5p[i][6]/(mm5p[i][4]+mm5p[i][5]+mm5p[i][6]+mm5p[i][7]);
            C_G[i+mismatchcyclelength] = C_T[i+mismatchcyclelength] + mm3p[i][6]/(mm3p[i][4]+mm3p[i][5]+mm3p[i][6]+mm3p[i][7]);

            C_C[i] = C_G[i] + mm5p[i][5]/(mm5p[i][4]+mm5p[i][5]+mm5p[i][6]+mm5p[i][7]);
            C_C[i+mismatchcyclelength] = C_G[i+mismatchcyclelength] + mm3p[i][5]/(mm3p[i][4]+mm3p[i][5]+mm3p[i][6]+mm3p[i][7]);
           
            /*
               CA	CC	CG	CT
            5' 2394	487148	885	89181
            3' 2506	591710	743	9419
            */
            /*
            fprintf(stderr,"position 5p %d \t CA %f\t CT %f\t CG %f\t CC %f\n",i,C_A[i],C_T[i],C_G[i],C_C[i]);
            fprintf(stderr,"position 5p %d \t CA %d\t CT %d\t CG %d\t CC %d\n",i,(int)mm5p[i][4],(int)mm5p[i][7],(int)mm5p[i][6],(int)mm5p[i][5]);
            fprintf(stderr,"position 3p %d \t CA %f\t CT %f\t CG %f\t CC %f\n",i,C_A[i+mismatchcyclelength],C_T[i+mismatchcyclelength],C_G[i+mismatchcyclelength],C_C[i+mismatchcyclelength]);
            fprintf(stderr,"position 3p %d \t CA %d\t CT %d\t CG %d\t CC %d\n",i,(int)mm3p[i][4],(int)mm3p[i][7],(int)mm3p[i][6],(int)mm3p[i][5]);
            */
        }

        // allocating memory for my mm5p and mm3p such that i can incorporate the tables        
        for (int i = 0; i < mismatchcyclelength; i++){
            free(mm5p[i]);
            free(mm3p[i]);
        }
        free(mm5p);
        free(mm3p);

        if(fileoutname!=NULL){
          FILE *MisIncorpFile = fopen(fileoutname, "w");
          for (int elem=0; elem<mismatchcyclelength; ++elem){
            fprintf(MisIncorpFile,"%f\t%f\t%f\t%f\n",A_A[elem],A_T[elem],A_G[elem],A_C[elem]);
          }
          for (int elem=0; elem<mismatchcyclelength; ++elem){
            fprintf(MisIncorpFile,"%f\t%f\t%f\t%f\n",T_A[elem],T_T[elem],T_G[elem],T_C[elem]);
          }
          for (int elem=0; elem<mismatchcyclelength; ++elem){
            fprintf(MisIncorpFile,"%f\t%f\t%f\t%f\n",G_A[elem],G_T[elem],G_G[elem],G_C[elem]);
          }
          for (int elem=0; elem<mismatchcyclelength; ++elem){
            fprintf(MisIncorpFile,"%f\t%f\t%f\t%f\n",C_A[elem],C_T[elem],C_G[elem],C_C[elem]);
          }
          for (int elem=0; elem<mismatchcyclelength; ++elem){
            fprintf(MisIncorpFile,"%f\t%f\t%f\t%f\n",A_A[elem+mismatchcyclelength],A_T[elem+mismatchcyclelength],A_G[elem+mismatchcyclelength],A_C[elem+mismatchcyclelength]);
          }
          for (int elem=0; elem<mismatchcyclelength; ++elem){
            fprintf(MisIncorpFile,"%f\t%f\t%f\t%f\n",T_A[elem+mismatchcyclelength],T_T[elem+mismatchcyclelength],T_G[elem+mismatchcyclelength],T_C[elem+mismatchcyclelength]);
          }
          for (int elem=0; elem<mismatchcyclelength; ++elem){
            fprintf(MisIncorpFile,"%f\t%f\t%f\t%f\n",G_A[elem+mismatchcyclelength],G_T[elem+mismatchcyclelength],G_G[elem+mismatchcyclelength],G_C[elem+mismatchcyclelength]);
          }
          for (int elem=0; elem<mismatchcyclelength; ++elem){
            fprintf(MisIncorpFile,"%f\t%f\t%f\t%f\n",C_A[elem+mismatchcyclelength],C_T[elem+mismatchcyclelength],C_G[elem+mismatchcyclelength],C_C[elem+mismatchcyclelength]);
          }
          
          fclose(MisIncorpFile);
        }
        

        // 5 prime end
        int offset = 0;
        for (int i = 0; i < mismatchcyclelength * 8; i++) {
            int section = i / mismatchcyclelength;
            int index_in_section = i % mismatchcyclelength;
            //fprintf(stderr,"section %d\t idx %d\n",section,index_in_section);
            switch (section){
                case 0:
                    freqval[i * 4] = A_A[index_in_section]; 
                    freqval[i * 4 + 1] = A_T[index_in_section]; 
                    freqval[i * 4 + 2] = A_G[index_in_section]; 
                    freqval[i * 4 + 3] = A_C[index_in_section];
                    //fprintf(stderr, "A index value %d \t %d\n", index_in_section, i * 4);
                    num_elem += 4;
                    break;
                case 1:
                    freqval[i * 4] = T_A[index_in_section]; 
                    freqval[i * 4 + 1] = T_T[index_in_section]; 
                    freqval[i * 4 + 2] = T_G[index_in_section]; 
                    freqval[i * 4 + 3] = T_C[index_in_section];
                    //fprintf(stderr, "T index value %d \t %d\t%d\n", index_in_section, i * 4);
                    num_elem += 4;
                    break;
                case 2:
                    freqval[i * 4] = G_A[index_in_section]; 
                    freqval[i * 4 + 1] = G_T[index_in_section]; 
                    freqval[i * 4 + 2] = G_G[index_in_section]; 
                    freqval[i * 4 + 3] = G_C[index_in_section];
                    //fprintf(stderr, "G index value %d \t %d\n", index_in_section, i * 4);
                    num_elem += 4;
                    break;
                case 3:
                    freqval[i * 4] = C_A[index_in_section]; 
                    freqval[i * 4 + 1] = C_T[index_in_section]; 
                    freqval[i * 4 + 2] = C_G[index_in_section]; 
                    freqval[i * 4 + 3] = C_C[index_in_section];
                    //fprintf(stderr, "C index value %d \t %d\n", index_in_section, i * 4);
                    num_elem += 4;                    
                    break;
                case 4:
                    freqval[i * 4] = A_A[index_in_section+mismatchcyclelength]; 
                    freqval[i * 4 + 1] = A_T[index_in_section+mismatchcyclelength]; 
                    freqval[i * 4 + 2] = A_G[index_in_section+mismatchcyclelength]; 
                    freqval[i * 4 + 3] = A_C[index_in_section+mismatchcyclelength];
                    //fprintf(stderr, "A index value %d \t %d\n", index_in_section, i * 4);
                    num_elem += 4;
                    break;
                case 5:
                    freqval[i * 4] = T_A[index_in_section+mismatchcyclelength]; 
                    freqval[i * 4 + 1] = T_T[index_in_section+mismatchcyclelength]; 
                    freqval[i * 4 + 2] = T_G[index_in_section+mismatchcyclelength]; 
                    freqval[i * 4 + 3] = T_C[index_in_section+mismatchcyclelength];
                    //fprintf(stderr, "T index value %d \t %d\t%d\n", index_in_section, i * 4);
                    num_elem += 4;
                    break;
                case 6:
                    freqval[i * 4] = G_A[index_in_section+mismatchcyclelength]; 
                    freqval[i * 4 + 1] = G_T[index_in_section+mismatchcyclelength]; 
                    freqval[i * 4 + 2] = G_G[index_in_section+mismatchcyclelength]; 
                    freqval[i * 4 + 3] = G_C[index_in_section+mismatchcyclelength];
                    //fprintf(stderr, "G index value %d \t %d\n", index_in_section, i * 4);
                    num_elem += 4;
                    break;
                case 7:
                    freqval[i * 4] = C_A[index_in_section+mismatchcyclelength]; 
                    freqval[i * 4 + 1] = C_T[index_in_section+mismatchcyclelength]; 
                    freqval[i * 4 + 2] = C_G[index_in_section+mismatchcyclelength]; 
                    freqval[i * 4 + 3] = C_C[index_in_section+mismatchcyclelength];
                    //fprintf(stderr, "C index value %d \t %d\n", index_in_section, i * 4);
                    num_elem += 4;
                    break;
            }
        }

        // Free the memory for all arrays
        free(A_A); free(A_T); free(A_G); free(A_C);
        free(T_A); free(T_T); free(T_G); free(T_C);
        free(G_A); free(G_T); free(G_G); free(G_C);
        free(C_A); free(C_T); free(C_G); free(C_C);
    }
}