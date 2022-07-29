#include <zlib.h>
#include "NGSNGS_func.h"
#include "NtSubModels.h"
#include "mrand.h"

#define LENS 4096
#define MAXBINS 100

void ErrorSub(double randval,char seqchar[], int pos){
  // Generates nucleotide substitutions
  //fprintf(stderr,"SUBERROR\n"); X, W, Z, Y
  if (seqchar[pos] == 'A' || seqchar[pos] == 'a'){
    if (0 < randval && randval <= 1.0/3.0){seqchar[pos] = 'C';} //X 
    else if (1.0/3.0 < randval && randval <= 2.0/3.0){seqchar[pos] = 'G';} //T
    else if (2.0/3.0 < randval && randval <= 1){seqchar[pos]  = 'T';} //C
  }
  else if (seqchar[pos] == 'C'|| seqchar[pos] == 'c'){ //'Z'
    if (0 < randval && randval <= 1.0/3.0){seqchar[pos] = 'G';} //Z
    else if (1.0/3.0 < randval && randval <= 2.0/3.0){seqchar[pos]  = 'T';} //G
    else if (2.0/3.0 < randval && randval <= 1){seqchar[pos]  = 'A';} //G
  }
  else if (seqchar[pos] == 'G'|| seqchar[pos] == 'g'){
    if (0 < randval && randval <= 1.0/3.0){seqchar[pos] = 'T';} //A
    else if (1.0/3.0 < randval && randval <= 2.0/3.0){seqchar[pos]  = 'A';} //T
    else if (2.0/3.0 < randval && randval <= 1){seqchar[pos]  = 'C';} //C
  }
  else if (seqchar[pos] == 'T'|| seqchar[pos] == 't'){
    if (0 < randval && randval <= 1.0/3.0){seqchar[pos] = 'A';} //G
    else if (1.0/3.0 < randval && randval <= 2.0/3.0){seqchar[pos]  = 'C';} //A
    else if (2.0/3.0 < randval && randval <= 1){seqchar[pos]  = 'G';} //C
  }
}

double* MisMatchFileArray(double* freqval,const char* filename,int &mismatchcyclelength){
    char buf[LENS];
    int i = 0;
    gzFile gz = Z_NULL;
    gz = gzopen(filename,"r");
    assert(gz!=Z_NULL);
    while(gzgets(gz,buf,LENS)){
        double val1;double val2;double val3;double val4;
        sscanf(buf,"%lf\t%lf\t%lf\t%lf\n",&val1,&val2,&val3,&val4);
        //std::cout << "iter " << i << std::endl;// " " << buf << std::endl;
        freqval[i*4] = val1; freqval[i*4+1] = val2; freqval[i*4+2] = val3; freqval[i*4+3] = val4;
        i++;
    }
    gzclose(gz);
    mismatchcyclelength = i/8;
  return freqval;
}

void MisMatchFile(char seq[],mrand_t *mr,double* freqval,int LEN){
  char ntdeam[4] = {'A', 'T', 'G', 'C'};//{'R', 'Q', 'S', 'U'};//{'X', 'Y', 'Z', 'W'}; //{'A', 'T', 'G', 'C'};
  double dtemp1;
  // 5' moving downwards from the 1 postion in the sequence 
  int Astart = 0;
  int Tstart = LEN*4; //4*15
  int Gstart = LEN*8; //15*8
  int Cstart = LEN*12; //15*12

  // moving upstream from the last position
  int Aend3 = LEN*20;
  int Tend3 = LEN*24; //4*15
  int Gend3 = LEN*28; //15*8
  int Cend3 = LEN*32; //15*12

  int seqlen = strlen(seq);
  //5'
  for (int row_idx = 0; row_idx < LEN;row_idx++){
    dtemp1 = mrand_pop(mr);//0.99;// mrand_pop(mr);
    //fprintf(stderr,"RANDOM VALUE %lf\n",dtemp1);
    if (seq[row_idx] == 'A' || seq[row_idx] == 'a'){
      if (dtemp1 <= freqval[Astart+(row_idx*4)]){seq[row_idx] = ntdeam[0];}
      else if (freqval[Astart+(row_idx*4)] < dtemp1 && dtemp1 <= freqval[Astart+(row_idx*4)+1]){seq[row_idx] = ntdeam[1];}
      else if (freqval[Astart+(row_idx*4)+1] < dtemp1 && dtemp1 <= freqval[Astart+(row_idx*4)+2]){seq[row_idx] = ntdeam[2];}
      else if (freqval[Astart+(row_idx*4)+2] < dtemp1 && dtemp1 <= freqval[Astart+(row_idx*4)+3]){seq[row_idx] = ntdeam[3];}
    }
    else if (seq[row_idx] == 'T' || seq[row_idx] == 't'){
      if (dtemp1 <= freqval[Tstart+(row_idx*4)]){seq[row_idx] = ntdeam[0];}
      else if (freqval[Tstart+(row_idx*4)] < dtemp1 && dtemp1 <= freqval[Tstart+(row_idx*4)+1]){seq[row_idx] = ntdeam[1];}
      else if (freqval[Tstart+(row_idx*4)+1] < dtemp1 && dtemp1 <= freqval[Tstart+(row_idx*4)+2]){seq[row_idx] = ntdeam[2];}
      else if (freqval[Tstart+(row_idx*4)+2] < dtemp1 && dtemp1 <= freqval[Tstart+(row_idx*4)+3]){seq[row_idx] = ntdeam[3];}
    }
    else if (seq[row_idx] == 'G' || seq[row_idx] == 'g'){
      if (dtemp1 <= freqval[Gstart+(row_idx*4)]){seq[row_idx] = ntdeam[0];}
      else if (freqval[Gstart+(row_idx*4)] < dtemp1 && dtemp1 <= freqval[Gstart+(row_idx*4)+1]){seq[row_idx] = ntdeam[1];}
      else if (freqval[Gstart+(row_idx*4)+1] < dtemp1 && dtemp1 <= freqval[Gstart+(row_idx*4)+2]){seq[row_idx] = ntdeam[2];}
      else if (freqval[Gstart+(row_idx*4)+2] < dtemp1 && dtemp1 <= freqval[Gstart+(row_idx*4)+3]){seq[row_idx] = ntdeam[3];}
    }
    else if (seq[row_idx] == 'C' || seq[row_idx] == 'c'){
      //fprintf(stderr,"FREQ VAL INDEX %d\n",Cstart+(row_idx));
      if (dtemp1 <= freqval[Cstart+(row_idx*4)]){seq[row_idx] = ntdeam[0];}
      else if (freqval[Cstart+(row_idx*4)] < dtemp1 && dtemp1 <= freqval[Cstart+(row_idx*4)+1]){seq[row_idx] = ntdeam[1];}
      else if (freqval[Cstart+(row_idx*4)+1] < dtemp1 && dtemp1 <= freqval[Cstart+(row_idx*4)+2]){seq[row_idx] = ntdeam[2];}
      else if (freqval[Cstart+(row_idx*4)+2] < dtemp1 && dtemp1 <= freqval[Cstart+(row_idx*4)+3]){seq[row_idx] = ntdeam[3];}
    }
  }
  //3'
  for (int row_idx = 0; row_idx < LEN;row_idx++){
    int row_idx_3p = seqlen-(LEN-row_idx);
    if (seq[row_idx_3p] == 'A' || seq[row_idx_3p] == 'a'){
      if (dtemp1 <= freqval[Aend3-((row_idx)*4)-4]){seq[row_idx_3p] = ntdeam[0];}
      else if (freqval[Aend3-((row_idx)*4)-4] < dtemp1 && dtemp1 <= freqval[Aend3-((row_idx)*4)-3]){seq[row_idx_3p] = ntdeam[1];}
      else if (freqval[Aend3-((row_idx)*4)-3] < dtemp1 && dtemp1 <= freqval[Aend3-((row_idx)*4)-2]){seq[row_idx_3p] = ntdeam[2];}
      else if (freqval[Aend3-((row_idx)*4)-2] < dtemp1 && dtemp1 <= freqval[Aend3-((row_idx)*4)-1]){seq[row_idx_3p] = ntdeam[3];}
    }
    else if (seq[row_idx_3p] == 'T' || seq[row_idx_3p] == 't'){
      if (dtemp1 <= freqval[Tend3-((row_idx)*4)-4]){seq[row_idx_3p] = ntdeam[0];}
      else if (freqval[Tend3-((row_idx)*4)-4] < dtemp1 && dtemp1 <= freqval[Tend3-((row_idx)*4)-3]){seq[row_idx_3p] = ntdeam[1];}
      else if (freqval[Tend3-((row_idx)*4)-3] < dtemp1 && dtemp1 <= freqval[Tend3-((row_idx)*4)-2]){seq[row_idx_3p] = ntdeam[2];}
      else if (freqval[Tend3-((row_idx)*4)-2] < dtemp1 && dtemp1 <= freqval[Tend3-((row_idx)*4)-1]){seq[row_idx_3p] = ntdeam[3];}
    }
    else if (seq[row_idx_3p] == 'G' || seq[row_idx_3p] == 'g'){
      if (dtemp1 <= freqval[Gend3-((row_idx)*4)-4]){seq[row_idx_3p] = ntdeam[0];}
      else if (freqval[Gend3-((row_idx)*4)-4] < dtemp1 && dtemp1 <= freqval[Gend3-((row_idx)*4)-3]){seq[row_idx_3p] = ntdeam[1];}
      else if (freqval[Gend3-((row_idx)*4)-3] < dtemp1 && dtemp1 <= freqval[Gend3-((row_idx)*4)-2]){seq[row_idx_3p] = ntdeam[2];}
      else if (freqval[Gend3-((row_idx)*4)-2] < dtemp1 && dtemp1 <= freqval[Gend3-((row_idx)*4)-1]){seq[row_idx_3p] = ntdeam[3];}
    }
    else if (seq[row_idx_3p] == 'C' || seq[row_idx_3p] == 'c'){
      if (dtemp1 <= freqval[Cend3-((row_idx)*4)-4]){seq[row_idx_3p] = ntdeam[0];}
      else if (freqval[Cend3-((row_idx)*4)-4] < dtemp1 && dtemp1 <= freqval[Cend3-((row_idx)*4)-3]){seq[row_idx_3p] = ntdeam[1];}
      else if (freqval[Cend3-((row_idx)*4)-3] < dtemp1 && dtemp1 <= freqval[Cend3-((row_idx)*4)-2]){seq[row_idx_3p] = ntdeam[2];}
      else if (freqval[Cend3-((row_idx)*4)-2] < dtemp1 && dtemp1 <= freqval[Cend3-((row_idx)*4)-1]){seq[row_idx_3p] = ntdeam[3];}
    }
  }
}
