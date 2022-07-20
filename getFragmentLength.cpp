#include <cstdlib>
#include <cstring>
#include <cstdio>

#include "getFragmentLength.h"

int getFragmentLength(sim_fragment *sf){
  int res;
  if(sf->type == 0)
    res = sf->fixlength;
  return res;
}

sim_fragment *sim_fragment_alloc(int type,double par1, double par2,const char *fname){
  sim_fragment *fp = new sim_fragment;
  fp->type = type;
  fp->fixlength = par1;
  return fp;
}


#ifdef __WITH_MAIN__

int main(int argc,char **argv){
  int nrep = 10;

  int type = 0;
  double par1 = 88;
  double par2 = 0.2;
  const char *fname = NULL;
  sim_fragment *sf = sim_fragment_alloc(type,par1,par2,fname);

  for(int i=0;i<nrep;i++)
    fprintf(stdout,"%d\n",getFragmentLength(sf));

  return 0;
}
#endif
