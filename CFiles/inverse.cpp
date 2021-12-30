#include <algorithm>
#include <cstdio>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>//for printing time

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>

#include <cstdlib>
#include <ctime>

#include <cstdio>
#include <cassert>
#include <cstdint>

#include <random>
#include <iterator>
#include <cmath>

#include <thread>         // std::thread
#include <mutex>        
#include <atomic>
#include <vector>

#include <getopt.h>
#define LENS 4096
double uniform()
{
    /*
     U(0,1): AS 183: Appl. Stat. 31:188-190
     Wichmann BA & Hill ID.  1982.  An efficient and portable
     pseudo-random number generator.  Appl. Stat. 31:188-190
     x, y, z are any numbers in the range 1-30000.  Integer operation up
     to 30323 required.
     
     Suggested to me by Ziheng Yang who also provided me with
     the source code used here.  I use it because it is both fast and portable.
     */
    static int x_rndu=11, y_rndu=23;int static z_rndu=137;

    double r;
    
    x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
    y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
    z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
    if (x_rndu<0) x_rndu+=30269;
    if (y_rndu<0) y_rndu+=30307;
    if (z_rndu<0) z_rndu+=30323;
    r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
    return (r-(int)r);
}

int BinarySearch_fraglength(double* SearchArray,int low, int high, double key)
{
    int ans = 0; 
    while (low <= high) {
        int mid = low + (high - low + 1) / 2;
        //fprintf(stderr,"test %lf\n",SearchArray[mid]);
        double midVal = SearchArray[mid];
 
        if (midVal < key) {
            ans = mid;
            low = mid + 1;
        }
        else if (midVal > key) {

            high = mid - 1;
        }
        else if (midVal == key) {
 
            high = mid - 1;
        }
    }
 
    return ans+1;
}

void FragArray(int& number,int*& Length, double*& Frequency,const char* filename){
  int* Frag_len = new int[LENS];
  double* Frag_freq = new double[LENS];
  int n =0;

  gzFile gz = Z_NULL;
  char buf[LENS];
  
  gz = gzopen(filename,"r");
  assert(gz!=Z_NULL);
  while(gzgets(gz,buf,LENS)){
    Frag_len[n] = atoi(strtok(buf,"\n\t "));
    Frag_freq[n] = atof(strtok(NULL,"\n\t "));
    n++;
  }
  gzclose(gz); 

  number = n;
  Length = Frag_len;
  Frequency = Frag_freq;
}

int Random_geometric_k(unsigned int  seed, const double p)
{
  double u = ((double) rand_r(&seed)/ RAND_MAX);
  int k;

  if (p == 1){k = 1;}
  else if(p == 0){k=0;}
  else{k = log (u) / log (1 - p);}

  return k;
}


double myrand(unsigned int persistent){
  return ((double) rand_r(&persistent)/ RAND_MAX);
}

// ------------------------------ //
/*
int main(int argc,char **argv){
  clock_t t=clock();
  time_t t2=time(NULL);
  
  unsigned int loc_seed = 4;
  
  char nt_qual[8] = {2+33,6+33,15+33,22+33,27+33,33+33,37+33,40+33};

  const char* nt_out;
  nt_out = &nt_qual[2];
  fprintf(stderr,"test %s lol\n",nt_out);
  return 0;
}
*/

  /*
  FILE *fp1;
  fp1 = fopen("geom_dist106v2.txt","wb");
  for(int i = 0; i<1e6; i++){
    int test = Random_geometric_k(loc_seed+i,0.36);
    fprintf(fp1,"%d\n",test);
  }
  fclose(fp1);*/
  /*unsigned int loc_seed = 4;
  double r;
  
 
    for(int i = 0; i<1e8; i++){
    double rand_val = myrand(loc_seed+i);
    fprintf(fp1,"Random key %lf\n ",rand_val);
  }

  }*/



  /*for(int i = 0; i<10; i++){
    double rand_val = ((double) rand_r(&loc_seed)/ RAND_MAX);
    fprintf(stderr,"Random key %lf\n ",rand_val);
  }*/

  /*
  int* Frag_len = new int[LENS];
  double* Frag_freq = new double[LENS];
  int number;

  FragArray(number,Frag_len,Frag_freq,"Size_dist/Size_dist_sampling.txt");

  for(int i = 0; i<10; i++){
    fprintf(stderr," %d\n",Frag_len[i]);
    fprintf(stderr," %lf\n",Frag_freq[i]);
  }

  for(int i = 0; i<4; i++){
    double rand_val = ((double) rand_r(&loc_seed)/ RAND_MAX);
    int guess = BinarySearch_fraglength(Frag_freq,0, number - 1, rand_val);
    fprintf(stderr,"Random key %lf, integer val %d with lookup %d\n ",rand_val,guess,Frag_len[guess]);
  }

  double keytest = 0.00;//((double) rand_r(&loc_seed)/ RAND_MAX);;
  int guess = BinarySearch_fraglength(Frag_freq,0, number - 1, keytest);
  fprintf(stderr,"Random key %lf, integer val %d with lookup %d\n ",keytest,guess,Frag_len[guess]);
  
  keytest =1;//((double) rand_r(&loc_seed)/ RAND_MAX);;
  guess = BinarySearch_fraglength(Frag_freq,0, number - 1, keytest);
  fprintf(stderr,"Random key %lf, integer val %d with lookup %d\n ",keytest,guess,Frag_len[guess]);
  
  std::cout << "-----------" << std::endl;
  */


/*
gzFile gz = Z_NULL;
  char buf[LENS];
  
  gz = gzopen("Size_dist/Size_dist_sampling.txt","r");
  assert(gz!=Z_NULL);

  int* Frag_len = new int[LENS];
  double* Frag_freq = new double[LENS];

  int n =0;
  while(gzgets(gz,buf,LENS)){
    if(strlen(buf)>4096-5){
      fprintf(stderr,"\t-> Problem parsing inputfile, increase buffer size\n");
      return 1;
    }
    Frag_len[n] = atoi(strtok(buf,"\n\t "));
    Frag_freq[n] = atof(strtok(NULL,"\n\t "));
    n++;
  }
  gzclose(gz); 


  double keytest = 0.00;//((double) rand_r(&loc_seed)/ RAND_MAX);;
  int guess = BinarySearch_fraglength(Frag_freq,0, n - 1, keytest);
  fprintf(stderr,"Random key %lf, integer val %d with lookup %d\n ",keytest,guess,Frag_len[guess]);
  
  double rand_val = 0;
  for(int i = 0; i<4; i++){
    double rand_val = ((double) rand_r(&loc_seed)/ RAND_MAX);
    int guess = BinarySearch_fraglength(Frag_freq,0, n - 1, rand_val);
    fprintf(stderr,"Random key %lf, integer val %d with lookup %d\n ",rand_val,guess,Frag_len[guess]);
  }*/

typedef struct{
  int reads;
  int threads;
  char *out;
}pars;

int HelpPage(FILE *fp){
  fprintf(fp,"./ngsngs [-type1,type2] [-typeX int,-out prefix]\n");
  return 0;
}

int print_pars(pars *p,FILE *fp){
  fprintf(fp,"out: %s typeis: %d\n",p->out,p->reads);
  return 0;
}

int print_type(int type){
  fprintf(stderr,"print_type is %d\n ",type * 2);
  return 0;
}

//returtype navn (par1,par2,par3)
pars *getpars(int argc,char ** argv){
  pars *mypars = new pars;
  while(*argv){
    //fprintf(stderr,"thisarg: %s \n",*argv);
    if(strcasecmp("-nreads",*argv)==0)
      mypars->reads = atoi(*(++argv));
      print_type(mypars->reads);
    if(strcasecmp("-nthreads",*argv)==0)
      mypars->threads = atoi(*(++argv));
      fprintf(stderr,"%d",mypars->threads);
    if(strcasecmp("-out",*argv)==0)
      mypars->out = strdup(*(++argv));
    ++argv;
  }
  return mypars;
}


int main(int argc,char**argv){
  pars *mypars = NULL;
  //if haven't been supplied with arguments, load default,print, and exit
  if(argc==1||(argc==2&&(strcasecmp(argv[1],"--version")==0||strcasecmp(argv[1],"--v")==0||
                        strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"--h")==0))){
    HelpPage(stderr);
    return 0;
  }
  else
    mypars = getpars(argc,argv);

  //print_pars(mypars,stdout);
}