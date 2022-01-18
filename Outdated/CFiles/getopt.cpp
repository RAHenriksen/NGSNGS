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
#include "getopt_func.h"

void Filecreate(const char* filename,int num){
  FILE *fp1;
  char file[80];
  strcpy(file,filename);
  strcat(file,".txt");
  fp1 = fopen(file,"wb");
  for(int i = 0; i<5; i++){
    fprintf(fp1,"Input test %d\t%d\n",i,num);
  }
  fclose(fp1);
}


int print_type(int type){
  fprintf(stderr,"print_type is %d\n ",type * 2);
  return 0;
}
