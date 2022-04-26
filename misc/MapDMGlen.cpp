#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm> //std::reverse
#include <iostream>
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>
#include <errno.h>
#include <stdint.h>
#include <random>

#include <pthread.h>

typedef struct{
  const char *MapDMG_len_in;
  const char *NGSNGS_len_out;
}argStruct;

int HelpPage(FILE *fp){
  fprintf(fp,"Fragment Length converter - converting MapDamage 2.0 lgdistribution files into NGSNGS length distribution format\n\n");
  fprintf(fp,"Usage\n./LenConvert -i <MapDamage Length Distribution> -o <NGSNGS length distribution file>\n");
  fprintf(fp,"\nExample\n./LenConvert -i lgdistribution.txt -o Ancient_hg19_CDF.txt\n");
  fprintf(fp,"./LenConvert --MapDMG_in lgdistribution.txt --NGSNGS_out Ancient_hg19_CDF.txt\n");
  fprintf(fp,"\nOptions: \n");
  fprintf(fp,"-h   | --help: \t\t\t Print help page.\n");
  fprintf(fp,"-v   | --version: \t\t Print help page.\n\n");
  fprintf(fp,"-i   | --MapDMG_in: \t\t ART platform profile intput in .txt format or .txt.gz\n");
  fprintf(fp,"-o   | --NGSNGS_out: \t\t NGSNGS nucletide quality profile in .txt format\n");
  exit(1);
  return 0;
}

argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->MapDMG_len_in = NULL;
  mypars->NGSNGS_len_out = NULL;
  ++argv;
  while(*argv){
    //fprintf(stderr,"ARGV %s\n",*argv);
    if(strcasecmp("-i",*argv)==0 || strcasecmp("--MapDMG_in",*argv)==0){
      mypars->MapDMG_len_in = strdup(*(++argv));
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--NGSNGS_out",*argv)==0){
      mypars->NGSNGS_len_out = strdup(*(++argv));
    }
    else{
      fprintf(stderr,"unrecognized input option %s, see NGSNGS help page\n\n",*(argv));
      exit(0);
    }
    ++argv;
  }
  return mypars;
}


int main_Length(int argc,char **argv){
    argStruct *mypars = NULL;
    if(argc==1||(argc==2&&(strcasecmp(argv[1],"--version")==0||strcasecmp(argv[1],"-v")==0||
                            strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"-h")==0))){
        HelpPage(stderr);
        return 0;
    }
    else{
        mypars = getpars(argc,argv);
        const char* MapDMG_input = mypars->MapDMG_len_in;
        const char* NGSNGS_output = mypars->NGSNGS_len_out;

        int length = 6000;
        char buf[length];
        gzFile gz = Z_NULL;

        int row_no = 0;
        int* Len_Val_tmp = new int[256];
        int* Len_Obs_tmp = new int[256];
        int Len_sum = 0;

        int Len_Val[256];
        int Len_Obs[256];

        gz = gzopen(MapDMG_input,"r");
        assert(gz!=Z_NULL);
        while(gzgets(gz,buf,length)){
            if (buf[0] == '+' || buf[0] == '-'){
                // buf :  +	50	10074
                char* tok_len; char* tok_count;
                strtok(buf,"\t"); // +
                tok_len = strtok(NULL,"\t"); // 50
                Len_Val_tmp[atoi(tok_len)] = atoi(tok_len);

                tok_count = strtok(NULL,"\t"); // 10074
                Len_sum += atoi(tok_count);
                Len_Obs_tmp[atoi(tok_len)] += atoi(tok_count);
            }
            else{
                continue;
            }   
        }
        gzclose(gz);

        double CumCount = 0;

        //Identify the lowest fragment length
        int low_len = 0;
        for (int elem=0; elem<256; ++elem){if (Len_Val_tmp[elem] != 0){low_len = Len_Val_tmp[elem];break;}};
        low_len--; // substract with one 

        FILE *ReadSizeFile = fopen(NGSNGS_output, "w");
        fprintf(ReadSizeFile,"%d \t %f\n",low_len,CumCount);
        for (int elem=0; elem<256; ++elem){
            if (Len_Val_tmp[elem] != 0){
                CumCount += (double)Len_Obs_tmp[elem]/(double)Len_sum;
                fprintf(ReadSizeFile,"%d \t %f\n",Len_Val_tmp[elem],CumCount);
            }
        }
        fclose(ReadSizeFile);
        
        delete[] Len_Val_tmp;
        delete[] Len_Obs_tmp;

    }
    
    return 0;
}

#ifdef __WITH_MAIN__
int main(int argc,char **argv){
    main_Length(argc,argv); 
    return 0;
}
#endif
//g++ MapDMGlen.cpp -std=c++11 -lm -lz -D __WITH_MAIN__ -o LenConvert