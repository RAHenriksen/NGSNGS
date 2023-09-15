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
  const char *MetaDMG_mis_in;
  const char *NGSNGS_mis_out;
}argStruct;

int HelpPage(FILE *fp){
  fprintf(fp,"Misincorporation converter - converting MetaDamage misincorporation files into NGSNGS misincorporation format\n");
  fprintf(fp,"Creates a global genome misincorporation file - not chromosome specific.\n\n");
  fprintf(fp,"Usage\n./MetaMisConvert -i <MetaDamage Misincorporation file> -o <NGSNGS misincorporation file>\n");
  fprintf(fp,"\nExample\n./MetaMisConvert -i MetaNtSub.txt -o NGSNGS_mis.txt\n");
  fprintf(fp,"./MetaMisConvert --input MetaNtSub.txt --output NGSNGS_mis.txt\n");
  fprintf(fp,"\nOptions: \n");
  fprintf(fp,"-h | --help: \t Print help page.\n");
  fprintf(fp,"-v | --version:  Print version.\n\n");
  fprintf(fp,"-i | --input: \t MetaDamage misincorporation input file in .txt or .txt.gz format.\n");
  fprintf(fp,"-o | --output: \t NGSNGS misincorporation output file in .txt format.\n");
  exit(1);
  return 0;
}

argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->MetaDMG_mis_in = NULL;
  mypars->NGSNGS_mis_out = NULL;
  ++argv;
  while(*argv){
    //fprintf(stderr,"ARGV %s\n",*argv);
    if(strcasecmp("-i",*argv)==0 || strcasecmp("--input",*argv)==0){
      mypars->MetaDMG_mis_in = strdup(*(++argv));
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--output",*argv)==0){
      mypars->NGSNGS_mis_out = strdup(*(++argv));
    }
    else{
      fprintf(stderr,"unrecognized input option %s, see NGSNGS help page\n\n",*(argv));
      exit(0);
    }
    ++argv;
  }
  return mypars;
}

int main_Mis(int argc,char **argv){
    argStruct *mypars = NULL;
    if(argc==1||(argc==2&&(strcasecmp(argv[1],"--version")==0||strcasecmp(argv[1],"-v")==0||
                            strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"-h")==0))){
        HelpPage(stderr);
        return 0;
    }
    else{
        mypars = getpars(argc,argv);
        const char* MetaDMG_mis_in = mypars->MetaDMG_mis_in;
        const char* NGSNGS_output = mypars->NGSNGS_mis_out;

        int length = 6000;
        char buf[length];

        gzFile gz = Z_NULL;
        char A5prime[4096];char T5prime[4096];char G5prime[4096];char C5prime[4096];
        char A3prime[4096];char T3prime[4096];char G3prime[4096];char C3prime[4096];
        gz = gzopen(MetaDMG_mis_in,"r");
        assert(gz!=Z_NULL);
        while(gzgets(gz,buf,length)){
          if (buf[0] != '#'){
            strtok(buf,"\t"); // chr1
            strtok(NULL,"\t");
            char* dir = strtok(NULL,"\t"); // 3p
            int token_count = 0;
            int pos = 0;
            if(dir[0]=='5'){
              char* tok1;
              double AA5, AT5, AG5, AC5;
              double TA5, TT5, TG5, TC5;
              double GA5, GT5, GG5, GC5;
              double CA5, CT5, CG5, CC5;
              while(((tok1 = strtok(NULL,"\t")))){
                if(token_count == 0){pos = atoi(tok1);}
                else if(token_count == 1){AA5 = atof(tok1);}
                else if(token_count == 2){AC5 = atof(tok1);}
                else if(token_count == 3){AG5 = atof(tok1);}
                else if(token_count == 4){AT5 = atof(tok1);}

                else if(token_count == 5){CA5 = atof(tok1);}
                else if(token_count == 6){CC5 = atof(tok1);}
                else if(token_count == 7){CG5 = atof(tok1);}
                else if(token_count == 8){CT5 = atof(tok1);}

                else if(token_count == 9){GA5 = atof(tok1);}
                else if(token_count == 10){GC5 = atof(tok1);}
                else if(token_count == 11){GG5 = atof(tok1);}
                else if(token_count == 12){GT5 = atof(tok1);}

                else if(token_count == 13){TA5 = atof(tok1);}
                else if(token_count == 14){TC5 = atof(tok1);}
                else if(token_count == 15){TG5 = atof(tok1);}
                else if(token_count == 16){TT5 = atof(tok1);}
                token_count++;
              }
              sprintf(A5prime+strlen(A5prime),"%lf \t%lf \t%lf \t%lf\n",AA5,(AA5+AT5),(AA5+AT5+AG5),1.000);
              sprintf(T5prime+strlen(T5prime),"%lf \t%lf \t%lf \t%lf\n",TA5,(TA5+TT5),(TA5+TT5+TG5),1.000);
              sprintf(G5prime+strlen(G5prime),"%lf \t%lf \t%lf \t%lf\n",GA5,(GA5+GT5),(GA5+GT5+GG5),1.000);
              sprintf(C5prime+strlen(C5prime),"%lf \t%lf \t%lf \t%lf\n",CA5,(CA5+CT5),(CA5+CT5+CG5),1.000);
            }
            else if(dir[0]=='3'){
              char* tok1;
              double AA3, AT3, AG3, AC3;
              double TA3, TT3, TG3, TC3;
              double GA3, GT3, GG3, GC3;
              double CA3, CT3, CG3, CC3;
              while(((tok1 = strtok(NULL,"\t")))){
                if(token_count == 0){pos = atoi(tok1);}
                else if(token_count == 1){AA3 = atof(tok1);}
                else if(token_count == 2){AT3 = atof(tok1);}
                else if(token_count == 3){AG3 = atof(tok1);}
                else if(token_count == 4){AC3 = atof(tok1);}

                else if(token_count == 5){CA3 = atof(tok1);}
                else if(token_count == 6){CC3 = atof(tok1);}
                else if(token_count == 7){CG3 = atof(tok1);}
                else if(token_count == 8){CT3 = atof(tok1);}

                else if(token_count == 9){GA3 = atof(tok1);}
                else if(token_count == 10){GC3 = atof(tok1);}
                else if(token_count == 11){GG3 = atof(tok1);}
                else if(token_count == 12){GT3 = atof(tok1);}

                else if(token_count == 13){TA3 = atof(tok1);}
                else if(token_count == 14){TC3 = atof(tok1);}
                else if(token_count == 15){TG3 = atof(tok1);}
                else if(token_count == 16){TT3 = atof(tok1);}
                token_count++;
              }
              sprintf(A3prime+strlen(A3prime),"%lf \t%lf \t%lf \t%lf\n",AA3,(AA3+AT3),(AA3+AT3+AG3),1.000);
              sprintf(T3prime+strlen(T3prime),"%lf \t%lf \t%lf \t%lf\n",TA3,(TA3+TT3),(TA3+TT3+TG3),1.000);
              sprintf(G3prime+strlen(G3prime),"%lf \t%lf \t%lf \t%lf\n",GA3,(GA3+GT3),(GA3+GT3+GG3),1.000);
              sprintf(C3prime+strlen(C3prime),"%lf \t%lf \t%lf \t%lf\n",CA3,(CA3+CT3),(CA3+CT3+CG3),1.000);
            }
          }
        }
        gzclose(gz);
        
        FILE *MetaDMGMis = fopen(NGSNGS_output, "w");
        fprintf(MetaDMGMis,"%s",A5prime);
        fprintf(MetaDMGMis,"%s",T5prime);
        fprintf(MetaDMGMis,"%s",G5prime);
        fprintf(MetaDMGMis,"%s",C5prime);
        fprintf(MetaDMGMis,"%s",A3prime);
        fprintf(MetaDMGMis,"%s",T3prime);
        fprintf(MetaDMGMis,"%s",G3prime);
        fprintf(MetaDMGMis,"%s",C3prime);
        fclose(MetaDMGMis);
    }
    return 0;
}

#ifdef __WITH_MAIN__
int main(int argc,char **argv){
    main_Mis(argc,argv); 
    return 0;
}
#endif
//g++ MapDMGmis.cpp -std=c++11 -lm -lz -D __WITH_MAIN__ -o MetaMisConvert
