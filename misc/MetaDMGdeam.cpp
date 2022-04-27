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
  fprintf(fp,"Misincorporation converter - converting MapDamage 2.0 misincorporation files into NGSNGS misincorporation format\n");
  fprintf(fp,"Creates a global genome misincorporation file - not chromosome specific.\n\n");
  fprintf(fp,"Usage\n./MetaMisConvert -i <MapDamage Misincorporation file> -o <NGSNGS misincorporation file>\n");
  fprintf(fp,"\nExample\n./MetaMisConvert -i misincorporation.txt -o NGSNGS_mis.txt\n");
  fprintf(fp,"./MetaMisConvert --MapDMG_in misincorporation.txt --NGSNGS_out NGSNGS_mis.txt\n");
  fprintf(fp,"\nOptions: \n");
  fprintf(fp,"-h   | --help: \t\t\t Print help page.\n");
  fprintf(fp,"-v   | --version: \t\t Print help page.\n\n");
  fprintf(fp,"-i   | --MapDMG_in: \t\t MapDamage misincorporation input file in .txt format or .txt.gz\n");
  fprintf(fp,"-o   | --NGSNGS_out: \t\t NGSNGS misincorporation output file in .txt format\n");
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
    if(strcasecmp("-i",*argv)==0 || strcasecmp("--MapDMG_in",*argv)==0){
      mypars->MetaDMG_mis_in = strdup(*(++argv));
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--NGSNGS_out",*argv)==0){
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
                    //fprintf(stderr,"----------------\nBUFFER1 %s\n",buf);
                    char* tok1;
                    double AA, AT, AG, AC;
                    double TA, TT, TG, TC;
                    double GA, GT, GG, GC;
                    double CA, CT, CG, CC;
                    while(((tok1 = strtok(NULL,"\t")))){
                        if(token_count == 0){pos = atoi(tok1);}
                        else if(token_count == 1){AA = atof(tok1);}
                        else if(token_count == 2){AT = atof(tok1);}
                        else if(token_count == 3){AG = atof(tok1);}
                        else if(token_count == 4){AC = atof(tok1);}

                        else if(token_count == 5){CA = atof(tok1);}
                        else if(token_count == 6){CC = atof(tok1);}
                        else if(token_count == 7){CG = atof(tok1);}
                        else if(token_count == 8){CT = atof(tok1);}

                        else if(token_count == 9){GA = atof(tok1);}
                        else if(token_count == 10){GC = atof(tok1);}
                        else if(token_count == 11){GG = atof(tok1);}
                        else if(token_count == 12){GT = atof(tok1);}

                        else if(token_count == 13){TA = atof(tok1);}
                        else if(token_count == 14){TC = atof(tok1);}
                        else if(token_count == 15){TG = atof(tok1);}
                        else if(token_count == 16){TT = atof(tok1);}
                        token_count++;
                    }
                    //fprintf(stderr,"Test %f \t%f \t%f \t%f\n",AA,(AA+AT),(AA+AT+AG),1.000);
                    sprintf(A5prime+strlen(A5prime),"%f \t%f \t%f \t%f\n",AA,(AA+AT),(AA+AT+AG),1.000);
                    sprintf(T5prime+strlen(T5prime),"%f \t%f \t%f \t%f\n",TA,(TA+TT),(TA+TT+TG),1.000);
                    sprintf(G5prime+strlen(G5prime),"%f \t%f \t%f \t%f\n",GA,(GA+GT),(GA+GT+GG),1.000);
                    sprintf(C5prime+strlen(C5prime),"%f \t%f \t%f \t%f\n",CA,(CA+CT),(CA+CT+CG),1.000);
                }
            }
        }
        gzclose(gz);
        
        FILE *MetaDMGMis = fopen(NGSNGS_output, "w");
        fprintf(MetaDMGMis,"%s",A5prime);
        fprintf(MetaDMGMis,"%s",T5prime);
        fprintf(MetaDMGMis,"%s",G5prime);
        fprintf(MetaDMGMis,"%s",C5prime);
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