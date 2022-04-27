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
  const char *MapDMG_mis_in;
  const char *NGSNGS_mis_out;
}argStruct;

int HelpPage(FILE *fp){
  fprintf(fp,"Misincorporation converter - converting MapDamage 2.0 misincorporation files into NGSNGS misincorporation format\n");
  fprintf(fp,"Creates a global genome misincorporation file - not chromosome specific.\n\n");
  fprintf(fp,"Usage\n./MisConvert -i <MapDamage Misincorporation file> -o <NGSNGS misincorporation file>\n");
  fprintf(fp,"\nExample\n./MisConvert -i misincorporation.txt -o NGSNGS_mis.txt\n");
  fprintf(fp,"./MisConvert --MapDMG_in misincorporation.txt --NGSNGS_out NGSNGS_mis.txt\n");
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
  mypars->MapDMG_mis_in = NULL;
  mypars->NGSNGS_mis_out = NULL;
  ++argv;
  while(*argv){
    //fprintf(stderr,"ARGV %s\n",*argv);
    if(strcasecmp("-i",*argv)==0 || strcasecmp("--MapDMG_in",*argv)==0){
      mypars->MapDMG_mis_in = strdup(*(++argv));
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
        const char* MapDMG_input = mypars->MapDMG_mis_in;
        const char* NGSNGS_output = mypars->NGSNGS_mis_out;

        int length = 6000;
        char buf[length];

        gzFile gz = Z_NULL;

        /*
        # table produced by mapDamage version 2.1.0
        # using mapped file SampleResult/align/Sim_2_s_se.bam and /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/hg19canon.fa as reference file
        # Chr: reference from sam/bam header, End: from which termini of DNA sequences, Std: strand of reads
        Chr	End	Std	Pos	A	C	G	T	Total	G>A	C>T	A>G	T>C	A>C	A>T	C>G	C>A	T>G	T>A	G>C	G>T	A>-	T>-	C>-	G>-	->A	->T	->C	->G	S
        */
        gz = gzopen(MapDMG_input,"r");
        assert(gz!=Z_NULL);
        std::vector<char *> all_chr_lines;   
        while(gzgets(gz,buf,length)){
            if (buf[0] != '#' && buf[0] != 'C'){//disregard the first couple of lines staring with '#' (description) or 'C' (column names)
                all_chr_lines.push_back(strdup(strtok(buf,"\t"))); //assign all chromosomes to a vector, this includes chr1	3p	+, chr1	3p	-, chr1	5p	+, chr1	5p	-
            }
        }
        gzclose(gz); 
        //occurence of the first chromosome -> i.e. number of positions for 5p+,5p-,3p+,3p-
        int pos_all = 0;
        for (unsigned long i = 0; i < all_chr_lines.size();i++){
            if(strcasecmp(all_chr_lines[0],all_chr_lines[i])==0){pos_all++;}
            else{break;}
        }
        int pos_val = pos_all/4;

        int* A_Count = new int[256]();int* A_T_Count = new int[256]();int* A_G_Count = new int[256]();int* A_C_Count = new int[256]();int* A_sum = new int[256]();
        
        int* T_A_Count = new int[256]();int* T_Count = new int[256]();int* T_G_Count = new int[256]();int* T_C_Count = new int[256]();int* T_sum = new int[256]();
        
        int* G_A_Count = new int[256]();int* G_T_Count = new int[256]();int* G_Count = new int[256]();int* G_C_Count = new int[256]();int* G_sum = new int[256]();

        int* C_A_Count = new int[256]();int* C_T_Count = new int[256]();int* C_G_Count = new int[256]();int* C_Count = new int[256]();int* C_sum = new int[256]();

        gz = gzopen(MapDMG_input,"r");
        assert(gz!=Z_NULL);
        int i = 0;
        while(gzgets(gz,buf,length)){
            if (buf[0] != '#' && buf[0] != 'C'){//disregard the first couple of lines staring with '#' (description) or 'C' (column names)
                //fprintf(stderr,"buffer 1 %s\n",buf);
                //chr1	3p	+	1	11689	8083	7781	11648	39201	1767	16	35	8	4	12	291	31	117	3	7	0	0	0	0	0	0	0	0	2077
                //chr1	5p	+	1	11854	7524	8101	11722	39201	13	2720	9	10	4	14	6	7	3	9	6	1495
                char* Frag_end; char* tok_count;
                strtok(buf,"\t"); // chr1
                Frag_end = strtok(NULL,"\t"); // 3p
                if(strcasecmp("5p",Frag_end)==0){
                    char* tok1;
                    tok1 = strtok(NULL,"\t"); // +
                    int token_count = 0;
                    int pos = 0;
                    while(((tok1 = strtok(NULL,"\t")))){
                        //Chr	End	Std	Pos	A	C	G	T	Total	G>A	C>T	A>G	T>C	A>C	A>T	C>G	C>A	T>G	T>A	G>C	G>T	A>-	T>-	C>-	G>-	->A	->T	->C	->G	S
                        if(token_count == 0){pos = atoi(tok1);}
                        else if(token_count == 1){A_Count[pos] += atoi(tok1);A_sum[pos] += atoi(tok1);}// nt A
                        else if(token_count == 2){C_Count[pos] += atoi(tok1);C_sum[pos] += atoi(tok1);}// nt C
                        else if(token_count == 3){G_Count[pos] += atoi(tok1);G_sum[pos] += atoi(tok1);}// nt G
                        else if(token_count == 4){T_Count[pos] += atoi(tok1);T_sum[pos] += atoi(tok1);}// nt T

                        else if(token_count == 6){G_A_Count[pos] += atoi(tok1);G_sum[pos] += atoi(tok1);} //G>A
                        else if(token_count == 7){C_T_Count[pos] += atoi(tok1);C_sum[pos] += atoi(tok1);} //C>T
                        else if(token_count == 8){A_G_Count[pos] += atoi(tok1);A_sum[pos] += atoi(tok1);} //A>G
                        else if(token_count == 9){T_G_Count[pos] += atoi(tok1);T_sum[pos] += atoi(tok1);} //A>G
                        
                        else if(token_count == 10){A_C_Count[pos] += atoi(tok1);A_sum[pos] += atoi(tok1);} //A>C
                        else if(token_count == 11){A_T_Count[pos] += atoi(tok1);A_sum[pos] += atoi(tok1);} //A>T
                        
                        else if(token_count == 12){C_G_Count[pos] += atoi(tok1);C_sum[pos] += atoi(tok1);} //C>G
                        else if(token_count == 13){C_A_Count[pos] += atoi(tok1);C_sum[pos] += atoi(tok1);} //C>A
                        
                        else if(token_count == 14){T_G_Count[pos] += atoi(tok1);T_sum[pos] += atoi(tok1);} //C>G
                        else if(token_count == 15){T_A_Count[pos] += atoi(tok1);T_sum[pos] += atoi(tok1);} //C>A

                        else if(token_count == 16){G_C_Count[pos] += atoi(tok1);G_sum[pos] += atoi(tok1);} //C>G
                        else if(token_count == 17){G_T_Count[pos] += atoi(tok1);G_sum[pos] += atoi(tok1);} //C>A
                        token_count++;
                    }
                }
                else if(strcasecmp("3p",Frag_end)==0){
                    char* tok1;
                    tok1 = strtok(NULL,"\t"); // +
                    int token_count = 0;
                    int pos = 0;
                    while(((tok1 = strtok(NULL,"\t")))){
                        //Chr	End	Std	Pos	A	C	G	T	Total	G>A	C>T	A>G	T>C	A>C	A>T	C>G	C>A	T>G	T>A	G>C	G>T	A>-	T>-	C>-	G>-	->A	->T	->C	->G	S
                        if(token_count == 0){pos = (pos_val*2)-atoi(tok1)+1;}

                        else if(token_count == 1){A_Count[pos] += atoi(tok1);A_sum[pos] += atoi(tok1);}// nt A
                        else if(token_count == 2){C_Count[pos] += atoi(tok1);C_sum[pos] += atoi(tok1);}// nt C
                        else if(token_count == 3){G_Count[pos] += atoi(tok1);G_sum[pos] += atoi(tok1);}// nt G
                        else if(token_count == 4){T_Count[pos] += atoi(tok1);T_sum[pos] += atoi(tok1);}// nt T

                        else if(token_count == 6){G_A_Count[pos] += atoi(tok1);G_sum[pos] += atoi(tok1);} //G>A
                        else if(token_count == 7){C_T_Count[pos] += atoi(tok1);C_sum[pos] += atoi(tok1);} //C>T
                        else if(token_count == 8){A_G_Count[pos] += atoi(tok1);A_sum[pos] += atoi(tok1);} //A>G
                        else if(token_count == 9){T_G_Count[pos] += atoi(tok1);T_sum[pos] += atoi(tok1);} //A>G
                        
                        else if(token_count == 10){A_C_Count[pos] += atoi(tok1);A_sum[pos] += atoi(tok1);} //A>C
                        else if(token_count == 11){A_T_Count[pos] += atoi(tok1);A_sum[pos] += atoi(tok1);} //A>T
                        
                        else if(token_count == 12){C_G_Count[pos] += atoi(tok1);C_sum[pos] += atoi(tok1);} //C>G
                        else if(token_count == 13){C_A_Count[pos] += atoi(tok1);C_sum[pos] += atoi(tok1);} //C>A
                        
                        else if(token_count == 14){T_G_Count[pos] += atoi(tok1);T_sum[pos] += atoi(tok1);} //C>G
                        else if(token_count == 15){T_A_Count[pos] += atoi(tok1);T_sum[pos] += atoi(tok1);} //C>A

                        else if(token_count == 16){G_C_Count[pos] += atoi(tok1);G_sum[pos] += atoi(tok1);} //C>G
                        else if(token_count == 17){G_T_Count[pos] += atoi(tok1);G_sum[pos] += atoi(tok1);} //C>A
                        token_count++;
                    }
                }
            }
        }
        gzclose(gz);

        FILE *MisIncorpFile = fopen(NGSNGS_output, "w");
        for (int elem=0; elem<256; ++elem){
            if (A_Count[elem] != 0){
                double A_freq = (double) A_Count[elem]/A_sum[elem]; double AT_freq = A_freq+(double)A_T_Count[elem]/A_sum[elem];
                double AG_freq = AT_freq+(double) A_G_Count[elem]/A_sum[elem]; double AC_freq = AG_freq+(double)A_C_Count[elem]/A_sum[elem];
                fprintf(MisIncorpFile,"%f\t%f\t%f\t%f\n",A_freq,AT_freq,AG_freq,AC_freq);
            }
        }
        for (int elem=0; elem<256; ++elem){
            if (T_Count[elem] != 0){
                double TA_freq = (double) T_A_Count[elem]/T_sum[elem]; double T_freq = TA_freq+(double) T_Count[elem]/T_sum[elem];
                double TG_freq = T_freq+(double) T_G_Count[elem]/T_sum[elem];double TC_freq = TG_freq+(double) T_C_Count[elem]/T_sum[elem];
                fprintf(MisIncorpFile,"%f\t%f\t%f\t%f\n",TA_freq,T_freq,TG_freq,TC_freq);
            }
        }
        for (int elem=0; elem<256; ++elem){
            if (G_Count[elem] != 0){
                double GA_freq = (double) G_A_Count[elem]/G_sum[elem]; double GT_freq = GA_freq+(double) G_T_Count[elem]/G_sum[elem];
                double G_freq = GT_freq+(double) G_Count[elem]/G_sum[elem];double GC_freq = G_freq+(double) G_C_Count[elem]/G_sum[elem];
                fprintf(MisIncorpFile,"%f\t%f\t%f\t%f\n",GA_freq,GT_freq,G_freq,GC_freq);
            }
        }
        for (int elem=0; elem<256; ++elem){
            if (C_Count[elem] != 0){
                double CA_freq = (double) C_A_Count[elem]/C_sum[elem]; double CT_freq = CA_freq+(double) C_T_Count[elem]/C_sum[elem];
                double CG_freq = CT_freq+(double) C_G_Count[elem]/C_sum[elem];double C_freq = CG_freq+(double) C_Count[elem]/C_sum[elem];
                fprintf(MisIncorpFile,"%f\t%f\t%f\t%f\n",CA_freq,CT_freq,CG_freq,C_freq);
            }
        }
        fclose(MisIncorpFile);
    }  
    return 0;
}

#ifdef __WITH_MAIN__
int main(int argc,char **argv){
    main_Mis(argc,argv); 
    return 0;
}
#endif
//g++ MapDMGmis.cpp -std=c++11 -lm -lz -D __WITH_MAIN__ -o MisConvert