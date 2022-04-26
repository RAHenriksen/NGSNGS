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
  const char *ART_qual_nt_profile_in;
  const char *NGSNGS_qual_nt_profile_out;
}argStruct;

int HelpPage(FILE *fp){
  fprintf(fp,"Read Quality Profile converter - converting ART profiles to NGSNGS format\n\n");
  fprintf(fp,"Usage\n./QualConvert -i <ART platform profile> -o <NGSNGS quality profile>\n");
  fprintf(fp,"\nExample\n./QualConvert -i HiSeq2500L150R1filter.txt -o HiSeq2500L150R1_CDF.txt\n");
  fprintf(fp,"./QualConvert --Art_in HiSeq2500L150R1filter.txt --NGSNGS_out HiSeq2500L150R1_CDF.txt\n");
  fprintf(fp,"\nOptions: \n");
  fprintf(fp,"-h   | --help: \t\t\t Print help page.\n");
  fprintf(fp,"-v   | --version: \t\t Print help page.\n\n");
  fprintf(fp,"-i   | --Art_in: \t\t ART platform profile intput in .txt format or .txt.gz\n");
  fprintf(fp,"-o   | --NGSNGS_out: \t\t NGSNGS nucletide quality profile in .txt format\n");
  exit(1);
  return 0;
}

argStruct *getpars(int argc,char ** argv){
  argStruct *mypars = new argStruct;
  mypars->ART_qual_nt_profile_in = NULL;
  mypars->NGSNGS_qual_nt_profile_out = NULL;
  ++argv;
  while(*argv){
    //fprintf(stderr,"ARGV %s\n",*argv);
    if(strcasecmp("-i",*argv)==0 || strcasecmp("--Art_in",*argv)==0){
      mypars->ART_qual_nt_profile_in = strdup(*(++argv));
    }
    else if(strcasecmp("-o",*argv)==0 || strcasecmp("--NGSNGS_out",*argv)==0){
      mypars->NGSNGS_qual_nt_profile_out = strdup(*(++argv));
    }
    else{
      fprintf(stderr,"unrecognized input option %s, see NGSNGS help page\n\n",*(argv));
      exit(0);
    }
    ++argv;
  }
  return mypars;
}

int main_qual(int argc,char **argv){
    argStruct *mypars = NULL;
    if(argc==1||(argc==2&&(strcasecmp(argv[1],"--version")==0||strcasecmp(argv[1],"-v")==0||
                            strcasecmp(argv[1],"--help")==0||strcasecmp(argv[1],"-h")==0))){
        HelpPage(stderr);
        return 0;
    }
    else{
        mypars = getpars(argc,argv);
        const char* ART_input = mypars->ART_qual_nt_profile_in;
        const char* NGSNGS_output = mypars->NGSNGS_qual_nt_profile_out; //"chr22_out";

        int count_line = 0;
        int qual_line = 1;
        int* ObsQS = new int[256];
        int col_no = 0;
        int FullQual[256];
        int length = 6000;
        char buf[length];
        gzFile gz = Z_NULL;

        // IDENTIFY NUMBER OF QUALITIES
        gz = gzopen(ART_input,"r");
        assert(gz!=Z_NULL);
        while(gzgets(gz,buf,length)){
            qual_line++;
            if (qual_line % 2){continue;}
            if (buf[0] == '.'){continue;}
            if (buf[0] == 'N'){continue;}
            else{
                strtok(buf,"\t");
                strtok(NULL,"\t");
                char* tok;
                while(((tok = strtok(NULL,"\t")))){ObsQS[atoi(tok)] = atoi(tok);}
            }
        }
        gzclose(gz);

        //Create the read quality file
        FILE *ReadQual_2 = fopen(NGSNGS_output, "w");
        for (int elem=0; elem<256; ++elem){
            if (ObsQS[elem] != 0){
                FullQual[col_no]=ObsQS[elem]; // creating an integer array only with qualities
                col_no++; //count the total number of qualities
                // save the first line of read qualities
                fprintf(ReadQual_2,"%d \t",elem);
            }
        }
        fprintf(ReadQual_2,"\n");

        // save the second line of probabilities
        for (int elem=0; elem<256; ++elem){
            if (ObsQS[elem] != 0){
                double prob = pow(10.0,(-ObsQS[elem]/10.0));
                fprintf(ReadQual_2,"%f \t",prob);
            }
        }
        fprintf(ReadQual_2,"\n");

        // open the file with quality counts again to get a new buffer
        gz = gzopen(ART_input,"r");
        assert(gz!=Z_NULL);

        while(gzgets(gz,buf,length)){
            count_line++;
            if (buf[0] == '.'){continue;}
            else if (buf[0] == 'N'){
                if (count_line % 2){
                    for (int i = 0; i < col_no; i++){
                        fprintf(ReadQual_2,"%e \t",1.00);
                    }
                    fprintf(ReadQual_2,"\n");
                }
                strtok(buf,"\t");
                strtok(NULL,"\t");
                char* tok;
                while(((tok = strtok(NULL,"\t")))){continue;}
            }
            else{
                /* creates three char arrays to copy both lines per position,e.g.
                A	18	3	7	16	23	28	34	38	41
                A	18	1	5554	47086	61240	125231	319599	997106	3188279
                */
                char Qual_line[256];char Count_line1[256];char Count_line2[256];

                // creates one copy of the nucleotide quality scores before moving on to the quality counts (next line)
                if (count_line % 2){snprintf(Qual_line, sizeof Qual_line, "%s", buf);} 
                else if (qual_line % 2){
                    double Freq_array[col_no];
                    double CDF_array[col_no];
                    int Qual_no[col_no];
                    int linesum = 0;
                    int line_elem = 0;

                    // create two copies of the line w. nt qualities
                    snprintf(Count_line1, sizeof Count_line1, "%s", buf);
                    snprintf(Count_line2, sizeof Count_line2, "%s", buf);
                    //fprintf(stderr,"Buffer 1: %s \nBuffer 2: %s \nBuffer 3: %s",Qual_line,Count_line1,Count_line2);
                    
                    // create first tokens for the count of qualities
                    strtok(Count_line1,"\t");
                    char* tok1;
                    tok1 = strtok(NULL,"\t"); //first token is position
                    //extract all the next tokens, i.e. the nucleotide quality counts to create the total sum of counts
                    while(((tok1 = strtok(NULL,"\t")))){linesum += atoi(tok1);line_elem++;}

                    // create new tokens since the Count_line1 tokens are destroyed
                    strtok(Count_line2,"\t");
                    strtok(NULL,"\t");
                    char* tok2;
                    tok2 = strtok(NULL,"\t");
                    int i = 1;
                    // first element of nucleotide quality frequency for that position
                    Freq_array[0] = (double) atof(tok2)/linesum;
                    // first element of nucleotide quality CDF for that position
                    CDF_array[0] = Freq_array[0];
                    // Since its the CDF the first values should be identical in both.
                    //fprintf(stderr,"Nucleotide quality frequency : %e and CDF value %e \n",Freq_array[0],CDF_array[0]);

                    while(((tok2 = strtok(NULL,"\t")))){
                        Freq_array[i] = (double) atof(tok2)/linesum;
                        CDF_array[i] = Freq_array[i]+CDF_array[i-1];
                        //fprintf(stderr,"Nucleotide quality frequency : %e and CDF value %e for nucelotide index %d\n",Freq_array[i],CDF_array[i],i);
                        i++;
                    }

                    // creates final token of the nucleotide qualities for that position
                    strtok(Qual_line,"\t");
                    strtok(NULL,"\t");
                    char* tok3;
                    tok3 = strtok(NULL,"\t");
                    int j = 1;
                    Qual_no[0] = atoi(tok3);
                    //fprintf(stderr,"Nucleotide quality %d for column index %d\n",Qual_no[0],0);
                    while(((tok3 = strtok(NULL,"\t")))){
                        Qual_no[j] = atoi(tok3);
                        //fprintf(stderr,"Nucleotide quality %d for column index %d\n",Qual_no[j],j);
                        j++;
                    }

                    //checking if the qualities from each line actual match the qualities from the array with total qualities
                    for (int col_i = 0; col_i < col_no; col_i++){
                        // if the index of nucleotide qualities match with the total nuceltide quality array, then save the CDF of that index
                        if(FullQual[col_i] == Qual_no[col_i]){
                            if (CDF_array[col_i] < 1.000000e-100){
                                fprintf(ReadQual_2,"%e \t",1.00);
                            }
                            else
                            {
                                fprintf(ReadQual_2,"%e \t",CDF_array[col_i]);
                            }
                            
                            //fprintf(stderr,"Full quality element : %d line quality %d CDF value %e for column index %d\n",FullQual[col_i],Qual_no[col_i],CDF_array[col_i],col_i);
                        }
                        /* if the nucleotide qualities for the full and for that line doesn't match, then it means that given line doesn't contain a given nt
                        Full nucleotide for total profile:
                            3	7   8	16	23	28	34  38  41
                        
                        Nucleotide values for a line missing quality value 7
                        A   46  3	8	16	23	28	34	38	41
                        A	46	1	20233	50182	99175	185278	493476	1569387	3280175
                        */
                        if(FullQual[col_i] != Qual_no[col_i]){
                            // if the last column values are below 1.0e-100 after moving several columns into the
                            // quality profile then its an artefact of skipping columns.
                            if (CDF_array[col_i-1] < 1.000000e-100){
                                fprintf(ReadQual_2,"%e \t",1.00);
                            }
                            else if (CDF_array[0] == 1.000000e+00 )
                            {   
                                fprintf(ReadQual_2,"%e \t",1.00);
                            }
                            else{fprintf(ReadQual_2,"%e \t",CDF_array[col_i-1]);} // so the column index which are skipped received the same value as previous to keep the CDF
                        }
                    }
                    fprintf(ReadQual_2,"\n");
                    // create the next buffer for next line to continue the iterations
                    strtok(buf,"\t");
                    strtok(NULL,"\t");
                    char* tok;
                    while(((tok = strtok(NULL,"\t")))){continue;}
                }
            }
        }
        fclose(ReadQual_2);
        gzclose(gz);
        delete[] ObsQS;
    }
    return 0;
}

#ifdef __WITH_MAIN__
int main(int argc,char **argv){
    main_qual(argc,argv); 
    return 0;
}
#endif

//g++ ReadQualConverter.cpp -std=c++11 -lm -lz -I /home/wql443/scratch/htslib/ /home/wql443/scratch/htslib/libhts.a -lpthread -lhts -lgsl -lgslcblas -lm -o QualConvert
//g++ ReadQualConverter.cpp -std=c++11 -lm -lz -D__WITH_MAIN__ -o QualConvert