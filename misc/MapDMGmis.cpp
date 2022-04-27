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


int main(int argc,char **argv){
    int length = 6000;
    char buf[length];

    gzFile gz = Z_NULL;

    /*
    # table produced by mapDamage version 2.1.0
    # using mapped file SampleResult/align/Sim_2_s_se.bam and /willerslev/users-shared/science-snm-willerslev-wql443/scratch/reference_files/Human/hg19canon.fa as reference file
    # Chr: reference from sam/bam header, End: from which termini of DNA sequences, Std: strand of reads
    Chr	End	Std	Pos	A	C	G	T	Total	G>A	C>T	A>G	T>C	A>C	A>T	C>G	C>A	T>G	T>A	G>C	G>T	A>-	T>-	C>-	G>-	->A	->T	->C	->G	S
    */
    gz = gzopen("misincorporation.txt","r");
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
    fprintf(stderr,"THE NUMBER OF POSITIONS ARE %d\n",pos_val);


    int* A_Count = new int[256]();int* A_T_Count = new int[256]();int* A_G_Count = new int[256]();int* A_C_Count = new int[256]();int* A_sum = new int[256]();
    
    int* T_A_Count = new int[256]();int* T_Count = new int[256]();int* T_G_Count = new int[256]();int* T_C_Count = new int[256]();int* T_sum = new int[256]();
    
    int* G_A_Count = new int[256]();int* G_T_Count = new int[256]();int* G_Count = new int[256]();int* G_C_Count = new int[256]();int* G_sum = new int[256]();

    int* C_A_Count = new int[256]();int* C_T_Count = new int[256]();int* C_G_Count = new int[256]();int* C_Count = new int[256]();int* C_sum = new int[256]();

    gz = gzopen("mistest.txt","r");
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
                    if(token_count == 0){pos = atoi(tok1);
                        //fprintf(stderr,"------------\nPOSITION IS %d \n",pos); //11722
                    }
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
                    //fprintf(stderr,"token count %d\n",token_count);
                }
                i++;
                //if(i>4){exit(0);}
            }
            else if(strcasecmp("3p",Frag_end)==0){
                char* tok1;
                tok1 = strtok(NULL,"\t"); // +
                int token_count = 0;
                int pos = 0;
                while(((tok1 = strtok(NULL,"\t")))){
                    //Chr	End	Std	Pos	A	C	G	T	Total	G>A	C>T	A>G	T>C	A>C	A>T	C>G	C>A	T>G	T>A	G>C	G>T	A>-	T>-	C>-	G>-	->A	->T	->C	->G	S
                    if(token_count == 0){
                        pos = (pos_val*2)-atoi(tok1)+1;
                        //fprintf(stderr,"------------\nPOSITION IN 3p IS %d \t %d\n",pos,atoi(tok1)); //11722
                    }
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
                    //fprintf(stderr,"token count %d\n",token_count);
                }
                i++;
                //if(i>4){exit(0);}
            }
        }
    }
    gzclose(gz);

    FILE *MisIncorpFile = fopen("MisOut2.txt", "w");
    for (int elem=0; elem<256; ++elem){
        if (A_Count[elem] != 0){
            double A_freq = (double) A_Count[elem]/A_sum[elem]; double AT_freq = A_freq+(double)A_T_Count[elem]/A_sum[elem];
            double AG_freq = AT_freq+(double) A_G_Count[elem]/A_sum[elem]; double AC_freq = AG_freq+(double)A_C_Count[elem]/A_sum[elem];
            //fprintf(MisIncorpFile,"%d\t%d\t%d\t%d\t%d\n",A_Count[elem],A_T_Count[elem],A_G_Count[elem],A_C_Count[elem],A_sum[elem]);
            //fprintf(MisIncorpFile,"%f\t%f\t%f\t%f\t%d\n",A_freq,AT_freq,AG_freq,AC_freq,A_sum[elem]);
            fprintf(MisIncorpFile,"%f\t%f\t%f\t%f\n",A_freq,AT_freq,AG_freq,AC_freq);
        }
    }
    for (int elem=0; elem<256; ++elem){
        if (T_Count[elem] != 0){
            double TA_freq = (double) T_A_Count[elem]/T_sum[elem]; double T_freq = TA_freq+(double) T_Count[elem]/T_sum[elem];
            double TG_freq = T_freq+(double) T_G_Count[elem]/T_sum[elem];double TC_freq = TG_freq+(double) T_C_Count[elem]/T_sum[elem];
            //fprintf(MisIncorpFile,"%d\t%d\t%d\t%d\t%d\n",A_Count[elem],A_T_Count[elem],A_G_Count[elem],A_C_Count[elem],A_sum[elem]);
            //fprintf(MisIncorpFile,"%f\t%f\t%f\t%f\t%d\n",A_freq,AT_freq,AG_freq,AC_freq,A_sum[elem]);
            fprintf(MisIncorpFile,"%f\t%f\t%f\t%f\n",TA_freq,T_freq,TG_freq,TC_freq);
        }
    }
    fclose(MisIncorpFile);
    return 0;
}

//g++ MapDMGlen.cpp -std=c++11 -lm -lz
/*
Chr	    End	Std	Pos	 A   	 C	 G	 T	    Total	G>A	C>T	A>G	T>C	A>C	A>T	C>G	C>A	T>G	T>A	G>C	G>T	A>-	T>-	C>-	G>-	->A	->T	->C	->G	S
chrY	5p	+	1	1261	780	805	1259	4105	0	280	0	0	1	0	0	2	0	0	0	137
chrY	5p	+	2	1319	670	831	1285	4105	0	125	0	0	1	0	0	0	2	0	0	137
chrY	5p	+	3	1250	772	870	1213	4105	1	80	1	1	1	0	0	2	0	0	0	110
chrY	5p	+	4	1260	763	816	1266	4105	2	58	0	1	0	0	0	0	1	0	0	76
chrY	5p	+	5	1222	804	802	1277	4105	0	34	0	0	0	3	2	0	1	0	0	39
chrY	5p	+	6	1255	806	835	1209	4105	3	38	1	0	1	1	0	0	0	0	0	19
*/

//chrY	5p	+	1	1261	780	805	1259	4105	0	280	0	0	1	0	0	2	0	0	0	137
//chrY	5p	-	1	1318	753	783	1250	4104	1	275	2	0	3	1	0	0	0	0	0	144
