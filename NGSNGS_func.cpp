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

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <zlib.h>
#include <errno.h>

#include <pthread.h>

#include "NGSNGS_func.h"

#define LENS 4096
#define MAXBINS 100

void DNA_complement(char seq[]){
  while (*seq) {
    switch(*seq) {
      case 'A':
      case 'a':
        *seq = 'T';
        break;
      case 'G':
      case 'g':
        *seq = 'C';
        break;
      case 'C':
      case 'c':
        *seq = 'G';
        break;
      case 'T':
      case 't':
        *seq = 'A';
        break;
      case 'N':
      case 'n':
        *seq = 'N';
        break;  
    }
    ++seq;
  }
}

void reverseChar(char* str) {
    std::reverse(str, str + strlen(str));
}

int Random_geometric_k(unsigned int  seed, const double p)
{
  double u = ((double) rand_r(&seed)/ RAND_MAX);
  double k;

  if (p == 1){k = 1;}
  else if(p == 0){k=0;}
  else{k = log (u) / log (1 - p);}

  return floor(k);
}

double myrand(unsigned int persistent){
  return ((double) rand_r(&persistent)/ RAND_MAX);
}

void SimBriggsModel(char* reffrag, char* frag, int L, double nv, double lambda, double delta_s, double delta, unsigned int seed){
    int l = 0;
    int r = L-1;
    while (l+r > L-2){
        l = 0;
        r = 0;
        double u_l = ((double) rand_r(&seed)/ RAND_MAX);//uniform();
        double u_r = ((double) rand_r(&seed)/ RAND_MAX);//uniform();
        
        if (u_l > 0.5){
            l = Random_geometric_k(seed+0,0.36);//distribution1(generator1);
            //fprintf(stderr,"U_L LEI %d\n",l);
            //cout << l << "\n";
        }
        if (u_r > 0.5){
          //fprintf(stderr,"U_R");
          r = Random_geometric_k(seed+1,0.36); //distribution1(generator2);
          //fprintf(stderr,"U_R LEI %d\n",r);
          //fprintf(stderr,"U_R %lf", r);
          //cout << r << "\n";
        }
    }
    //cout << l << " " << r <<"\n";
    for (int i = 0; i<l; i++){
        if (reffrag[i] == 'C' || reffrag[i] == 'c' ){
            double u = ((double) rand_r(&seed)/ RAND_MAX);//uniform();
            if (u < delta_s){
                frag[i] = 'T';
            }else{
                frag[i] = 'C';
            }
        }else{
            frag[i] = reffrag[i];
        }
    }
    for (int i = 0; i < r; i++){
        if (reffrag[L-i-1] == 'G' || reffrag[L-i-1] == 'g'){
            double u = ((double) rand_r(&seed)/ RAND_MAX);//uniform();
            if (u < delta_s){
                frag[L-i-1] = 'A';
            }else{
                frag[L-i-1] = 'G';
            }
        }else{
            frag[L-i-1] = reffrag[L-i-1];
        }
    }
    double u_nick = ((double) rand_r(&seed)/ RAND_MAX);//uniform();
    double d = nv/((L-l-r-1)*nv+1-nv);
    int p_nick = l;
    double cumd = 0;
    while (u_nick > cumd & p_nick < L-r-1){
        cumd += d;
        p_nick +=1;
    }
    for (int i = l; i < L-r; i++){
        if ((reffrag[i] == 'C' || reffrag[i] == 'c') && i<=p_nick){
            double u = ((double) rand_r(&seed)/ RAND_MAX);//uniform();
            if (u < delta){
                frag[i] = 'T';
            }else{
                frag[i] = 'C';
            }
        }else if ((reffrag[i] == 'G' || reffrag[i] == 'g') && i>p_nick){
            double u = ((double) rand_r(&seed)/ RAND_MAX);//uniform();
            if (u < delta){
                frag[i] = 'A';
            }else{
                frag[i] = 'G';
            }
        }else{
            frag[i] = reffrag[i];
        }
    }
}

const char* Error_lookup(double a,double err[6000],int nt_offset, int read_pos,int outputoffset){

  //const char* nt_qual[8] = {"#", "\'", "0", "7" ,"<", "B", "F","I"}; // 35,39,48,55,66,70,73 //then withouth & in nt_qual[0]
  //const char* nt_qual[8] = {"#", "\'", "0", "7" ,"<", "B", "F","I"};//{"#", "\'", "0", "7" ,"<", "B", "F","I"};

  char nt_qual[8] = {(char) (2+outputoffset),
                     (char) (6+outputoffset),
                     (char) (15+outputoffset),
                     (char) (22+outputoffset),
                     (char) (27+outputoffset),
                     (char) (33+outputoffset),
                     (char) (37+outputoffset),
                     (char) (40+outputoffset)};
  //const char* nt_qual[8] = (const char*) nt_qual_int[9];
  int offset = ((nt_offset+read_pos)*8);
  //printf("offset %d \n", offset);
  const char* nt_out;
  if (a <= err[offset]){
    nt_out = &nt_qual[0];
  }
  else if (err[offset] < a && a <= err[offset+1]){
    nt_out = &nt_qual[1];
    }
  else if (err[offset+1] < a && a <= err[offset+2]){
    nt_out = &nt_qual[2];
  }
  else if (err[offset+2] < a && a <= err[offset+3]){
    nt_out = &nt_qual[3];
  }
  else if (err[offset+3] < a && a<= err[offset+4]){
    nt_out = &nt_qual[4];
  }
  else if (err[offset+4] < a && a<= err[offset+5]){
    nt_out = &nt_qual[5];
  }
  else if (err[offset+5] < a && a<= err[offset+6]){
    nt_out = &nt_qual[6];
  }
  else if (err[offset+6] < a && a <= err[offset+7]){
    nt_out = &nt_qual[7];
  }
  return nt_out;
}

double* Qual_array(double* freqval,const char* filename){
  int length = 6000;
  char buf[length];
  gzFile gz = Z_NULL;
  gz = gzopen(filename,"r");
  assert(gz!=Z_NULL);
  int i = 0;
  while(gzgets(gz,buf,length)){
    double val1;double val2;double val3;double val4;double val5;double val6;double val7;double val8;
    sscanf(buf,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",&val1,&val2,&val3,&val4,&val5,&val6,&val7,&val8);
    //std::cout << "iter " << i << std::endl;// " " << buf << std::endl;
    freqval[i*8] = val1; freqval[i*8+1] = val2; freqval[i*8+2] = val3; freqval[i*8+3] = val4;
    freqval[i*8+4] = val5; freqval[i*8+5] = val6; freqval[i*8+6] = val7; freqval[i*8+7] = val8;
    i++;
  }
  gzclose(gz);
  return freqval;
}

void Read_Qual_new(char *seq,char *qual,unsigned int seed,double* freqval,int outputoffset){
  //memset(qual, '\0', sizeof(qual));
  //fprintf(stderr,"INSIDE FUNCITON NOW \n");
  int Tstart = 150;
  int Gstart = 300;
  int Cstart = 450;
  int Nstart = 600;

  //srand(time( NULL ));
  //srand(seed);

  //srand48(seed);
  int seqlen = strlen(seq);

  for (int row_idx = 0; row_idx < seqlen; row_idx++){
    //fprintf(stderr,"index %d \n",row_idx);
    double r = ((double) rand_r(&seed)/ RAND_MAX);
    //double r = 0.43; // drand48();//0.43;//(double) rand()/RAND_MAX;
    //fprintf(stderr,"random value %lf \n",r);
    //fprintf(stderr,"Seq idx %c \n",seq[row_idx]);
    switch(seq[row_idx]){
      //fprintf(stderr,"qual %s \n",qual);
      case 'A':
      case 'a':
        //fprintf(stderr,"err %s \n",Error_lookup(r,freqval,0,row_idx,outputoffset));
        qual[row_idx] = 'F'; //Error_lookup(r,freqval,0,row_idx,outputoffset)[0];
        //qual[row_idx+1] = '\0';
        //strncat(qual, Error_lookup(r,freqval,0,row_idx,outputoffset), 1);
        //puts(qual);
        break;
      case 'T':
      case 't':
        qual[row_idx] =  'I'; // Error_lookup(r,freqval,Tstart,row_idx,outputoffset)[0];
        //qual[row_idx+1] = '\0';
        //strncat(qual, Error_lookup(r,freqval,Tstart,row_idx,outputoffset), 1);
        break;  
      case 'G':
      case 'g':
        qual[row_idx] = '!'; // Error_lookup(r,freqval,Gstart,row_idx,outputoffset)[0];
        //qual[row_idx+1] = '\0';
        // fprintf(stderr,"string error %c \n",Error_lookup(r,freqval,Gstart,row_idx,outputoffset)[0]);
        //strncat(qual, Error_lookup(r,freqval,Gstart,row_idx,outputoffset), 1);
        break;
      case 'C':
      case 'c':
        qual[row_idx] = 'B'; // Error_lookup(r,freqval,Cstart,row_idx,outputoffset)[0];
        //qual[row_idx+1] = '\0';
        //strncat(qual, Error_lookup(r,freqval,Cstart,row_idx,outputoffset), 1);
        break;
      case 'N':
      case 'n':
        qual[row_idx] = '<'; // Error_lookup(r,freqval,Nstart,row_idx,outputoffset)[0];
        //qual[row_idx+1] = '\0';
        //strncat(qual, Error_lookup(r,freqval,Nstart,row_idx,outputoffset), 1);
        break;
    }
  }
  qual[seqlen] = '\0';
  //memset(qual, 0, sizeof(qual)); //remove this to create the full nucleotide string but also the 4 mistakes?
  //fprintf(stderr,"EXITING FUNCITON NOW \n");
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
  //int LENS = 4096;
  //int* Frag_len = new int[LENS];
  //double* Frag_freq = new double[LENS];
  int n =0;

  gzFile gz = Z_NULL;
  char buf[LENS];
  
  gz = gzopen(filename,"r");
  assert(gz!=Z_NULL);
  while(gzgets(gz,buf,LENS)){
    Length[n] = atoi(strtok(buf,"\n\t ")); //before it was Frag_len[n]
    Frequency[n] = atof(strtok(NULL,"\n\t "));
    n++;
  }
  gzclose(gz); 

  number = n;
  //Length = Frag_len;
  //Frequency = Frag_freq;
  //delete[] Frag_len;
  //delete[] Frag_freq;
}

void printTime(FILE *fp){
  time_t rawtime;
  struct tm * timeinfo; 
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  fprintf (fp, "\t-> %s", asctime (timeinfo) );
}

void Header_func(htsFormat *fmt_hts,const char *outfile_nam,samFile *outfile,sam_hdr_t *header,faidx_t *seq_ref,int chr_total, size_t genome_len){
  // Creates a header for the bamfile. The header is initialized before the function is called //

  if (header == NULL) { fprintf(stderr, "sam_hdr_init");}
    
  // Creating header information
  
  char genome_len_buf[1024];
  for(int i=0;i<chr_total;i++){
    const char *name = faidx_iseq(seq_ref,i);

    int name_len =  faidx_seq_len(seq_ref,name);
    snprintf(genome_len_buf,1024,"%d", name_len);
    
    // reference part of the header, int r variable ensures the header is added
    int r = sam_hdr_add_line(header, "SQ", "SN", name, "LN", genome_len_buf, NULL);
    if (r < 0) { fprintf(stderr,"sam_hdr_add_line");}
    memset(genome_len_buf,0, sizeof(genome_len_buf));
  }
  // saving the header to the file
  if (sam_hdr_write(outfile, header) < 0) fprintf(stderr,"writing headers to %s", outfile_nam); //outfile
}

char* full_genome_create(faidx_t *seq_ref,int chr_total,int chr_sizes[],const char *chr_names[],size_t chr_size_cumm[]){
  
  size_t genome_size = 0;
  chr_size_cumm[0] = 0;
  for (int i = 0; i < chr_total; i++){
    const char *chr_name = faidx_iseq(seq_ref,i);
    int chr_len = faidx_seq_len(seq_ref,chr_name);
    chr_sizes[i] = chr_len;
    chr_names[i] = chr_name;
    genome_size += chr_len;
    chr_size_cumm[i+1] = genome_size;
  }

  char* genome = (char*) malloc(sizeof(char) * (genome_size+chr_total));
  genome[0] = 0; //Init to create proper C string before strcat
  //chr_total
  for (int i = 0; i < chr_total; i++){

    const char *data = fai_fetch(seq_ref,chr_names[i],&chr_sizes[i]);
    //sprintf(&genome[strlen(genome)],data);
    //strcat(genome,data);  //Both gives conditional jump or move error
    if (data != NULL){
      sprintf(genome+strlen(genome),data); 
    }
    // several of the build in functions allocates memory without freeing it again.
    free((char*)data); //Free works on const pointers, so we have to cast into a const char pointer
  }
  return genome;
}


//! Allocate workspace for random-number sampling.
ransampl_ws* ransampl_alloc( int n )
{
    ransampl_ws *ws;
    if ( !(ws = (ransampl_ws*) malloc( sizeof(ransampl_ws) )) ||
         !(ws->alias = (int*) malloc( n*sizeof(int) )) ||
         !(ws->prob = (double*) malloc( n*sizeof(double) )) ) {
        fprintf( stderr, "ransampl: workspace allocation failed\n" );
        exit(ENOMEM); // ENOMEM
    }
    ws->n = n;
    return ws;
}

//! Initialize workspace by precompute alias tables from given probabilities.
void ransampl_set( ransampl_ws *ws, const double* p )
{
    const int n = ws->n;
    int i, a, g;

    // Local workspace:
    double *P;
    int *S, *L;
    if ( !(P = (double*) malloc( n*sizeof(double) ) ) ||
         !(S = (int*) malloc( n*sizeof(int) ) ) ||
         !(L = (int*) malloc( n*sizeof(int) ) ) ) {
        fprintf( stderr, "ransampl: temporary allocation failed\n" );
        exit(ENOMEM);
    }

    // Normalise given probabilities:
    double sum=0;
    for ( i=0; i<n; ++i ) {
        if( p[i]<0 ) {
            fprintf( stderr, "ransampl: invalid probability p[%i]<0\n", i );
            exit(EINVAL);
        }
        sum += p[i];
    }
    if ( !sum ) {
        fprintf( stderr, "ransampl: no nonzero probability\n" );
        exit(EINVAL);
    }
    for ( i=0; i<n; ++i )
        P[i] = p[i] * n / sum;

    // Set separate index lists for small and large probabilities:
    int nS = 0, nL = 0;
    for ( i=n-1; i>=0; --i ) {
        // at variance from Schwarz, we revert the index order
        if ( P[i]<1 )
            S[nS++] = i;
        else
            L[nL++] = i;
    }

    // Work through index lists
    while ( nS && nL ) {
        a = S[--nS]; // Schwarz's l
        g = L[--nL]; // Schwarz's g
        ws->prob[a] = P[a];
        ws->alias[a] = g;
        P[g] = P[g] + P[a] - 1;
        if ( P[g] < 1 )
            S[nS++] = g;
        else
            L[nL++] = g;
    }

    while ( nL )
        ws->prob[ L[--nL] ] = 1;

    while ( nS )
        // can only happen through numeric instability
        ws->prob[ S[--nS] ] = 1;

    // Cleanup:
    free( P );
    free( S );
    free( L );
}

int ransampl_draw2( ransampl_ws *ws,double r1, double r2)
{   
    //fprintf(stderr,"%lf\t%lf\n",r1,r2);
    const int i = (int) (ws->n * r1);
    return r2 < ws->prob[i] ? i : ws->alias[i];
}

//! Free the random-number sampling workspace.
void ransampl_free( ransampl_ws *ws )
{
    free( ws->alias );
    free( ws->prob );
    free( ws );
}

ransampl_ws ***ReadQuality(char *ntqual, int ntcharoffset,const char *freqfile,unsigned long &readcycle){
  ransampl_ws ***dists = new ransampl_ws**[5];

  std::vector<char *> all_lines;
  gzFile gz = Z_NULL;
  assert(((gz = gzopen(freqfile,"rb")))!=Z_NULL);
  char buf[LENS];
  while(gzgets(gz,buf,LENS))
    all_lines.push_back(strdup(buf));
  gzclose(gz);
  
  fprintf(stderr,"All lines: %lu\n",all_lines.size());
  unsigned long readcyclelength = (all_lines.size()-1)/5; //all_lines.size()-1

  fprintf(stderr,"Inferred read cycle lengths: %lu\n",readcyclelength);
  readcycle = readcyclelength; 
  //loop over inputdata
  int nbins = -1;
  double probs[MAXBINS];
  for(int b=0;b<5;b++){
    dists[b] = new ransampl_ws *[readcyclelength];
    for(int pos = 0 ; pos<readcyclelength;pos++){
      int at = 0;
      probs[at++] = atof(strtok(all_lines[1+b*readcyclelength+pos],"\n\t ")); //1+b*readcyclelength+pos
      char *tok = NULL;
      while(((tok=strtok(NULL,"\n\t ")))){
	      probs[at++] = atof(tok);
	      assert(at<MAXBINS);
      }
      if(nbins==-1){
	      nbins = at;
	      fprintf(stderr,"Number of qualities/bins in inputfile: %d\n",nbins);
      }
      if(nbins!=at){
	      fprintf(stderr,"Problems, number of columns is different nbins: %d at: %d\n",nbins,at);
	      exit(0);
      }
      dists[b][pos] =  ransampl_alloc( nbins );
      ransampl_set(dists[b][pos],probs);
    }
  }

  //printf(all_lines[0]);
  int idx = 1;
  //char nt_qual[nbins]; //{(char) (2+outputoffset)}; char nt_qual[8]
  //nt_qual[idx] = //(char) atoi(strtok(all_lines[0],"\n\t ")); //1+b*readcyclelength+pos
  fprintf(stderr,"Number of qualities/bins in inputfile: %d\n",nbins);
 
  //extract the first token
  ntqual[0] = (char) (atoi(strtok(all_lines[0],"\n\t "))+ntcharoffset);
  char *tok = NULL;

  //extract the next
  while(((tok=strtok(NULL,"\n\t ")))){
	  ntqual[idx++] = (char) (atoi(tok)+ntcharoffset);
  }
  
  //strdup function allocate necessary memory to store the sourcing string implicitly, i need to free the returned string
  for (int i = 0; i < all_lines.size(); i++){free(all_lines[i]);}
  
  return dists;
}