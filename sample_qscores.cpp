#include <zlib.h>
#include <cstring>
#include <cassert>
#include <htslib/kstring.h>

#include "RandSampling.h"
#include "sample_qscores.h"
#include "mrand.h"
#include "NGSNGS_misc.h"

#define LENS 10000
#define MAXBINS 100

extern int refToInt[256];
extern char intToRef[4];

double phred2Prob[256];

ransampl_ws ***ReadQuality(char *ntqual, double *ErrProb, int ntcharoffset,const char *freqfile,int &inferredreadcycle){
  /*
  ReadQuality - reads quality scores and error probabilities from the input frequency file and sets up alias sampling distributions for cycle specific sequencing error depending on nucleotide
  The function also initializes a lookup table for Phred scores to probability conversions.
 
  @param ntqual            A character array where nucleotide quality scores will be stored.
  @param ErrProb           An array where error probabilities corresponding to each read position will be stored.
  @param ntcharoffset      An integer offset to adjust quality scores.
  @param freqfile          Path to the input frequency file containing quality score data.
  @param inferredreadcycle Reference to an integer where the inferred read cycle count will be stored.
  
  @return                  A 3D array of pointers (`ransampl_ws ***`) representing the alias method 
                           workspaces set up for each of the 5 bins representing possible nucleotide substitution errors (A,G,C,T,N) given read positions.
  */

  // Initialize Phred score to probability conversion table
  for(int qscore =0 ;qscore<256;qscore++){
    double d = qscore;
    phred2Prob[qscore] = pow(10,((-d)/10));
  }

  // Allocate memory for the 3D array of distribution
  ransampl_ws ***dists = new ransampl_ws**[5];

  // Read lines from the file and store them in the vector

  std::vector<char *> all_lines;
  gzFile gz = Z_NULL;
  assert(((gz = gzopen(freqfile,"rb")))!=Z_NULL);
  char buf[LENS];
  while(gzgets(gz,buf,LENS))
    all_lines.push_back(strdup(buf));
  gzclose(gz);

  inferredreadcycle = (all_lines.size()-2)/5; // Infer the number of sequencing cycles from the input data

  // Loop over each of the 5 bins
  int nbins = -1;
  double probs[MAXBINS];

  for(int b=0;b<5;b++){
    dists[b] = new ransampl_ws *[inferredreadcycle];

    // Loop over each read position
    for(int pos = 0 ; pos<inferredreadcycle;pos++){
      int at = 0;
      probs[at++] = atof(strtok(all_lines[2+b*inferredreadcycle+pos],"\n\t "));
      char *tok = NULL;

      // Parse probabilities from the line
      while(((tok=strtok(NULL,"\n\t ")))){
	      probs[at++] = atof(tok);
	      assert(at<MAXBINS);
      }

      // Validate that all bins have the same number of columns (the substitutions)
      if(nbins==-1){
	      nbins = at;
      }
      if(nbins!=at){
	      fprintf(stderr,"Problems, number of columns is different nbins: %d at: %d\n",nbins,at);
	      exit(1);
      }

      // Allocate and set up the alias method for the current position
      dists[b][pos] =  ransampl_alloc( nbins );
      ransampl_set(dists[b][pos],probs);
    }
  }

  // Parse the nucleotide quality scores from the first line
  int qualidx = 1;
  ntqual[0] = (char) (atoi(strtok(all_lines[0],"\n\t "))+ntcharoffset);
  char *qualtok = NULL;
  while(((qualtok=strtok(NULL,"\n\t ")))){
    ntqual[qualidx++] = (char) (atoi(qualtok)+ntcharoffset);
  }
  int Err_idx = 1;

  // Parse the error probabilities from the second line
  ErrProb[0] = atof(strtok(all_lines[1],"\n\t"));
  char *Errtok = NULL;
  while ((Errtok=strtok (NULL,"\n\t"))){
    ErrProb[Err_idx++] = (double) atof(Errtok);
  }

  // Free the memory allocated by strdup for each line
  for (unsigned long i = 0; i < all_lines.size(); i++){
    free(all_lines[i]);
  }
  
  return dists;
}

int sample_qscores(char *bases, char *qscores,int len,ransampl_ws ***ws,char *NtQuals,mrand_t *mr,int simError, int ntcharoffset){
  /*
  sample_qscores -  Samples quality scores for a given sequence of bases using a precomputed alias structure and generate sequencing error (by replacing nucleotide) based on quality scores

  @param bases         Character array representing the sequence of bases within original DNA molecule.
  @param qscores       Character array where sampled quality scores will be stored.
  @param len           Length of the sequence (`bases` and `qscores` arrays).
  @param ws            3D array of alias sampling workspaces for each base and position.
  @param NtQuals       Array mapping index to nucleotide quality scores.
  @param mr            Pointer to the random number generator state (type `mrand_t`).
  @param simError      Boolean flag indicating whether sequencing errors should be simulated.
  @param ntcharoffset  Offset used to adjust quality scores to Phred scale.
  @return              Returns 1 if an error was simulated in the sequence, otherwise 0.
 */

  int seq_err = 0;

  // Iterate over each base in the input sequence
  for(int i = 0;i<len;i++){
    // Generate two random numbers for sampling
    double dtemp1 = mrand_pop(mr);
    double dtemp2 = mrand_pop(mr);

    // Convert the input base to its integer representation
    char inbase = refToInt[bases[i]];

    // Draw a quality score index for the current base at position `i`
    int qscore_idx =  ransampl_draw2(ws[inbase][i],dtemp1,dtemp2);

    // Assign the corresponding quality score from NtQuals
    qscores[i] = NtQuals[qscore_idx];

    // If simulating sequencing errors, determine if an error should be introduced
    if (simError){
      double tmprand = mrand_pop(mr);
      // Check if the generated random number is less than the probability of error
      if ( tmprand < phred2Prob[qscores[i]-ntcharoffset]){
        // Calculate a random output base, ensuring it differs from the input base and set the error flag
        int outbase=(int)floor(4.0*phred2Prob[qscores[i]-ntcharoffset]*tmprand);
        while (((outbase=((int)floor(4*mrand_pop(mr))))) == inbase);
	      bases[i] = intToRef[outbase];
        seq_err = 1;
      }
    }
  }

  return seq_err;
}

int sample_qscores_amplicon(kstring_t* seq,kstring_t* qual,ransampl_ws ***ws,char *NtQuals,mrand_t *mr,int simError, int ntcharoffset){
  /*
  sample_qscores_amplicon -  Samples quality scores for an amplicon sequence and optionally simulates sequencing errors

  @param seq           Pointer to the kstring_t sequence representing the empirical sequencing reads.
  @param qual          Pointer to the kstring_t sequence representing the empirical quality strings reads.
  @param ws            3D array of alias sampling workspaces for each base and position.
  @param NtQuals       Array mapping index to nucleotide quality scores.
  @param mr            Pointer to the random number generator state (type `mrand_t`).
  @param simError      Boolean flag indicating whether sequencing errors should be simulated.
  @param ntcharoffset  Offset used to adjust quality scores to Phred scale.
  @return              Returns 1 if an error was simulated in the sequence, otherwise 0.
  */

  int seq_err = 0;

  // Iterate over each base in the empirical sequence read
  for(int i = 0;i<seq->l;i++){
    // Generate two random numbers for sampling
    double dtemp1 = mrand_pop(mr);
    double dtemp2 = mrand_pop(mr);

    // Convert the input base to its integer representation
    char inbase = refToInt[seq->s[i]];
    
    // Draw a quality score index for the current base at position `i`
    int qscore_idx =  ransampl_draw2(ws[refToInt[seq->s[i]]][i],dtemp1,dtemp2);

    // Assign the corresponding quality score from NtQuals
    qual->s[i] = NtQuals[qscore_idx];

    // If simulating sequencing errors, determine if an error should be introduced
    if (simError){
      double tmprand = mrand_pop(mr);
      // Check if the generated random number is less than the probability of error
      if ( tmprand < phred2Prob[qual->s[i]-ntcharoffset]){
        // Calculate a random output base, ensuring it differs from the input base and set the error flag
        int outbase=(int)floor(4.0*phred2Prob[qual->s[i]-ntcharoffset]*tmprand);//DRAGON
        while (((outbase=((int)floor(4*mrand_pop(mr))))) == inbase);
        seq->s[i] = intToRef[outbase];
        seq_err = 1;
      }
    }  
  }

  return seq_err;
}

int sample_qscores_fix(char *bases, char *qscores, int qscoresval,int len,mrand_t *mr,int simError, int ntcharoffset){
  /*
  sample_qscores_fix -  Samples quality scores for sequence and optionally simulates sequencing errors using a fixed quality score.

  @param bases         Pointer to the array of sequence bases.
  @param qscores       Pointer to the array where quality scores will be stored.
  @param qscoresval    The fixed integer quality score value to be applied to all bases.
  @param len           The length of the sequence.
  @param mr            Pointer to the random number generator state (type `mrand_t`).
  @param simError      Boolean flag indicating whether sequencing errors should be simulated.
  @param ntcharoffset  Offset used to adjust quality scores to Phred scale.
  @return              Returns 1 if an error was simulated in the sequence, otherwise 0.

  */
  int seq_err = 0;
      
  if (simError){
    // Precompute Phred scale probabilities for all possible quality scores
    for(int qscore =0 ;qscore<256;qscore++){
      double d = qscore;
      phred2Prob[qscore] = pow(10,((-d)/10));
    }
  }

  // Convert fized integer score into the corresponding character considering the offset (different from fq or bam)
  char fixedscore = (char)(qscoresval + ntcharoffset);
  
  for(int i = 0;i<len;i++){
    qscores[i] = fixedscore;
    char inbase = refToInt[bases[i]];

    if (simError){
      double tmprand = mrand_pop(mr);
      // Check if the generated random number is less than the probability of error
      if ( tmprand < phred2Prob[fixedscore-ntcharoffset]){
        // Calculate a random output base, ensuring it differs from the input base and set the error flag
        int outbase=(int)floor(4.0*phred2Prob[fixedscore-ntcharoffset]*tmprand);
        while (((outbase=((int)floor(4*mrand_pop(mr))))) == inbase);
	      bases[i] = intToRef[outbase];
        seq_err = 1;
      }
    }
  }
  return seq_err;
}

//void add_indel_amplicon_fa(mrand_t *mr,kstring_t* seq,double *pars,int* ops){
//add_indel_amplicon_fqbam(mr,&FQseq->seq,&FQseq->qual,pars,ops,ErrProbTypeOffset);

int sample_qscores_fix_amplicon(mrand_t *mr,kstring_t* seq,int qscoresval,int ntcharoffset){
  /*
  sample_qscores_fix_amplicon - Simulates sequencing errors for a given sequence based on a fixed quality score for empirical sequencing data (so with fasta input and fastq output).

  @param mr            Pointer to the random number generator state (type `mrand_t`).
  @param seq           Pointer to the sequence object (`kstring_t`) to be updated with simulated errors.
  @param qscoresval    The fixed integer quality score value used to determine error probabilities.
  @param ntcharoffset  Offset used to adjust quality scores to the Phred scale.
  @return              Returns 1 if an error was simulated in the sequence, otherwise 0.

  */
  int seq_err = 0;
  
  // temporary sequences
  kstring_t seq_intermediate;
  seq_intermediate.l = seq->l;
  seq_intermediate.m = seq->l;
  seq_intermediate.s = (char *)malloc((seq->l + 1) * sizeof(char)); // Allocate memory
  strcpy(seq_intermediate.s, seq->s); // Copy seq->s to seq_intermediate.s
  
  // Precompute Phred scale probabilities for all possible quality scores
  for(int qscore =0 ;qscore<256;qscore++){
    double d = qscore;
    phred2Prob[qscore] = pow(10,((-d)/10));
  }

  for(int i = 0;i<seq->l;i++){
    char inbase = refToInt[seq->s[i]];
    double tmprand = mrand_pop(mr);
    // Check if the generated random number is less than the probability of error
    if ( tmprand < phred2Prob[qscoresval-ntcharoffset]){
      // Calculate a random output base, ensuring it differs from the input base and set the error flag
      int outbase=(int)floor(4.0*phred2Prob[qscoresval]*tmprand);
      while (((outbase=((int)floor(4*mrand_pop(mr))))) == inbase);
	    seq_intermediate.s[i] = intToRef[outbase];
      seq->s[i] = seq_intermediate.s[i];
      seq_err = 1;
    }
  }

  seq->s[seq->l] = '\0';

  free(seq_intermediate.s);
  seq_intermediate.s = NULL;
  seq_intermediate.l = seq_intermediate.m = 0;

  return seq_err;
}

#ifdef __WITH_MAIN__
int main(int argc, char **argv){
  const char *profile_fname = "Test_Examples/AccFreqL150R1.txt";
  char ntquals[1024];
  double errorArray[1024];
  int maxreadcycles = 150;
  //ransampl_ws ***ws = ReadQuality(ntquals,errorArray,33,profile_fname);
  mrand_t *mr = mrand_alloc(3,88);
  char bases[30];
  char qscores[30];
  memset(bases,'\0',30);
  memset(qscores,'\0',30);

  char specific_char = 'I';
  char* char_array = (char*)malloc((30 + 1) * sizeof(char));
  for(int i=0;i<30;i++){
    bases[i] = intToRef[(int)floor(drand48()*4)];
    char_array[i] = specific_char;
    std::cout << " i " << i << " bases " << bases[i] << std::endl;
  }
  sample_qscores_fix(bases,qscores,10,30,mr,1,33);
  //sample_qscores(bases,qscores,30,ws,ntquals,mr,0,0);
  fprintf(stderr,"@readname\n%s\n+\n%s\n",bases,qscores);
  //gtggTAGAGATAAAGCACATTCTTTAGGAGTGAATATGGNNTNNCTGCNCGCANANTGNNATTGNNTTGCNNNTNNANCGNNNCNNTNNNGNTTNGCNACAGCNANGNNA
  return 0;
}
#endif

//g++ sample_qscores.cpp RandSampling.o mrand.o NGSNGS_misc.o -std=c++11 -lm -lz -D__WITH_MAIN__ -o Scores