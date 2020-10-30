#include <cstdio>
#include <cassert>
#include <htslib/faidx.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <vector>
#include <stdio.h>
#include <typeinfo>
#include <random>
#include <iterator>
#include <ctime>
#include <cmath>

// I would like to create a function with TK's code since its optimal in case we wish to 
//simulate a given number of fragments

struct Fastq_structure
{
    std::string ID;
    std::string sequence;
    std::string strand;
    std::string quality_value;
};

struct Fasta_structure
{
    std::string ID;
    std::string sequence;
};


void random_seq(faidx_t *seq_ref){
  // choose a random sequence -> still ned to change it so it saves the output to a single file.
  int readlength=35;
  int nreads = 8;
  for(int i=0;i<nreads;i++){
    char buf[96];//assume maxlength for readid is 96bytes
    int whichref = lrand48() % faidx_nseq(seq_ref);
    fprintf(stderr,"\t-> Whichref: %d\n",whichref);
    const char *name = faidx_iseq(seq_ref,whichref);
    int name_len =  faidx_seq_len(seq_ref,name);
    fprintf(stderr,"\t-> name: \'%s\' name_len: %d\n",name,name_len);

    int start = lrand48() % name_len;
    int stop = start+readlength;
    if(stop>name_len)
      stop = name_len;
    snprintf(buf,96,"%s:%d-%d",name,start,stop);
    fprintf(stderr,"buf: %s\n",buf);
    fprintf(stdout,"%s\n+\n",buf);
    char *data = fai_fetch(seq_ref,name,&name_len);

    for(int i=start;i<stop;i++)
      fprintf(stdout,"%c",data[i]);
    fprintf(stdout,"\n");

    for(int i=start;i<stop;i++)
      fprintf(stdout,"F");
    fprintf(stdout,"\n");
  }
}

void Deamin_2(std::string& Seq, std::string nt_find, int seed,
              double alpha=1.0,double beta=2.0,int start=0,int end=25){
    //declaring random seed to build random generator for gamma distribution
    //std::srand(Seq.size());
    std::default_random_engine generator(seed);
    std::gamma_distribution<double> distr(alpha,beta);

    // searching for T nt within the first 25 nt
    std::vector<int> Index_vec;
    std::string Sub = Seq.substr(start,end);

    int pos = Sub.find(nt_find);
    while(pos < Seq.size()) //std::string::npos
    {
        Index_vec.push_back(pos);
        pos =Sub.find(nt_find, pos + nt_find.size());
    }

    for (int i = 0; i < Index_vec.size(); i++){
        if (Index_vec.at(i) == int(distr(generator))) {
            Seq.replace(Index_vec[i],1,"U");

        }
		else {
            continue;
		}
    }

}

std::string Illumina_Qual(std::string Seq){
  int Q = 40;
  char ASCII = char(Q+33);
  double prob = 10^(-Q/10);
  return std::string(Seq.size(),ASCII);
}

int main(int argc,char **argv){

  const char *fastafile = "chr22.fa";
  //we use structure faidx_t from htslib to load in a fasta
  faidx_t *ref = NULL;
  ref  = fai_load(fastafile);
  assert(ref!=NULL);//check that we could load the file

  fprintf(stderr,"\t-> Number of contigs/scaffolds/chromosomes in file: \'%s\': %d\n",fastafile,faidx_nseq(ref));
  
  // choosing random sequences using -> random_seq(ref);

  // is lrand48() in order to pick a random sequence if containing more?
  int whichref = lrand48() % faidx_nseq(ref);
  std::cout << "reference number " << whichref << std::endl;
  const char *name = faidx_iseq(ref,whichref);
  int name_len =  faidx_seq_len(ref,name);
  std::cout << name << std::endl;
  std::cout << name_len << std::endl;
  
  int start_pos = 1;
  int end_pos = name_len; //30001000
  
  std::ofstream outfa("output.fa");
  std::ofstream outfq("output.fq");
  while(start_pos <= end_pos){
        // creates random number in the range of the fragment size rand() % ( high - low + 1 ) + low
        int rand_len = (std::rand() % (80 - 30 + 1)) + 30;
        // std::cout << "random number " << rand_len << std::endl;
        
        char *sequence = faidx_fetch_seq(ref,name,start_pos,start_pos+rand_len,&name_len);
          if (sequence[1] == 'N'){
            start_pos += rand_len;
          }
          else {
            std::string Damage = std::string(sequence); // Casting the character to string
            Deamin_2(Damage,"T",rand_len,1.0,2.0);
            if (Damage.size() < 150){
              char adapter = 'X';
              double No = (150-Damage.size())/2.0;

              outfa << ">" << name << ":" << start_pos << "-" << start_pos+name_len << "_Length:" << Damage.size() << std::endl;
              outfa << std::string(floor(No),adapter) << Damage << std::string(ceil(No),adapter) << std::endl;

              outfq << "@" << name << ":" << start_pos << "-" << start_pos+name_len << "_Length:" << Damage.size() << std::endl;
              outfq << std::string(floor(No),adapter) << Damage << std::string(ceil(No),adapter) << std::endl;
              outfq << "+" << std::endl;
              outfq << Illumina_Qual(Damage) << std::endl;
            }
            start_pos += rand_len;
          }
          //resets random seed each time with a new random number
          //std::srand(rand_len);
  }
}


  //std::cout << sequence << std::endl;
  /*
          if (sequence[1] == 'N'){
          //std::cout << sequence[1] << std::endl;
          //std::cout << sequence << std::endl;
          continue;
        }
        else{

  std::cout << sequence << std::endl;
  std::cout << "-----------" << std::endl;
  // Casting the character to string
  std::string deaminated = std::string(sequence);
  Deamin_2(deaminated,"T",1.0,2.0);
  std::cout << deaminated << std::endl;

  std::ofstream out("output.fa");
  out << ">" << name << ":" << start_pos << "-" << start_pos+name_len << std::endl;
  out << sequence;
  return 0;
  */

        /*
        if (Damage.find_last_of("N")+1 == Damage.size()){
          std::cout << "HURA";
          continue;
        }
        else{
          std::cout << "start " << start_pos << " and end " << end_pos << std::endl;
          Deamin_2(Damage,"T",rand_len,1.0,2.0);
          std::cout << Damage << std::endl;
          if (Damage.size() < 150){
            char adapter = 'N';
            double No = (150-Damage.size())/2.0;

            std::cout << Illumina_Qual(sequence) << std::endl;
            outfa << ">" << name << ":" << start_pos << "-" << start_pos+name_len << std::endl;
            outfa << std::string(floor(No),adapter) << Damage << std::string(ceil(No),adapter) << std::endl;
            std::cout << "------------" << std::endl;
          }
        } */
