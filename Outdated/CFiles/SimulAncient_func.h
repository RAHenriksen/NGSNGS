#ifndef SIMULANCIENT_FUNC_H
#define SIMULANCIENT_FUNC_H

void Deamin_char(char* str,char nt[],int seed,double alpha=1.0,double beta=2.0,int start=0,int end=25);

double** create2DArray(const char* filename,int width,int height);

void Read_Qual2(char *seq,char *qual,std::discrete_distribution<>*Dist,std::default_random_engine &gen);

void Bam_baseQ(char *seq,char *qual,std::discrete_distribution<>*Dist,std::default_random_engine &gen);

void DNA_complement(char seq[]);

void reverseChar(char* str);

void Qual_dist(double** Array2d,std::discrete_distribution<> dist[],int size);

void Seq_err(double** Array2d,std::discrete_distribution<> nt_sub[],int size);

void Ill_err(char *seq,std::discrete_distribution<>*Dist,std::default_random_engine &gen);

void Size_freq_dist(std::ifstream &infile, std::discrete_distribution<> dist[], int seed);

int *Size_select_dist(std::ifstream &infile);

double uniform();

void SimBriggsModel(char* reffrag, char* frag, int L, double nv, double lambda, double delta_s, double delta, unsigned int seed);

const char* Error_lookup(double a,double err[6000],int nt_offset, int read_pos);

double* Qual_array(double* freqval,const char* filename);

void Read_Qual_new(char *seq,char *qual,unsigned int seed,double* freqval);

int BinarySearch_fraglength(double* SearchArray,int low, int high, double key);

void FragArray(int& number,int*& Length, double*& Frequency,const char* filename);

double myrand(unsigned int persistent);

#endif