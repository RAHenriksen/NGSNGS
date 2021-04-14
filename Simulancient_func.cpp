#include "SimulAncient_func.h"

void Deamin_char(char* str,char nt[],int seed,double alpha,double beta,int start,int end){   
  // Deamination of nucleotides
  std::vector<int> Index_vec;
  std::srand(seed+std::time(nullptr));
  std::random_device rd;
  std::default_random_engine generator(rd());
  std::gamma_distribution<double> distr(alpha,beta);

  int i = strcspn(str,nt);
  Index_vec.push_back(i);

  while(i < end) {
    int tmp = strcspn(str+i+1,nt);
    i += tmp + 1;
      Index_vec.push_back(i);
  }

  for (int i = 0; i < Index_vec.size(); i++){
    if (Index_vec.at(i) == int(distr(generator))) {
      //remember to create an input for the nt so it works for both 5' and 3' 
      str[Index_vec.at(i)]='T';
      }
		else {
      continue;
		}
  }
}