#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <iterator>
#include <ctime>


void nt_index(std::vector<int> & Idx_vec, std::string Seq, std::string nt_find)
{   // use & to pass by reference to the vector
    // Get the first occurrence
    int pos = Seq.find(nt_find);

    // continue while the position index is smaller than sequence
    while( pos < Seq.size()) //std::string::npos
    {
        // append position to the vector
        Idx_vec.push_back(pos);

        // Get the next occurrence from the current position
        pos =Seq.find(nt_find, pos + nt_find.size());
    }
}

std::string DNA_deamination(std::string Seq, std::vector<int> & Idx_vec, int type){
    // initially i though i would create a function with type to define the damage type

    double deamin_rate = 0.40;
    int vec_size = Idx_vec.size();
    int Nt_no = vec_size * deamin_rate;

    std::srand(std::time(nullptr)); // initialize random with current time as seed generator

    std::string Deamin = Seq;
    for (int i = 0; i < Nt_no;i++){
        int rand_index = std::rand() % vec_size; // pick a random index of the T-index vector
        Deamin.replace(Idx_vec[rand_index],1,"U"); // changes the randomly picked T
    }
    //std::cout << Deamin << std::endl;
    return Deamin;
}


void Deamin_2(std::string& Seq, std::string nt_find,
              double alpha=1.0,double beta=2.0,int start=0,int end=10){
    //declaring random seed to build random generator for gamma distribution
    std::srand(std::time(nullptr));
    std::default_random_engine generator(std::rand());
    std::gamma_distribution<double> distr(alpha,beta);

    std::vector<int> Index_vec;

    std::string Sub = Seq.substr(start,end);
    // std::cout << "string pos " << Sub << std::endl;

    int pos = Sub.find(nt_find);
    while( pos < Seq.size()) //std::string::npos
    {
        Index_vec.push_back(pos);
        pos =Sub.find(nt_find, pos + nt_find.size());
    }

    for (int i = 0; i < Index_vec.size(); i++){
        std::cout << Index_vec.at(i) << std::endl;
        if (Index_vec.at(i) == int(distr(generator))) {
            //std::cout << Index_vec.at(i) << std::endl;
            std::cout << "rand number " << int(distr(generator)) << std::endl;
            //std::cout << "test1  " << int(distr(generator)) << std::endl;
            Seq.replace(Index_vec[i],1,"U");
        }
		else {
            continue;
		}
    }
}

int main(){
    std::string Seq1 = "TTTTTTTTTTTTTTTTTT";
    std::string Seq2 = "TTCCATATCAGCA";
    std::string Seq3 = "TATCAGCATATC";
    std::string Seq4 = "TCGCGAATACTGAT";
    std::cout << "BEFORE "<< Seq2 << std::endl;
    //Deamin_2(Seq1,"T",0.8,4.0);
    Deamin_2(Seq2,"T");
    std::cout << "AFTER " << Seq2;
    /*

    std::vector<int> Idx_vec; // declaration of unsigned integer vector index

    nt_index(Idx_vec, Seq , "T"); //the nt_index function assigns position values to the declared Idx_vec
    std::cout << "Size of sequence are " << Seq.size() << " nucleotides" << std::endl;

    for(int pos : Idx_vec)
        std::cout << pos << "\t";

    std::string test = DNA_deamination(Seq,Idx_vec,2);
    */
    return 0;

}

