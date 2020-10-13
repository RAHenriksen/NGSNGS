#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <iterator>
#include <ctime>


void nt_index(std::vector<int> & Idx_vec, std::string Seq, std::string nt_find)
{
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

std::string DNA_fragmentation(){
    return 0;
}

std::string DNA_deamination(std::string Seq, std::vector<int> & Idx_vec, int type){
    // initially i though i would create a function with type to define the damage type

    double deamin_rate = 0.40;
    int vec_size = Idx_vec.size();
    int Nt_no = vec_size * deamin_rate;
    std::cout << "vecsize " << vec_size << " and number of changes " << Nt_no << std::endl;
    std::srand(std::time(nullptr)); // initialize random with current time as seed generator
    std::cout << Seq << std::endl;

    std::string Deamin = Seq;
    for (int i = 0; i < Nt_no;i++){
        int rand_index = std::rand() % vec_size; // pick a random index
        int value = Idx_vec[rand_index]; // selecting random position of nt
        Deamin.replace(value,1,"U");
    }
    std::cout << Deamin << std::endl;
    return Deamin;
}

int main()
{
    std::string Seq = "ATGGTACCATAGTGGACCATATCAGCATATCATAATCGTTCGCGAATACTGAT";
    std::vector<int> Idx_vec; // declaration of unsigned integer vector index

    nt_index(Idx_vec, Seq , "T");
    std::cout << "Size of sequence are " << Seq.size() << " nucleotides" << std::endl;

    for(int pos : Idx_vec)
        std::cout << pos << "\t";

    std::string test = DNA_deamination(Seq,Idx_vec,2);
    std::cout << test ;
    return 0;

}

