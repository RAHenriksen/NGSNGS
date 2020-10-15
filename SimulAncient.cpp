#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <iterator>
#include <ctime>
#include<stdio.h>


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

std::string DNA_fragmentation(std::string Seq,int min_len, int max_len){
    std::srand(std::time(nullptr));
    int rand = std::rand() % max_len + min_len;
    std::cout << "random number" << rand << std::endl;
    std::cout << Seq.size() << std::endl;
    std::string DNA_frag;
    int start = 0;
    if(Seq.size() > max_len) {
            DNA_frag = Seq.substr(start,max_len);
    }
    return DNA_frag;
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
    return Deamin;
}

std::string fafq_seq(std::string in_name, std::string out_name) {
    std::ifstream myfile(in_name);
    std::ofstream out_file(out_name);
    if (myfile.is_open() && out_file.is_open())
    {
        int i = 1;
        std::string line;
        while( std::getline(myfile,line) )
        {
            if (line.rfind("@",0) == 0) {
                out_file << "@Read_" << i << std::endl;
                i += 1;
            }
            else if (line.rfind("A", 0) == 0 ||
                line.rfind("T", 0) == 0 ||
                line.rfind("G", 0) == 0 ||
                line.rfind("C", 0) == 0){

                std::string Seq = line;
                std::vector<int> Idx_vec; // declaration of unsigned integer vector index

                nt_index(Idx_vec, Seq , "T"); //the nt_index function assigns position values to the declared Idx_vec
                std::string test = DNA_deamination(Seq,Idx_vec,2);

                out_file << test << std::endl;
            }
            else {
                out_file << line << std::endl;
            }
        }
        out_file.close();
        myfile.close();
    }
}

void File_type(std::string in_name, std::string out_name) {
    // Checks file type for the input files
    std::string fileext = in_name;
    int ext_pos = fileext.rfind(".");
    std::string ext = fileext.substr(ext_pos+1);

    //enum extensions {fq,fastq,fa,fasta,sam,bam,cram,vcf}
    if (ext == "fq" || ext == "fastq") {
        std::cout << ext << std::endl;
        fafq_seq(in_name,out_name);
    }
    else if (ext == "fa" || ext == "fasta")
        {
        std::cout << ext << std::endl;
        fafq_seq(in_name,out_name);
    }
    else if (ext == "sam")
        {
        std::cout << ext << std::endl;
    }
    else if (ext == "bam")
        {
            std::cout << ext << std::endl;
    }
    else if (ext == "cram")
        {
        std::cout << ext << std::endl;
    }
    else if (ext == "vcf")
        {
        std::cout << ext << std::endl;
    }
    else {
        std::cout << "Unable to open file, wrong format, the files needs to be fa,fasta...";
    }
}

int main() {

    //std::string File = "BC06_7.fastq";
    //File_type(File,"read222.fq");

    //std::string Seq = "ATGGTACCATAGTGGACCATATCAGCATATCATAATCGTTCGCGAATACTGAT";
    //std::cout << DNA_fragmentation(Seq,10,20);

    int frag_min = 30;
    int frag_max = 80;
    int seq_max = 170;

    std::srand(std::time(nullptr));
    int rand = std::rand() % (frag_max - frag_min + 1) + frag_min; // creating range: rand() % ( high - low + 1 ) + low
    std::cout << "random number " << rand << std::endl;
    int i = rand;
    while(i <= seq_max){
        std::cout << i << std::endl;
        i+=rand;
    }

    /*
        std::string Seq = "ATGGATGGTACCATAGTGGACCATATCAGCATATCATAATCGTTCGCGAATACTGATTATGGTACCATAGTGGACCATATCAGCATATCATAATCGTTCGCGAATACTGATACCATAGTGGACCATATCAGCATATCATAATCGTTCGCGAATACTGAT";
    int frag_min = 30;
    int frag_max = 80;
    int seq_max = 170;

    std::srand(std::time(nullptr));
    int rand = std::rand() % (frag_max - frag_min + 1) + frag_min; // creating range: rand() % ( high - low + 1 ) + low
    std::cout << "random number " << rand << std::endl;
    int i = rand;
    int j = 0;
    while(i <= Seq.size()){
        std::cout << i << std::endl;
        i+=rand;
        std::cout << Seq.substr(j,i) << std::endl;
        j+=i;
    }
    */
    return 0;
}
