#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <iterator>
#include <ctime>
#include<stdio.h>

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

std::vector<std::string> DNA_fragmentation(std::string Seq,std::string nt_qual,int min_len, int max_len){
    int seq_max = Seq.size();
    int sub_str_start = 0;

    std::vector<std::string> Frag_Deamin_seq;

    std::srand(std::time(nullptr)); // initialize with computer time
    while(sub_str_start <= seq_max){
        // creates random number in the range of the fragment size rand() % ( high - low + 1 ) + low
        int rand = (std::rand() % (max_len - min_len + 1)) + min_len;

        std::string Frag = Seq.substr(sub_str_start,rand);
        std::string Qual = nt_qual.substr(sub_str_start,rand);
        sub_str_start += rand;

        Frag_Deamin_seq.push_back(Frag);
        Frag_Deamin_seq.push_back(Qual);
        //std::cout << Frag << std::endl;
        std::srand(rand); // restarts seed everytime with a new random value
    }
    return Frag_Deamin_seq;
}

void Deamin_2(std::string& Seq, std::string nt_find,
              double alpha=1.0,double beta=2.0,int start=0,int end=25){
    //declaring random seed to build random generator for gamma distribution
    //std::srand(Seq.size());
    std::srand(std::time(nullptr));
    std::default_random_engine generator(std::rand());
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

std::istream& operator>>(std::istream& in, Fastq_structure& frag){
    //extraction operator

    std::getline(in, frag.ID);

    if (frag.ID.size() == 0 || frag.ID[0] != '@') {
        // checks the read format is correct by setting the stream as faulty
        in.setstate(std::ios_base::failbit);
        return in;
    }

    std::getline(in, frag.sequence);

    std::getline(in, frag.strand);
    if (frag.strand.size() == 0) {
        in.setstate(std::ios_base::failbit);
        return in;
    }

    std::getline(in, frag.quality_value);
    return in;
}


std::ostream& operator<<(std::ostream& out, const Fastq_structure& frag)
{

    std::vector<std::string> Damage = DNA_fragmentation(frag.sequence,frag.quality_value,30,80);

    int j = 1;
    for (int i=0; i<Damage.size();i+=2){
        //out << frag.ID << '\n';
        std::cout << "----------------------------" << '\n';
        out << "read" << '\n';
        Deamin_2(Damage[i],"T",1.0,2.0);
        std::cout << "test " << Damage[i] << std::endl;
        out << Damage[i] << '\n';
        out << frag.strand << '\n';
        out << Damage[i+1] << '\n';
        j += 1;

    }

    return out;
}

int main()
{
    std::ifstream in("Input.fastq");

    std::vector<Fastq_structure> frags;

    for (Fastq_structure tmp; in >> tmp;) {
        frags.push_back(tmp);
    }

    /*
    for (int i = 0; i < frags.size(); i++) {
		std::cout << frags.at(i) << ' ';
	}*/

    // Insert code for mutating the fragments

    std::ofstream out("output.fq");
    for (const auto& f : frags)
        out << f;
}

