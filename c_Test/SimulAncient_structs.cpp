#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <iterator>
#include <ctime>
#include<stdio.h>


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
        // std::cout << Frag << std::endl;
        std::srand(rand); // restarts seed everytime with a new random value
    }
    return Frag_Deamin_seq;
}

struct FastqFragment
{
    std::string ID;
    std::string sequence;
    std::string strand;
    std::string quality_value;
};


std::istream& operator>>(std::istream& in, FastqFragment& frag){
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

std::ostream& operator<<(std::ostream& out, const FastqFragment& frag)
{

    std::vector<std::string> Damage = DNA_fragmentation(frag.sequence,frag.quality_value,30,80);

    int j = 1;
    for (int i=0; i<Damage.size();i+=2){
        out << frag.ID << "_" << j << '\n';
        out << Damage[i] << '\n';
        out << frag.strand << '\n';
        out << Damage[i+1] << '\n';
        j += 1;
    }
    return out;
}

int main()
{
    std::ifstream in("example.fq");

    std::vector<FastqFragment> frags;

    for (FastqFragment tmp; in >> tmp;) {
        frags.push_back(tmp);
    }

    /*
    for (int i = 0; i < frags.size(); i++) {
		std::cout << frags.at(i) << ' ';
	}*/

    // Insert code for mutating the fragments

    std::ofstream out("output.txt");
    for (const auto& f : frags)
        out << f;
}
