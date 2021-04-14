#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include<stdio.h>

std::string fafq_seq(std::string file_name) {
    std::string line;
    std::ifstream myfile (file_name);

    if (myfile.is_open())
    {
        while( std::getline(myfile,line) )
        {
            if (line.rfind("A", 0) == 0 ||
                line.rfind("T", 0) == 0 ||
                line.rfind("G", 0) == 0 ||
                line.rfind("C", 0) == 0){
                std::cout << "hurra" << line << std::endl;
            }
            else {
                std::cout << "lol" << '\n';
            }
        }
        myfile.close();
    }
}

void File_type(std::string file_name) {
    std::string fileext = file_name;
    int ext_pos = fileext.rfind(".");
    std::string ext = fileext.substr(ext_pos+1);

    //enum extensions {fq,fastq,fa,fasta,sam,bam,cram,vcf}
    if (ext == "fq" || ext == "fastq") {
        std::cout << ext << std::endl;
        fafq_seq(file_name);
    }
    else if (ext == "fa" || ext == "fasta")
        {
        std::cout << ext << std::endl;
        fafq_seq(file_name);
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
    //std::string File = "Test.fastq";
    std::string File = "Read1.bam";

    File_type(File);

}
