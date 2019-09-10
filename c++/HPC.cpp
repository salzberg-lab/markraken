//
// Object to compress homopolymers of fasta file
//

#include "HPC.h"

HPC::HPC(const std::string &filepath_in, const std::string & filepath_out) {
    this->filepath_in = filepath_in; // private copies of file paths
    this->filepath_out = filepath_out;
}

void HPC::compress() {
    //    read stuff here
    const std::string stringfoo = "ABCDDDR";

    int i = stringfoo.size();
    std::string letters;
    letters.reserve(i); // HPC sequence will always be less than or equal to size of original seq

    homopolymer_compress(stringfoo, letters);

    // write stuff here
}

void HPC::homopolymer_compress(const std::string &str, std::string &compressed) {
    int i = str.size();
    for (int j=0; j<i; ++j){
        while (str[j] == str[j+1]){
            j++;
        }
        compressed.push_back(str[j]);
    }
}


