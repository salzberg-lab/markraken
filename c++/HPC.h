//
// Created by marku on 9/10/2019.
//

#include <iostream>
#include <fstream>

#ifndef MARKRAKEN_HPC_H
#define MARKRAKEN_HPC_H


class HPC {
public:
//    HPC() = default;
    HPC(const std::string & filepath_in, const std::string & filepath_out);
    ~HPC() = default;
    void compress();

private:
    std::string filepath_in; // private copies of file paths
    std::string filepath_out;

    void homopolymer_compress(const std::string & str, std::string & compressed);

};


#endif //MARKRAKEN_HPC_H
