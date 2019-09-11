//
//  Object to compress homopolymers of fasta file
//

#include <iostream>
#include <fstream>
#include <vector>
//#include <zlib.h>
//#include <stdio.h>
//#include "include/kseq/kseq.hpp"
#include "include/FastaTools.h"
#include "include/NcbiTaxonomy.h"


#ifndef MARKRAKEN_HPC_H
#define MARKRAKEN_HPC_H


class HPC {
public:
    HPC() = default;
    ~HPC() = default;
    void read(const std::string &filepath_in);
    void compress();
    void write(const std::string &filepath_out);

    void extract_miniseq_taxids();
    std::vector<TaxID> taxids;


private:
    std::vector <std::string> seqs;
    std::vector <std::string> infolines;
    std::vector <std::string> seqs_compressed;
    int n_seqs;

    void homopolymer_compress(const std::string &s, std::string &s_compressed);


};


#endif //MARKRAKEN_HPC_H
