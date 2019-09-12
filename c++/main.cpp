#include <iostream>
#include <algorithm>
#include "HPC.h"
#include "include/NcbiTaxonomy.h"
#include "index.h"
#include "markerizer.h"
//#include "include/smhasher/MurmurHash3.h"
//#include "include/kraken2/compact_hash.h"
//#include <vector>

std::string random_string( size_t length ){
    // https://stackoverflow.com/questions/440133/how-do-i-create-a-random-alpha-numeric-string-in-c
    auto randchar = []() -> char{
        const char charset[] = "ATCG";
        const size_t max_index = (sizeof(charset) - 1);
        return charset[ rand() % max_index ];
    };
    std::string str(length,0);
    std::generate_n( str.begin(), length, randchar );
    return str;
}


int main() {
//    std::string filepath_in = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/DB_small.fa";
//    std::string filepath_in = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/DB.fna";
//    std::string filepath_in = "/ccb/salz4-4/markus/genome_data/ecoli_MG1655/GCF_000005845.2_ASM584v2_genomic.fna";
    std::string DB_HPC_path = "/ccb/salz4-4/markus/markraken/data/compressed/DB_HPC_small.fa";
    std::string index_path = "/ccb/salz4-4/markus/markraken/data/compressed/index.markmap";
    std::string marker_path = "/ccb/salz4-4/markus/markraken/data/compressed/m.marker";



    // homopolymer compress fasta and save as new fasta
//    HPC compressor;
//    compressor.read(DB_path);
//    compressor.compress();
//    compressor.write(DB_HPC_path);

    // read in HPC compressed seqs and extract taxids from miniseq format
    HPC reader;
    reader.read(DB_HPC_path);
    reader.extract_miniseq_taxids();
    std::vector<int> &taxids = reader.taxids;
    std::vector<std::string> &seqs = reader.seqs;

    // generate markers and save to file for all future use
    uint32_t n_markers = 8000;
    uint32_t marker_length = 10;
    uint32_t seed = 2019;

    markerizer M;
    M.generate_markers(n_markers,marker_length, seed);
    M.save_markers(marker_path);
    M.load_markers(marker_path);

//    std::vector<std::string> markers;
//    markers.reserve(n_markers);
//    std::cout << random_string(marker_length) << std::endl;


    // extract all markers from the compressed sequence
    std::vector<std::vector<uint32_t>> marker_seqs;

//    marker_seqs  =
//    std::cout << seqs[0] << std::endl;


//    // load taxonomy for LCA calculation
//    std::string names_path = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/taxonomy/names.dmp";
//    std::string nodes_path = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/taxonomy/nodes.dmp";
//    std::string merged_path = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/taxonomy/merged.dmp";
////    std::string names_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/names.dmp";
////    std::string nodes_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/nodes.dmp";
////    std::string merged_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/merged.dmp";
//
//
//
//    // build index example
//    class index marker_index;
//    marker_index.read_taxonomy(names_path, nodes_path, merged_path);
//
//    uint32_t foo_1 = 8000;
//    uint32_t foo_2 = 1;
//    uint32_t foo_3 = 99;
//
//    TaxID &id_1 = taxids[0];
//    TaxID &id_2 = taxids[1];
//
//    std::vector<uint32_t> marker = {foo_1, foo_2, foo_3};
//
//    marker_index.add_pair(marker, TaxID(id_1));
//    marker_index.add_pair(marker, TaxID(id_2));
//
//    std::size_t foo = marker_index.predict_taxid(marker);
//    std::cout << foo << std::endl;
//    marker_index.save_index(index_path);
//
//    class index marker_index2;
//    marker_index2.load_index(index_path);
//    uint32_t bar = marker_index2.predict_taxid(marker);
//    std::cout << bar << std::endl;
//    return 0;
}