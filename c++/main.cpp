#include <iostream>
#include "HPC.h"
#include "include/NcbiTaxonomy.h"
#include "include/smhasher/MurmurHash3.h"
//#include "include/kraken2/compact_hash.h"
//#include <vector>

std::size_t combine(std::vector<uint16_t> const& vec){
    //https://en.wikipedia.org/wiki/Tiny_Encryption_Algorithm
    //https://github.com/HowardHinnant/hash_append/issues/7
    std::size_t seed = vec.size();
    for(auto& i : vec) {
        seed ^= i + 0x9E37 + (seed << 6) + (seed >> 2);
    }
    return seed;
}

int main() {
//    std::string filepath_in = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/DB_small.fa";
//    std::string filepath_in = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/DB.fna";
//    std::string filepath_in = "/ccb/salz4-4/markus/genome_data/ecoli_MG1655/GCF_000005845.2_ASM584v2_genomic.fna";
    std::string DB_HPC_path = "/ccb/salz4-4/markus/markraken/data/compressed/DB_HPC_small.fa";

    // homopolymer compress fasta and save as new fasta
//    HPC compressor;
//    compressor.read(DB_path);
//    compressor.compress();
//    compressor.write(DB_HPC_path);

//    // read in HPC compressed seqs and extract taxids from miniseq format
//    HPC reader;
//    reader.read(DB_HPC_path);
//    reader.extract_miniseq_taxids();
//    std::vector<int> &taxids = reader.taxids;
//    std::vector<std::string> &seqs = reader.seqs;

    // markerize the compressed sequence
//    std::cout << seqs[0] << std::endl;


//    // load taxonomy for LCA calculation
//    std::string names_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/names.dmp";
//    std::string nodes_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/nodes.dmp";
//    std::string merged_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/merged.dmp";
//    NcbiTaxonomy tax(names_path, nodes_path, merged_path);
//
//    // example taxid calculation
//    TaxID &id_1 = taxids[0];
//    TaxID &id_2 = taxids[1];
//    std::vector<TaxID> taxids2;
//    taxids2.push_back(id_1);
//    taxids2.push_back(id_2);
//    std::cout << tax.LCA(taxids2)->id << std::endl;


    uint16_t foo_1 = 8000;
    uint16_t foo_2 = 1;
    uint16_t foo_3 = 99;

    std::vector<uint16_t> marker = {foo_1, foo_2, foo_3};

    std::size_t foo = combine(marker);
    std::cout << foo << std::endl;

    std::unordered_map<uint32_t, uint64_t> markmap;
//    markmap.reserve();
    markmap.insert({foo, 1});
    std::cout << markmap[foo] << std::endl;


    return 0;

}