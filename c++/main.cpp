#include <iostream>
#include "HPC.h"
#include "include/NcbiTaxonomy.h"
//#include <vector>

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

    // read in HPC compressed seqs and extract taxids from miniseq format
    HPC reader;
    reader.read(DB_HPC_path);
    reader.extract_miniseq_taxids();
    std::vector<int> taxids = reader.taxids;

    // load taxonomy for LCA calculation
    std::string names_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/names.dmp";
    std::string nodes_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/nodes.dmp";
    std::string merged_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/merged.dmp";
    NcbiTaxonomy tax(names_path, nodes_path, merged_path);

    // example taxid calculation
    TaxID &id_1 = taxids[0];
    TaxID &id_2 = taxids[1];
    std::vector<TaxID> taxids2;
    taxids2.push_back(id_1);
    taxids2.push_back(id_2);
    std::cout << tax.LCA(taxids2)->id << std::endl;



    return 0;

}