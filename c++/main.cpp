#include <iostream>
#include "HPC.h"
#include "include/NcbiTaxonomy.h"

int main() {
//    std::string filepath_in = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/DB_small.fa";
//    std::string filepath_in = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/DB.fna";
//    std::string filepath_in = "/ccb/salz4-4/markus/genome_data/ecoli_MG1655/GCF_000005845.2_ASM584v2_genomic.fna";
//    std::string DB_HPC_path = "/ccb/salz4-4/markus/markraken/data/compressed/DB_HPC_small.fa";

//    HPC compressor;
//    compressor.read(DB_path);
//    compressor.compress();
//    compressor.write(DB_HPC_path);

//    HPC reader;
//    reader.read(DB_HPC_path);

    std::string names_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/names.dmp";
    std::string nodes_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/nodes.dmp";
    std::string merged_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/merged.dmp";
    NcbiTaxonomy tax(names_path, nodes_path, merged_path);

    TaxID id_1 = 562;
    TaxID id_2 = 632;

    TaxID id_LCA = tax.LCA(id_1, id_2);
    std::cout << id_LCA << std::endl;



    return 0;

}