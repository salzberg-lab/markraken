#include <iostream>
#include "HPC.h"
#include "include/NcbiTaxonomy.h"

int main() {
//    std::string filepath_in = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/DB_small.fa";
    std::string filepath_in = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/DB.fna";
//    std::string filepath_in = "/ccb/salz4-4/markus/genome_data/ecoli_MG1655/GCF_000005845.2_ASM584v2_genomic.fna";
    std::string filepath_out = "/ccb/salz4-4/markus/markraken/data/compressed/DB_HPC_small.fa";

    HPC compressor;
    compressor.read(filepath_in);
    compressor.compress();
    compressor.write(filepath_out);

    return 0;

}