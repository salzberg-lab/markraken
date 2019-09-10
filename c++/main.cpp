#include <iostream>
#include "HPC.h"

int main() {
//    std::cout << "Hello, World!" << std::endl;
    std::string filepath_in = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq/DB.fa";
    std::string filepath_out = "/ccb/salz4-4/markus/markraken/data/compressed/DB_HPC.fa";

    HPC compressor(filepath_in, filepath_out);

    return 0;

}