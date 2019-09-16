#include <iostream>
#include <algorithm>
#include "HPC.h"
#include "include/NcbiTaxonomy.h"
#include "index.h"
#include "markerizer.h"
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

//#include "include/smhasher/MurmurHash3.h"
//#include "include/kraken2/compact_hash.h"
//#include <vector>

int main() {
//    std::string filepath_in = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/DB_small.fa";
//    std::string DB_path = "/ccb/salz4-4/markus/markraken/data/refseq_singlestrain/refseq.fa";
    std::string DB_HPC_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_HPC.fa";

//    std::string index_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_index.markmap";
//    std::string marker_path = "/ccb/salz4-4/markus/markraken/data/compressed/mn1000.marker";
//    std::string markerized_seq_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_HPC_ml10_mn_1000.markerized";
//    std::string taxid_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_HPC_ml10_mn_1000.taxid";

    std::string index_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_index_mn8000.markmap";
    std::string marker_path = "/ccb/salz4-4/markus/markraken/data/compressed/mn8000.marker";
    std::string markerized_seq_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_HPC_ml10_mn_8000.markerized";
    std::string taxid_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_HPC_ml10_mn_8000.taxid";

//    // homopolymer compress fasta and save as new fasta
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


    // marker parameters
    uint32_t n_markers = 8000;
    uint32_t marker_length = 10;
    uint32_t seed = 2019;


    // generate markers and markerize HPC sequence
    markerizer M;
    M.generate_markers(n_markers,marker_length, seed);
    M.save_markers(marker_path);
//    M.load_markers(marker_path);


    // extract all markers from the compressed sequence
    std::vector<std::vector<uint32_t>> marker_seqs;
    std::vector<uint32_t> mseq;
    for (int i=0; i<seqs.size(); i++){
        if (i%100==0){
            std::cout << i  << " of " << seqs.size() << std::endl;
        }
        mseq = M.markerize(seqs[i]);
        marker_seqs.push_back(mseq);
    }


    // save markerized sequence and taxids
    std::ofstream f_out(markerized_seq_path, std::ios::binary);
    cereal::BinaryOutputArchive oarchive(f_out);
    oarchive(marker_seqs);
    f_out.close();

    std::ofstream f_out2(taxid_path, std::ios::binary);
    cereal::BinaryOutputArchive oarchive2(f_out2);
    oarchive2(marker_seqs);
    f_out2.close();


//    // load markerized sequence and taxids
//    std::vector<std::vector<uint32_t>> marker_seqs;
//
//    std::ifstream f_in(markerized_seq_path, std::ios::binary);
//    cereal::BinaryInputArchive iarchive(f_in);
//    iarchive(marker_seqs);
//    f_in.close();
//
//    std::vector<int> taxids;
//    std::ifstream f_in2(taxid_path, std::ios::binary);
//    cereal::BinaryInputArchive iarchive2(f_in2);
//    iarchive2(taxids);
//    f_in2.close();
//
//    for (int i=0; i<10; i++){
//        for (int j=0; j<marker_seqs[i].size();j++){
//            std::cout << marker_seqs[i][j] << std::endl;
//        }
//    }


//    // load taxonomy for LCA calculation
//    std::string names_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/names.dmp";
//    std::string nodes_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/nodes.dmp";
//    std::string merged_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/merged.dmp";


    // marKmer parameters
    int marK = 3; // length of marKmers, like k for kmers

//    // build index from marKmers
//    class index marker_index;
//    marker_index.read_taxonomy(names_path, nodes_path, merged_path);
//
//    for (int i=0; i<marker_seqs.size(); i++){
//        std::vector<int> FOOSKIP; // TODO figure out why some make it stall, this is a temporary solution
//        FOOSKIP.push_back(1763);
//        FOOSKIP.push_back(2054);
//        FOOSKIP.push_back(2857);
//        FOOSKIP.push_back(2858);
//        FOOSKIP.push_back(2859);
//        FOOSKIP.push_back(2867);
//        FOOSKIP.push_back(3018);
//        FOOSKIP.push_back(3380);
//        FOOSKIP.push_back(4607);
//        FOOSKIP.push_back(5029);
//        FOOSKIP.push_back(5040);
//        FOOSKIP.push_back(5858);
//
//        if (std::find(FOOSKIP.begin(), FOOSKIP.end(), i) == FOOSKIP.end()) {
//
//            std::cout << taxids[i] << std::endl;
//
//            if (i % 100 == 0) {
//                std::cout << i << " of " << marker_seqs.size() << std::endl;
//            }
//            std::cout << i << " of " << marker_seqs.size() << std::endl;
//
//            std::vector<uint32_t> &singleseq = marker_seqs[i];
//            TaxID id = taxids[i];
//            for (int j = 0; j < singleseq.size() - marK; j++) {
//                std::vector<uint32_t> marker(&singleseq[j], &singleseq[j + marK]); // TODO do this without copying data
//                marker_index.add_pair(marker, TaxID(id));
//            }
//        }
//    }
//    marker_index.save_index(index_path);


//    // load index
//    // TODO index requires loading a taxonomy to not crash, probably due to way it is started, fix this
//    std::string names_path = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/taxonomy/names.dmp";
//    std::string nodes_path = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/taxonomy/nodes.dmp";
//    std::string merged_path = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/taxonomy/merged.dmp";
//
//    class index marker_index;
//    marker_index.read_taxonomy(names_path, nodes_path, merged_path);
//    marker_index.load_index(index_path);
//
//    int idx = 112;
//    for (int i=idx; i<idx+1; i++){
//        for (int j=0; j<marker_seqs[i].size()-marK; j++){
//            std::vector<uint32_t> testmarKmer(&marker_seqs[i][j], &marker_seqs[i][j+marK]);
//            uint32_t pred = marker_index.predict_taxid(testmarKmer);
//            std::cout << pred << " "; //<< std::endl;
//        }
//        std::cout << "\ncorrect taxid: " << taxids[idx] << std::endl;
//    }













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