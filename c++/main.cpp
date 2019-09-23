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


//template<typename K, typename V>
//void print_map(std::unordered_map<K,V> const &m)
//{
//    for (auto const& pair: m) {
//        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
//    }
//}

int main() {
//    std::string filepath_in = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/DB_small.fa";
//    std::string DB_path = "/ccb/salz4-4/markus/markraken/data/refseq_singlestrain/refseq.fa";
    std::string DB_HPC_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_HPC.fa";

//    std::string index_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_index_ml10_mn1000_mk4.markmap";
////    std::string index_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_index_ml10_mn1000_mk5.markmap";
//    std::string marker_path = "/ccb/salz4-4/markus/markraken/data/compressed/mn1000.marker";
//    std::string markerized_seq_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_HPC_ml10_mn_1000.markerized";
//    std::string taxid_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_HPC_ml10_mn_1000.taxid";

//    std::string index_path = "/scratch0/markus/refseq_index_mn8000.markmap";
//    std::string index_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_index_mn8000.markmap";
//    std::string marker_path = "/ccb/salz4-4/markus/markraken/data/compressed/mn8000.marker";
//    std::string markerized_seq_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_HPC_ml10_mn_8000.markerized";
//    std::string taxid_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_HPC_ml10_mn_8000.taxid";

    std::string marker_path = "/ccb/salz4-4/markus/markraken/data/compressed/mn80000.marker";
    std::string markerized_seq_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_HPC_ml10_mn_80000.markerized";
    std::string taxid_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_HPC_ml10_mn_80000.taxid";
    std::string index_path = "/ccb/salz4-4/markus/markraken/data/compressed/refseq_index_ml10_mn80000_mk4.markmap";





//    // homopolymer compress fasta and save as new fasta
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


    // marker parameters
    uint32_t n_markers = 80000;
    uint32_t marker_length = 10;
    uint32_t seed = 2019;


    // generate markers and markerize HPC sequence
    markerizer M;
//    M.generate_markers(n_markers,marker_length, seed);
//    M.save_markers(marker_path);
    M.load_markers(marker_path);


//    // extract all markers from the compressed sequence
//    std::vector<std::vector<uint32_t>> marker_seqs;
//    std::vector<uint32_t> mseq;
//    for (int i=0; i<seqs.size(); i++){
//        mseq = M.markerize(seqs[i]);
//        marker_seqs.push_back(mseq);
//        if (i%100==0){
//            std::cout << i  << " of " << seqs.size() << std::endl;
//            std::cout << mseq[0] << std::endl;
//        }
//    }


//    // save markerized sequence and taxids
//    std::ofstream f_out(markerized_seq_path, std::ios::binary);
//    cereal::BinaryOutputArchive oarchive(f_out);
//    oarchive(marker_seqs);
//    f_out.close();
//
//    std::ofstream f_out2(taxid_path, std::ios::binary);
//    cereal::BinaryOutputArchive oarchive2(f_out2);
//    oarchive2(taxids);
//    f_out2.close();



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

//    for (int i=0; i<10; i++){
//        for (int j=0; j<marker_seqs[i].size();j++){
//            std::cout << marker_seqs[i][j] << std::endl;
//        }
//    }




    // marKmer parameters
    int marK = 4; // length of marKmers, like k for kmers
//
//
//    // load taxonomy for LCA calculation
//    std::string names_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/names.dmp";
//    std::string nodes_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/nodes.dmp";
//    std::string merged_path = "/ccb/salz4-4/markus/markraken/data/ncbi_taxonomy/merged.dmp";
//
//    // build index from marKmers
//    class index marker_index;
//    marker_index.read_taxonomy(names_path, nodes_path, merged_path);
//
//    for (int i=0; i<marker_seqs.size(); i++){
//        if (i % 100 == 0) {
//            std::cout << i << " of " << marker_seqs.size() << std::endl;
//        }
//        std::vector<uint32_t> &singleseq = marker_seqs[i];
//
//        std::cout << i << " of " << marker_seqs.size() << "\t" << singleseq.size()  << std::endl;
//
//        TaxID id = taxids[i];
//        size_t len = (marK > singleseq.size()) ? 0 : singleseq.size() - marK;
//        for (int j = 0; j < len; j++) {
//            std::vector<uint32_t> marker(&singleseq[j], &singleseq[j + marK]); // TODO do this without copying data
//            marker_index.add_pair(marker, id);
//
////            for (int q=0; q<marker.size(); q++){
////                std::cout << marker[q] << " ";
////            }
////            std::cout << std::endl;
//
//        }
//    }
//    marker_index.save_index(index_path);



    // load index
    // TODO index requires loading a taxonomy to not crash, probably due to way it is started, fix this
    std::string names_path = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/taxonomy/names.dmp";
    std::string nodes_path = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/taxonomy/nodes.dmp";
    std::string merged_path = "/ccb/salz4-4/markus/markraken/data/databases/miniSeq+H/taxonomy/merged.dmp";

    class index marker_index;
    marker_index.read_taxonomy(names_path, nodes_path, merged_path);
    marker_index.load_index(index_path);
    std::cout << "index loaded" << std::endl;


//    // example markmer index LCA retrieval
//    std::vector<int> idx_v = {112, 1, 2, 3, 555, 2881, 10, 11};
//    for (int q=0; q<idx_v.size(); q++) {
//        int idx = idx_v[q];
//        for (int i = idx; i < idx + 1; i++) {
//            for (int j = 0; j < marker_seqs[i].size() - marK; j++) {
//                std::vector<uint32_t> testmarKmer(&marker_seqs[i][j], &marker_seqs[i][j + marK]);
//                uint32_t pred = marker_index.predict_taxid(testmarKmer);
//                std::cout << pred << " "; //<< std::endl;
//            }
//            std::cout << "\ncorrect taxid: " << taxids[idx] << "\n\n\n";
//        }
//    }


//    // load simulated reads
//    std::string reads_path = "/ccb/salz4-4/markus/markraken/data/simulated/ecoli/simlord_40k_defaults/reads_small.fasta";
//    // homopolymer compress fasta and save as new fasta
//    HPC compressor;
//    compressor.read(reads_path);
//    compressor.compress();
//    std::vector<std::string> &seqs = compressor.seqs_compressed;
////    std::cout << seqs[0].size();
//
////    markerizer M;
////    for (int i=0; i<seqs.size(); i++){
//    for (int i=0; i<1; i++){
//        std::string &s = seqs[i];
//        std::vector<uint32_t> markers = M.markerize(s);
//        std::vector<uint32_t> taxid_predict;
////        for (int q=0; q<markers.size(); q++){
////            std::cout << markers[q] << " ";
////        }
////        std::cout << "\n";
//        for (int j=0; j<markers.size()-marK; j++){
//            std::vector<uint32_t> marKmer(&markers[j], &markers[j+marK]);
//            uint32_t t = marker_index.predict_taxid(marKmer);
//            taxid_predict.push_back(t);
//        }
//        for (int j=0; j<taxid_predict.size(); j++){
//            std::cout << taxid_predict[j] << " ";
//        }
//
//    }

    // assign reads
//    std::string reads_path = "/ccb/salz4-4/markus/markraken/data/Zymo/Zymo-GridION-EVEN-BB-SN-PCR-R10HC-flipflop.fasta";
    std::string reads_path = "/ccb/salz4-4/markus/markraken/data/Zymo/small.fasta";
//    std::string reads_path = "/ccb/salz4-4/markus/markraken/data/simulated/ecoli/art/072119.fa";
    // homopolymer compress fasta and save as new fasta
    HPC compressor;
    compressor.read(reads_path);
    compressor.compress();
    std::vector<std::string> &seqs = compressor.seqs_compressed;
//    std::cout << seqs[0].size();

//    markerizer M;
//    for (int i=0; i<seqs.size(); i++){
    for (int i=0; i<50; i++){
        std::string &s = seqs[i];
        std::vector<uint32_t> markers = M.markerize(s);
        std::vector<uint32_t> taxid_predict;

        size_t len = (marK > markers.size()) ? 0 : markers.size() - marK;
        for (int j=0; j<len; j++){
            std::vector<uint32_t> marKmer(&markers[j], &markers[j+marK]);
            uint32_t t = marker_index.predict_taxid(marKmer);
            taxid_predict.push_back(t);
        }
        for (int j=0; j<taxid_predict.size(); j++){
            std::cout << taxid_predict[j] << " ";
        }
        std::cout << "\n\n";

    }




//    print_map(marker_index.markmap);




}