//
// Created by marku on 9/11/2019.
//

#include <set>
#include <cstring>
#include <unordered_map>
#include "markerizer.h"

void markerizer::generate_markers(const uint32_t &n_markers, const uint32_t &marker_length, const uint32_t seed) {
    this->marker_length = marker_length;
    // generates a set of ~n unique homopolymer compressed markers and their reverse complements
    std::set<std::string> markerset;
    char alphabet[] = "ATCG";
    srand(seed);
    std::string marker;
    while (markerset.size() <= n_markers){
        marker.reserve(marker_length);
        char m = alphabet[rand()%strlen(alphabet)];
        marker += m;
        for (int i=0; i<marker_length; i++){
            std::set<char> reduced_alphabet;
            for (int j=0; j<strlen(alphabet); j++){
                reduced_alphabet.insert(alphabet[j]);
            }
            reduced_alphabet.erase(marker[marker.size()-1]);
            std::vector<char> v(reduced_alphabet.begin(), reduced_alphabet.end());
            char s = v[rand()%v.size()];
            marker += s;
        }
        markerset.insert(marker);
        markerset.insert(reverse_complement(marker));
        marker.clear();
        marker.reserve(marker_length);
    }
    markers.assign(markerset.begin(), markerset.end());
}

std::string markerizer::reverse_complement(std::string seq) {
    std::string revcomp;
    std::unordered_map<char, char> mapper;
    mapper.insert(std::make_pair('A', 'T'));
    mapper.insert(std::make_pair('T', 'A'));
    mapper.insert(std::make_pair('G', 'C'));
    mapper.insert(std::make_pair('C', 'G'));
    for (int i=seq.size()-1; i>=0; i--){
        revcomp += mapper[seq[i]];
    }
    return revcomp;
}

void markerizer::save_markers(const std::string &filepath) {
    std::ofstream f_out(filepath, std::ios::binary);
    cereal::BinaryOutputArchive oarchive(f_out);
    oarchive(markers);
    f_out.close();
}

void markerizer::load_markers(const std::string &filepath) {
    std::ifstream f_in(filepath, std::ios::binary);
    cereal::BinaryInputArchive iarchive(f_in);
    iarchive(markers);
    f_in.close();
    this->marker_length = markers[0].size();
}

std::vector<uint32_t> markerizer::markerize(std::string nuc_seq) {
    std::vector<uint32_t> marker_seq;
    marker_seq.reserve(nuc_seq.size());
    std::string subseq;
    subseq.reserve(marker_length);

    //https://stackoverflow.com/questions/19481662/c-index-of-element-in-sorted-stdvector
    for (int i=0; i<nuc_seq.size()-marker_length; i++){
        subseq = nuc_seq.substr(i, marker_length); // copying is not great, c++17 is required for good string viewing, probably better to use char*
        auto lower = std::lower_bound(markers.begin(), markers.end(), subseq);
        bool found = lower != markers.end() && *lower == subseq;
        if (found){
            auto idx = std::distance(markers.begin(), lower);
            marker_seq.push_back(idx);
        }
    }
    return marker_seq;
}
