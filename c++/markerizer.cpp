//
// Created by marku on 9/11/2019.
//

#include <set>
#include <cstring>
#include <unordered_map>
#include "markerizer.h"

void markerizer::generate_markers(const uint32_t &n_markers, const uint32_t &marker_length, const uint32_t seed) {
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
