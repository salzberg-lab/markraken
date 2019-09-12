//
// Created by marku on 9/11/2019.
//

#ifndef MARKRAKEN_MARKERIZER_H
#define MARKRAKEN_MARKERIZER_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <fstream>

class markerizer {
public:
    markerizer() = default;
    ~markerizer() = default;
    void generate_markers(const uint32_t &n_markers, const uint32_t &marker_length, const uint32_t seed);
    void save_markers(const std::string &filepath);
    void load_markers(const std::string &filepath);

    std::vector<uint16_t> markerize(std::vector<std::string>);
private:
    std::vector<std::string> markers;
    std::string reverse_complement(std::string seq);
};


#endif //MARKRAKEN_MARKERIZER_H
