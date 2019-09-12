//
// Created by marku on 9/11/2019.
//

#ifndef MARKRAKEN_MARKERIZER_H
#define MARKRAKEN_MARKERIZER_H

#include <iostream>
#include <vector>

class markerizer {
public:
    markerizer() = default;
    ~markerizer() = default;
    void generate_markers(const int &n_markers, const int &marker_length, const int seed);
    void save_markers(const std::string &filepath);
    void load_markers(const std::string &filepath);

    std::vector<uint16_t> markerize(std::vector<std::string>);
};


#endif //MARKRAKEN_MARKERIZER_H
