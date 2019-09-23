//
// Created by marku on 9/11/2019.
//

#ifndef MARKRAKEN_INDEX_H
#define MARKRAKEN_INDEX_H


#include <vector>
#include <iostream>
#include <unordered_map>
#include "include/NcbiTaxonomy.h"
#include <cereal/archives/binary.hpp>
#include <cereal/types/unordered_map.hpp>
#include <fstream>

class index {
public:
    index() = default;
    ~index();

    void read_taxonomy(const std::string &names_path, const std::string &nodes_path, const std::string &merged_path);

    void add_pair(const std::vector<uint32_t> &markers, const uint32_t &taxid);
    void reserve_space(int n_markmers);

    void save_index(const std::string &filepath);
    void load_index(const std::string &filepath);

    uint32_t predict_taxid(const std::vector<uint32_t> &markers);

    std::unordered_map<uint32_t, uint32_t> markmap;

private:
    NcbiTaxonomy * tax;
    void update(const uint32_t &singleint, const uint32_t &taxid);
    uint32_t combine(std::vector<uint32_t> const &vec);
};


#endif //MARKRAKEN_INDEX_H
