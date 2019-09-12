//
// Created by marku on 9/11/2019.
//

#ifndef MARKRAKEN_INDEX_H
#define MARKRAKEN_INDEX_H


#include <vector>
#include <iostream>
#include <unordered_map>
#include "include/NcbiTaxonomy.h"

class index {
public:
    index() = default;
    ~index();

    void read_taxonomy(const std::string &names_path, const std::string &nodes_path, const std::string &merged_path);


    void add_pair(const std::vector<uint16_t> &markers, const std::size_t &taxid);
//    uint32_t
    void save_index(const std::string &filepath);
    void reserve_space(int n_markmers);

    std::size_t predict_taxid(const std::vector<uint16_t> &markers);
private:
//    NcbiTaxonomy tax(names_path, nodes_path, merged_path);
    NcbiTaxonomy * tax;

    std::unordered_map<std::size_t, std::size_t> markmap;

    void update(const std::size_t &singleint, const std::size_t &taxid);

    std::size_t combine(std::vector<uint16_t> const &vec);
};


#endif //MARKRAKEN_INDEX_H
