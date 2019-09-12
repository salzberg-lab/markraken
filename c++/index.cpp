//
// Created by marku on 9/11/2019.
//

#include "index.h"

std::size_t index::combine(std::vector<uint16_t> const &vec) {
    //https://en.wikipedia.org/wiki/Tiny_Encryption_Algorithm
    //https://github.com/HowardHinnant/hash_append/issues/7
    std::size_t seed = vec.size();
    for(auto& i : vec) {
        seed ^= i + 0x9E37 + (seed << 6) + (seed >> 2);
    }
    return seed;
}

void index::add_pair(const std::vector<uint16_t> &markers, const std::size_t &taxid) {
    const size_t singleint = combine(markers);
    update(singleint, taxid);
}

void index::reserve_space(int n_markmers) {
    markmap.reserve(n_markmers);
}

void index::update(const std::size_t &singleint, const std::size_t &taxid) {
    // insert if not already present
    // calculate LCA and insert if insert returns false
    std::pair<std::unordered_map<std::size_t, std::size_t>::iterator, bool> inserter;
    inserter = markmap.insert(std::make_pair(singleint, taxid));

    if (!inserter.second){
        uint32_t taxid_present = inserter.first->second;
        std::vector<TaxID> taxids;
        taxids.push_back(taxid);
        taxids.push_back(taxid_present);
        uint32_t ancestor = tax->LCA(taxids)->taxId;
        inserter.first->second = ancestor;
    }
}

std::size_t index::predict_taxid(const std::vector<uint16_t> &markers) {
    const size_t singleint = combine(markers);
    std::size_t t_p = markmap[singleint];
    return t_p;
}

void index::read_taxonomy(const std::string &names_path, const std::string &nodes_path, const std::string &merged_path) {
    this->tax = new NcbiTaxonomy(names_path, nodes_path, merged_path);
}

index::~index() {
    delete this->tax;
}
