//
// Created by marku on 9/11/2019.
//

#include "index.h"

uint32_t index::combine(std::vector<uint32_t> const &vec) {
    //https://en.wikipedia.org/wiki/Tiny_Encryption_Algorithm
    //https://github.com/HowardHinnant/hash_append/issues/7
    std::size_t seed = vec.size();
    for(auto& i : vec) {
        seed ^= i + 0x9E3779B9 + (seed << 6) + (seed >> 2);
    }
    return seed;
}

void index::add_pair(const std::vector<uint32_t> &markers, const uint32_t &taxid) {
    const uint32_t singleint = combine(markers);
    update(singleint, taxid);
}

void index::reserve_space(int n_markmers) {
    markmap.reserve(n_markmers);
}

void index::update(const uint32_t &singleint, const uint32_t &taxid) {
    // insert if not already present
    // calculate LCA and insert if insert returns false
    std::pair<std::unordered_map<uint32_t, uint32_t>::iterator, bool> inserter;
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

uint32_t index::predict_taxid(const std::vector<uint32_t> &markers) {
    const uint32_t singleint = combine(markers);
    uint32_t t_p = markmap[singleint];
    return t_p;
}

void index::read_taxonomy(const std::string &names_path, const std::string &nodes_path, const std::string &merged_path) {
    this->tax = new NcbiTaxonomy(names_path, nodes_path, merged_path);
}

index::~index() {
    delete this->tax;
}

void index::save_index(const std::string &filepath) {
    std::ofstream f_out(filepath, std::ios::binary);
    cereal::BinaryOutputArchive oarchive(f_out);
    oarchive(markmap);
    f_out.close();
}

void index::load_index(const std::string &filepath) {
    std::ifstream f_in(filepath, std::ios::binary);
    cereal::BinaryInputArchive iarchive(f_in);
    iarchive(markmap);
    f_in.close();
}
