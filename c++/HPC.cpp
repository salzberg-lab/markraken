//
// Object to compress homopolymers of fasta file
//

#include "HPC.h"
#include "include/NcbiTaxonomy.h"


void HPC::read(const std::string &filepath_in) {
    FastaReader fasta_reader(filepath_in);
    FastaRecord current_contig;
    while (fasta_reader.good()){
        fasta_reader.next(current_contig);
        seqs.push_back(current_contig.seq_);
        infolines.push_back(current_contig.id_);
    }
    n_seqs = seqs.size();
}

void HPC::compress() {
    seqs_compressed.reserve(seqs.size());
    for (int i=0; i<n_seqs; i++){
        int s = seqs[i].size();
        std::string seq_c;
        seq_c.reserve(s); // HPC sequence will always be less than or equal to size of original seq
        homopolymer_compress(seqs[i], seq_c);
        seqs_compressed.push_back(seq_c);
    }
}

void HPC::homopolymer_compress(const std::string &str, std::string &compressed) {
    int i = str.size();
    for (int j=0; j<i; ++j){
        while (str[j] == str[j+1]){
            j++;
        }
        compressed.push_back(str[j]);
    }
}

void HPC::write(const std::string &filepath_out) {
    FastaWriter fasta_writer(filepath_out);
    FastaRecord current_contig;
    for (int i=0; i<n_seqs; i++){
        current_contig.seq_ = seqs_compressed[i];
        current_contig.id_ = infolines[i];
        fasta_writer.write(current_contig);
    }
}

const std::vector<std::string> explode(const std::string& s, const char& c)
//http://www.cplusplus.com/articles/2wA0RXSz/
{
    std::string buff{""};
    std::vector<std::string> v;

    for(auto n:s)
    {
        if(n != c) buff+=n; else
        if(n == c && buff != "") { v.push_back(buff); buff = ""; }
    }
    if(buff != "") v.push_back(buff);

    return v;
}

void HPC::extract_miniseq_taxids() {
    const char delimiter = '|';
    for (int i=0; i<n_seqs; i++){
        std::string &line = infolines[i];
        std::vector<std::string> split = explode(line, delimiter);
        std::string &id = split[2];
        if (id.at(0) == 'x'){
            id.erase(0, 1);
        }

        TaxID id_num = std::stoi(id);
        taxids.push_back(id_num);
    }
}








