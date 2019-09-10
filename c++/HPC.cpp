//
// Object to compress homopolymers of fasta file
//

#include "HPC.h"

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








