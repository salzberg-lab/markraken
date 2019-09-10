#include <iostream>
using namespace std;


string homopolymer_compress(const string & str){
    int i = str.size();
    string letters;
    letters.reserve(i); // HPC sequence will always be less than or equal to size of original seq
    for (int j=0; j<i; ++j){
        while (str[j] == str[j+1]){
            j++;
        }
        letters.push_back(str[j]);
    }
    return letters;
}

int main() {
//    std::cout << "Hello, World!" << std::endl;
    string foo = "AAATTCGGGAAAA";
    string HPC = homopolymer_compress(foo);
//    std::cout << i
    return 0;
}