//
// Created by marku on 9/11/2019.
//

#ifndef MARKRAKEN_HASHER_H
#define MARKRAKEN_HASHER_H


#include <cstdint>
#include <cstdlib>

class hasher {
public:
    hasher() = default;
    ~hasher() = default;
    uint32_t murmur3_32(const uint8_t* key, size_t len, uint32_t seed);
};


#endif //MARKRAKEN_HASHER_H
