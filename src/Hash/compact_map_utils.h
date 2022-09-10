//
// Created by fritsche on 10/03/2021.
//

#pragma once

#include <cstdint>
#include <sys/stat.h>

namespace CompactMapUtils {
    static inline uint32_t reduce(uint32_t x, uint32_t N) {
        return x % N;
//        return ((uint64_t) x * (uint64_t) N) >> 32 ;
    }


    static inline uint32_t MurmurHash32(uint32_t k) {
        k *= 0xff51afd7;
        k ^= k >> 17;
        k *= 0xc4ceb9fe;
        k ^= k >> 17;
        return k;
    }

    static inline uint32_t Hash(uint64_t key, uint32_t capacity) {
        return reduce(MurmurHash32(key), capacity);
    }

    static inline uint64_t MurmurHash3(uint64_t k) {
        k ^= k >> 33;
        k *= 0xff51afd7ed558ccd;
        k ^= k >> 33;
        k *= 0xc4ceb9fe1a85ec53;
//        k ^= k >> 33;
        return k;
    }

    static inline uint64_t Hash(uint64_t key, size_t value_bits, size_t max) {
        auto hash = MurmurHash3(key);
        if (hash > max)
            return hash >> value_bits;
        return hash;
    }
}