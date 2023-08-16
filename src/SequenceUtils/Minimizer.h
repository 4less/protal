//
// Created by fritsche on 17/02/2022.
//

#pragma once

#include <cstdint>
#include <cstddef>
#include <xxhash64.h>
#include <iostream>

namespace protal {

    template <typename T>
    concept IsConceptFreeMinimizer =
    requires(T t, size_t s) {
        { t.IsMinimizer(s) } -> std::same_as<bool>;
    };


    class None {};

    class Minimizer {
    public:
        virtual bool IsMinimizer (uint64_t& key) = 0;
    };

    class Syncmer {
    private:
        const uint8_t m_k;
        const uint8_t m_kbits;
        const uint8_t m_s;
        const uint8_t m_sbits;
        const uint8_t m_t;

        const uint32_t m_shift_size;
        const uint8_t m_smask;

        uint32_t m_smers[16] { 0, 0, 0, 0, 0, 0, 0, 0 };
        uint16_t m_shifts[16] { 0, 0, 0, 0, 0, 0, 0, 0 };

    public:
        Syncmer(uint8_t k, uint8_t s, uint8_t t) :
                m_k(k), m_kbits(k*2), m_s(s), m_sbits(s*2), m_smask((1u << (2*s)) - 1), m_shift_size(k - s + 1), m_t(t) {

            double shift_step = (((double) m_kbits - m_sbits + 1) / 16);

//            std::cout << "Syncmer submer: " << m_shift_size << std::endl;
            for (int i = 0; i < m_shift_size; i++) {
                m_shifts[i] = (m_shift_size*2) - ((i+1)*2);
            }
        }

        static uint8_t mindex(const uint32_t *ary, size_t len) {
            uint32_t min = 255, index = 0;
            for (int i = 0; i < len; i++) {
                if (ary[i] < min) {
                    min = ary[i];
                    index = i;
                }
            }
            return index;
        }

        inline bool IsMinimizer(uint64_t& key) {
            size_t min_index = 0;
            size_t min = UINT64_MAX;
#pragma omp simd
            for (int i = 0; i < m_shift_size; i++) {
                m_smers[i] = (key >> m_shifts[i]) & m_smask;
                min_index = (m_smers[i] < min) * i + (m_smers[i] >= min) * min_index;
                min =  (m_smers[i] < min) * m_smers[i] + (m_smers[i] >= min) * min;
            }
            return min_index == m_t || min_index == (m_shift_size - 1 - m_t);
        }
    };


    // Adheres to ContextFreeMinimalizer
    class ClosedSyncmer {
    private:
        const uint8_t m_k;
        const uint8_t m_kbits;
        const uint8_t m_s;
        const uint8_t m_sbits;
        const uint8_t m_t;

        const uint32_t m_shift_size;
        const uint8_t m_smask;

        uint32_t m_smers[16] { 0, 0, 0, 0, 0, 0, 0, 0 };
        uint16_t m_shifts[16] { 0, 0, 0, 0, 0, 0, 0, 0 };

    public:
        ClosedSyncmer(uint8_t k, uint8_t s, uint8_t t) :
                m_k(k), m_kbits(k*2), m_s(s), m_sbits(s*2), m_smask((1u << (2*s)) - 1), m_shift_size(k - s + 1), m_t(t) {
//            std::cout << "ClosedSyncmer" << std::endl;
            for (int i = 0; i < m_shift_size; i++) {
                m_shifts[i] = (m_shift_size*2) - ((i+1)*2);
            }
        }

        inline bool operator () (uint64_t& key) {
            size_t min_index = 0;
            size_t min = UINT64_MAX;
#pragma omp simd
            for (int i = 0; i < m_shift_size; i++) {
                m_smers[i] = (key >> m_shifts[i]) & m_smask;
                min_index = (m_smers[i] < min) * i + (m_smers[i] >= min) * min_index;
                min =  (m_smers[i] < min) * m_smers[i] + (m_smers[i] >= min) * min;
            }
            return min_index == m_t || min_index == (m_shift_size - 1 - m_t);
        }
    };




    class UniversalMin : Minimizer  {
    private:
        double m_threshold = 0;
        XXHash64 hashgen{ 0 };

    public:
        UniversalMin(double threshold) : m_threshold(threshold) {};

        static inline uint64_t MurmurHash3(uint64_t k) {
            k ^= k >> 33;
            k *= 0xff51afd7ed558ccd;
            k ^= k >> 33;
            k *= 0xc4ceb9fe1a85ec53;
            return k;
        }

        static uint64_t inline Hash2(uint64_t key) {
            return XXHash64::hash(&key, 8, 0);
        }

        inline bool IsMinimizer(uint64_t& key) {

            return (double)Hash2(key)/UINT64_MAX < m_threshold;
        }

        static bool IsMinimizer(uint64_t& key, double threshold) {
            return (double)Hash2(key)/UINT64_MAX < threshold;
        }

        static bool IsMinimizer2(uint64_t& key) {
            const uint64_t THRESHOLD = 0.15 * UINT64_MAX;
            return Hash2(key) < THRESHOLD;
        }

        static bool IsMinimizerMM(uint64_t& key, double threshold) {
            return (double)MurmurHash3(key)/UINT64_MAX < threshold;
        }

        static bool IsMinimizerMM2(uint64_t& key) {
            static const uint64_t THRESHOLD = 0.15 * UINT64_MAX;
            return MurmurHash3(key) < THRESHOLD;
//            return (double)MurmurHash3(key)/UINT64_MAX < threshold;
        }
    };


}