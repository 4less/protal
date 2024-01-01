//
// Created by fritsche on 31/10/22.
//

#pragma once

#include "wfa2-lib/bindings/cpp/WFAligner.hpp"
#include <immintrin.h>

using namespace::wfa;

namespace test {
    static void TestEndsFree() {
        // ACTAGCTAGCGACTTACTTTCTGCTCTTCTCTGCTGCTGCTGATGCTAGCTGCTAGCTGATCGATCGTAGCTATC
//        std::string query = "ACTAGCTAGCGACTTACTTTCTGCTCTTCTCTGCTGCTGCTGATGCTA";
//        std::string reference = "CTTACTATCTGCTCTTCTCTGCTGCTGCTGATGCTAGCTGCTAGCTGATCGATCGTAGCTATC";
        std::string query =     "ACTAGCTAGCGACTTACTTTCTGCTCTTCTCTG";     // PATTERN
        std::string reference =    "AGCTAGCGACTTACTTTCTGCTCTTCTCTGCTG";  // TEXT
//        std::string query =     "AATTAATTTAAGTCTAGGCTACTTTCGGTACTTTGTTCTT";
//        std::string reference =    "AATTTAAGTCTAGGCTACTTTCGGTACTTTCTT";


        std::cout << "QUERY:  " << query << std::endl;
        std::cout << "REFERENCE: " << reference << std::endl;
        WFAlignerGapAffine aligner{ 4, 6, 2, WFAligner::Alignment, WFAligner::MemoryHigh };
        aligner.setVerbose(2);



        aligner.alignEndsFree(query, 3, 0, reference, 0, 3);


        std::cout << "Status: " << aligner.getAlignmentStatus() << std::endl;
        std::cout << "Score: " << aligner.getAlignmentScore() << std::endl;
        aligner.cigarPrintPretty(stdout, query, reference);
    }

    static void RegularSim(uint32_t query, uint32_t* flex_start, size_t flex_length) {
        std::cout << "RegularSim" << std::endl;
        for (auto i = 0; i < flex_length; i += 8) {
            for (auto j = 0; j < 8; j++) {
                std::cout << (query & flex_start[i+j]) << " ";
            }
            std::cout << std::endl;
        }
    }

    template<class T> inline void Log(const __m256i & value)
    {
        const size_t n = sizeof(__m256i) / sizeof(T);
        T buffer[n];
        _mm256_storeu_si256((__m256i*)buffer, value);
        std::cout << "SIMD [";
        for (int i = 0; i < n; i++)
            std::cout << buffer[i] << " ";
        std::cout << "]" << std::endl;
    }

    inline __m256i popcount_pshufb32(const __m256i v) {

        __m256i lookup = _mm256_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2,
                       2, 3, 2, 3, 3, 4, 0, 1, 1, 2, 1, 2, 2, 3,
                       1, 2, 2, 3, 2, 3, 3, 4);
        __m256i low_mask = _mm256_set1_epi8(0x0f);
        __m256i lo = _mm256_and_si256(v, low_mask);
        __m256i hi = _mm256_and_si256(_mm256_srli_epi16(v, 4), low_mask);
        __m256i popcnt1 = _mm256_shuffle_epi8(lookup, lo);
        __m256i popcnt2 = _mm256_shuffle_epi8(lookup, hi);
        __m256i sum8 = _mm256_add_epi8(popcnt1, popcnt2);
        return _mm256_srli_epi32(
            _mm256_mullo_epi32(sum8, _mm256_set1_epi32(0x01010101)), 24);
        // vpmulld is slowish (2 uops) on most recent Intel CPUs
        // but still single-uop on AMD
    }


    static void SIMDSim(uint32_t query, uint32_t* flex_start, size_t flex_length) {

        constexpr uint32_t mask1 = 0b01010101010101010101010101010101;
        constexpr uint32_t mask2 = 0b10101010101010101010101010101010;
        constexpr uint32_t invert = 0b11111111111111111111111111111111;
        const __m256i MASK1_256_VEC = _mm256_set1_epi32(static_cast<int>(mask1));
        const __m256i MASK2_256_VEC = _mm256_set1_epi32(static_cast<int>(mask2));
        const __m256i INVERT = _mm256_set1_epi32(static_cast<int>(invert));

        const __m256i QUERY_256_VEC = _mm256_set1_epi32(static_cast<int>(query));
        test::Log<uint32_t>(QUERY_256_VEC);

        uint32_t* sequence_header;
        size_t sequence_header_length;
        // Store sequence and similarity
        __m256i SEQ_256_VEC, SIM_256_VEC, M1_256_VEC, M2_256_VEC, XOR_256_VEC, XOR_256_VEC2, POPCOUNT;
        // Initialize simd query vector
        // First do similarity then do popcount
        for (auto i = 0; i < flex_length; i += 8) {
            SEQ_256_VEC = _mm256_loadu_si256((__m256i *) (flex_start+i));
            XOR_256_VEC = _mm256_xor_si256(_mm256_xor_si256(SEQ_256_VEC, QUERY_256_VEC),  INVERT);
            M1_256_VEC = _mm256_and_si256(XOR_256_VEC, MASK1_256_VEC);
            M2_256_VEC = _mm256_srli_epi32(_mm256_and_si256(XOR_256_VEC, MASK2_256_VEC), 1);
            SIM_256_VEC = _mm256_and_si256(M1_256_VEC, M2_256_VEC);
            POPCOUNT  = popcount_pshufb32(SIM_256_VEC);
        }
        test::Log<uint32_t>(SEQ_256_VEC);
        test::Log<uint32_t>(SIM_256_VEC);
        test::Log<uint32_t>(POPCOUNT);
    }

    static void SIMDSimilarity() {
        const uint32_t query = 0b11110010011111111111100111111111;
        std::vector<uint32_t> flex {
            0b11001111111111111111001111111111,
            0b11110010011111111111100111111111,
            0b11111111111111111111111111111111,
            0b11111111111000111111111111111111,
            0b11111111111111111111111111111111,
            0b11111111111111111111111111111111,
            0b11111111111100001111111111111111,
            0b11110001111111110001111111111111
        };
        test::SIMDSim(query, flex.data(), flex.size());
        test::RegularSim(query, flex.data(), flex.size());
    };
}
