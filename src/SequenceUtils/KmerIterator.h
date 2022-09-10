//
// Created by fritsche on 14/08/22.
//

#pragma once

#include <cstring>
#include <cassert>
#include "FastxReader.h"
#include "KmerUtils.h"
#include "Minimizer.h"
#include "Constants.h"

namespace protal {



    // make templated
//    template<typename ContextFreeMinimizer=Syncmer15>
    class CanonicalKmerIterator {
    private:
//        FastxRecord m_record;
        uint32_t m_key_bytes;
        uint32_t m_last_pos_shift;
        uint32_t m_k;
        size_t m_mask;
        uint32_t m_pos = 0;
//        std::string m_sequence;
        char* m_seq = nullptr;
        size_t m_seq_length;
        bool m_first = true;

        size_t m_canonical_kmer_fwd = 0;
        size_t m_canonical_kmer_rev = 0;

        void GetKey(size_t& key) {
            for (int i = 0; i < m_k; i++) {
                key |= KmerUtils::BaseToInt(m_seq[i], 0) << (2 * (m_k - i - 1));
            }
        }

        void GetKeyC(size_t& key) {
            for (int i = 0; i < m_k; i++) {
                key |= KmerUtils::BaseToIntC(m_seq[i], 0) << (2llu * i);
            }
        }

        void RollKey(size_t &key) {
            key <<= 2;
            key &= m_mask;
            key |= KmerUtils::BaseToInt(m_seq[m_pos + m_k - 1], 0);
        }

        void RollKeyC(size_t &key) {
            key >>= 2;
            key |= KmerUtils::BaseToIntC(m_seq[m_pos + m_k - 1], 0) << 2llu * (m_k - 1);
        }

    public:
        int64_t GetPos() const {
            return m_pos - 1;
        }

        void Reset() {
            m_pos = 0;
            m_first = true;
            m_canonical_kmer_fwd = 0;
            m_canonical_kmer_rev = 0;
        }
//
//        void SetRecord(FastxRecord& record) {
//            m_record = record;
//            m_sequence = record.sequence;
//            Reset();
//        }

        // TODO this is not ideal -  should not have to convert const char* to char*
        void SetSequence(char* seq, size_t len) {
            m_seq = seq;
            m_seq_length = len;
            Reset();
        }
//
//        void SetSequence(std::string const& sequence) {
//            m_sequence = sequence;
//            Reset();
//        }

        bool HasNext() {
            return (m_pos < m_seq_length - m_k + 1);
        }


        void Init(size_t k) {
            m_k = k;
            m_key_bytes = (2 * m_k + 8 - 1) / 8;
            m_last_pos_shift = (8 - ((2 * m_k) % 8)) % 8;

            m_mask = (1llu << (m_k*2)) -1;
        }

        void Init(CanonicalKmerIterator const& other) {
            m_k = other.m_k;
            Init(m_k);
        }

        std::string ToString() const {
            std::string str = "";
            str += "k: " + std::to_string(m_k);
            return str;
        }

        CanonicalKmerIterator() {};

        CanonicalKmerIterator(size_t k) : m_k(k) {
            Init(m_k);
        }


        ~CanonicalKmerIterator() {}

        inline void KmerList(KmerList& list) {
            size_t kmer;
            list.clear();
            while ((*this)(kmer)) {
                list.emplace_back( KmerElement(kmer, GetPos()) );
            }
        }

        inline bool operator () (size_t& key) {
            if (!HasNext()) return false;
            if (m_first) {
                GetKey(m_canonical_kmer_fwd);
                GetKeyC(m_canonical_kmer_rev);
                m_first = false;
            } else {
                RollKey(m_canonical_kmer_fwd);
                RollKeyC(m_canonical_kmer_rev);
            }
            bool take_normal = m_canonical_kmer_fwd < m_canonical_kmer_rev;
            key = (take_normal * m_canonical_kmer_fwd) + (!take_normal * m_canonical_kmer_rev);

            m_pos++;
            return true;
        }
    };

    template<typename CFMinimizer>
    requires ContextFreeMinimizer<CFMinimizer>
    class SimpleKmerHandler {
    private:
        CFMinimizer m_minimizer{};

        uint32_t m_k{};
        size_t m_kmer{};
        size_t m_mask{};
        uint32_t m_pos = 0;

//        char* m_seq = nullptr;
        std::string_view m_seq{};
        bool m_first = true;

        size_t m_total_kmers = 0;
        size_t m_total_minimizers = 0;

        size_t m_canonical_kmer_fwd = 0;
        size_t m_canonical_kmer_rev = 0;

        void GetKey(size_t& key) {
            for (int i = 0; i < m_k; i++) {
                key |= KmerUtils::BaseToInt(m_seq[i], 0) << (2 * (m_k - i - 1));
            }
        }

        void GetKeyC(size_t& key) {
            for (int i = 0; i < m_k; i++) {
                key |= KmerUtils::BaseToIntC(m_seq[i], 0) << (2llu * i);
            }
        }

        void RollKey(size_t &key) {
            key <<= 2;
            key &= m_mask;
            key |= KmerUtils::BaseToInt(m_seq[m_pos + m_k - 1], 0);
        }

        void RollKeyC(size_t &key) {
            key >>= 2;
            key |= KmerUtils::BaseToIntC(m_seq[m_pos + m_k - 1], 0) << 2llu * (m_k - 1);
        }

    public:
        /**
         * Implementing KmerStatisticsConcept
         * @returns Total Kmers processed
         */
        size_t TotalKmers() {
            return m_total_kmers;
        }

        /**
         * Implementing MinimizerStatisticsConcept
         * @returns Total Kmers processed
         */
        size_t TotalMinimizers() {
            return m_total_minimizers;
        }

        void Init(size_t k) {
            m_k = k;
            m_mask = (1llu << (m_k*2)) -1;
        }

        SimpleKmerHandler(size_t k, CFMinimizer& minimizer) : m_k(k), m_minimizer(minimizer) {
            Init(m_k);
        }
        SimpleKmerHandler(SimpleKmerHandler const& other) :
                m_minimizer(other.m_minimizer) {
            Init(other.m_k);
        }
        Syncmer minimizer{15, 7, 2};


        int64_t GetPos() const {
            return m_pos - 1;
        }

        void Reset() {
            m_pos = 0;
            m_first = true;
            m_canonical_kmer_fwd = 0;
            m_canonical_kmer_rev = 0;

            m_total_kmers = 0;
            m_total_minimizers = 0;
        }

        void SetSequence(std::string_view const& sequence) {
            m_seq = sequence;
            Reset();
        }
        void SetSequence(std::string_view const&& sequence) {
            m_seq = sequence;
            Reset();
        }

        bool HasNext() {
            return (m_pos < m_seq.length() - m_k + 1);
        }


        [[nodiscard]] std::string ToString() const {
            std::string str = "";
            str += "k: " + std::to_string(m_k);
            return str;
        }

        inline void operator () (std::string_view const& sequence, KmerList& list) {
            list.clear();
            SetSequence(sequence);
            while (NextWrapper(m_kmer)) {
                if (m_minimizer(m_kmer)) {
                    m_total_minimizers++;
                    list.emplace_back(KmerElement(m_kmer, GetPos()));
                }
            }
        }
        inline void operator () (std::string_view const&& sequence, KmerList& list) {
            list.clear();
            SetSequence(sequence);
            while (NextWrapper(m_kmer)) {
                if (m_minimizer(m_kmer)) {
                    m_total_minimizers++;
                    list.emplace_back(KmerElement(m_kmer, GetPos()));
                }
            }
        }

        inline bool Next(size_t& key) {
            while (NextWrapper(m_kmer)) {
                if (m_minimizer(m_kmer)) {
                    key = m_kmer;
                    return true;
                }
            }
            return false;
        }

        inline bool NextWrapper(size_t& key) {
            if (!HasNext()) return false;
            m_total_kmers++;
            if (m_first) {
                GetKey(m_canonical_kmer_fwd);
                GetKeyC(m_canonical_kmer_rev);
                m_first = false;
            } else {
                RollKey(m_canonical_kmer_fwd);
                RollKeyC(m_canonical_kmer_rev);
            }
            bool take_normal = m_canonical_kmer_fwd < m_canonical_kmer_rev;
            key = (take_normal * m_canonical_kmer_fwd) + (!take_normal * m_canonical_kmer_rev);

            m_pos++;
            return true;
        }
    };
}