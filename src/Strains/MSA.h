//
// Created by fritsche on 04/09/23.
//

#pragma once

#include "Strain.h"
#include <format>

namespace protal {
    class VariantCaller {
        size_t m_min_var_qual_sum = 0;
        size_t m_min_var_cov = 0;
    public:
        inline Variant &Call(VariantBin &variants) const {
            std::sort(variants.begin(), variants.end(), [](Variant const &a, Variant const &b) {
                return a.QualitySum() > b.QualitySum();
            });
            return variants.front();
        }

        inline bool Pass(Variant &call, VariantBin &variants) const {
            return !(call.QualitySum() < m_min_var_qual_sum || call.Observations() < m_min_var_cov);
        }
    };

    class IUPACConverter {
        uint32_t m_key;

        static constexpr const char int_to_iupac[]{
                'N', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'
        };

    public:
        void Reset() {
            m_key = 0;
        }

        char Get() {
            return int_to_iupac[m_key];
        }

        void Add(char base) {
            switch (base) {
                case 'A':
                    m_key |= 1;
                    break;
                case 'C':
                    m_key |= 2;
                    break;
                case 'G':
                    m_key |= 4;
                    break;
                case 'T':
                    m_key |= 8;
                    break;
                default:
                    break;
            }
        }
    };

    class BaseCaller {
        static const char NO_COVERAGE = '-';
        static const char NO_BASE_PASSES = 'N';
        IUPACConverter m_iupac;

    public:
        inline char operator()(VariantBin &variants, uint16_t const &coverage, VariantCaller const &caller) {
            if (coverage == 0) {
                return NO_COVERAGE;
            }

            auto &call = caller.Call(variants);
            auto pass = caller.Pass(call, variants);

            uint32_t iupac_code;
            m_iupac.Reset();
            for (auto &var: variants) {
                if (caller.Pass(var, variants)) m_iupac.Add(var.GetVariant());
            }
            return m_iupac.Get();
        }
    };

    class MSABuilder {
        struct PartitionSequence {
            std::string name = "";
            std::string sequence = "";

        };

        struct PartitionEntry {
            std::string name = "";
            size_t start = 0;
            size_t end = 0;

            std::string ToString() {
                //DNA, gene1 = 0-395
                return "DNA, " + name + " = " + std::to_string(start) + "-" + std::to_string(end);
            }
        };

    private:
//        BaseCaller m_bcaller;
        VariantCaller m_vcaller;

        std::vector<PartitionSequence> m_sequences;
        std::vector<PartitionEntry> m_partitions;

    public:
        MSABuilder(VariantCaller& variant_caller) : m_vcaller(variant_caller) {
        }

        void Clear() {
            m_sequences.clear();
            m_partitions.clear();
        }

        void AddPartition(std::string sequence_name, std::string partition_name, VariantHandler& variant_handler, SequenceRangeHandlerRef sequence_range_handler) {
            auto& varbins = variant_handler.GetVariants();
            auto covs = sequence_range_handler.get().CalculateCoverageVector2();
            std::vector<Variant> variants(10, Variant());


        }
    };

//    std::string Consensus(std::string reference, std::vector<>) {
//
//    }
}
