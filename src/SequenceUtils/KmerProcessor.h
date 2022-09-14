//
// Created by fritsche on 24/08/22.
//

#pragma once

#include <utility>
#include <cstddef>
#include <vector>
#include "KmerLookup.h"
#include "Constants.h"
#include "GenomeLoader.h"
#include "WFA2Wrapper.h"
#include "BufferedOutput.h"
#include "SeedingStrategy.h"
//#include "fmt/core.h"

namespace protal {
    template<typename KmerLookup>
    requires KmerLookupConcept<KmerLookup>
    class KmerProcessor {
    private:
        WFA2Wrapper m_aligner = WFA2Wrapper(4, 6, 2);
        AlignmentResult m_alignment_result;

        // Kmer lookup must be thread-safe
        GenomeLoader& m_genome_loader;
//        std::ostream& m_os;

        // Each KmerLookup must be copyable and for thread safety
        // it cannot be a reference
        KmerLookup m_kmer_lookup;
        LookupList m_anchor_first;
        LookupList m_anchor_second;

        // Debug anchor
        LookupList m_anchor_third;
        LookupList m_anchor_fourth;
        LookupList m_anchor_fifth;

        AnchorPairList m_anchor_pairs;

        size_t m_idx_from_left = 0;
        size_t m_idx_from_right = 0;
        size_t m_kmer_size = 0;

//        BufferedStringOutput m_output_buffer;

        Benchmark anchor_lookup_time = Benchmark("AlignmentAnchor lookup time");
        Benchmark anchor_pair_time = Benchmark("AlignmentAnchor pair time");
        Benchmark alignment_time = Benchmark("Alignment time");
        Benchmark output_time = Benchmark("Output time");

        size_t m_num_alignments = 0;

        Seeding m_seeding{2};

        inline void SearchAnchor(LookupList &anchor_list, KmerList &kmer_list, bool from_left) {

            if (from_left) {
                for (; m_idx_from_left < kmer_list.size(); m_idx_from_left++) {
                    auto [mmer, pos] = kmer_list[m_idx_from_left];
                    m_kmer_lookup.Get(anchor_list, mmer, pos);

                    if (!anchor_list.empty()) {
                        m_idx_from_left++;
                        break;
                    }
                }
            } else {
                for (; m_idx_from_right < kmer_list.size(); m_idx_from_right++) {
                    auto [mmer, pos] = kmer_list[kmer_list.size() - m_idx_from_right - 1];

                    m_kmer_lookup.Get(anchor_list, mmer, pos);
                    if (!anchor_list.empty()) {
                        m_idx_from_right++;
                        break;
                    }
                }
            }
        }

        inline size_t LookUpAllDummy(KmerList &kmer_list, LookupList &anchor_list) {
            size_t total = 0;
            for (auto idx = 0; idx < kmer_list.size(); idx++) {
                auto [mmer, pos] = kmer_list[idx];
                    m_kmer_lookup.Get(anchor_list, mmer, pos);
            }
            return anchor_list.size();
        }

        void GetAnchorPairsFromAnchors(LookupList &reference_anchor_list, LookupList &other_anchor_list) {
            for (auto& ref_anchor : reference_anchor_list) {
                for (auto& other_anchor : other_anchor_list) {
                    if (ref_anchor.FromSameSequence(other_anchor) && !ref_anchor.Equals(other_anchor)) {
                        if (ref_anchor.genepos < other_anchor.genepos) {
                            m_anchor_pairs.emplace_back( AnchorPair(ref_anchor, other_anchor) );
                        } else {
                            m_anchor_pairs.emplace_back( AnchorPair(other_anchor, ref_anchor) );
                        }
                    }
                }
            }
        }

    public:
        size_t m_dummy = 0;
        size_t total_alignments = 0;
        KmerProcessor(KmerLookup& kmer_lookup, GenomeLoader& genome, size_t kmer_size) :
                m_kmer_lookup(kmer_lookup),
                m_genome_loader(genome),
                m_kmer_size(kmer_size) {
//            std::cout << "Regular constructor KmerProcessor" << std::endl;
        };

        KmerProcessor(KmerProcessor const& other) :
                m_kmer_lookup(other.m_kmer_lookup),
                m_genome_loader(other.m_genome_loader),
                m_kmer_size(other.m_kmer_size) {
//            std::cout << "Copy constructor KmerProcessor" << std::endl;
        };

        void Reset() {
            m_idx_from_left = 0;
            m_idx_from_right = 0;
            m_anchor_first.clear();
            m_anchor_second.clear();
            m_anchor_third.clear();
            m_anchor_fourth.clear();
            m_anchor_pairs.clear();
        }

        static bool IsReverse(LookupResult const& first_anchor, LookupResult const& second_anchor) {
            return first_anchor.readpos > second_anchor.readpos;
        }

        inline static bool IsReverse(AnchorPair const& anchor_pair) {
            return anchor_pair.first.readpos > anchor_pair.second.readpos;
        }

        inline void ReverseAnchorReadPos(LookupResult& anchor, size_t&& read_length) {
            anchor.readpos = read_length - anchor.readpos - m_kmer_size;
        }
        inline void ReverseAnchorPairReadPos(AnchorPair& anchor_pair, size_t&& read_length) {
            ReverseAnchorReadPos(anchor_pair.first, read_length);
            ReverseAnchorReadPos(anchor_pair.second, read_length);
        }

        inline void ReverseAnchorReadPos(LookupResult& anchor, size_t& read_length) {
            anchor.readpos = read_length - anchor.readpos - m_kmer_size;
        }
        inline void ReverseAnchorPairReadPos(AnchorPair& anchor_pair, size_t& read_length) {
            ReverseAnchorReadPos(anchor_pair.first, read_length);
            ReverseAnchorReadPos(anchor_pair.second, read_length);
        }

        size_t Dummy() {
            return m_dummy;
        }
//
//        inline void WriteOutput() {
//#pragma omp critical(m_os_lock)
//            {
//                m_output_buffer.Write(m_os);
//            }
//        }

        inline void PrintSummary() {
            std::cout << "#Alignments: " << m_num_alignments << std::endl;
        }

        inline void PrintBMs() {
            alignment_time.PrintResults();
            anchor_lookup_time.PrintResults();
            anchor_pair_time.PrintResults();
            output_time.PrintResults();
        }

        inline void operator () (AlignmentResultList& results, KmerList &kmer_list, FastxRecord record) {
//            Utils::Input();

            constexpr bool heavy_debug = false;
//            std::cout << "Reference to record " << &record << std::endl;
            constexpr bool more_anchor = true;
            constexpr bool naive_seeding = true;

            Reset();

            anchor_lookup_time.Start();
            SearchAnchor(m_anchor_first, kmer_list, true);
            SearchAnchor(m_anchor_second, kmer_list, false);
            SearchAnchor(m_anchor_third, kmer_list, true);
            SearchAnchor(m_anchor_fourth, kmer_list, false);
            anchor_lookup_time.Stop();


            anchor_pair_time.Start();

//            if (m_anchor_left_1.empty() || m_anchor_right_1.empty()) return;
//            std::cout << "Compute pairs" << std::endl;

            if constexpr (naive_seeding) {
                if (m_anchor_first.empty() || m_anchor_second.empty()) return;
                GetAnchorPairsFromAnchors(m_anchor_first, m_anchor_second);
            } else {
                m_seeding.Add(0, m_anchor_first.begin(), m_anchor_first.end());
                m_seeding.Add(1, m_anchor_second.begin(), m_anchor_second.end());
                m_seeding.ComputeSeedPairs();
            }
            auto& pairs = naive_seeding ? m_anchor_pairs : m_seeding.GetAnchorList();


            anchor_pair_time.Stop();

            auto& fwd = record.sequence;
            auto rev = KmerUtils::ReverseComplement(fwd);

            char* read = nullptr;
            size_t read_len = record.sequence.length();
            bool reversed = false;

            m_dummy += pairs.size();

            if constexpr(heavy_debug) {
                std::cout << "\n\n\n\n" << record.header << std::endl;
                std::cout << "\n LIST ANCHORS: " << std::endl;
                for (auto &anchor: m_anchor_pairs) {
                    std::cout << anchor.first.ToString() << " ---- " << anchor.second.ToString() << std::endl;
                }
                std::cout << "\n" << std::endl;
            }

            size_t pair_count = 0;
//           for (auto& pair : m_anchor_pairs) {
            for (auto& pair : pairs) {

                if (IsReverse(pair)) {
                    reversed = true;
                    ReverseAnchorPairReadPos(pair, read_len);
                    read = const_cast<char *>(rev.c_str());
                } else {
                    read = const_cast<char *>(record.sequence.c_str());
                }

                total_alignments++;
//                auto& anchor1 = pair.first;
//                auto& anchor2 = pair.second;

                int rpos_diff = static_cast<int>(pair.second.readpos) - static_cast<int>(pair.first.readpos);
                int gpos_diff = static_cast<int>(pair.second.genepos) - static_cast<int>(pair.first.genepos);

                int64_t anchor_distance_deviation = rpos_diff - gpos_diff;


                int abs_pos = pair.first.genepos - pair.first.readpos;

                auto& genome = m_genome_loader.GetGenome(pair.first.taxid);

                if (!genome.IsLoaded()) {
                    genome.LoadGenomeOMP();
                }

                auto &gene = genome.GetGene(pair.first.geneid);

                size_t read_start = abs_pos < 0 ? -abs_pos : 0;
                size_t gene_start = std::max(abs_pos, 0);

                size_t overlap = std::min(gene.Sequence().length() - gene_start, read_len - read_start);
                m_dummy += genome.GetKey();

                if ((gene_start + overlap) > gene.Sequence().length()) {
//#pragma omp critical(fail)
//                    {
                        std::cout << gene_start + overlap << std::endl;
                        std::cout << gene.Sequence().length() << std::endl;
                        std::cout << "pair count: " << pair_count << std::endl;
                        std::cout << "read_len:   " << read_len << std::endl;
                        std::cout << "read_start: " << read_start << std::endl;
                        std::cout << "gene_start: " << gene_start << std::endl;
                        std::cout << "gene_len:   " << gene.Sequence().length() << std::endl;
                        exit(9);
//                    }
//                    return;
                }

                assert((gene_start + overlap) <= gene.Sequence().length());



                std::string query = std::string(read + read_start, overlap);
                std::string reference = gene.Sequence().substr(gene_start, overlap);

                alignment_time.Start();
                m_aligner.Alignment(query,
                                    reference, m_alignment_result);
                m_dummy += m_aligner.GetAligner().getAlignmentScore();
                m_num_alignments++;
                alignment_time.Stop();



                double ani_proxy = WFA2Wrapper::CigarANI(m_aligner.GetAligner().getAlignmentCigar());
                if (ani_proxy > 0.93) {
                    results.emplace_back(AlignmentResult(m_aligner.GetAligner().getAlignmentScore(),
                                                         m_aligner.GetAligner().getAlignmentCigar(),
                                                         pair.first.taxid,
                                                         pair.first.geneid,
                                                         abs_pos,
                                                         !reversed));
                }

//                auto info = WFA2Wrapper::GetAlignmentInfo(m_aligner.Cigar());
//                std::cout << info.ToString() << std::endl;
//
//                if (info.alignment_ani < 1) {
//                    m_aligner.PrintAlignment();
//                    Utils::Input();
//                }


                if constexpr(heavy_debug) {
                    std::cout << "\n -> " << pair.first.ToString() << " " << pair.second.ToString() << std::endl;
                    std::cout << genome.GetKey() << std::endl;
                    std::cout << pair.first.taxid << " -> AlignmentScore: "
                              << m_aligner.GetAligner().getAlignmentScore() << std::endl;
                    m_aligner.PrintAlignment();
                    std::cout << "CigarANI: " << m_aligner.CigarANI(m_aligner.Cigar());
                }

                // TODO do not handle output here.
//                output_time.Start();
//                m_output_buffer << record.id << "\t" << std::to_string(pair.first.taxid) << "\t" << std::to_string(pair.first.geneid)<< "\t"
//                                << std::to_string(WFA2Wrapper::CigarANI(m_aligner.Cigar())) << "\n";
//                if (m_output_buffer.Full()) WriteOutput();
//
//                output_time.Stop();


                pair_count++;
            }
            if constexpr(heavy_debug) {
                Utils::Input();
            }
        }
    };
}