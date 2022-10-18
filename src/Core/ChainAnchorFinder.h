//
// Created by fritsche on 12/10/22.
//

#pragma once
#include "ChainingStrategy.h"
#include "Constants.h"
#include "Benchmark.h"


namespace protal {

    using SeedList = LookupList;
    using Anchor = ChainAlignmentAnchor;
    using ChainAnchorList = std::vector<Anchor>;

    template<typename KmerLookup>
    requires KmerLookupConcept<KmerLookup>
    class ChainAnchorFinder {
        KmerLookup m_kmer_lookup;

        SeedList m_all_seeds;
        SeedList m_seed_tmp;
        SeedList m_seed_tmp_take;
        SeedList m_seed_tmp_other;

        size_t m_k = 15;

        inline void FindSeeds(KmerList &kmer_list, SeedList& seeds) {
            size_t prev_size = 0;
            for (auto pair : kmer_list) {
                auto [mmer, pos] = pair;
                m_kmer_lookup.Get(seeds, mmer, pos);
            }
        }


        using OffsetPair = std::pair<int, int>;


        static bool OffsetPairMatch(OffsetPair const& a, OffsetPair const& b) {
            return a.first == b.first || a.second == b.second;
        }

        static std::string OffsetPairToString(OffsetPair const& a) {
            return "[" + std::to_string(a.first) + "," + std::to_string(a.second) + "]";
        }

        static OffsetPair Offset(Seed const& s) {
            return OffsetPair { static_cast<int>(s.genepos) + s.readpos, static_cast<int>(s.genepos) - s.readpos };
        }

        static bool IsReverse(Anchor& a) {
            return a.Front().readpos > a.Back().readpos;
        }

        static inline void ReverseSeed(Seed& seed, size_t& read_length, size_t k) {
            seed.readpos = read_length - seed.readpos - k;
        }
        inline void ReverseSeedList(SeedList &seeds, size_t& read_length) {
            for (auto& seed : seeds) {
                ReverseSeed(seed, read_length, m_k);
            }
            std::reverse(seeds.begin(), seeds.end());
        }

        Anchor ExtractAnchor(SeedList &seeds, size_t read_length) {
            bool forward = seeds.front().genepos < seeds.back().genepos;

//            std::cout << "Init " << std::string(50, '-') << "fwd " << forward << std::endl;
//            for (auto &seed: seeds) {
//                std::cout << seed.ToString() << std::endl;
//            }

            if (!forward) {
                ReverseSeedList(seeds, read_length);
//                std::cout << "Reverse " << std::string(50, '-') << std::endl;
//                for (auto &seed: seeds) {
//                    std::cout << seed.ToString() << std::endl;
//                }
            }

            Anchor anchor(seeds.front().taxid, seeds.front().geneid, forward);
            bool stop = false;
            bool skip_first = true;
            for (auto& seed : seeds) {
                if (!skip_first && (seeds.front().genepos >= seed.genepos)) {
//                    std::cout << "Skip " << seed.ToString() << " " << (seeds.front().genepos < seed.genepos != forward) << "  detected: " << forward << std::endl;
                    stop = true;
                    skip_first = false;
                    continue;
                }
                skip_first = false;
                anchor.AddSeed(seed, m_k);
            }

            /////////////////////////////////////////////
//            if (stop) {
//                std::cout << "Debug " << std::string(50, '-') << std::endl;
//                for (auto &seed: seeds) {
//                    std::cout << seed.ToString() << std::endl;
//                }
//                std::cout << " -> " << anchor.ToString() << std::endl;
//                Utils::Input();
//            }
            /////////////////////////////////////////////

//            std::cout << " -> " << anchor.ToString() << std::endl;
            return anchor;
        }

        void FindAnchorsSingleRef(SeedList &seeds, ChainAnchorList &anchors, size_t read_length) {
            while (seeds.size() > 1) {
                auto& init_seed = seeds[0];
                auto init_offset = Offset(init_seed);

                for (auto i = 1; i < seeds.size(); i++) {
                    auto& seed = seeds[i];
                    auto seed_offset = Offset(seed);
                    if (OffsetPairMatch(init_offset, seed_offset)) {
                        m_seed_tmp_take.emplace_back(seed);
                    } else {
                        m_seed_tmp_other.emplace_back(seed);
                    }
                }
                if (m_seed_tmp_take.size() > 1) {
                    anchors.emplace_back(ExtractAnchor(m_seed_tmp_take, read_length));
                }

                seeds.clear();
                m_seed_tmp_take.clear();
                std::swap(m_seed_tmp_other, seeds);
            }
        }

        inline void FindPairs(SeedList &seeds, ChainAnchorList &anchors, size_t read_length) {
            if (seeds.empty()) return;

            for (auto i = 0; i < seeds.size()-1; i++) {
                auto& seed = seeds[i];
                auto& next_seed = seeds[i+1];

                if (seed == next_seed) {
                    size_t group_size = 2;
                    while (i+group_size < seeds.size() && seeds[i+group_size] == seed) group_size++;

                    auto start_it = seeds.begin() + i;
                    auto end_it = start_it + group_size;
                    m_seed_tmp.clear();
                    m_seed_tmp.insert(m_seed_tmp.end(),
                                      std::make_move_iterator(start_it),
                                      std::make_move_iterator(end_it));
                    FindAnchorsSingleRef(m_seed_tmp, anchors, read_length);
                    i += group_size - 1;
                }
            }
        }

        void Sort(SeedList &list) {
            std::sort(list.begin(), list.end(), Seed::SortByReadComparator);
        }

        void Reset() {
            // Assess what time seeding takes.
            m_all_seeds.clear();
        }

    public:
        size_t dummy = 0;
        Benchmark m_bm_seeding{"Seeding"};
        Benchmark m_bm_processing{"Sorting Seeds"};
        Benchmark m_bm_pairing{"Pairing"};
        Benchmark m_bm_sorting_anchors{"Sorting Anchors"};

        ChainAnchorFinder(KmerLookup& lookup) :
                m_kmer_lookup(lookup) {
        };

        ChainAnchorFinder(ChainAnchorFinder const& other) :
                m_kmer_lookup(other.m_kmer_lookup) {};

        void operator () (KmerList& kmer_list, SeedList& seeds, ChainAnchorList& anchors, size_t read_length) {
            m_bm_seeding.Start();
            FindSeeds(kmer_list, seeds);
            m_bm_seeding.Stop();

            m_bm_processing.Start();
            Sort(seeds);
            m_bm_processing.Stop();

            m_bm_pairing.Start();
            FindPairs(seeds, anchors, read_length);
            m_bm_pairing.Stop();

            m_bm_sorting_anchors.Start();
            std::sort(anchors.begin(), anchors.end() ,
                      [](Anchor const& a, Anchor const& b) {
                return a.total_length > b.total_length;
            });
            m_bm_sorting_anchors.Stop();
        }
    };

}