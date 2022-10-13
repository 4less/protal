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

        Anchor ExtractAnchor(SeedList &seeds) {
            bool forward = seeds.front().genepos < seeds.back().genepos ? true : false;

            Anchor anchor(seeds.front().taxid, seeds.front().geneid, forward);
            for (auto& seed : seeds) {
                anchor.AddSeed(seed, m_k);
            }
            return anchor;
        }

        void FindAnchorsSingleRef(SeedList &seeds, ChainAnchorList &anchors) {
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
                    anchors.emplace_back(ExtractAnchor(m_seed_tmp_take));
                }

                seeds.clear();
                m_seed_tmp_take.clear();
                std::swap(m_seed_tmp_other, seeds);
            }
        }

        inline void FindPairs(SeedList &seeds, ChainAnchorList &anchors) {
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
                    FindAnchorsSingleRef(m_seed_tmp, anchors);
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
        Benchmark m_bm_processing{"Process Seeds"};
        Benchmark m_bm_pairing{"Pairing"};

        ChainAnchorFinder(KmerLookup& lookup) :
                m_kmer_lookup(lookup) {
        };

        ChainAnchorFinder(ChainAnchorFinder const& other) :
                m_kmer_lookup(other.m_kmer_lookup) {};

        void operator () (KmerList& kmer_list, SeedList& seeds, ChainAnchorList& anchors) {
            m_bm_seeding.Start();
            FindSeeds(kmer_list, seeds);
            m_bm_seeding.Stop();

            m_bm_processing.Start();
            Sort(seeds);
            m_bm_processing.Stop();

            m_bm_pairing.Start();
            FindPairs(seeds, anchors);
            m_bm_pairing.Stop();
        }
    };

}