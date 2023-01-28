//
// Created by fritsche on 12/10/22.
//

#pragma once
#include "ChainingStrategy.h"
#include "Constants.h"
#include "Benchmark.h"
#include "robin_map.h"


namespace protal {

    using SeedList = LookupList;
    using Anchor = ChainAlignmentAnchor;
    using ChainAnchorList = std::vector<Anchor>;
    using TaxaCounts = tsl::robin_map<uint32_t, uint16_t>;
    using TaxonCountPair = std::pair<uint32_t, uint16_t>;
    using TaxonCountPairs = std::vector<TaxonCountPair>;
    using TaxonSet = tsl::robin_set<uint32_t>;
    using LookupResultList = std::vector<LookupPointer>;

    template<typename KmerLookup>
    requires KmerLookupConcept<KmerLookup>
    class ChainAnchorFinder {
        KmerLookup m_kmer_lookup;
        LookupResultList m_lookups;

        SeedList m_all_seeds;
        SeedList m_seed_tmp;
        SeedList m_seed_tmp_take;
        SeedList m_seed_tmp_other;

        size_t m_k = 15;// Update and put hat in as parameter

        TaxonCountPairs m_count_pairs;
        TaxaCounts m_taxa_counts;
        TaxonSet m_selected_taxa;
        SeedList m_subset;

        size_t m_min_successful_lookups = 4;
        size_t m_max_seed_size = 10;

        using RecoverySet = tsl::robin_set<std::pair<uint32_t, uint32_t>>;

        inline void FindSeeds(KmerList &kmer_list, SeedList& seeds) {
            m_lookups.clear();
            for (auto [mmer, pos] : kmer_list) {
                m_kmer_lookup.Get(m_lookups, mmer, pos);
            }
            std::sort(m_lookups.begin(), m_lookups.end(), [](LookupPointer const& a, LookupPointer const& b) {
                return a.size < b.size;
            });

            uint32_t previous_size = 0;
            uint32_t successful_lookups = 0;
            uint32_t total_lookups = 0;
            for (auto i = 0; i < m_lookups.size(); i++) {
                auto& lookup = m_lookups[i];
                m_kmer_lookup.GetFromLookup(seeds, lookup);
                total_lookups++;
                successful_lookups += (seeds.size() > previous_size);
                if (seeds.size() > m_max_seed_size && successful_lookups > m_min_successful_lookups) {
                    break;
                }
                previous_size = seeds.size();
            }
//            std::cout << successful_lookups << "/" << total_lookups << ": " << seeds.size() << std::endl;
        }

        void RecoverAnchors(ChainAnchorList& anchors, bool expect_fwd, RecoverySet& recover) {

        }


        inline void FindSeeds2(KmerList &kmer_list, SeedList& seeds) {
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

        static inline void ReverseSeed(Seed& seed, size_t& read_length, size_t k, size_t subk = 16) {
//            std::cout << "Reverse: " << read_length << ", " << seed.readpos << ", " << k << " = " << (read_length - seed.readpos - k) << std::endl;
            seed.readpos = read_length - seed.readpos - k;
//            seed.readpos = read_length - seed.readpos + 1;
        }
        inline void ReverseSeedList(SeedList &seeds, size_t& read_length) {
            for (auto& seed : seeds) {
                ReverseSeed(seed, read_length, m_k);
            }
            std::reverse(seeds.begin(), seeds.end());
        }

        Anchor ExtractAnchor(SeedList &seeds, size_t read_length) {
            bool forward = seeds.front().genepos < seeds.back().genepos;
            if (!forward) {
                ReverseSeedList(seeds, read_length);
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

            return anchor;
        }

        void FindAnchorsSingleRef(SeedList &seeds, ChainAnchorList &anchors, size_t read_length) {
//            std::cout << "Find Anchors Single Ref" << std::endl;
//            for (auto& seed : seeds) {
//                std::cout << seed.ToString() << std::endl;
//            }
            while (seeds.size() > 1) {
                auto& init_seed = seeds[0];
                auto init_offset = Offset(init_seed);

                m_seed_tmp_take.emplace_back(init_seed);

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

        void Subset(SeedList &list, SeedList &subset, size_t take_top=10) {
            // Collect taxa counts
            size_t max_taxon = 0;
            size_t max_count = 0;
            for (auto seed : list) {
                m_taxa_counts[seed.taxid]++;
                if (m_taxa_counts[seed.taxid] > max_count) {
                    max_taxon = seed.taxid;
                    max_count = m_taxa_counts[seed.taxid];
                }
            }
            m_count_pairs.clear();
            int threshold = max_count >= 4 ? 1 : 0;
            for (auto pair : m_taxa_counts) {
                if (pair.second > threshold)
                    m_count_pairs.emplace_back(pair);
            }
            std::sort(m_count_pairs.begin(), m_count_pairs.end(), [](TaxonCountPair const& a, TaxonCountPair const& b) {
                return a.second > b.second;
            });
            m_selected_taxa.clear();

            for (auto i = 0; i < m_count_pairs.size() && i < take_top; i++) {
                m_selected_taxa.insert(m_count_pairs[i].first);
            }
            m_subset.clear();
            for (auto seed : list) {
//                if (m_selected_taxa.contains(seed.taxid)) {
                    m_subset.emplace_back(seed);
//                }
            }
        }

        void Reset() {
            // Assess what time seeding takes.
            m_all_seeds.clear();
        }

    public:
        size_t dummy = 0;
        Benchmark m_bm_seeding{"Seeding"};
        Benchmark m_bm_seed_subsetting{"Subset Seeds"};
        Benchmark m_bm_processing{"Sorting Seeds"};
        Benchmark m_bm_pairing{"Pairing"};
        Benchmark m_bm_sorting_anchors{"Sorting Anchors"};

        ChainAnchorFinder(KmerLookup& lookup, size_t k, size_t min_successful_lookups, size_t max_seed_size) :
                m_kmer_lookup(lookup), m_max_seed_size(max_seed_size), m_min_successful_lookups(min_successful_lookups), m_k(k) {
        };

        ChainAnchorFinder(ChainAnchorFinder const& other) :
                m_kmer_lookup(other.m_kmer_lookup), m_max_seed_size(other.m_max_seed_size), m_min_successful_lookups(other.m_min_successful_lookups), m_k(other.m_k) {};

        void operator () (KmerList& kmer_list, SeedList& seeds, ChainAnchorList& anchors, size_t read_length) {
            m_bm_seeding.Start();
            FindSeeds(kmer_list, seeds);
//            FindSeeds2(kmer_list, seeds);
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