//
// Created by fritsche on 03/09/22.
//

#pragma once

#include "Seedmap.h"
#include "Constants.h"
#include "SeedingStrategy.h"

namespace protal {

    using MapKey = std::pair<int32_t,int32_t>;
    struct MapKeyComp {
        bool operator() (MapKey const& a, MapKey const& b) {
            return (a.first != b.first && a.first < b.first) || (a.first == b.first && a.second < b.second);
        }
    };
    struct MapKeyHash {
        size_t operator() (MapKey const& a) const {
            return static_cast<size_t>(a.first) | (static_cast<size_t>(a.second) << 32lu);
        }
    };
    MapKeyComp mkc;

    template<typename KmerLookup>
    requires KmerLookupConcept<KmerLookup>
    class HashMapAnchorFinder {
        KmerLookup m_kmer_lookup;


        using SeedSeedList = std::vector<SeedList>;
        SeedList m_all_seeds;

        SeedList m_seed_tmp_take;
        SeedList m_seed_tmp_other;

        SeedList m_best_seed;
        std::vector<std::pair<int32_t,int32_t>> m_seed_size;

        size_t m_k = 15;


        tsl::sparse_map<MapKey,int32_t, MapKeyHash> m_map;

        inline void FindSeeds(KmerList &kmer_list) {
            size_t prev_size = 0;
            for (auto pair : kmer_list) {
                auto [mmer, pos] = pair;
                m_kmer_lookup.Get(m_all_seeds, mmer, pos);
                if (m_all_seeds.size() != prev_size) {
                    m_seed_size.emplace_back(std::pair<int32_t,int32_t>( pos, m_all_seeds.size() - prev_size ));
                    prev_size = m_all_seeds.size();
                }
            }
        }

        static MapKey GetKey(Seed const& seed) {
            return MapKey(static_cast<int32_t>(seed.taxid), static_cast<int32_t>(seed.geneid));
        }

        inline void ProcessSeeds(SeedList& seed_list) {
            m_map.clear();
            for (auto& seed : seed_list) {
                auto key = GetKey(seed);
                if (!m_map.contains(key)) {
                    m_map.insert( { key, 0 } );
                }
                m_map[key]++;
            }
        }

        void Reset() {

            // Assess what time seeding takes.
            m_all_seeds.clear();
            m_seed_size.clear();
            m_best_seed.clear();
        }

        using OffsetPair = std::pair<int, int>;

        size_t AnchorBaseCoverage(size_t k, SeedList &seeds) {
            size_t base_cov = 0;
            size_t last_seed_end = 0;
            size_t seed_start = 0;
            size_t seed_end = 0;
            for (auto& seed : seeds) {
                seed_start = seed.readpos;
                seed_end = seed.readpos + k;
                base_cov += seed_end - std::max(seed_start, last_seed_end);
                last_seed_end = seed_end;
            }
            return base_cov;
        }

        static bool OffsetPairMatch(OffsetPair const& a, OffsetPair const& b) {
            return a.first == b.first || a.second == b.second;
        }

//        static bool OneOutOfTwoInRange(OffsetPair const& a, OffsetPair const& b) {
//            return a.first == b.first || a.second == b.second;
//        }

        static std::string OffsetPairToString(OffsetPair const& a) {
            return "[" + std::to_string(a.first) + "," + std::to_string(a.second) + "]";
        }

        static OffsetPair Offset(Seed const& s) {
            return OffsetPair { static_cast<int>(s.genepos) + s.readpos, static_cast<int>(s.genepos) - s.readpos };
        }

        Anchor ExtractAnchor(SeedList &seeds) {
            auto& first = seeds.front();
            auto& second = seeds.back();
            size_t base_cov = AnchorBaseCoverage(m_k, seeds);
            if (first.genepos < second.genepos) {
                return AlignmentAnchor(first, second, base_cov);
            } else {
                return AlignmentAnchor(second, first, base_cov);
            }
        }

        void FindAnchorsSingleRef(SeedList &seeds, AlignmentAnchorList &anchors) {
            anchors.clear();


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
                    anchors.template emplace_back(ExtractAnchor(m_seed_tmp_take));
                }


//                if (!m_seed_tmp_other.empty()) {
//                    std::cout << std::string(40, '-') << " first" << std::endl;
//                    std::cout << init_seed.ToString() << " " << OffsetPairToString(Offset(init_seed)) << std::endl;
//                    std::cout << std::string(40, '-') << " take" << std::endl;
//                    for (auto &seed: m_seed_tmp_take) std::cout << seed.ToString() << " " << OffsetPairToString(Offset(seed)) << std::endl;
//                    std::cout << std::string(40, '-') << " other" << std::endl;
//                    for (auto &seed: m_seed_tmp_other) std::cout << seed.ToString() << " " << OffsetPairToString(Offset(seed)) << std::endl;
//                    Utils::Input();
//                }

                seeds.clear();
                m_seed_tmp_take.clear();
                std::swap(m_seed_tmp_other, seeds);
            }
        }

    public:
        size_t dummy = 0;
        size_t utilized_anchors = 0;
        Benchmark m_bm_seeding{"Seeding"};
        Benchmark m_bm_processing{"Process Seeds"};
        Benchmark m_bm_pairing{"Pairing"};

        HashMapAnchorFinder(KmerLookup& lookup) :
                m_kmer_lookup(lookup) {};

        HashMapAnchorFinder(HashMapAnchorFinder const& other) :
                m_kmer_lookup(other.m_kmer_lookup) {};

        void operator () (KmerList& kmer_list, AlignmentAnchorList& anchors) {
            anchors.clear();
            Reset();

            m_bm_seeding.Start();
            FindSeeds(kmer_list);
            m_bm_seeding.Stop();
            m_bm_processing.Start();
            ProcessSeeds(m_all_seeds);
            m_bm_processing.Stop();


            dummy += m_map.size();

            m_bm_pairing.Start();
            MapKey best;
            size_t best_count = 0;
            for (auto& [ckey, occ] : m_map) {
                if (occ > best_count) {
                    best = ckey;
                    best_count = occ;
                }
            }

//            std::cout << std::string(69, '-') << std::endl;
            for (auto& seed : m_all_seeds) {
                auto key = GetKey(seed);
                if (best.first == key.first && best.second == key.second) {
                    m_best_seed.emplace_back(seed);
//                    std::cout << seed.ToString() << std::endl;
                }
            }
            m_bm_pairing.Stop();

            // Insert best seed set as anchor

            FindAnchorsSingleRef(m_best_seed, anchors);
        }
    };

    template<typename KmerLookup>
    requires KmerLookupConcept<KmerLookup>
    class SimpleAnchorFinder {
        KmerLookup m_kmer_lookup;
        Seeding m_seeding{4};

        size_t m_idx_from_left = 0;
        size_t m_idx_from_right = 0;

//        size_t m_advance_after_find = 1;
        size_t m_advance_after_find = 3;


        // Anchors. Find better solution?!
        LookupList m_anchor_left_1;
        LookupList m_anchor_right_1;
        LookupList m_anchor_left_2;
        LookupList m_anchor_right_2;

        SeedList m_all_seeds;
        std::vector<std::pair<int32_t,int32_t>> m_seed_size;

        inline void FindSeeds(KmerList &kmer_list) {
            size_t prev_size = 0;
            for (auto pair : kmer_list) {
                auto [mmer, pos] = pair;
                m_kmer_lookup.Get(m_all_seeds, mmer, pos);
                if (m_all_seeds.size() != prev_size) {
                    m_seed_size.emplace_back(std::pair<int32_t,int32_t>( pos, m_all_seeds.size() - prev_size ));
                    prev_size = m_all_seeds.size();
                }
            }
        }

        inline void SearchAnchor(LookupList &anchor_list, KmerList &kmer_list, bool from_left) {
            if (from_left) {
                for (; m_idx_from_left < kmer_list.size(); m_idx_from_left++) {
                    auto [mmer, pos] = kmer_list[m_idx_from_left];
                    m_kmer_lookup.Get(anchor_list, mmer, pos);

                    if (!anchor_list.empty()) {
                        m_idx_from_left += m_advance_after_find;
                        break;
                    }
                }
            } else {
                for (; m_idx_from_right < kmer_list.size(); m_idx_from_right++) {
                    auto [mmer, pos] = kmer_list[kmer_list.size() - m_idx_from_right - 1];

                    m_kmer_lookup.Get(anchor_list, mmer, pos);
                    if (!anchor_list.empty()) {
                        m_idx_from_right += m_advance_after_find;
                        break;
                    }
                }
            }
        }
        void Reset() {
            m_idx_from_left = 0;
            m_idx_from_right = 0;

            m_anchor_left_1.clear();
            m_anchor_right_1.clear();
            m_anchor_left_2.clear();
            m_anchor_right_2.clear();

            // Assess what time seeding takes.
            m_all_seeds.clear();
            m_seed_size.clear();
        }

    public:
        size_t dummy = 0;
        size_t utilized_anchors = 0;

        SimpleAnchorFinder(KmerLookup& lookup) :
                m_kmer_lookup(lookup) {};

        SimpleAnchorFinder(SimpleAnchorFinder const& other) :
                m_kmer_lookup(other.m_kmer_lookup) {};

        void operator () (KmerList& kmer_list, AlignmentAnchorList& anchors) {
            anchors.clear();
            Reset();

            SearchAnchor(m_anchor_left_1, kmer_list, true);
            SearchAnchor(m_anchor_right_1, kmer_list, false);
            SearchAnchor(m_anchor_left_2, kmer_list, true);
            SearchAnchor(m_anchor_right_2, kmer_list, false);

            m_seeding.Add(0, m_anchor_left_1.begin(), m_anchor_left_1.end());
            m_seeding.Add(1, m_anchor_left_2.begin(), m_anchor_left_2.end());
            m_seeding.Add(2, m_anchor_right_2.begin(), m_anchor_right_2.end());
            m_seeding.Add(3, m_anchor_right_1.begin(), m_anchor_right_1.end());

            utilized_anchors += m_anchor_left_1.size();
            utilized_anchors += m_anchor_left_2.size();
            utilized_anchors += m_anchor_right_2.size();
            utilized_anchors += m_anchor_right_1.size();

            FindSeeds(kmer_list);
            dummy += m_all_seeds.size();

            std::sort(m_seed_size.begin(), m_seed_size.end(), [](std::pair<int32_t, int32_t> const& a, std::pair<int32_t, int32_t> const& b) {
                return a.second < b.second;
            });

            if (m_seed_size.size() > 0 ) {
                dummy += m_seed_size.front().first + m_seed_size.back().first;
            }

//            for (auto s : m_seed_size) {
//                std::cout << s.first << "," << s.second << "\t";
//            }
//            std::cout << std::endl;


            m_seeding.ComputeSeedPairs();

//            m_seeding.PrintSeeds();
//            m_seeding.PrintAnchors();

            for (auto& anchor_pair :  m_seeding.GetAnchorList()) {
                if (anchor_pair.first.genepos < anchor_pair.second.genepos) {
                    anchors.emplace_back( AlignmentAnchor(anchor_pair.first, anchor_pair.second, 2) );
                } else {
                    anchors.emplace_back( AlignmentAnchor(anchor_pair.second, anchor_pair.first, 2) );
                }
            }

//            m_seeding.PrintAnchors();
        }
    };

    template<typename KmerLookup>
    requires KmerLookupConcept<KmerLookup>
    class NaiveAnchorFinder {
        KmerLookup m_kmer_lookup;

        size_t m_idx_from_left = 0;
        size_t m_idx_from_right = 0;

        size_t m_advance_after_find = 10;

        // Anchors. Find better solution?!
        LookupList m_anchor_first;
        LookupList m_anchor_second;
        std::vector<std::pair<LookupResult, LookupResult>> m_anchor_pairs;

        inline void SearchAnchor(LookupList &anchor_list, KmerList &kmer_list, bool from_left) {
            if (from_left) {
                for (; m_idx_from_left < kmer_list.size(); m_idx_from_left++) {
                    auto [mmer, pos] = kmer_list[m_idx_from_left];
                    m_kmer_lookup.Get(anchor_list, mmer, pos);

                    if (!anchor_list.empty()) {
                        m_idx_from_left += m_advance_after_find;
                        break;
                    }
                }
            } else {
                for (; m_idx_from_right < kmer_list.size(); m_idx_from_right++) {
                    auto [mmer, pos] = kmer_list[kmer_list.size() - m_idx_from_right - 1];

                    m_kmer_lookup.Get(anchor_list, mmer, pos);
                    if (!anchor_list.empty()) {
                        m_idx_from_right += m_advance_after_find;
                        break;
                    }
                }
            }
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

        void Reset() {
            m_idx_from_left = 0;
            m_idx_from_right = 0;

            m_anchor_first.clear();
            m_anchor_second.clear();
            m_anchor_pairs.clear();
        }

    public:
        NaiveAnchorFinder(KmerLookup& lookup) :
                m_kmer_lookup(lookup) {};

        NaiveAnchorFinder(NaiveAnchorFinder const& other) :
                m_kmer_lookup(other.m_kmer_lookup) {};

        void operator () (KmerList& kmer_list, AlignmentAnchorList& anchors) {
            anchors.clear();
            Reset();
            SearchAnchor(m_anchor_first, kmer_list, true);
            SearchAnchor(m_anchor_second, kmer_list, false);

            GetAnchorPairsFromAnchors(m_anchor_first, m_anchor_second);

            for (auto& anchor_pair :  m_anchor_pairs) {
                anchors.emplace_back( AlignmentAnchor(anchor_pair.first, anchor_pair.second, 2) );
            }
        }
    };
}
