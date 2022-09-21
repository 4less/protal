//
// Created by fritsche on 03/09/22.
//

#pragma once

#include "Seedmap.h"
#include "Constants.h"
#include "SeedingStrategy.h"

namespace protal {
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

        SeedList m_seeds_tmp;
        SeedList m_all_seeds;

        inline bool FindSeeds(KmerList &kmer_list) {
            for (auto pair : kmer_list) {
                auto [mmer, pos] = pair;
                m_kmer_lookup.Get(m_seeds_tmp, mmer, pos);
                if (!m_seeds_tmp.empty()) {
                    m_
                }
            }
            return false;
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
        }

    public:
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
