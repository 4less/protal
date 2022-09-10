//
// Created by fritsche on 27/08/22.
//

#pragma once

#include "KmerLookup.h"

namespace protal {
    struct AlignmentAnchor {
        LookupResult a;
        LookupResult b;
        size_t hit_anchor_count = 0;

        AlignmentAnchor(LookupResult& a, LookupResult &b, size_t hit_anchor_count) :
                a(a), b(b), hit_anchor_count(hit_anchor_count) {};
        AlignmentAnchor(LookupResult&& a, LookupResult &&b, size_t hit_anchor_count) :
                a(a), b(b), hit_anchor_count(hit_anchor_count) {};
    };

    using Seed = LookupResult;
    using SeedList = std::vector<Seed>;
    using SeedIterator = SeedList::iterator;
    using SeedListSet = std::vector<SeedList>;
    using SeedListSetIterators = std::vector<SeedIterator>;
    using Anchor = AlignmentAnchor;
//    using AnchorList = std::vector<Anchor>;
    using AnchorList = std::vector<std::pair<LookupResult, LookupResult>>;

    class Seeding {
        size_t m_seeding_position_count;
        SeedListSet m_seedlist_set;

        SeedListSetIterators m_seedlist_set_iterators;
        AnchorList m_anchor_list;

    public:
        Seeding(size_t seeding_position_count) : m_seeding_position_count(seeding_position_count) {
            m_seedlist_set.resize(seeding_position_count);
            m_seedlist_set_iterators.resize(seeding_position_count);
        }

        void Add(size_t idx, SeedIterator begin, SeedIterator end) {
            assert(idx < m_seedlist_set.size());
            m_seedlist_set[idx].clear();
            while (begin < end) {
                m_seedlist_set[idx].emplace_back(*begin);
                begin++;
            }
            m_seedlist_set_iterators[idx] = m_seedlist_set[idx].begin();
        }

        bool LoopCondition() const {
            size_t not_end_count = 0;
            for (auto i = 0; i < m_seedlist_set.size(); i++) {
                if (m_seedlist_set_iterators[i] != m_seedlist_set[i].end()) not_end_count++;
                if (not_end_count >= 2) return true;
            }
            return false;
        }

        [[nodiscard]] int MinAnchorIdx() const {

            int current_min_idx = 0;
            auto current_min = m_seedlist_set_iterators[current_min_idx];


            // Loop over "iterators"
            for (int i = 1; i < m_seedlist_set_iterators.size(); i++) {

                auto& current_element = m_seedlist_set_iterators[i];
                auto current_end = m_seedlist_set[i].end();

                if (current_element == current_end) continue;

                if (*current_element < *current_min) {
                    current_min = current_element;
                    current_min_idx = i;
                }
            };
            return current_min_idx;
        }

        [[nodiscard]] inline bool IsSmallerThanOtherIndices(int idx) const {
            auto& current_iterator = m_seedlist_set_iterators[idx];
            return std::all_of(m_seedlist_set_iterators.begin(), m_seedlist_set_iterators.end(), [&current_iterator](SeedIterator const& it) {
               return current_iterator == it || *current_iterator < *it;
            });
        }

        inline bool HasIdenticalPair(int idx) {
            auto& current_iterator = m_seedlist_set_iterators[idx];
            return std::any_of(m_seedlist_set_iterators.begin(), m_seedlist_set_iterators.end(), [&current_iterator](SeedIterator const& it) {
                return !current_iterator->Equals(*it) && *current_iterator == *it;
            });
        }

        inline void AdvanceSeedIterator(int idx) {
            assert(idx < m_seedlist_set.size());
            auto end = m_seedlist_set[idx].end();
            do {
                std::advance(m_seedlist_set_iterators[idx], 1);
            } while (m_seedlist_set_iterators[idx] != end && IsSmallerThanOtherIndices(idx));
        }

        inline void PrintSeeds() const {
            std::cout << "Print Seedlists " << m_seedlist_set.size() << std::endl;
            for (auto& seedlist : m_seedlist_set) {
                std::cout << "Seedlist: " << seedlist.size() << std::endl;
                for (auto& seed : seedlist) {
                    std::cout << seed.ToString() << std::endl;
                }
            }
        };

        inline void PrintAnchors() const {
            std::cout << "Print Anchors " << m_anchor_list.size() << std::endl;
            for (auto& pair : m_anchor_list) {
                std::cout << pair.first.ToString() << "\t" << pair.second.ToString() << std::endl;
            }
        };

        inline void SortSeedListSet() {
            std::cout << "Before sort" << std::endl;
            for (auto& seedlist : m_seedlist_set) {
                std::cout << "Sort list" << std::endl;
                for (auto& seed : seedlist) {
                    std::cout << seed.ToString() << std::endl;
                }
            }

            std::cout << "Sort lists" << std::endl;
            for (auto& seedlist : m_seedlist_set) {
                std::cout << "Sort list" << std::endl;
                for (auto& seed : seedlist) {
                    std::cout << seed.ToString() << std::endl;
                }
                std::sort(seedlist.begin(), seedlist.end(), [](Seed const& a, Seed const& b) {
                    return a < b;
                });
            }
        }

        AnchorList& GetAnchorList() {
            return m_anchor_list;
        }

        void ComputeSeedPairs() {
            m_anchor_list.clear();

            std::cout << "Anchorize" << std::endl;
            std::cout << "m_anchor_list: " << m_anchor_list.size() << std::endl;
            while (LoopCondition()) {
                auto min_anchor_idx = MinAnchorIdx();

                for (auto &it: m_seedlist_set_iterators) {
                    std::cout << it->ToString() << "\t\t";
                }
                std::cout << std::endl;

                Utils::Input();

                if (HasIdenticalPair(min_anchor_idx)) {
                    m_anchor_list.emplace_back(std::pair<Seed, Seed>(*(m_seedlist_set_iterators[0]), *(m_seedlist_set_iterators[1])));
                }

                AdvanceSeedIterator(min_anchor_idx);
            }
        }
    };
}