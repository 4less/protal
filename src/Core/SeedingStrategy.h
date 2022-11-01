//
// Created by fritsche on 27/08/22.
//

#pragma once

#include "KmerLookup.h"

namespace protal {
    using Seed = LookupResult;
    using SeedList = std::vector<Seed>;
    using SeedIterator = SeedList::iterator;
    using SeedListSet = std::vector<SeedList>;
    using SeedListSetIterators = std::vector<SeedIterator>;

    struct AlignmentAnchor {
        LookupResult a;
        LookupResult b;
        size_t hit_anchor_count = 0;

        AlignmentAnchor(LookupResult& a, LookupResult &b, size_t hit_anchor_count) :
                a(a), b(b), hit_anchor_count(hit_anchor_count) {};
        AlignmentAnchor(LookupResult&& a, LookupResult &&b, size_t hit_anchor_count) :
                a(a), b(b), hit_anchor_count(hit_anchor_count) {};
    };


    using Anchor = CAlignmentAnchor;
//    using AnchorList = std::vector<Anchor>;
    using AnchorList = std::vector<std::pair<LookupResult, LookupResult>>;

    class Seeding {
        size_t m_seeding_position_count;
        SeedListSet m_seedlist_set;

        HTSLIB_SAM_H

        SeedListSetIterators m_seedlist_set_iterators;
        AnchorList m_anchor_list;
        std::vector<uint32_t> m_paired_idces;

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
            while (!IteratorHasNext(current_min_idx) && current_min_idx < m_seedlist_set_iterators.size()) {
                current_min_idx++;
            }
            if (current_min_idx == m_seedlist_set_iterators.size()) {
                std::cout << "Error in MinAnchorIdx" << std::endl;
                exit(9);
            }

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

        inline bool IteratorHasNext(int i) const {
            return m_seedlist_set_iterators[i] != m_seedlist_set[i].end();
        }
        inline bool IteratorHasIdenticalNext(int i) const {
            auto curr = m_seedlist_set_iterators[i];
            auto next = std::next(curr);
            auto end = m_seedlist_set[i].end();
            return curr != end &&
                   next != end &&
                   *next == *curr;
        }

        [[nodiscard]] inline bool IsSmallerThanOtherIndices(int idx) const {
            assert(IteratorHasNext(idx));
            auto& current_iterator = m_seedlist_set_iterators[idx];
            bool result = true;
            for (auto i = 0; i < m_seedlist_set_iterators.size(); i++) {
                auto& it = m_seedlist_set_iterators[i];
                result &= (!IteratorHasNext(i) || current_iterator == it || *current_iterator < *it);
            }
            return result;
        }

        inline bool HasIdenticalPair(int idx) {
            auto& current_iterator = m_seedlist_set_iterators[idx];
            return std::any_of(m_seedlist_set_iterators.begin(), m_seedlist_set_iterators.end(), [&current_iterator](SeedIterator const& it) {
                return !current_iterator->Equals(*it) && *current_iterator == *it;
            });
        }


        inline void IdenticalPairIdces(int idx, std::vector<uint32_t>& idces) {
            assert(IteratorHasNext(idx));
            idces.clear();
            auto& current_iterator = m_seedlist_set_iterators[idx];
            for (auto i = 0; i < m_seedlist_set_iterators.size(); i++) {
                auto& it = m_seedlist_set_iterators[i];
                if (IteratorHasNext(i) && *current_iterator == *it) idces.emplace_back(i);
            }
        }

        inline void AdvanceSeedIterator(int idx) {
            assert(idx < m_seedlist_set.size());
            auto end = m_seedlist_set[idx].end();

            if (!IteratorHasNext(idx)) return;

            do {
                std::advance(m_seedlist_set_iterators[idx], 1);
            } while (IteratorHasNext(idx) && IsSmallerThanOtherIndices(idx));
        }


        inline void AdvanceSeedIterators(std::vector<uint32_t> indices) {
            bool not_unique = std::any_of(indices.begin(), indices.end(), [this](int i) {
                auto curr = m_seedlist_set_iterators[i];
                auto next = std::next(curr);
                return (IteratorHasNext(i) && next != m_seedlist_set[i].end() && *curr == *next);
            });
            if (not_unique) {
//                std::cout << "Not unique. ";
                for (auto idx : indices) {
                    if (IteratorHasIdenticalNext(idx)) {
                        AdvanceSeedIterator(indices.at(0));
//                        std::cout << idx << std::endl;
                        break;
                    }
                }
            } else {
                for (auto idx : indices) {
                    assert(idx < m_seedlist_set.size());
                    AdvanceSeedIterator(idx);
                }
            }
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

        void PrintIterators() const {
            std::cout << "Print iterators: ";
            bool stop = false;
            for (auto i = 0; i < m_seedlist_set_iterators.size(); i++) {
                auto& it = m_seedlist_set_iterators[i];
                bool end = it == m_seedlist_set[i].end();

                if (!end && it->taxid > 1'000'000) {
                    stop = true;
                }

                std::cout << (end ? "End" : it->ToString()) << "\t\t";
            }
            std::cout << "(stop: " << stop << ")" << std::endl;

            if (stop) {
                std::cout << "__________________" << m_seedlist_set_iterators.size() << std::endl;
                for (auto i = 0; i < m_seedlist_set_iterators.size(); i++) {
                    auto &it = m_seedlist_set_iterators[i];
                    std::cout << i << " Predecessor: " << std::prev(it)->ToString() << std::endl;
                    bool end = it == m_seedlist_set[i].end();
                }
                std::cout << "__________________" << std::endl;

                PrintSeeds();
                Utils::Input();
                exit(5);
            }
        }

        void ComputeSeedPairs() {
            m_anchor_list.clear();

            int ident_pair_num = 0;

            int debug_counter = 0;
            while (LoopCondition()) {
                auto min_anchor_idx = MinAnchorIdx();


                IdenticalPairIdces(min_anchor_idx, m_paired_idces);

                if (m_paired_idces.size() > 1) {
//                    std::cout << "______________Indices: " <<  m_paired_idces.size() << std::endl;
//                    for (auto& idx : m_paired_idces) std::cout << idx << ", "; std::cout << endl;
                    ident_pair_num++;

                    // Todo fix - dont store same anchors.
                    auto& first = *(m_seedlist_set_iterators[m_paired_idces.front()]);
                    auto& second = *(m_seedlist_set_iterators[m_paired_idces.back()]);

                    if (!first.Equals(second)) {
                        m_anchor_list.emplace_back(
                                std::pair<Seed, Seed>(*(m_seedlist_set_iterators[m_paired_idces.front()]),
                                                      *(m_seedlist_set_iterators[m_paired_idces.back()])));
                    }
                    AdvanceSeedIterators(m_paired_idces);
                } else {
                    AdvanceSeedIterator(min_anchor_idx);
                }

                if (debug_counter++ > 1000) {
                    std::cout << "debug_counter:  " << debug_counter << std::endl;
                    std::cout << "Seeds" << std::endl;
                    PrintSeeds();

                    std::cout << "   \n" << std::endl;
                    std::cout << "Iterators" << std::endl;
                    PrintIterators();
                    exit(2);

                }

            }
//            Utils::Input();
        }
    };
}