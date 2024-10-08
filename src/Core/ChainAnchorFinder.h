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
        int m_indel = 6;

        TaxonCountPairs m_count_pairs;
        TaxaCounts m_taxa_counts;
        TaxonSet m_selected_taxa;
        SeedList m_subset;

        size_t m_min_successful_lookups = 4;
        size_t m_max_seed_size = 10;
        size_t m_lookup_index = 0;

        GenomeLoader& m_genome_loader;
        std::string* m_fwd = nullptr;
        std::string m_rev = "";

        ChainAnchorList m_anchors;

        bool m_error_in_read = false;

        size_t skip_count = 0;


        using RecoverySet = tsl::robin_set<uint64_t>;
        RecoverySet m_recovery;

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
            for (m_lookup_index = 0; m_lookup_index < m_lookups.size(); m_lookup_index++) {
                auto& lookup = m_lookups[m_lookup_index];
                m_kmer_lookup.GetFromLookup(seeds, lookup);
                total_lookups++;
                successful_lookups += (seeds.size() > previous_size);
                if (seeds.size() > m_max_seed_size && successful_lookups > m_min_successful_lookups) {
                    break;
                }
                previous_size = seeds.size();
            }
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

        static bool OffsetPairMatch(OffsetPair const& a, OffsetPair const& b, int indel = 0) {
            return min(abs(a.first - b.first), abs(a.second - b.second)) <= indel;
        }

        static bool OffsetPairAmbiguous(OffsetPair const& a, OffsetPair const& b) {
            return a.first == b.first && a.second == b.second;
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

        static inline void ReverseSeed(Seed& seed, size_t read_length, size_t k, size_t subk = 16) {
            seed.readpos = read_length - seed.readpos - k;
        }
        inline void ReverseSeedList(SeedList &seeds, size_t& read_length) {
            for (auto& seed : seeds) {
                ReverseSeed(seed, read_length, m_k);
            }
            std::reverse(seeds.begin(), seeds.end());
        }

        static bool IdenticalIgnoreN(std::string_view a, std::string_view b) {
            for (auto i = 0; i < a.length(); i++) {
                if (a[i] != b[i] && a[i] != 'N' && b[i] != 'N') {
                    return false;
                }
            }
            return true;
        }

        static bool IdenticalIgnoreAmbig(std::string_view a, std::string_view b) {
            for (auto i = 0; i < a.length(); i++) {
                if (a[i] != b[i] && (a[i] == 'A' || a[i] == 'C' || a[i] == 'G' || a[i] == 'T') && (b[i] == 'A' || b[i] == 'C' || b[i] == 'G' || b[i] == 'T')) {
                    return false;
                }
            }
            return true;
        }

        Anchor ExtractAndExtendAnchorFromSeed(Seed& seed, std::string& fwd, std::string& rev) {
//            if (seed.taxid == 0) exit(23); // remove
            auto& genome = m_genome_loader.GetGenome(seed.taxid);
            auto& gene = genome.GetGeneOMP(seed.geneid);

            ChainLink fwd_link = ChainLink(seed.genepos, seed.readpos, m_k);
            std::string_view qseedf(fwd.c_str() + seed.readpos, fwd_link.length);//query.substr(s.readpos, s.length);

            ReverseSeed(seed, fwd.length(), m_k);
            ChainLink rev_link = ChainLink(seed.genepos, seed.readpos, m_k);
            std::string_view qseedr(rev.c_str() + seed.readpos, rev_link.length);//query.substr(s.readpos, s.length);
            std::string_view rseed(gene.Sequence().c_str() + seed.genepos, fwd_link.length);

            bool forward = qseedf == rseed;
            if (qseedf != rseed && qseedr != rseed) {
                bool fwd_valid = IdenticalIgnoreN(qseedf, rseed);
                bool rev_valid = IdenticalIgnoreN(qseedr, rseed);
                forward = fwd_valid;
                // std::cerr << "fwd " << qseedf  << " " << qseedr << ": " << rseed << std::endl;
                // std::cerr << seed.ToString() << std::endl;
                // if (!fwd_valid && !rev_valid) {
                //     std::cerr << "completely wrong" << std::endl;
                // }
            }

            if (forward) {
                ExtendSeed(fwd_link, fwd, gene.Sequence());
            } else {
                ExtendSeed(rev_link, rev, gene.Sequence());
            }

            ChainAlignmentAnchor anchor{ seed.taxid, seed.geneid, forward };
            anchor.chain.emplace_back(forward ? fwd_link : rev_link);
            anchor.total_length = forward ? fwd_link.length : rev_link.length;

            return anchor;
        }

        Anchor ExtractAnchor(SeedList &seeds, size_t read_length) {
            if (seeds.size() == 1) {
                std::cout << "Seeds==1 "<< std::endl;
                std::cout << seeds.back().ToString() << std::endl;
                Utils::Input();
            }

            if (seeds.front().genepos == seeds.back().genepos) {
                seeds.pop_back();
            }

            bool forward = seeds.front().genepos < seeds.back().genepos;
            if (seeds.front().genepos)
            if (!forward) {
                ReverseSeedList(seeds, read_length);
            }

            auto uniques = std::accumulate(seeds.begin(), seeds.end(), 0, [](auto acc, Seed& seed){ return acc + seed.unique; });
            auto uniques_two = std::accumulate(seeds.begin(), seeds.end(), 0, [](auto acc, Seed& seed){ return acc + seed.unique_dist_two; });


            Anchor anchor(seeds.front().taxid, seeds.front().geneid, forward, uniques, uniques_two);
            bool stop = false;
            bool is_first = true;
            auto& first = seeds.front();
            auto& prev = first;
            for (auto& seed : seeds) {
//                if (!(seed.readpos > first.readpos && seed.genepos > first.genepos)) {
//                    std::cout << first.ToString() << std::endl;
//                    std::cout << seed.ToString() << std::endl;
//                    Utils::Input();
//                }

                if (!is_first && (prev.genepos >= seed.genepos || prev.readpos >= seed.readpos)) {
//                    std::cout << "Skip seeds because gene positions do not work out" << std::endl;
//                    std::cout << prev.ToString() << std::endl;
//                    std::cout << seed.ToString() << std::endl;
//                    std::cout << "Skip" << std::endl;
//                    skip_count += 1;
//                    std::cout << skip_count << std::endl;
                    continue;
                }
                is_first = false;
                anchor.AddSeed(seed, m_k);
                prev = seed;
            }

            return anchor;
        }

        void FindAnchorsSingleRef(SeedList &seeds, ChainAnchorList &anchors, size_t read_length) {

            while (seeds.size() > 1) {
                auto& init_seed = seeds[0];
                auto init_offset = Offset(init_seed);

                m_seed_tmp_take.emplace_back(init_seed);

                auto& prev = init_seed;
                for (auto i = 1; i < seeds.size(); i++) {
                    auto& seed = seeds[i];
                    auto seed_offset = Offset(seed);

                    if (seed.genepos == prev.genepos || seed.readpos == prev.readpos) {
//                        std::cout << "________________Skip" << std::endl;
//                        std::cout << prev.ToString() << std::endl;
//                        std::cout << seed.ToString() << std::endl;
                        skip_count += 1;
//                        std::cout << skip_count << std::endl;
                        continue;
                    }

                    if (OffsetPairMatch(init_offset, seed_offset, m_indel)) {
                        m_seed_tmp_take.emplace_back(seed);
                        prev = seed;
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
        Benchmark m_bm_operator{"Seed-finding operator"};
        Benchmark m_bm_reverse_complement{"Reverse complementing read"};
        Benchmark m_bm_seeding{"Seeding"};
        Benchmark m_bm_seed_subsetting{"Subset Seeds"};
        Benchmark m_bm_processing{"Sorting Seeds"};
        Benchmark m_bm_pairing{"Pairing"};
        Benchmark m_bm_extend_anchors{"Extending Anchors"};
        Benchmark m_bm_sorting_anchors{"Sorting Anchors"};
        Benchmark m_bm_recovering_anchors{"Recovering Anchors"};

        size_t recovered_count = 0;
        size_t total_count = 0;

        ChainAnchorFinder(KmerLookup& lookup, size_t k, size_t min_successful_lookups, size_t max_seed_size, GenomeLoader& genome_loader) :
                m_kmer_lookup(lookup), m_max_seed_size(max_seed_size), m_min_successful_lookups(min_successful_lookups),
                m_k(k), m_genome_loader(genome_loader) {
        };

        ChainAnchorFinder(ChainAnchorFinder const& other) :
                m_kmer_lookup(other.m_kmer_lookup), m_max_seed_size(other.m_max_seed_size),
                m_min_successful_lookups(other.m_min_successful_lookups), m_k(other.m_k),
                m_genome_loader(other.m_genome_loader) {};

        bool Success() {
            return !m_error_in_read;
        }

        void NextObservationBM() {
            m_bm_seeding.AddObservation();
        }

        size_t RecoverAnchors(ChainAnchorList& own_anchors, RecoverySet& recover, size_t align_top) {
            m_bm_recovering_anchors.Start();
            m_seed_tmp.clear();
            m_anchors.clear();
            size_t recovered = 0;
            for (auto& ele : m_recovery) {
                auto find = recover.find(ele);
                if (find != recover.end()) {
                    recover.erase(ele);
                }
            }
            for (auto i = 0; i < own_anchors.size(); i++) {
                auto& own_anchor = own_anchors[i];
                auto key = (static_cast<uint64_t>(own_anchor.taxid) << 40llu) | (static_cast<uint64_t>(own_anchor.geneid) << 20llu);
                if (recover.contains(key)) {
                    if (i >= align_top) {
                        m_anchors.emplace_back(own_anchor);
                        own_anchors.erase(own_anchors.begin() + i);
                        recovered++;
                    }
                    recover.erase(key);
                    i--;
                }
            }


            total_count++;
            if (recover.size() == 0) {
                if (!m_anchors.empty()) {
                    m_anchors.insert(m_anchors.end(), own_anchors.begin(), own_anchors.end());
                    std::swap(m_anchors, own_anchors);
                }
                return recovered;
            }
            recovered_count++;

            size_t to_recover = recover.size();
            for (; m_lookup_index < m_lookups.size(); m_lookup_index++) {
                auto& pointers = m_lookups[m_lookup_index];
                m_kmer_lookup.RecoverFromLookup(m_seed_tmp, pointers, recover);
                if (to_recover != recover.size()) {
                    recovered++;
                    to_recover = recover.size();
                }
                if (!m_seed_tmp.empty() && m_recovery.empty()) break;
            }

            for (auto& seed : m_seed_tmp) {
                // std::cerr << "From single seed " << seed.ToString() << std::endl;
                auto anchor = ExtractAndExtendAnchorFromSeed(seed, *m_fwd, m_rev);
                m_anchors.emplace_back(anchor);
            }

            if (!m_anchors.empty()) {
                m_anchors.insert(m_anchors.end(), own_anchors.begin(), own_anchors.end());
                std::swap(m_anchors, own_anchors);
            }


            m_bm_recovering_anchors.Stop();
            return recovered;
        }

        RecoverySet& BestAnchors() {
            return m_recovery;
        }

        //static
        std::pair<size_t, size_t> ExtendSeed(ChainLink& s, std::string const& query, std::string const& gene, uint16_t query_left_limit=0, uint16_t query_right_limit=0) {
            size_t extension_left = 0;
            size_t extension_right = 0;
            if (query_right_limit == 0) query_right_limit = query.length();

            constexpr bool check_middle = true;

            std::string_view qseed(query.c_str() + s.readpos, s.length);//query.substr(s.readpos, s.length);
            std::string_view rseed(gene.c_str() + s.genepos, s.length);

            if constexpr (check_middle) {
                std::string_view qseed(query.c_str() + s.readpos, s.length);
                std::string_view rseed(gene.c_str() + s.genepos, s.length);

                bool faulty = false;
                if (qseed != rseed) {
                    auto n_count_q = 0;
                    auto n_count_r = 0;
                    for (int qpos = s.readpos, rpos = s.genepos; qpos < s.readpos + s.length && rpos >= 0 && s.genepos + s.length; qpos++, rpos++) {
                        n_count_q += query[qpos] == 'N';
                        n_count_q += gene[rpos] == 'N';
                        if (query[qpos] != gene[rpos] && query[qpos] != 'N' && gene[rpos] != 'N') {
                            faulty = true;
                        }
                    }
                    // std::cerr << "N in query: " << n_count_q << ", N in gene: " << n_count_r << std::endl;
                }
                if (faulty) {
                    // auto seed = query.substr(s.readpos, s.length);
                    // std::cerr << "2 Seed is not exact match? " << s.length << std::endl;
                    m_error_in_read = true;
                    return { 0, 0 };
                }
            }

            for (int qpos = s.readpos - 1, rpos = s.genepos - 1;
                 qpos >= query_left_limit && rpos >= 0 && query[qpos] == gene[rpos];
                 qpos--, rpos--) {
                extension_left++;
            }

            for (int qpos = s.readpos + s.length, rpos = s.genepos + s.length;
                 qpos < query_right_limit && rpos < gene.length() && query[qpos] == gene[rpos];
                 qpos++, rpos++) {
                extension_right++;
            }

            assert(s.ReadStart() >= extension_left);
            assert(s.GeneStart() >= extension_left);
            assert(s.ReadEnd() + extension_right <= query.length());
            assert(s.GeneEnd() + extension_right <= gene.length());

            s.ExtendLength(extension_left, extension_right);
            return { extension_left, extension_right };
        }

        void ExtendSeed1(ChainLink& seed, ChainLink* prev, ChainLink* next, std::string const& query, std::string const& ref) {
            auto [lefta, righta] = ExtendSeed(seed, query, ref,
                                              (prev != nullptr ? prev->readpos + prev->length : 0),
                                              (next != nullptr ? next->readpos : ref.length()));

        }

        void ExtendSeed2(ChainLink& seed, std::string const& query, std::string const& ref) {
            ExtendSeed1(seed, nullptr, nullptr, query, ref);
        }

        bool CheckSeedsInAnchor(ChainAlignmentAnchor& anchor, std::string const& query) {
            auto& genome = m_genome_loader.GetGenome(anchor.taxid);
            auto& gene = genome.GetGeneOMP(anchor.geneid);

            bool faulty = false;
            for (auto i = 0; i < anchor.chain.size(); i++) {
                auto& seed = anchor.chain[i];

                auto seed_q = query.substr(seed.readpos, seed.length);
                auto seed_r = gene.Sequence().substr(seed.genepos, seed.length);

                if (!IdenticalIgnoreN(seed_q, seed_r) && IdenticalIgnoreN(KmerUtils::ReverseComplement(seed_q), seed_r)) {
                    faulty = true;
                    break;
                }
            }
            if (faulty) {
//#pragma omp critical(print_error)
//{
//                std::cerr << "Error with seed -------------------------" << anchor.forward << std::endl;
//                std::cerr << "Original orientation: " << anchor.forward << std::endl;
//                std::cerr << query << std::endl;
//                std::cerr << anchor.ToString() << std::endl;
//                std::cerr << anchor.ToVisualString2() << std::endl;
//                for (auto i = 0; i < anchor.chain.size(); i++) {
//                    auto& seed = anchor.chain[i];
//
//                    auto seed_q = query.substr(seed.readpos, seed.length);
//                    auto seed_r = gene.Sequence().substr(seed.genepos, seed.length);
//
//                    std::cerr << seed.ToString() << std::endl;
//                    std::cerr << std::string(seed.readpos, ' ') << seed_q << " fwd? " << anchor.forward << std::endl;
//                    std::cerr << std::string(seed.readpos, ' ') << KmerUtils::ReverseComplement(seed_q) << " fwd? " << !anchor.forward << std::endl;
//                    std::cerr << std::string(seed.readpos, ' ') << seed_r << " Reference" << std::endl;
//                }
//                std:cerr << " -------------------------" << anchor.forward << std::endl;
//                m_error_in_read = true;
//}
            }
            return faulty;
        }

        void ExtendAnchor(ChainAlignmentAnchor& anchor, std::string const& query) {
            bool faulty = CheckSeedsInAnchor(anchor, query);

//            if (faulty) {
//                std::cerr << "Faulty___________" << anchor.chain.size() << std::endl;
//                std::cerr << anchor.ToString() << std::endl;
//            }

            auto& genome = m_genome_loader.GetGenome(anchor.taxid);
            auto& gene = genome.GetGeneOMP(anchor.geneid);

            // std::cerr << anchor.ToVisualString2() << std::endl;
            // std::cerr << anchor.ToString() << std::endl;
            bool changed_once = false;
            for (auto i = 0; i < anchor.chain.size(); i++) {
                auto& seed = anchor.chain[i];

                auto seed_q = query.substr(seed.readpos, seed.length);
                auto seed_r = gene.Sequence().substr(seed.genepos, seed.length);
                bool validseed = IdenticalIgnoreAmbig(seed_q, seed_r);

                if (!validseed) {
//                    std::cerr << "Invalid seed" << std::endl;
                    skip_count += 1;
                    return;
                }

                if (!IdenticalIgnoreAmbig(seed_q, seed_r) && IdenticalIgnoreAmbig(KmerUtils::ReverseComplement(seed_q), seed_r)) {

                    std::cout << "Is forward? " << anchor.forward << std::endl;
                    std::cout << anchor.ToString() << std::endl;
                    std::cout << "Query: " << query << std::endl;

                    m_error_in_read = true;
                    anchor.forward = !anchor.forward;
                    changed_once = true;
                }

                if (m_error_in_read) {
                    std::cerr << seed.ToString() << std::endl;
                    std::cerr << seed_q << std::endl;
                    std::cerr << seed_r << std::endl;
                }

                if (i > 0) {
                    auto offa1 = (seed.genepos - seed.readpos);
                    auto offa2 = (anchor.chain[i-1].genepos - anchor.chain[i-1].readpos);
                    int indel = abs(int(offa1) - int(offa2));
                    if (indel > 6 || !validseed) {
                        std::cout << "Index: " << i << std::endl;
                        std::cout << anchor.ToString() << std::endl;
                        std::cout << anchor.ToVisualString2() << std::endl;
                        std::cout << "Offa1: " << offa1 << " Offa2: " << offa2 << ", Valid seed? " << validseed << std::endl;
                        std::cout << "Indels: " << indel << std::endl;
                        std::cout << "Previous: " << anchor.chain[i-1].ToString() << std::endl;
                        std::cout << "Current:  " << seed.ToString() << std::endl;

                        std::cout << seed_q << std::endl;
                        std::cout << seed_r << std::endl;
                    }
                }



                ExtendSeed1(seed,
                           (i > 0) ? &anchor.chain[i-1] : nullptr,
                           (i+1) < anchor.chain.size() ? &anchor.chain[i+1] : nullptr,
                           query, gene.Sequence());

                if (i > 0 && seed.OverlapsWithLeft(anchor.chain[i-1])) {
                    anchor.chain[i-1].Merge(seed.readpos, seed.length);
                    anchor.chain.erase(anchor.chain.begin() + i);
                    i--;
                }
            }

            if (anchor.total_length == query.size()) {

            }

            anchor.UpdateLength();
        }

        void operator () (KmerList& kmer_list, SeedList& seeds, ChainAnchorList& anchors, std::string& query) {
            m_error_in_read = false;

            m_bm_operator.Start();

            m_bm_reverse_complement.Start();
            auto read_length = query.length();
            m_fwd = &query;
            m_rev = KmerUtils::ReverseComplement(query);
            m_bm_reverse_complement.Stop();

            m_recovery.clear();

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

            m_bm_extend_anchors.Start();
            for (auto& anchor : anchors) {
                if (anchor.chain.size() == 1 && anchor.total_length < 20) {
                    auto& seed = anchor.chain.back();
                    auto& genome = m_genome_loader.GetGenome(anchor.taxid);
                    auto& gene = genome.GetGeneOMP(anchor.geneid);
                    auto seed_q_fwd = m_fwd->substr(seed.readpos, seed.length);
                    auto seed_q_rev = m_rev.substr(seed.readpos, seed.length);
                    auto seed_r = gene.Sequence().substr(seed.genepos, seed.length);

//                    std::cout << "Single anchor " << anchor.forward << std::endl;
//                    std::cout << seed_q_fwd << " FWD" << std::endl;
//                    std::cout << seed_q_rev << " REV" << std::endl;
//                    std::cout << seed_r << " REF" << std::endl;

                    if (IdenticalIgnoreAmbig(seed_q_fwd, seed_r)) {
                        anchor.forward = true;
                    } else if (IdenticalIgnoreAmbig(seed_q_rev, seed_r)) {
                        anchor.forward = false;
                    } else {
                        std::cout << "Skip" << std::endl;
                        continue;
                    }

//                    Utils::Input();
                }

                ExtendAnchor(anchor, anchor.forward ? *m_fwd : m_rev);
            }
            m_bm_extend_anchors.Stop();

//            if (m_error_in_read) {
//#pragma omp critical(print_error)
//                {
//                    std::cerr << "\nMinimizers (" << kmer_list.size() << ")" << std::endl;
//                    for (auto& [kmer, pos] : kmer_list) {
//                        std::cerr << std::string(pos, ' ') << string_view(query.c_str() + pos, m_k) << "(" << KmerUtils::ToString(kmer, m_k*2) << ")" << std::endl;
//                    }
//                    std::cerr << "\nSeeds (" << seeds.size() << ")" << std::endl;
//                    for (auto& seed : seeds) {
//                        auto& gene = m_genome_loader.GetGenome(seed.taxid).GetGene(seed.geneid).Sequence();
//                        std::cerr << std::string(seed.readpos, ' ') << string_view(query.c_str() + seed.readpos, m_k) << "  " << seed.ToString() << std::endl;
//                        std::cerr << std::string(seed.readpos, ' ') << string_view(gene.c_str() + seed.genepos, m_k) << std::endl;
//                    }
//                    std::cerr << "\nAnchors (" << anchors.size() << ")" << std::endl;
//                    for (auto& anchor : anchors) {
//                        std::cerr << anchor.ToString() << " Forward? " << anchor.forward << std::endl;
//                        std::cerr << query << std::endl;
//                        std::cerr << KmerUtils::ReverseComplement(query) << std::endl;
//                        std::cerr << anchor.ToString() << std::endl;
//                        std::cerr << anchor.ToVisualString2() << std::endl;
//                        auto& gene = m_genome_loader.GetGenome(anchor.taxid).GetGene(anchor.geneid).Sequence();
//                        for (auto& s : anchor.chain) {
//                            std::cerr << std::string(s.readpos, ' ') << string_view(query.c_str() + s.readpos, m_k) << "  " << s.ToString() << std::endl;
//                            std::cerr << std::string(s.readpos, ' ') << string_view(gene.c_str() + s.genepos, m_k) << std::endl;
//                        }
//                        std::cerr << std::endl;
//                    }
//                    std::cerr << "-----------------------------------------------------------" << std::endl;
//                }
//            }

            m_bm_sorting_anchors.Start();
            std::sort(anchors.begin(), anchors.end() ,
                      [](Anchor const& a, Anchor const& b) {
                return a.total_length > b.total_length;
            });
            m_bm_sorting_anchors.Stop();

            m_bm_operator.Stop();
            for (int i = 0; i < 3 && i < anchors.size(); i++) {
                m_recovery.insert((static_cast<uint64_t>(anchors[i].taxid) << 40llu) | (static_cast<uint64_t>(anchors[i].geneid) << 20llu));
            }
        }
    };

}