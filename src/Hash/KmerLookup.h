//
// Created by fritsche on 18/08/22.
//

#pragma once

#include "Seedmap.h"
#include "Constants.h"

namespace protal {
    // KmerPutter is a Metatemplate programming "Interface"  and needs to satisfy the following functions.
    // KmerPutter needs to handle OMP collisions or have the map datastructure deal with it internally
    // -> Put function:  Put(size_t key, size_t taxid, size_t geneid, size_t genepos)
    // -> Boolean that says if it needs pre_processing
    //    template <typename T, typename S=size_t>
    //    concept KmerPutter = requires (T t, S s) {
    //        { t.Put(s, s, s, s) } -> std::same_as<void>;
    //        { t.FirstPut(s, s, s, s) } -> std::same_as<void>;
    //    };

    struct LookupResult {
        uint32_t taxid = UINT32_MAX;
        uint32_t geneid = UINT32_MAX;
        uint32_t genepos = UINT32_MAX;
        uint32_t readpos = UINT32_MAX;

        LookupResult() {};

        LookupResult(uint32_t taxid, uint32_t geneid, uint32_t genepos, uint32_t readpos) :
                taxid(taxid), geneid(geneid), genepos(genepos), readpos(readpos) {};

        inline bool SameTaxon(LookupResult const& other) const {
            return taxid == other.taxid;
        }

        inline bool SameGene(LookupResult const& other) const {
            return geneid == other.geneid;
        }

        inline bool FromSameSequence(LookupResult const& other) const {
            return geneid == other.geneid && taxid == other.taxid;
        }

        inline bool Equals(LookupResult const& other) const {
            return geneid == other.geneid && taxid == other.taxid && genepos == other.genepos;
        }

        static inline int IntComp(uint32_t const& a, uint32_t const& b) {
            return (a < b) * 1 + (a > b) * -1;
        }

//        inline int operator <=> (LookupResult const& other) const {
//            bool taxid_comp = IntComp(this->taxid, other.taxid);
//            if (taxid_comp != 0) {
//                return taxid_comp;
//            }
//            return IntComp(this->geneid, other.geneid);
//        }

        inline bool operator > (LookupResult const& other) const {
            if (this->taxid != other.taxid) {
                return this->taxid > other.taxid;
            }
            if (this->geneid != other.geneid) {
                return this->geneid > other.geneid;
            }
            return false;
        };

        template<uint64_t taxshift, uint64_t geneshift>
        inline uint64_t ToUINT64_t() const {
            return (static_cast<uint64_t>(taxid) << taxshift) | (static_cast<uint64_t>(geneid) << geneshift) | static_cast<uint64_t>(readpos);
        }

        inline uint64_t ToUINT64_t2() const {
            return (static_cast<uint64_t>(taxid) << 40llu) | (static_cast<uint64_t>(geneid) << 20llu) | static_cast<uint64_t>(readpos);
        }

        inline bool operator < (LookupResult const& other) const {
            if (this->taxid != other.taxid) {
                return this->taxid < other.taxid;
            }
            if (this->geneid != other.geneid) {
                return this->geneid < other.geneid;
            }
            return false;
        };
//        inline bool operator > (LookupResult const& other) const = default;
        inline int operator == (LookupResult const& other) const {
            return this->taxid == other.taxid && this->geneid == other.geneid;
        }

        static bool SortByReadComparator(LookupResult const& a, LookupResult const& b) {
            if (a.taxid != b.taxid) {
                return a.taxid < b.taxid;
            }
            if (a.geneid != b.geneid) {
                return a.geneid < b.geneid;
            }
            if (a.readpos != b.readpos) {
                return a.readpos < b.readpos;
            }
            return false;
        }

        static bool SortByReadComparator2(LookupResult const& a, LookupResult const& b) {
            return a.ToUINT64_t<40,20>() < b.ToUINT64_t<40,20>();
        }

        std::string ToString() const {
            std::string str;
            str += "[LUR: ";
            str += std::to_string(taxid) + ", ";
            str += std::to_string(geneid) + ", ";
            str += std::to_string(genepos) + ", ";
            str += std::to_string(readpos) + "]";
            return str;
        }
    };

    struct LookupPointer {
        Entry<20,20,20>* values_begin = nullptr;
        Entry<20,20,20>* values_end = nullptr;
        uint32_t* flex_begin = nullptr;
        uint32_t* flex_end = nullptr;
        uint32_t flex_key = 0;
        uint32_t read_pos = 0;
        uint32_t size = 0;
        bool retrieved = false;

        std::string ToString() {
            std::string str;
            str += "[" + std::to_string(read_pos) + ", " + std::to_string(size) + ", " + KmerUtils::ToString(flex_key, 32) + "]";
            return str;
        }
    };



    class KmerLookupSM {
        std::vector<LookupResult> m_results;
        std::vector<LookupPointer> m_lookups;
//        std::shared_ptr<Seedmap> m_sm;
        Seedmap& m_sm;
        Entry<20,20,20>* m_entry_begin = nullptr;
        Entry<20,20,20>* m_entry_end = nullptr;
        uint32_t* m_flex_begin = nullptr;
        uint32_t* m_flex_end = nullptr;

        std::vector<uint16_t> flex_vector;

        size_t m_flex_k = 12;
        size_t m_flex_k_half = m_flex_k/2;
        size_t m_flex_k_bits = m_flex_k*2;
        size_t m_max_ubiquity = 16;

        size_t m_taxid{};
        size_t m_geneid{};
        size_t m_genepos{};

        LookupPointer m_lookup_tmp;
    public:
//        inline void Load(std::istream& ifs) {
////            m_sm = std::make_shared<Seedmap>(Seedmap());
////            m_sm->Load(ifs);
//        }
        using RecoverySet = tsl::robin_set<uint64_t>;

        KmerLookupSM(Seedmap& map, size_t max_ubiquity=128) :
                m_sm(map), m_flex_k(map.m_flex_k), m_flex_k_bits(map.m_flex_k_bits), m_flex_k_half(map.m_flex_k/2),
                m_max_ubiquity(max_ubiquity) {
        };

        KmerLookupSM(KmerLookupSM const& other) :
                m_sm(other.m_sm), m_flex_k(other.m_sm.m_flex_k), m_flex_k_bits(other.m_sm.m_flex_k_bits),
                m_max_ubiquity(other.m_max_ubiquity), m_flex_k_half(other.m_sm.m_flex_k/2) {
        };

        void Clear() {
            m_lookups.clear();
        }

        inline void Get(std::vector<LookupPointer>& result, size_t &kmer, uint32_t readpos) {
            size_t flex_key = m_sm.FlexKey(kmer);
            m_sm.Get(kmer, m_lookup_tmp.values_begin, m_lookup_tmp.values_end, m_lookup_tmp.flex_begin, m_lookup_tmp.flex_end);
            if (m_lookup_tmp.values_begin == nullptr) return;
            m_lookup_tmp.flex_key = flex_key;
            m_lookup_tmp.read_pos = readpos;
            m_lookup_tmp.size = m_lookup_tmp.values_end - m_lookup_tmp.values_begin;
            result.emplace_back(std::move(m_lookup_tmp));
        }

        inline void GetFromLoopupSIMD(LookupList& result, LookupPointer& pointers) {

            if (pointers.flex_begin != nullptr) {
                flex_vector.clear();
                auto max = 0;
                auto max_count = 0;
                for (auto begin = pointers.flex_begin; begin < pointers.flex_end; begin++) {
                    auto sim = Seedmap::Similarity(*begin, pointers.flex_key);
                    flex_vector.emplace_back(sim);
                    if (sim > max)  {
                        max = sim;
                        max_count = 0;
                    }
                    max_count += (sim == max);
                }

                if (max_count > m_max_ubiquity) {
                    //                    std::cout << "No recovery " << max_count << "/" << m_entry_end - m_entry_begin << std::endl;
                    return;
                }

                //                std::cout << "Recovery" << std::endl;
                for (auto i = 0; i < flex_vector.size(); i++) {
                    if (flex_vector[i] == max) {
                        pointers.values_begin[i].Get(m_taxid, m_geneid, m_genepos);

                        if (m_taxid == 0) {
                            // Debug
                            std::cout << pointers.values_begin[i].ToString() << std::endl;
                            for (auto& e : result) {
                                std::cout << e.ToString() << std::endl;
                            }
                            continue;
                        }

                        result.emplace_back(LookupResult( m_taxid, m_geneid, m_genepos, pointers.read_pos + m_flex_k_half ));
                    }
                }
            }
        }

        inline void GetFromLookup(LookupList& result, LookupPointer& pointers) {
            if (pointers.flex_begin != nullptr) {
                flex_vector.clear();
                auto max = 0;
                auto max_count = 0;
                for (auto begin = pointers.flex_begin; begin < pointers.flex_end; begin++) {
                    auto sim = Seedmap::Similarity(*begin, pointers.flex_key);
                    flex_vector.emplace_back(sim);
                    if (sim > max)  {
                        max = sim;
                        max_count = 0;
                    }
                    max_count += (sim == max);
                }

                if (max_count > m_max_ubiquity) {
//                    std::cout << "No recovery " << max_count << "/" << m_entry_end - m_entry_begin << std::endl;
                    return;
                }

//                std::cout << "Recovery" << std::endl;
                for (auto i = 0; i < flex_vector.size(); i++) {
                    if (flex_vector[i] == max) {
                        pointers.values_begin[i].Get(m_taxid, m_geneid, m_genepos);

                        if (m_taxid == 0) {
                            // Debug
                            std::cout << pointers.values_begin[i].ToString() << std::endl;
                            for (auto& e : result) {
                                std::cout << e.ToString() << std::endl;
                            }
                            continue;
                        }

                        result.emplace_back(LookupResult( m_taxid, m_geneid, m_genepos, pointers.read_pos + m_flex_k_half ));
                    }
                }
            } else {
                m_entry_begin = pointers.values_begin;
                for (; m_entry_begin < pointers.values_end; m_entry_begin++) {
                    m_entry_begin->Get(m_taxid, m_geneid, m_genepos);
                    result.emplace_back(LookupResult( m_taxid, m_geneid, m_genepos, pointers.read_pos + m_flex_k_half ));
                }
            }
        }

        inline void RecoverFromLookup(LookupList& result, LookupPointer& pointers, RecoverySet& recovery) {
            for (auto iter = pointers.values_begin; iter < pointers.values_end; iter++) {
                if (recovery.contains(iter->MaskPosition())) {
//                    std::cout << "Recovery: " << iter->ToString() << std::endl;
                    auto [taxid, geneid, genepos] = iter->Get();
                    result.emplace_back(LookupResult(taxid, geneid, genepos, pointers.read_pos + m_flex_k_half));
                    recovery.erase(recovery.find(iter->MaskPosition()));
                    return;
                }
            }
        }

        inline void Get(LookupList& result, size_t &kmer, uint32_t readpos, RecoverySet* choose=nullptr, LookupList* recovered_results=nullptr) {
            m_sm.Get(kmer, m_entry_begin, m_entry_end, m_flex_begin, m_flex_end);


            if (m_entry_begin == nullptr || m_entry_end == nullptr) {
                return;
            }

            if (m_flex_begin != nullptr) {
                size_t flex_key = m_sm.FlexKey(kmer);
                flex_vector.clear();
                auto max = 0;
                auto max_count = 0;
                for (auto begin = m_flex_begin; begin < m_flex_end; begin++) {
                    auto sim = Seedmap::Similarity(*begin, flex_key);
                    flex_vector.emplace_back(sim);
                    if (sim > max)  {
                        max = sim;
                        max_count = 0;
                    }
                    max_count += (sim == max);
                }

                if (max_count > m_max_ubiquity) {
//                    std::cout << "No recovery " << max_count << "/" << m_entry_end - m_entry_begin << std::endl;
                    return;
                }

//                std::cout << "Recovery" << std::endl;
                for (auto i = 0; i < flex_vector.size(); i++) {
                    if (flex_vector[i] == max) {
                        m_entry_begin[i].Get(m_taxid, m_geneid, m_genepos);
                        result.emplace_back(LookupResult( m_taxid, m_geneid, m_genepos, readpos + m_flex_k_half ));
                    }
//                    // Recovering step.
//                    if (choose != nullptr && choose->contains(std::pair<uint32_t, uint32_t>(static_cast<uint32_t>(m_taxid), static_cast<uint32_t>(m_geneid)))) {
//                        recovered_results->emplace_back(LookupResult( m_taxid, m_geneid, m_genepos, readpos + m_flex_k_half ));
//                    }
                }
            } else {
                for (; m_entry_begin < m_entry_end; m_entry_begin++) {
                    m_entry_begin->Get(m_taxid, m_geneid, m_genepos);
                    result.emplace_back(LookupResult( m_taxid, m_geneid, m_genepos, readpos + m_flex_k_half ));
                }
            }

        }
    };
}