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

    class KmerLookupSM {
        std::vector<LookupResult> m_results;
//        std::shared_ptr<Seedmap> m_sm;
        Seedmap& m_sm;
        Entry<20,20,20>* m_entry_begin = nullptr;
        Entry<20,20,20>* m_entry_end = nullptr;
        size_t m_taxid{};
        size_t m_geneid{};
        size_t m_genepos{};
    public:
//        inline void Load(std::istream& ifs) {
////            m_sm = std::make_shared<Seedmap>(Seedmap());
////            m_sm->Load(ifs);
//        }

        KmerLookupSM(Seedmap& map) : m_sm(map) {
//            std::cout << "Call regular constructor of KmerLookupSM" << std::endl;
        };

        KmerLookupSM(KmerLookupSM const& other) :
                m_sm(other.m_sm) {
//            std::cout << "Call copy constructor of KmerLookupSM" << std::endl;
        };

        inline void Get(LookupList& result, size_t &kmer, uint32_t readpos) {
            m_sm.Get(kmer, m_entry_begin, m_entry_end);

            if (m_entry_begin == nullptr || m_entry_end == nullptr) {
//                std::cerr << "shoops " << kmer << std::endl;
                return;
            }

            for (; m_entry_begin < m_entry_end; m_entry_begin++) {
                m_entry_begin->Get(m_taxid, m_geneid, m_genepos);
                result.emplace_back(LookupResult( m_taxid, m_geneid, m_genepos, readpos ));
            }
        }
    };
}