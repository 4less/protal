//
// Created by fritsche on 18/08/22.
//

#pragma once

#include "Seedmap.h"

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


    template <typename T>
    concept HasPut =
            requires(T t, size_t s) {
                { t.Put(s, s, s, s) } -> std::same_as<void>;
            };

    template <typename T>
    concept HasFirstPut =
            requires(T t, size_t s) {
                { t.FirstPut(s) } -> std::same_as<void>;
                { t.InitializeForPut() } -> std::same_as<void>;
            };

    class KmerPutterSM {
        Seedmap m_sm;
    public:
        inline void FirstPut(size_t key) {
            m_sm.CountUpKey(key);
        }

        inline void Put(size_t& key, size_t& taxid, size_t& geneid, size_t&& genepos) {
            m_sm.PutOMP(key, taxid, geneid, genepos);
        }

        inline void Put(size_t& key, size_t& taxid, size_t& geneid, size_t& genepos) {
            m_sm.PutOMP(key, taxid, geneid, genepos);
        }

        inline void InitializeForPut() {
            m_sm.BuildValuePointers();
        }

        inline void Save(std::ostream& ofs) {
            m_sm.Save(ofs);
        }
    };
}