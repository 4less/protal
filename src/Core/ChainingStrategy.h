//
// Created by fritsche on 12/10/22.
//

#pragma once


#include <cstdint>
#include <vector>
#include "KmerLookup.h"

namespace protal {
    using Seed = LookupResult;

    struct ChainLink {
        uint32_t genepos;
        uint16_t readpos;
        uint16_t length;

        uint32_t ReadStart() const {
            return readpos + length;
        }
        uint32_t ReadEnd() const {
            return readpos + length;
        }
        uint32_t GeneStart() const {
            return genepos + length;
        }
        uint32_t GeneEnd() const {
            return genepos + length;
        }
        uint32_t ExtendLength(uint32_t left, uint32_t right) {
            genepos -= left;
            readpos -= left;
            length += left;
            length += right;
        }

        ChainLink(uint32_t genepos, uint16_t readpos, uint16_t length) :
                genepos(genepos),
                readpos(readpos),
                length(length) {}

        void Merge(uint32_t other_readpos, size_t other_length) {
            length = other_readpos + other_length - readpos;
        }

        std::string ToString() {
            return "(RGL " + std::to_string(readpos) + "," + std::to_string(genepos) + ',' + std::to_string(length) + ')';
        }
    };

    using ChainList = std::vector<ChainLink>;
    struct ChainAlignmentAnchor {
        uint32_t taxid;
        uint32_t geneid;
        uint16_t total_length;
        bool forward;
        ChainList chain;

        ChainAlignmentAnchor(uint32_t taxid, uint32_t geneid, bool forward) :
                taxid(taxid),
                geneid(geneid),
                forward(forward) {};

        ChainLink& Back() noexcept {
            return chain.back();
        }

        ChainLink& Front() noexcept {
            return chain.front();
        }

        bool OverlapWithLeft(Seed const& seed) {
            return seed.readpos <= Back().readpos + Back().length;
        }
        void MergeWithLeft(Seed const& seed, size_t length) {
            Back().Merge(seed.readpos, length);
        }

        void AddSeed(Seed const& seed, size_t length) {
            if (!chain.empty() && OverlapWithLeft(seed)) {
                MergeWithLeft(seed, length);
            }
            chain.emplace_back(ChainLink(seed.genepos, seed.readpos, length));
        }

        std::string ToString() {
            std::string str = "";
            str += "[ChainAnchor " + std::to_string(taxid) + "," + std::to_string(geneid) + "  ";
            for (auto& link : chain) {
                str += link.ToString() + " ";
            }
            std::cout << "]";
        }
    };

}