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
            return readpos;
        }
        uint32_t ReadEnd() const {
            return readpos + length;
        }
        uint32_t GeneStart() const {
            return genepos;
        }
        uint32_t GeneEnd() const {
            return genepos + length;
        }
        void ExtendLength(uint32_t left, uint32_t right) {
            assert(genepos >= left);
            assert(readpos >= left);
            genepos -= left;
            readpos -= left;
            length += left;
            length += right;
        }


        ChainLink(uint32_t genepos, uint16_t readpos, uint16_t length) :
                genepos(genepos),
                readpos(readpos),
                length(length) {}

        bool OverlapsWithLeft(ChainLink const& left) const {
            return ReadStart() <= left.ReadEnd();
        }

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
        uint16_t total_length = 0;
        bool forward = true;
        ChainList chain;

        ChainAlignmentAnchor(uint32_t taxid, uint32_t geneid, bool forward) :
                taxid(taxid),
                geneid(geneid),
                forward(forward) {};

        ChainLink& Back() noexcept {
            if (chain.empty()) exit(27);
            return chain.back();
        }

        ChainLink& Front() noexcept {
            if (chain.empty()) exit(27);
            return chain.front();
        }

        bool OverlapWithLeft(Seed const& seed) {
            return seed.readpos <= Back().readpos + Back().length;
        }
        void MergeWithLeft(Seed const& seed, size_t length) {
            total_length -= Back().length;
            Back().Merge(seed.readpos, length);
            total_length += Back().length;
        }

        void AddSeed(Seed const& seed, size_t length) {
            if (!chain.empty() && OverlapWithLeft(seed)) {
                MergeWithLeft(seed, length);
            } else {
                chain.emplace_back(ChainLink(seed.genepos, seed.readpos, length));
                total_length += length;
            }
        }
        size_t UpdateLength() {
            total_length = std::accumulate(chain.begin(), chain.end(), 0, [](size_t acc, ChainLink const& link) {
                return acc + link.length;
            });
            return total_length;
        }

        std::string ToString() {
            std::string str = "";
            str += "[ChainAnchor " + std::to_string(taxid) + "," + std::to_string(geneid) + "  ";
            size_t prev_end = 0;
            for (auto& link : chain) {
//                str += std::string(link.readpos,' ');
//                str += std::string(link.length,'X');
                str += link.ToString() + " ";
//                if (prev_end) {
//                    str += " -" + std::to_string(link.readpos - prev_end) + "- ";
//                }
//                str += '(' + std::to_string(link.readpos) + ',' + std::to_string(link.readpos + link.length) + ')';
//                prev_end = link.readpos + link.length;
            }
            str += "... ";
            str += std::to_string(total_length);
            str += ']';
            return str;
        }
    };

}