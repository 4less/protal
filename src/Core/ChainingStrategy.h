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
            return ReadStart() <= left.ReadEnd() && (readpos - genepos) == (left.readpos - left.genepos);
        }

        void Merge(uint32_t other_readpos, size_t other_length) {
            length = other_readpos + other_length - readpos;
        }

        void MergeWithRight(uint32_t other_readpos, size_t other_length) {
//            std::cout << "Length before: " << length << std::endl;
            if (other_readpos < (readpos+length)) {
                other_length -= readpos + length - other_readpos;
            }
            length += other_length;
//            std::cout << "Length after: " << length << std::endl;
        }

        std::string ToString() const {
            return "(RGL " + std::to_string(readpos) + "," + std::to_string(genepos) + ',' + std::to_string(length) + ')';
        }
    };

    using ChainList = std::vector<ChainLink>;
    struct ChainAlignmentAnchor {
        uint32_t taxid;
        uint32_t geneid;
        uint16_t total_length = 0;
        uint16_t unique = 0;
        uint16_t unique_best_two = 0;
        bool forward = true;
        ChainList chain;

        ChainAlignmentAnchor(uint32_t taxid, uint32_t geneid, bool forward, uint8_t unique = 0, uint8_t unique_best_two = 0) :
                taxid(taxid),
                geneid(geneid),
                forward(forward),
                unique(unique),
                unique_best_two(unique_best_two) {};

        ChainLink& Back() noexcept {
            if (chain.empty()) exit(27);
            return chain.back();
        }

        ChainLink& Front() noexcept {
            if (chain.empty()) exit(27);
            return chain.front();
        }

        bool OverlapWithLeft(Seed const& seed) {
            auto offa1 = (seed.genepos - seed.readpos);
            auto offa2 = (Back().genepos - Back().readpos);
            int indel = abs(int(offa1) - int(offa2));
//            if (seed.readpos <= Back().readpos + Back().length && indel != 0) {
//                std::cout << ToString() << std::endl;
//                std::cout << ToVisualString2() << std::endl;
//                std::cout << Back().ToString() << std::endl;
//                std::cout << seed.ToString() << std::endl;
//                std::cout << "indel: " << indel << std::endl;
////                Utils::Input();
//            }
            return seed.readpos <= Back().readpos + Back().length && indel == 0;
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

        std::string ToString() const {
            std::string str = "";
            str += "[ChainAnchor " + std::to_string(taxid) + "," + std::to_string(geneid) + "  ";
            size_t prev_end = 0;
            for (auto& link : chain) {
                str += link.ToString() + " ";
            }
            str += "... ";
            str += std::to_string(total_length);
            str += "\tunique:";
            str += std::to_string(unique);
            str += "\tunique2:";
            str += std::to_string(unique_best_two);
            str += ']';
            return str;
        }

        std::string ToVisualString() {
            std::string str = "";
            size_t prev_end = 0;
            size_t last_end = 0;
            for (auto& link : chain) {
                str += std::string(link.readpos, ' ');
                str += std::string(link.length, '#');
                str += '\n';
            }
            return str;
        }

        std::string ToVisualString2() {
            std::string str = "";
            size_t prev_end = 0;
            uint16_t last_end = 0;
            for (auto& link : chain) {
                if (link.readpos > last_end) {
                    str += std::string(link.readpos - last_end, ' ');
                }
                str += std::string((link.readpos + link.length) - std::max(link.readpos, last_end), '#');
                last_end = link.readpos + link.length;
            }
            return str;
        }
    };

}