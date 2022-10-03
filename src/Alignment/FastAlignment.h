//
// Created by fritsche on 28/09/22.
//

#pragma once

#include <string>

namespace protal {
    class FastAligner {
    public:
        static std::pair<int, std::string> FastAlign(std::string const& a, std::string const& b, size_t mismatch_penalty = 4) {
            std::string cigar = "";
            double matches = 0;
            int alignment_score = 0;
            for (auto i = 0; i < a.length(); i++) {
                if (a[i] == b[i]) {
                    cigar += 'M';
                    matches++;
                } else {
                    cigar += 'X';
                    alignment_score -= mismatch_penalty;
                }
            }
            return { alignment_score, cigar };
        }

        static std::pair<int, std::string> FastAlign(std::string_view a, std::string_view b, size_t mismatch_penalty = 4) {
            std::string cigar = "";
            double matches = 0;
            int alignment_score = 0;
            for (auto i = 0; i < a.length(); i++) {
                if (a[i] == b[i]) {
                    cigar += 'M';
                    matches++;
                } else {
                    cigar += 'X';
                    alignment_score -= mismatch_penalty;
                }
            }
            return { alignment_score, cigar };
        }
    };

}