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

        /**
         *
         * @param cigar uncompressed cigar string
         * @return
         */
        static double MismatchTail(std::string &cigar) {
            double a = 0;
            double b = 0;
            for (int i = 0; i < 12; i++) a += cigar[i] == 'M';
            for (int i = 0; i < 12; i++) b += cigar[cigar.length() - i - 1] == 'M';
            return std::min(a/12.0f, b/12.0f);
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