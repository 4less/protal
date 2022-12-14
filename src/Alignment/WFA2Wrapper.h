//
// Created by fritsche on 17/08/22.
//

#pragma once

#include <string>
#include <istream>
#include <iostream>
#include "wfa2-lib/bindings/cpp/WFAligner.hpp"


using namespace::wfa;

namespace protal {
    struct AlignmentInfo {
        uint32_t alignment_start = 0;
        uint32_t alignment_length = 0;
        float alignment_ani = 0.0f;
        std::string compressed_cigar;

        std::string ToString() const {
            std::string str = std::to_string(alignment_start);
            str += ", " + std::to_string(alignment_length);
            str += ", " + std::to_string(alignment_ani);
            str += ", " + compressed_cigar;
            return str;
        }
    };

    class WFA2Wrapper {
        WFAlignerGapAffine m_aligner;

        int m_mismatch;
        int m_gap_opening;
        int m_gap_extension;

        std::string seq1;
        std::string seq2;

    public:
        WFA2Wrapper(int mismatch, int gap_opening, int gap_extension) :
                m_mismatch(mismatch),
                m_gap_opening(gap_opening),
                m_gap_extension(gap_extension),
                m_aligner(mismatch, gap_opening, gap_extension, WFAligner::Alignment, WFAligner::MemoryHigh) {
        }

        WFA2Wrapper(const WFA2Wrapper &other) :
                m_mismatch(other.m_mismatch),
                m_gap_opening(other.m_gap_opening),
                m_gap_extension(other.m_gap_extension),
                m_aligner(other.m_mismatch, other.m_gap_opening, other.m_gap_extension, WFAligner::Alignment,
                          WFAligner::MemoryHigh) {
        };

        static double CigarANI(std::string cigar) {
            size_t matches = std::count_if(cigar.begin(), cigar.end(), [](char const &c) {
                return c == 'M';
            });
            size_t mismatches = std::count_if(cigar.begin(), cigar.end(), [](char const &c) {
                return c == 'X';
            });
            size_t indel = std::count_if(cigar.begin(), cigar.end(), [](char const &c) {
                return c == 'I' || c == 'D';
            });

            return (static_cast<double>(matches) / (matches + mismatches + indel));
        }


        static void GetAlignmentInfo(AlignmentInfo& info, std::string const& cigar) {
            size_t alignment_start = 0;
            size_t alignment_length = 0;
            size_t indel_count = 0;
            size_t match_count = 0;
            size_t mismatch_count = 0;

            info.compressed_cigar = "";

            char last_instruction = ' ';
            size_t instruction_counter = 0;

            // Jump to pos where alignment starts (ignore hard 'H' and softclipping 'S')
            while (cigar[alignment_start] == 'H' || cigar[alignment_start] == 'S')  {
                if (cigar[alignment_start] != last_instruction) {
                    if (instruction_counter > 0) {
                        info.compressed_cigar += std::to_string(instruction_counter) + last_instruction;
                    }
                    last_instruction = cigar[alignment_start];
                    instruction_counter = 1;
                } else {
                    instruction_counter++;
                }
                alignment_start++;
            }

            // Get remainder of operations
            if (instruction_counter > 0) {
                info.compressed_cigar += std::to_string(instruction_counter) + last_instruction;
            }

            // Get Alignment length
            for (auto i = alignment_start; i < cigar.length(); i++) {
                char c = cigar[i];

                // Compressed cigar
                if (cigar[i] != last_instruction) {
                    if (instruction_counter > 0) {
                        info.compressed_cigar += std::to_string(instruction_counter) + last_instruction;
                    }
                    last_instruction = cigar[i];
                    instruction_counter = 1;
                } else {
                    instruction_counter++;
                }
                // compressed cigar end

                bool insertion = c == 'I';
                bool deletion = c == 'D';

                indel_count += deletion || insertion;
                match_count += c == 'M';
                mismatch_count  += c == 'X';

                alignment_length += deletion;
                alignment_length -= (insertion ||  c == 'H' || c == 'S');
                alignment_length++;
            }

            // Get remainder of operations
            if (instruction_counter > 0) {
                info.compressed_cigar += std::to_string(instruction_counter) + last_instruction;
            }

            info.alignment_start = alignment_start;
            info.alignment_length = alignment_length;
            info.alignment_ani = (static_cast<double>(match_count) / (match_count + mismatch_count + indel_count));
        }

        static AlignmentInfo GetAlignmentInfo(std::string const& cigar) {
            AlignmentInfo info;
            GetAlignmentInfo(info, cigar);
            return info;
        }

        inline std::string Cigar() {
            return m_aligner.getAlignmentCigar();
        }

        void Alignment(std::string query, std::string ref, AlignmentResult &result) {
            this->seq1 = query;
            this->seq2 = ref;
            m_aligner.alignEnd2End(seq1, seq2);
        }

        void Alignment(std::string query, std::string ref) {
            this->seq1 = query;
            this->seq2 = ref;
            m_aligner.alignEnd2End(seq1, seq2);
        }

        void PrintCigar(std::ostream &os = std::cout) {
            std::string cigar = m_aligner.getAlignmentCigar();
            os << cigar << std::endl;
        }

        void PrintAlignmentScore(std::ostream &os = std::cout) {
            os << m_aligner.getAlignmentScore() << std::endl;
        }

        WFAlignerGapAffine &GetAligner() {
            return m_aligner;
        }

        void PrintAlignment(std::ostream &os = std::cout) {
            PrintCigar();
            PrintAlignmentScore();
            std::cout << seq1 << std::endl;
            std::cout << seq2 << std::endl;

            m_aligner.cigarPrintPretty(stdout, seq1, seq2);
        }
    };
}