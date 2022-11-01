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
        uint32_t read_start_offset = 0;
        uint32_t gene_start_offset = 0;
        uint32_t alignment_length = 0;
        float alignment_ani = 0.0f;

        uint16_t deletions = 0;
        uint16_t insertions = 0;

        std::string cigar;
        std::string compressed_cigar;

        std::string ToString() const {
            std::string str = std::to_string(read_start_offset);
            str += ", " + std::to_string(gene_start_offset);
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
        size_t m_x_drop;

        std::string_view seq1;
        std::string_view seq2;

        WFAligner::AlignmentStatus m_status;
    public:
        WFA2Wrapper(int mismatch, int gap_opening, int gap_extension, size_t x_drop) :
                m_mismatch(mismatch),
                m_gap_opening(gap_opening),
                m_gap_extension(gap_extension),
                m_x_drop(x_drop),
                m_aligner(mismatch, gap_opening, gap_extension, WFAligner::Alignment, WFAligner::MemoryHigh) {
        }

        WFA2Wrapper(const WFA2Wrapper &other) :
                m_mismatch(other.m_mismatch),
                m_gap_opening(other.m_gap_opening),
                m_gap_extension(other.m_gap_extension),
                m_x_drop(other.m_x_drop),
                m_aligner(other.m_mismatch, other.m_gap_opening, other.m_gap_extension, WFAligner::Alignment,
                          WFAligner::MemoryHigh) {
        };

        static double CigarANI(std::string cigar) {
            if (cigar.empty()) return 0;
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

        static double CompressedCigarANI(std::string cigar) {
            if (cigar.empty()) return 0;
            size_t count;
            char c;
            int digit_start = -1;

            size_t matches = 0;
            size_t mismatches = 0;
            size_t indel = 0;

            for (auto i = 0; i < cigar.length(); i++) {
                if (!std::isdigit(cigar[i])) {
                    count = stoi(cigar.substr(digit_start, i-digit_start));
                    c = cigar[i];
                    digit_start = -1;
                    if (c == 'M') matches += count;
                    if (c == 'X') mismatches += count;
                    if (c == 'I' || c == 'D') indel += count;
                } else if (digit_start == -1) {
                    digit_start = i;
                }
            }
            return (static_cast<double>(matches) / static_cast<double>(matches + mismatches + indel));
        }

        static int CompressedCigarScore(std::string cigar, int mismatch_pen = 4, int gapopen_pen = 6, int gapext_pen = 2) {
            if (cigar.empty()) return 0;

            int score = 0;
            int count;
            char c;
            int digit_start = -1;


            for (auto i = 0; i < cigar.length(); i++) {
                if (!std::isdigit(cigar[i])) {
                    count = stoi(cigar.substr(digit_start, i-digit_start));
                    c = cigar[i];
                    digit_start = -1;
                    if (c == 'X') score -= mismatch_pen * count;
                    if (c == 'I' || c == 'D') score -= (gapopen_pen + (count-1) * gapext_pen);
                } else if (digit_start == -1) {
                    digit_start = i;
                }
            }
            return score;
        }

        static int CigarScore(std::string cigar, int mismatch_pen = 4, int gapopen_pen = 6, int gapext_pen = 2) {
            int score = 0;
            char last = 'M';
            for (char c : cigar) {
                if (c == 'X') score -= 4;
                if (c == 'I' && last == 'I') score -= gapext_pen;
                if (c == 'I' && last != 'I') score -= gapopen_pen;
                if (c == 'D' && last == 'D') score -= gapext_pen;
                if (c == 'D' && last != 'D') score -= gapopen_pen;
                last = c;
            }
            return score;
        }


        static void GetAlignmentInfo(AlignmentInfo& info, std::string const& cigar, int dove_left_read=0, int dove_left_reference=0, int dove_right_read=0, int dove_right_reference=0) {
            size_t alignment_start = 0;
            size_t alignment_length = 0;
            size_t indel_count = 0;
            size_t match_count = 0;
            size_t mismatch_count = 0;

            info.compressed_cigar = "";
            info.deletions = 0;
            info.insertions = 0;

            int cigar_start_read = 0;
            while (cigar[cigar_start_read] == 'D' && cigar_start_read < dove_left_read)
                cigar_start_read++;

            int cigar_start_gene = 0;
            while (cigar[cigar_start_gene] == 'I' && cigar_start_gene < dove_left_reference)
                cigar_start_gene++;

            int cigar_end_read = cigar.length() - 1;
            while (cigar[cigar_end_read] == 'D' && cigar_end_read >= (cigar.length() - dove_right_read))
                cigar_end_read--;

            int cigar_end_gene = cigar.length() - 1;
            while (cigar[cigar_end_gene] == 'I' && cigar_end_gene >= (cigar.length() - dove_right_reference))
                cigar_end_gene--;

            cigar_end_read++;
            cigar_end_gene++;

            int cigar_start = std::max(cigar_start_read, cigar_start_gene);
            int cigar_end = std::max(cigar_end_read, cigar_end_gene);
            info.cigar = cigar.substr(cigar_start, cigar_end - cigar_start);

            // Get compressed cigar
            char last_instruction = ' ';
            size_t instruction_counter = 0;

            info.compressed_cigar = "";
            for (auto i = 0; i < info.cigar.length(); i++) {
                char c = info.cigar[i];
                // Compressed cigar
                if (info.cigar[i] != last_instruction) {
                    if (instruction_counter > 0) {
                        info.compressed_cigar += std::to_string(instruction_counter) + last_instruction;
                    }
                    last_instruction = info.cigar[i];
                    instruction_counter = 1;
                } else {
                    instruction_counter++;
                }
                // compressed cigar end

                bool insertion = c == 'I';
                bool deletion = c == 'D';
                info.insertions += insertion;
                info.deletions += deletion;

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

            info.read_start_offset = cigar_start_read;
            info.gene_start_offset = cigar_start_gene;

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

        bool Success() {
            return m_status == WFAligner::AlignmentStatus::StatusSuccessful;
        }

        void Alignment(std::string query, std::string ref, AlignmentResult &result) {
            m_aligner.setMaxAlignmentScore(INT32_MAX);
            m_status = m_aligner.alignEnd2End(query, ref);
        }

        void Alignment(std::string query, std::string ref) {
            m_aligner.setMaxAlignmentScore(INT32_MAX);
            m_status = m_aligner.alignEnd2End(query, ref);
        }

        void Alignment(std::string_view query, std::string_view ref, int max_score=INT32_MAX) {
            seq1 = query;
            seq2 = ref;
            m_aligner.setMaxAlignmentScore(max_score);
            m_status = m_aligner.alignEnd2End(query, ref);
        }

        void Alignment(std::string &query, std::string &ref, int q_begin_free, int q_end_free, int r_begin_free, int r_end_free, int max_score=INT32_MAX) {
            seq1 = query;
            seq2 = ref;
            m_aligner.setMaxAlignmentScore(max_score);
            m_status = m_aligner.alignEndsFree(query, q_begin_free, q_end_free, ref, r_begin_free, r_end_free);
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

        int GetAlignmentScore() const {
            return m_aligner.getAlignmentScore();
        }

        std::string GetAlignmentCigar() const {
            return m_aligner.getAlignmentCigar();
        }

        void PrintAlignment(std::ostream &os = std::cout) {
            PrintCigar();
            PrintAlignmentScore();
            std::cout << seq1 << std::endl;
            std::cout << seq2 << std::endl;
            if (seq1.empty() || seq2.empty()) {
                std::cerr << "seq1 or seq2 not initialized (string_view in WFA2Wrapper)";
            }

            m_aligner.cigarPrintPretty(stdout, seq1, seq2);
        }
    };
}