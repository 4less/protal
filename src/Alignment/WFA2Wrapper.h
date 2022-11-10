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