//
// Created by fritsche on 17/08/22.
//

#pragma once

#include <string>
#include <istream>
#include <iostream>
#include "wfa2-lib-2.3.4/bindings/cpp/WFAligner.hpp"
// #include "wfa2-lib-2.3.4/bindings/cpp/WFAligner.hpp"


using namespace::wfa;

namespace protal {
    class WFA2Wrapper2 {
        WFAlignerGapAffine m_aligner;

        int m_mismatch;
        int m_gap_opening;
        int m_gap_extension;
        size_t m_x_drop;

//        std::string_view seq1;
//        std::string_view seq2;


        WFAligner::AlignmentStatus m_status;
    public:
        std::string seq1;
        std::string seq2;
        WFA2Wrapper2(int mismatch, int gap_opening, int gap_extension, size_t x_drop) :
                m_mismatch(mismatch),
                m_gap_opening(gap_opening),
                m_gap_extension(gap_extension),
                m_x_drop(x_drop),
                m_aligner(mismatch, gap_opening, gap_extension, WFAligner::Alignment, WFAligner::MemoryHigh) {
        }

        WFA2Wrapper2(const WFA2Wrapper2 &other) :
                m_mismatch(other.m_mismatch),
                m_gap_opening(other.m_gap_opening),
                m_gap_extension(other.m_gap_extension),
                m_x_drop(other.m_x_drop),
                m_aligner(other.m_mismatch, other.m_gap_opening, other.m_gap_extension, WFAligner::Alignment,
                          WFAligner::MemoryHigh) {
        };

        inline std::string Cigar() {
            return m_aligner.getCIGAR(true);
        }

        void Reset() {
            m_status = WFAligner::AlignmentStatus::StatusMaxStepsReached;
        }

        bool Success() {
            return m_status == WFAligner::AlignmentStatus::StatusAlgCompleted;
        }

        void Alignment(std::string& query, std::string& ref, int max_score=INT32_MAX) {
//            std::cout << "---------------------- Alignment End2End Start" << std::endl;
//            std::cout << query << std::endl;
//            std::cout << ref << std::endl;
            seq1 = query;
            seq2 = ref;
//            m_aligner.setMaxAlignmentScore(INT32_MAX); //TODO
//            std::cout << "alignEnd2End()" << std::endl;
            m_status = m_aligner.alignEnd2End(ref, query);
            auto& cc = m_aligner.GetAligner()->cigar;
            std::cout << "CIG: " << std::string(cc->operations + cc->begin_offset, cc->end_offset - cc->begin_offset) << std::endl;

            if (Success()) {
#pragma omp critical(debug_out)
                {
                    std::string a = "TAGGTAAGCGCTACTCGGAGTTGGAGCCGATCAATAACGTGGACCGTGATTTGGTGTCTGCCCAGGAAGATCTTGAGGCACCCCGGGAAATGGCACATGAGGACCATGAGTTACAGGCCGAGGCTGAGCGCCTTGAGGTGGAAGTTGTGG";
                    std::string b = "TAGGTAAGCGCTACTCGGAGTTGCAGCCGATCATTAACGTGCACCGTGATTTGGTGTCTGCCCAGGAAGATCTTGAGGCAGCCCGTGAAATGGCACATGAGGACCATGAGTTCCAGGCCGAGGCTGAGCGCCTTGAGGTGGAAGTTGTGG";
                    std::cout << "Alignment End2End -> " << m_status << std::endl;
                    std::cout << a << std::endl;
                    std::cout << b << std::endl;
                    std::cout << m_aligner.getCIGAR(true) << std::endl;
                    PrintAlignment(std::string_view(query), std::string_view(ref));
                    std::cout << std::flush << std::endl;

                    WFAlignerGapAffine aligner(4, 2, 6, WFAligner::Alignment, WFAligner::MemoryHigh);
                    std::cout << "Align ends free" << std::endl;
                    aligner.alignEndsFree(a, 0, 0,
                                          b, 0, 0);
                    auto& cigar = aligner.GetAligner()->cigar;
                    std::cout << "CIG: " << std::string(cigar->operations + cigar->begin_offset, cigar->end_offset - cigar->begin_offset) << std::endl;
                    if (aligner.getAlignmentStatus() == 0) {
                        aligner.printPretty(stdout, a.c_str(), a.size(), b.c_str(), b.size());
                        fflush(stdout);
                    }
                    std::cout << "Align end2end" << std::endl;
                    aligner.alignEnd2End(query, ref);
                    std::cout << "CIG: " << std::string(cigar->operations + cigar->begin_offset, cigar->end_offset - cigar->begin_offset) << std::endl;
                    if (aligner.getAlignmentStatus() == 0) {
                        aligner.printPretty(stdout, a.c_str(), a.size(), b.c_str(), b.size());
                        fflush(stdout);
                    }
                    std::cout << "Old aligner but new strings" << std::endl;
                    m_aligner.alignEnd2End(a, b);
                    auto& cig = m_aligner.GetAligner()->cigar;
                    std::cout << "lastCIG: " << std::string(cig->operations + cig->begin_offset, cig->end_offset - cig->begin_offset) << std::endl;
                    if (m_aligner.getAlignmentStatus() == 0) {
                        aligner.printPretty(stdout, a.c_str(), a.size(), b.c_str(), b.size());
                        fflush(stdout);
                    }
                    std::cout << "-----" << std::endl;
                }
            }
            std::cout << "Alignment End2End End -----------------------" << std::endl;
        }

//        void Alignment(std::string_view query, std::string_view ref, int max_score=INT32_MAX) {
//            seq1 = query;
//            seq2 = ref;
//            m_aligner.setMaxAlignmentScore(max_score);
//            m_status = m_aligner.alignEnd2End(ref, query);
//        }

//        void Alignment(const char* query, size_t query_length, const char* ref, size_t ref_length, int max_score=INT32_MAX) {
//            m_aligner.setMaxAlignmentScore(max_score);
//            m_status = m_aligner.alignEnd2End(query, query_length, ref, ref_length);
//        }

        void Alignment(std::string query, std::string ref, int q_begin_free, int q_end_free, int r_begin_free, int r_end_free, int max_score=INT32_MAX) {
            seq1 = query;
            seq2 = ref;
//            m_aligner.setMaxAlignmentScore(max_score); //TODO
            m_status = m_aligner.alignEndsFree(ref, r_begin_free, r_end_free, query, q_begin_free, q_end_free);
//            if (Success()) {
//#pragma omp critical(debug_out)
//                std::cerr << "Print Alignment" << std::endl;
//                PrintAlignmentErr();
//                std::cerr << std::flush << std::endl;
//            }
        }

        void Alignment(std::string_view &query, std::string_view &ref, int q_begin_free, int q_end_free, int r_begin_free, int r_end_free, int max_score=INT32_MAX) {
            seq1 = query;
            seq2 = ref;
//            m_aligner.setMaxAlignmentScore(max_score); //TODO
            m_status = m_aligner.alignEndsFree(ref.cbegin(), ref.length(), r_begin_free, r_end_free, query.cbegin(), query.length(), q_begin_free, q_end_free);
        }

        void Alignment(const char* query, size_t query_len, const char* ref, size_t ref_length, int q_begin_free, int q_end_free, int r_begin_free, int r_end_free, int max_score=INT32_MAX) {
//            m_aligner.setMaxAlignmentScore(max_score); //TODO
            m_status = m_aligner.alignEndsFree(query, query_len, q_begin_free, q_end_free, ref, ref_length, r_begin_free, r_end_free);
//            m_status = m_aligner.alignEndsFree(query, q_begin_free, q_end_free, ref, ref_length, r_begin_free, r_end_free);
        }

        void PrintCigar(std::ostream &os = std::cout) {
            std::string cigar = m_aligner.getCIGAR(true);
            os << cigar << std::endl;
        }

        void PrintAlignmentScore(std::ostream &os = std::cout) {
            os << m_aligner.getAlignmentScore() << std::endl;
        }

        WFAlignerGapAffine &GetAligner() {
            return m_aligner;
        }

        int GetAlignmentScore() {
            return m_aligner.getAlignmentScore();
        }

        std::string GetAlignmentCigar() {
            return m_aligner.getCIGAR(true);
        }

        void PrintAlignment() {
            PrintCigar(std::cout);
            PrintAlignmentScore(std::cout);
            std::cout << "Seq1: " << seq1 << std::endl;
            std::cout << "Seq2: " << seq2 << std::endl;
            if (seq1.empty() || seq2.empty()) {
                std::cerr << "seq1 or seq2 not initialized (string_view in WFA2Wrapper)";
            }

//            m_aligner.cigarPrintPretty(stdout, seq1, seq2);
            m_aligner.printPretty(stdout, seq1.c_str(), seq1.size(), seq2.c_str(), seq2.size());
        }
        void PrintAlignment(std::string_view& query, std::string_view& ref) {
            PrintCigar(std::cout);
            PrintAlignmentScore(std::cout);
            std::cout << "Seq1: " << query << std::endl;
            std::cout << "Seq2: " << ref << std::endl;
            if (query.empty() || ref.empty()) {
                std::cerr << "seq1 or seq2 not initialized (string_view in WFA2Wrapper)";
            }

//            m_aligner.cigarPrintPretty(stdout, seq1, seq2);
            m_aligner.printPretty(stdout, seq1.c_str(), seq1.size(), seq2.c_str(), seq2.size());
        }
        void PrintAlignment(std::string_view&& query, std::string_view&& ref) {
            std::cout << "PRINT START -------------------_" << std::endl;
            PrintCigar(std::cout);
            PrintAlignmentScore(std::cout);
            std::cout << "Seq1: " << query << std::endl;
            std::cout << "Seq2: " << ref << std::endl;
            if (query.empty() || ref.empty()) {
                std::cerr << "seq1 or seq2 not initialized (string_view in WFA2Wrapper)";
            }

//            m_aligner.cigarPrintPretty(stdout, seq1, seq2);
            m_aligner.printPretty(stdout, ref.cbegin(), ref.size(), query.cbegin(), query.size());
            fflush(stdout);
            std::cout << "PRINT END -------------------_" << std::endl;
        }

        void PrintAlignmentErr() {
            PrintCigar(std::cerr);
            PrintAlignmentScore(std::cerr);
            std::cerr << "Seq1: " << seq1 << std::endl;
            std::cerr << "Seq2: " << seq2 << std::endl;
            if (seq1.empty() || seq2.empty()) {
                std::cerr << "seq1 or seq2 not initialized (string_view in WFA2Wrapper)";
            }

//            m_aligner.cigarPrintPretty(stdout, seq1, seq2);
            m_aligner.printPretty(stderr, seq2.c_str(), seq2.size(), seq1.c_str(), seq1.size());
        }
    };
}