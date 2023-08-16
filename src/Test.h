//
// Created by fritsche on 31/10/22.
//

#pragma once

#include "wfa2-lib/bindings/cpp/WFAligner.hpp"

using namespace::wfa;

namespace test {
    static void TestEndsFree() {
        // ACTAGCTAGCGACTTACTTTCTGCTCTTCTCTGCTGCTGCTGATGCTAGCTGCTAGCTGATCGATCGTAGCTATC
//        std::string query = "ACTAGCTAGCGACTTACTTTCTGCTCTTCTCTGCTGCTGCTGATGCTA";
//        std::string reference = "CTTACTATCTGCTCTTCTCTGCTGCTGCTGATGCTAGCTGCTAGCTGATCGATCGTAGCTATC";
        std::string query =     "ACTAGCTAGCGACTTACTTTCTGCTCTTCTCTG";     // PATTERN
        std::string reference =    "AGCTAGCGACTTACTTTCTGCTCTTCTCTGCTG";  // TEXT
//        std::string query =     "AATTAATTTAAGTCTAGGCTACTTTCGGTACTTTGTTCTT";
//        std::string reference =    "AATTTAAGTCTAGGCTACTTTCGGTACTTTCTT";


        std::cout << "QUERY:  " << query << std::endl;
        std::cout << "REFERENCE: " << reference << std::endl;
        WFAlignerGapAffine aligner{ 4, 6, 2, WFAligner::Alignment, WFAligner::MemoryHigh };
        aligner.setVerbose(2);



        aligner.alignEndsFree(query, 3, 0, reference, 0, 3);


        std::cout << "Status: " << aligner.getAlignmentStatus() << std::endl;
        std::cout << "Score: " << aligner.getAlignmentScore() << std::endl;
        aligner.cigarPrintPretty(stdout, query, reference);
    }
}