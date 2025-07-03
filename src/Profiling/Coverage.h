// #pragma once

// #include <vector>
// #include <utility>
// #include <string>
// #include "CigarIterator.h"
// #include "SamHandler.h"

// using SeqRange = std::pair<size_t, size_t>;
// using SeqRangeList = std::vector<SeqRange>;

// using SamEntry = protal::SamEntry;
// using POS_t = protal::POS_t;
// using CIGAR_t = protal::CIGAR_t;

// class Coverage {
// private:
//     size_t m_gene_length = 0;
//     SeqRangeList m_ranges;

//     static SeqRangeList GetRangesFromCigar(CIGAR_t const& cigar, POS_t pos);
//     static void MergeRanges(SeqRangeList& a, const SeqRangeList& b);

//     friend class CoverageTest_GetRangesFromCigar_Test; //Testing
//     friend class CoverageTest_MergeRanges_Test; //Testing
// public:

//     Coverage(size_t gene_length) : m_gene_length(gene_length) {}

//     void AddSam(SamEntry const& sam);

//     void AddRange(CIGAR_t const& cigar, POS_t pos);

//     const SeqRangeList& GetRanges() const {
//         return m_ranges;
//     }
    
// };