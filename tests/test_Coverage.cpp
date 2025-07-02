#include <gtest/gtest.h>
#include <string>
#include <sstream>
#include "../src/Profiling/Coverage.h"


// Define operator== for SeqRange to allow comparison in EXPECT_EQ
inline bool operator==(const SeqRange& lhs, const SeqRange& rhs) {
    return lhs.first == rhs.first && lhs.second == rhs.second;
}

std::string PrintRanges(SeqRangeList const& ranges) {
    std::ostringstream oss;
    for (size_t i = 0; i < ranges.size(); ++i) {
        oss << "[" << ranges[i].first << ", " << ranges[i].second << "]";
        if (i + 1 < ranges.size()) oss << ", ";
    }
    return oss.str();
}

TEST(CoverageTest, GetRangesFromCigar) {
    std::string cigar1 = "10X4D2I4D";
    POS_t pos1 = 10;

    // 
    

    auto ranges1 = Coverage::GetRangesFromCigar(cigar1, pos1);

    // Print m_ranges as a string formatted by "[start, end], [start, end]"

    // for (size_t i = 0; i < ranges.size(); ++i) {
    //     oss << "[" << ranges[i].first << ", " << ranges[i].second << "]";
    //     if (i + 1 < ranges.size()) oss << ", ";
    // }
    std::cout << "Ranges: " << PrintRanges(ranges1) << std::endl;

    // Expected: ranges before and after each 'D' block
    // 10X: [10,20), 4D: [20,20), 2I: [20,22), 4D: [22,22)
    // But the actual logic may vary depending on how m_ranges is used.
    // Let's check the number and values of the ranges:
    ASSERT_EQ(ranges1.size(), 1);
    EXPECT_EQ(ranges1[0].first, 10);
    EXPECT_EQ(ranges1[0].second, 20); // 10X: 10..19, end=20


    std::string cigar2 = "5M4D4M";
    auto ranges2 = Coverage::GetRangesFromCigar(cigar2, pos1);

    // 1         2
    // 0123456789012
    // ACGCTCGCTCGCT
    // ACGCT----CGCT

    std::cout << "Ranges: " << PrintRanges(ranges2) << std::endl;
    ASSERT_EQ(ranges2.size(), 2);
    EXPECT_EQ(ranges2[0].first, 10);
    EXPECT_EQ(ranges2[0].second, 15);
    EXPECT_EQ(ranges2[1].first, 19);
    EXPECT_EQ(ranges2[1].second, 23);


    std::string cigar3 = "2D2M2X2D2M3I5M";
    auto ranges3 = Coverage::GetRangesFromCigar(cigar3, 0);

    // 1         
    // 0123456789   012345 rpos
    // GCATTAATAC---GTAAT
    // --ATAT--ACTAGGTAAT
    // DDMMXXDDMMIIIMMMMM

    std::cout << "Ranges: " << PrintRanges(ranges3) << std::endl;
    ASSERT_EQ(ranges3.size(), 2);
    EXPECT_EQ(ranges3[0].first, 2);
    EXPECT_EQ(ranges3[0].second, 6);
    EXPECT_EQ(ranges3[1].first, 8);
    EXPECT_EQ(ranges3[1].second, 15);


    std::string cigar4 = "2D2M2X2D2M3I2D5M";
    auto ranges4 = Coverage::GetRangesFromCigar(cigar4, 0);

    // 1         
    // 0123456789   0123456 rpos
    // GCATTAATAC---CCGTAAT
    // --ATAT--ACTAG--GTAAT
    // DDMMXXDDMMIIIDDMMMMM

    std::cout << "Ranges: " << PrintRanges(ranges4) << std::endl;
    ASSERT_EQ(ranges4.size(), 3);
    EXPECT_EQ(ranges4[0].first, 2);
    EXPECT_EQ(ranges4[0].second, 6);
    EXPECT_EQ(ranges4[1].first, 8);
    EXPECT_EQ(ranges4[1].second, 10);
    EXPECT_EQ(ranges4[2].first, 12);
    EXPECT_EQ(ranges4[2].second, 17);
}


TEST(CoverageTest, MergeRanges) {
    auto range1 = SeqRangeList{ {2, 6}, {10, 100} };
    auto range2 = SeqRangeList{ {5, 9}, {12, 17}, {20, 23}, {50, 101} };
    auto range3 = SeqRangeList{ {9, 10}, {101, 103}, {1000,2000} };

    Coverage::MergeRanges(range1, range2);

    std::cout << "Merged Ranges: " << PrintRanges(range1) << std::endl;

    EXPECT_EQ(range1.size(), 2);
    EXPECT_EQ(range1[0], SeqRange(2, 9));
    EXPECT_EQ(range1[1], SeqRange(10, 101));


    Coverage::MergeRanges(range1, range3);
    std::cout << "Merged Ranges: " << PrintRanges(range1) << std::endl;

    EXPECT_EQ(range1.size(), 2);
    EXPECT_EQ(range1[0], SeqRange(2, 103));
    EXPECT_EQ(range1[1], SeqRange(1000, 2000));
}