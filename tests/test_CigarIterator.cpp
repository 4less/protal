#include <gtest/gtest.h>
#include "../src/Profiling/CigarIterator.h"

TEST(CigarIteratorTest, BasicIteration) {
    std::string cigar = "2M1D2I3X";
    CigarIterator it(cigar);

    // Expected: {cpos, qpos, rpos, op}
    std::vector<std::tuple<size_t, size_t, size_t, char>> expected = {
        {0, 0, 0, 'M'}, {1, 1, 1, 'M'},
        {2, 2, 2, 'D'},
        {3, 2, 3, 'I'}, {4, 3, 3, 'I'},
        {5, 4, 3, 'X'}, {6, 5, 4, 'X'}, {7, 6, 5, 'X'}
    };

    size_t idx = 0;

    for (auto [cpos, qpos, rpos, op] : it) {
        std::cout << "CIGAR pos: " << cpos << ", Query pos: " << qpos
                  << ", Ref pos: " << rpos << ", Operation: " << op
                  << " idx: " << idx << std::endl;
        ASSERT_LT(idx, expected.size());
        EXPECT_EQ(cpos, std::get<0>(expected[idx]));
        EXPECT_EQ(qpos, std::get<1>(expected[idx]));
        EXPECT_EQ(rpos, std::get<2>(expected[idx]));
        EXPECT_EQ(op,   std::get<3>(expected[idx]));
        ++idx;
    }
    EXPECT_EQ(idx, expected.size());
}