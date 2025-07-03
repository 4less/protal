// #include "Coverage.h"
// #include <algorithm> // for std::min, std::max

// SeqRangeList Coverage::GetRangesFromCigar(CIGAR_t const &cigar, POS_t pos)
// {
//     SeqRangeList ranges{};
//     CigarIterator it(cigar);

//     auto start = pos;
//     auto prev_rpos = 0;
//     bool in_del = false;

//     for (auto [cpos, qpos, rpos, op] : it) {
//         std::cout << "cpos: " << cpos << ", qpos:" << qpos << ", rpos:" << rpos << ", " << op << std::endl;
//         if ((op == 'D')) { // || op == 'I'
//             if (!in_del && start < pos + prev_rpos) {
//                 ranges.emplace_back(start, pos + rpos);
//             }
//             in_del = op == 'D';
//             start = pos + rpos;
//         } else if (in_del) {
//             start = pos + rpos;
//             in_del = false;
//         }
//         prev_rpos = rpos;
//     }
//     // Add the last range if needed
//     if (start < pos + prev_rpos) {
//         ranges.emplace_back(start, pos + prev_rpos + 1);
//     }

//     return ranges;
// }

// void Coverage::AddSam(SamEntry const &sam)
// {
//     Coverage::AddRange(sam.m_cigar, sam.m_pos);
// }

// void Coverage::AddRange(CIGAR_t const& cigar, POS_t pos) {
//     auto ranges = GetRangesFromCigar(cigar, pos);
//     MergeRanges(m_ranges, ranges);
// }


// void Coverage::MergeRanges(SeqRangeList& a, const SeqRangeList& b) {
//     // Append b to a
//     a.insert(a.end(), b.begin(), b.end());

//     // Sort by start position
//     std::sort(a.begin(), a.end());

//     // Merge overlapping ranges in a
//     size_t idx = 0;
//     for (size_t i = 1; i < a.size(); ++i) {
//         if (a[idx].second < a[i].first) {
//             // No overlap, move to next
//             ++idx;
//             a[idx] = a[i];
//         } else {
//             // Overlap, merge
//             a[idx].second = std::max(a[idx].second, a[i].second);
//         }
//     }
//     // Resize to keep only merged ranges
//     a.resize(idx + 1);
// }