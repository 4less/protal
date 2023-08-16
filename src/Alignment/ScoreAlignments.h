//
// Created by fritsche on 04/11/22.
//

#pragma once

#include "Constants.h"
#include "InternalReadAlignment.h"
#include "AlignmentUtils.h"

namespace protal {
    class ScoreAlignments {
    public:
        double Score(OptIRAPair const& pair, int mismatch_penalty=4) {
            return 0.0;
        }

        double Score(SamEntry& sam, int mismatch_penalty=4) {
            CigarInfo info;
            CompressedCigarInfo(sam.m_cigar, info);
            return CompressedCigarScore(sam.m_cigar) / (mismatch_penalty * info.clipped_alignment_length);
        }

        double Score(SamEntry& sam1, SamEntry& sam2, int mismatch_penalty=4) {
            CigarInfo info1;
            CigarInfo info2;
            CompressedCigarInfo(sam1.m_cigar, info1);
            CompressedCigarInfo(sam2.m_cigar, info2);

            return (CompressedCigarScore(sam1.m_cigar) + CompressedCigarScore(sam2.m_cigar)) /
                (mismatch_penalty * (info1.clipped_alignment_length * info2.clipped_alignment_length));
        }
    };
}