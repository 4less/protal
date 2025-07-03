//
// Created by fritsche on 10/10/23.
//

#pragma once

#include <cstdint>
#include "Constants.h"
#include <vector>
#include <string>
#include "SamHandler.h"
#include <optional>

namespace Profiler {
    struct AlignmentSpan;

    using ReadLength = uint32_t;

//    using GeneId = size_t;
    using TaxonID = uint32_t;
    using GeneID = uint32_t;
    using AlignmentID = size_t;
    using ReadID = size_t;
    using CoverageVector = std::vector<uint8_t>;
    using GeneIDVector = std::vector<GeneID>;
    using Alignment = protal::ReducedSam;
    using AlignmentPair = std::pair<std::optional<protal::ReducedSam>, std::optional<protal::ReducedSam>>;
    using AlignmentVector = std::vector<Alignment>;
    using AlignmentPairVector = std::vector<AlignmentPair>;
    using AlignmentIDVector = std::vector<AlignmentID>;
    using AlignmentSpanVector = std::vector<AlignmentSpan>;




    /***
     * AlignmentSpan for sorting alignments quickly.
     * Size 16 Byte (128 bits, 64 + 32 + 32)
     */
    struct AlignmentSpan {
        AlignmentID id = SIZE_MAX;
        protal::POS_t read_pos = UINT32_MAX;
        ReadLength alignment_length = UINT32_MAX;

        bool IsUninitialized() const {
            return id == SIZE_MAX;
        }
    };
}


