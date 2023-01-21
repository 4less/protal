//
// Created by fritsche on 24/11/22.
//

#pragma once

#include <cstddef>

namespace protal {
    class ProfilerOptions {
    public:
        // _hq is for high quality
        // Thresholds to rate
        double max_score_dist_to_second_hq = 0.02;
        double min_score_hq = 0.97;
        double max_score_dist_to_second_mq = 0.01;
        double min_score_mq = 0.95;
        double max_score_dist_to_second_lq = 0.005;
        double min_score_lq = 0.9;
        int min_mapq = 1;

        size_t min_alignment_length = 70;

        bool IsHQ(double score) {
            return score >= min_score_hq;
        }

        bool IsMQ(double score) {
            return score >= min_score_mq;
        }

        bool IsLQ(double score) {
            return score >= min_score_lq;
        }
    };
}