//
// Created by fritsche on 11/03/2021.
//

#pragma once

#include <cstdio>
#include <iostream>
#include <chrono>
#include "Utils.h"

using namespace std::chrono;


class ProgressBar {
    using TimePoints = std::vector<std::pair<time_point<high_resolution_clock>, size_t>>;

    size_t progress_size_ = 0;
    size_t progress_pos_ = 0;
    size_t steps_ = 100;
    size_t bar_width_ = 70;
    double progress_ = 0.0;
    int old_pos = -1;

    std::string bar = "";

    size_t last_length = 0;

    time_point<high_resolution_clock> last_time_point;

    size_t window_size = 10;
    TimePoints time_points;

    void print(int pos, double progress, double expected_runtime=0) {
        bar.clear();
        bar.reserve(150);
        bar = "[";

        for (int i = 0; i < bar_width_; ++i) {
            if (i < pos) bar += "=";
            else if (i == pos) bar += ">";
            else bar += " ";
        }
        std::string exp = Utils::FormatSeconds(std::max(0, int(expected_runtime)));
        bar += "] ";
        bar += std::to_string(int(progress * 100.0));
        bar += " %  ";
        bar += exp;
        bar += " left    \r";

        std::cout << bar;

        std::cout.flush();
    }

    size_t ExpectedRuntime() {
        size_t expected_runtime;

        auto& first = time_points[0];
        auto& last = time_points[time_points.size() - 1];

        auto duration = duration_cast<std::chrono::microseconds>(last.first - first.first);

        auto progress_diff = last.second - first.second;

        double progress_per_s = progress_diff/((duration.count()/1000000.f));
        auto remaining_progress = progress_size_ - last.second;

        if (progress_per_s == 0) return 0;
        expected_runtime = (remaining_progress/progress_per_s);

        return expected_runtime;
    }


public:
    ProgressBar (size_t progress_size) : progress_size_(progress_size) {
        bar.reserve(150);
        reset(progress_size_);
        time_points.reserve(window_size);
    };

    ProgressBar() = default;

    void reset(size_t new_progress_size) {
        progress_size_ = new_progress_size;
        progress_pos_ = 0;
        progress_ = 0.0;
        old_pos = -1;

        time_points.clear();
        time_points.push_back( { high_resolution_clock::now(), 0 } );

        print(bar_width_ * progress_, progress_, 0);
    }

    void UpdateAdd(size_t progress_pos) {
        Update(progress_pos_ + progress_pos);
    }

    void Update(size_t new_progress_pos) {

        time_point<high_resolution_clock> new_time_point = high_resolution_clock::now();

        // shift window
        if (time_points.size() == window_size)
            time_points.erase(time_points.begin());
        time_points.push_back( { high_resolution_clock::now(), new_progress_pos } );

        size_t expected_runtime = ExpectedRuntime();


        progress_pos_ = new_progress_pos;
        last_time_point = new_time_point;

        progress_ = (progress_pos_ < progress_size_) * ((double) progress_pos_ / progress_size_) + (progress_pos_ >= progress_size_);


        int pos = bar_width_ * progress_;

//        if (pos == old_pos) return;
//        if (pos == old_pos) return;
//        else {
            print(pos, progress_, expected_runtime);
            old_pos = pos;
//        }
    }
};