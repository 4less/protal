//
// Created by joachim on 05/06/2020.
//

#pragma once

#include <string>
#include <chrono>
#include <iostream>

using namespace std::chrono;
using namespace std;

namespace protal {
    enum Time {
        days, hours, minutes, seconds, milliseconds, microseconds
    };

    class Benchmark {
        string name;
        long time_sum = 0;
        size_t samplings = 0;
        time_point<high_resolution_clock> start_time;

    public:
        Benchmark(string name) : name(name) {
            Start();
        }

        void Start() {
            start_time = high_resolution_clock::now();
        }

        void Stop() {
//            if (time_sum != 0) return;
            samplings++;
            auto stop_time = high_resolution_clock::now();
            auto duration = duration_cast<std::chrono::microseconds>(stop_time - start_time);
//            std::cout << "Duration count : " << duration.count() << std::endl;
            time_sum += duration.count();
//            std::cout << "time_sum : " << time_sum << std::endl;
        }

        uint64_t GetDuration(Time format) {
            if (time_sum == 0) {
                Stop();
            }

            switch (format) {
                case Time::microseconds:
                    return time_sum;
                case Time::milliseconds:
                    return time_sum / 1'000;
                case Time::seconds:
                    return time_sum / 1'000'000;
                case Time::minutes:
                    return time_sum / 60'000'000;
                case Time::hours:
                    return time_sum / 3'600'000'000;
                case Time::days:
                    return time_sum / 86'400'000'000;
            }
            return 0L;
        }

        void PrintResults() {
            if (time_sum == 0) {
                Stop();
            }

            uint hours = time_sum / 3600000000;
            time_sum -= hours * 3600000000;
            uint mins = time_sum / 60000000;
            time_sum -= mins * 60000000;
            uint secs = time_sum / 1000000;
            time_sum -= secs * 1000000;
            uint msecs = time_sum / 1000;

            std::cout << name << " took ";
            if (hours) std::cout << hours << "h ";
            if (mins) std::cout << mins << "m ";
            if (secs) std::cout << secs << "s ";
            if (msecs) std::cout << msecs << "ms";
            if (samplings) std::cout << " (" << std::to_string(samplings) << ")";
            std::cout << std::endl;
        }
    };
}
