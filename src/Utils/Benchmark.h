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
        size_t threads_sampled = 1;
        time_point<high_resolution_clock> start_time;

    public:
        Benchmark(string name, size_t threads_sampled=1) : name(name), threads_sampled(threads_sampled) {
            //Start();
        }

        std::string GetName() const {
            return name;
        }

        void AddObservation() {
            samplings++;
        }

        void Start(bool new_sample=true) {
            samplings += new_sample;
            start_time = high_resolution_clock::now();
        }

        void Stop() {
//            if (time_sum != 0) return;
            auto stop_time = high_resolution_clock::now();
            auto duration = duration_cast<std::chrono::microseconds>(stop_time - start_time);
//            std::cout << "Duration count : " << duration.count() << std::endl;
            time_sum += duration.count();
        }

        uint64_t GetDuration(Time format) {
//            if (time_sum == 0) {
//                Stop();
//            }

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

        void Join(Benchmark const& other, bool add_thread=true) {
            if (add_thread && time_sum > 0) threads_sampled++;
            time_sum += other.time_sum;
        }

        void PrintResults() {
            if (time_sum == 0) {
                Stop();
            }

            if (threads_sampled > 1) {
                time_sum /= static_cast<double>(threads_sampled);
            }

            auto time_sum_local = time_sum;

            uint hours = time_sum_local / 3600000000;
            time_sum_local -= hours * 3600000000;
            uint mins = time_sum_local / 60000000;
            time_sum_local -= mins * 60000000;
            uint secs = time_sum_local / 1000000;
            time_sum_local -= secs * 1000000;
            uint msecs = time_sum_local / 1000;

            std::cout << name << " took ";
            if (hours) std::cout << hours << "h ";
            if (mins) std::cout << mins << "m ";
            if (secs) std::cout << secs << "s ";
            if (msecs) std::cout << msecs << "ms";
            if (!hours && !mins && !secs && !msecs) std::cout << " less than 0ms";
            if (samplings > 1) std::cout << " (" << std::to_string(samplings) << ")";
            if (threads_sampled > 1) std::cout << " mean over " << threads_sampled << " threads";
            std::cout << std::endl;
        }
    };
}
