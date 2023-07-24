//
// Created by fritsche on 10/03/23.
//

#pragma once

#include <vector>
#include <string>
#include <iostream>

namespace protal {
    class TaxonStatisticsOutput {
        using ColumnPair = std::pair<std::string, size_t>;
        static inline ColumnPair SAMPLEID = {"#SAMPLEID", 0};
        static inline ColumnPair VERTICALOV = {"VerticalCoverage", 2};
        static inline ColumnPair TOTALREADS = {"TotalReads", 3};
        static inline ColumnPair TOTALLENGTH = {"TotalLength", 4};
        static inline ColumnPair MEANANI = {"MeanAni", 5};
        static inline ColumnPair MEANMAPQ = {"MeanMAPQ", 6};
        static inline ColumnPair ACCEPTED = {"Accepted", 1};
        std::vector<ColumnPair> m_header { VERTICALOV, TOTALREADS, TOTALLENGTH, MEANANI, MEANMAPQ, ACCEPTED };
        std::string sep= "\t";

    public:
        void PrintHeader(std::ostream& os) {
            std::sort(m_header.begin(), m_header.end(), [](ColumnPair const& a, ColumnPair const& b) {
               return a.second < b.second;
            });

            os << SAMPLEID.first;
            for (auto& [name, pos] : m_header) {
                os << sep << name;
            }
        }

        void PrintLine(std::ostream& os, std::string sample_name, double vertical_coverage, size_t total_reads,
                       size_t total_length, double mean_ani, double mean_mapq, bool accepted) {
            os << '\n' << sample_name;
            for (auto& cpair : m_header) {
                auto& [name, pos] = cpair;

                if (cpair == VERTICALOV) {
                    os << sep << vertical_coverage;
                } else if (cpair == TOTALREADS) {
                    os << sep << total_reads;
                } else if (cpair == TOTALLENGTH) {
                    os << sep << total_length;
                } else if (cpair == MEANANI) {
                    os << sep << mean_ani;
                } else if (cpair == MEANMAPQ) {
                    os << sep << mean_mapq;
                } else if (cpair == ACCEPTED) {
                    os << sep << accepted;
                }
            }
        }
    };
}