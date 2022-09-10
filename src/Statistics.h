//
// Created by fritsche on 17/08/22.
//

#pragma once

#include <cstddef>
#include <fstream>
#include <iostream>

namespace protal {
    struct Statistics {
        size_t reads = 0;
        size_t kmers_accepted = 0;
        size_t kmers_total = 0;
        int thread_num = -1;

        void Join(Statistics& statistics) {
            reads += statistics.reads;
            kmers_total += statistics.kmers_total;
            kmers_accepted += statistics.kmers_accepted;
        }

        void WriteStats(std::ostream& os=std::cout) {
            if (thread_num != -1)
                os << "Thread num:      " << thread_num << std::endl;

            os << "Processed reads: \t" << reads << std::endl;
            os << "Total kmers:     \t" << kmers_total << std::endl;
            os << "Accepted kmers:  \t" << kmers_accepted << std::endl;
            double compression = (double) kmers_total/kmers_accepted;
            os << "Compression:     \t" << compression << std::endl;
        }
    };
}