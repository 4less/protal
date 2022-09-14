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
        size_t total_alignments = 0;
        size_t output_alignments = 0;
        int thread_num = -1;
        size_t num_threads = 0;

        void Join(Statistics& statistics) {
            reads += statistics.reads;
            kmers_total += statistics.kmers_total;
            kmers_accepted += statistics.kmers_accepted;
            total_alignments += statistics.total_alignments;
            output_alignments += statistics.output_alignments;
            num_threads++;
        }

        void WriteStats(std::ostream& os=std::cout) {
            if (thread_num != -1)
                os << "Thread num:      " << thread_num << std::endl;

            os << "Processed reads:    \t" << reads << std::endl;
            os << "Total kmers:        \t" << kmers_total << std::endl;
            os << "Minimizers kmers:   \t" << kmers_accepted << std::endl;
            os << "Total alignments:   \t" << total_alignments << std::endl;
            os << "Output alignments:  \t" << output_alignments << std::endl;
            os << "Threads:            \t" << num_threads << std::endl;
//            double compression = (double) kmers_total/kmers_accepted;
//            os << "Compression:     \t" << compression << std::endl;
        }
    };
}