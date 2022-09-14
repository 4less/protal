//
// Created by fritsche on 14/08/22.
//

#pragma once

#include "Options.h"
#include "SequenceUtils/SeqReader.h"
#include "SequenceUtils/KmerIterator.h"
#include "Statistics.h"
#include <iostream>
#include <fstream>
#include <omp.h>
#include "Constants.h"
#include "Hash/KmerLookup.h"
#include "GenomeLoader.h"
#include "AlignmentOutputHandler.h"

namespace protal::classify {

    template<typename KmerHandler, typename AnchorFinder, typename AlignmentHandler, typename OutputHandler, DebugLevel debug>
    requires KmerHandlerConcept<KmerHandler> && AnchorFinderConcept<AnchorFinder>  && AlignmentHandlerConcept<AlignmentHandler>
    static Statistics Run(std::string const& read_se_path, protal::Options const& options, AnchorFinder& anchor_finder_global, AlignmentHandler& alignment_handler_global, OutputHandler& output_handler_global, KmerHandler& kmer_handler_global) {//, GenomeLoader& loader, Seedmap& map) {
        std::cout << "Run classify " << read_se_path << std::endl;

        // Shared
        std::ifstream is(read_se_path, std::ios::in);
        size_t dummy = 0;

        // Set Thread Num
        omp_set_num_threads(options.GetThreads());

        Statistics statistics{};

        std::ofstream varkit_output(options.GetOutputFile());

#pragma omp parallel default(none) shared(std::cout, varkit_output, options, is, dummy, kmer_handler_global, statistics, anchor_finder_global, alignment_handler_global, output_handler_global)//, loader, map)
        {
            // Private variables
            FastxRecord record;

            OutputHandler output_handler(output_handler_global);

            // Extract variables from kmi_global
            KmerHandler kmer_handler(kmer_handler_global);
            AnchorFinder anchor_finder(anchor_finder_global);
            AlignmentHandler alignment_handler(alignment_handler_global);

            // IO
            SeqReader reader;
            Statistics thread_statistics;
            thread_statistics.thread_num = omp_get_thread_num();

            // Intermediate storage objects
            KmerList kmers;
            AlignmentAnchorList anchors;
            AlignmentResultList alignment_results;

            while (reader(is, record)) {
                // Implement logger
                if constexpr (debug == DEBUG_VERBOSE) {
                    std::cout << "Record: " << record.id << std::endl;
                }

                thread_statistics.reads++;

                // Clear intermediate storage objects
                kmers.clear();
                anchors.clear();
                alignment_results.clear();

                // Retrieve kmers
                kmer_handler(std::string_view(record.sequence), kmers);

                if constexpr(KmerStatisticsConcept<KmerHandler>) {
                    thread_statistics.kmers_total += kmer_handler.TotalKmers();
                }
                if constexpr(KmerStatisticsConcept<KmerHandler>) {
                    thread_statistics.kmers_accepted += kmer_handler.TotalMinimizers();
                }


                // Calculate Anchors
                anchor_finder(kmers, anchors);
                thread_statistics.total_anchors += anchors.size();

                // Do Alignment
                alignment_handler(anchors, alignment_results, record.sequence);

                thread_statistics.total_alignments += alignment_results.size();

                // Output alignments
                output_handler(alignment_results, record);
            }

#pragma omp critical(statistics)
            {
                thread_statistics.output_alignments = output_handler.alignments;
                statistics.Join(thread_statistics);
            }
        }
        is.close();

        return statistics;
    }




    template<typename KmerHandler, typename AnchorFinder, typename AlignmentHandler, typename OutputHandler, DebugLevel debug>
    requires KmerHandlerConcept<KmerHandler> && AnchorFinderConcept<AnchorFinder>  && AlignmentHandlerConcept<AlignmentHandler>
    static Statistics RunPairedEnd(SeqReaderPE& reader_global, protal::Options const& options, AnchorFinder& anchor_finder_global, AlignmentHandler& alignment_handler_global, OutputHandler& output_handler_global, KmerHandler& kmer_handler_global) {//, GenomeLoader& loader, Seedmap& map) {

        size_t dummy = 0;

        // Set Thread Num
        omp_set_num_threads(options.GetThreads());

        Statistics statistics{};

#pragma omp parallel default(none) shared(std::cout, reader_global, options, dummy, kmer_handler_global, statistics, anchor_finder_global, alignment_handler_global, output_handler_global)
        {
            // Private variables
            FastxRecord record1;
            FastxRecord record2;

            // TODO this is a fix to get varkit output as quickly as possible
            OutputHandler output_handler(output_handler_global);

            // Extract variables from kmi_global
            KmerHandler kmer_handler(kmer_handler_global);
            AnchorFinder anchor_finder(anchor_finder_global);
            AlignmentHandler alignment_handler(alignment_handler_global);

            // IO
            SeqReaderPE reader{reader_global};
            Statistics thread_statistics;
            thread_statistics.thread_num = omp_get_thread_num();

            // Intermediate storage objects
            KmerList kmers1;
            AlignmentAnchorList anchors1;
            AlignmentResultList alignment_results1;
            KmerList kmers2;
            AlignmentAnchorList anchors2;
            AlignmentResultList alignment_results2;

            while (reader(record1, record2)) {
                thread_statistics.reads++;

                // Clear intermediate storage objects
                kmers1.clear();
                kmers2.clear();
                anchors1.clear();
                anchors2.clear();
                alignment_results1.clear();
                alignment_results2.clear();

                // Retrieve kmers
                kmer_handler(std::string_view(record1.sequence), kmers1);

                if constexpr(KmerStatisticsConcept<KmerHandler>) {
                    thread_statistics.kmers_total += kmer_handler.TotalKmers();
                }
                if constexpr(KmerStatisticsConcept<KmerHandler>) {
                    thread_statistics.kmers_accepted += kmer_handler.TotalMinimizers();
                }

                kmer_handler(std::string_view(record2.sequence), kmers2);

                if constexpr(KmerStatisticsConcept<KmerHandler>) {
                    thread_statistics.kmers_total += kmer_handler.TotalKmers();
                }
                if constexpr(KmerStatisticsConcept<KmerHandler>) {
                    thread_statistics.kmers_accepted += kmer_handler.TotalMinimizers();
                }

                // Calculate Anchors
                anchor_finder(kmers1, anchors1);
                anchor_finder(kmers2, anchors2);

                thread_statistics.total_anchors += anchors1.size();
                thread_statistics.total_anchors += anchors2.size();

                // Do Alignment
                alignment_handler(anchors1, alignment_results1, record1.sequence);
                alignment_handler(anchors2, alignment_results2, record2.sequence);

                thread_statistics.total_alignments += alignment_results1.size();
                thread_statistics.total_alignments += alignment_results2.size();

                // Output alignments
                output_handler(alignment_results1, record1);
                output_handler(alignment_results2, record2);

                if constexpr(debug == DEBUG_VERBOSE) {
#pragma omp critical(write)
                    {
                        thread_statistics.WriteStats(std::cout);
                    }
                }
                if constexpr(debug == DEBUG_EXTRAVERBOSE) {

                }
            }

#pragma omp critical(statistics)
            {
                thread_statistics.output_alignments = output_handler.alignments;
                statistics.Join(thread_statistics);
            }
        }

        return statistics;
    }

}