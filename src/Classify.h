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
#include "CoreBenchmark.h"

namespace protal::classify {

    template<typename KmerHandler, typename AnchorFinder, typename AlignmentHandler, typename OutputHandler, DebugLevel debug, typename AlignmentBenchmark=NoBenchmark>
    requires KmerHandlerConcept<KmerHandler> && AnchorFinderConcept<AnchorFinder>  && AlignmentHandlerConcept<AlignmentHandler>
    static Statistics
    Run(SeqReader& reader_global, protal::Options const& options, AnchorFinder& anchor_finder_global, AlignmentHandler& alignment_handler_global, OutputHandler& output_handler_global, KmerHandler& kmer_handler_global, AlignmentBenchmark benchmark_global={}) {
        constexpr bool benchmark_active = !std::is_same<AlignmentBenchmark, NoBenchmark>();

        // Shared
        size_t dummy = 0;

        // Set Thread Num
        omp_set_num_threads(options.GetThreads());

        Statistics statistics{};
        Benchmark bm_alignment_global{"Alignment handler"};

#pragma omp parallel default(none) shared(std::cout, bm_alignment_global, benchmark_global, options, dummy, kmer_handler_global, statistics, anchor_finder_global, alignment_handler_global, output_handler_global, reader_global)//, loader, map)
        {
            // Private variables
            FastxRecord record;


            OutputHandler output_handler(output_handler_global);

            // Extract variables from kmi_global
            KmerHandler kmer_handler(kmer_handler_global);
            AnchorFinder anchor_finder(anchor_finder_global);
            AlignmentHandler alignment_handler(alignment_handler_global);

            // IO
            SeqReader reader{ reader_global };

            Statistics thread_statistics;
            thread_statistics.thread_num = omp_get_thread_num();
            AlignmentBenchmark thread_core_benchmark{ benchmark_global };

            // Intermediate storage objects
            KmerList kmers;
            SeedList seeds;
            AlignmentAnchorList anchors;
            AlignmentResultList alignment_results;

            Benchmark bm_alignment{"Alignment handler"};

            while (reader(record)) {
                // Implement logger
                thread_statistics.reads++;

                // Clear intermediate storage objects
                kmers.clear();
                seeds.clear();
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
//                std::cout << record.id << std::endl;
                anchor_finder(kmers, seeds, anchors);
                thread_statistics.total_anchors += anchors.size();

                // Do Alignment
                bm_alignment.Start();
                alignment_handler(anchors, alignment_results, record.sequence);
                bm_alignment.Stop();

                thread_statistics.total_alignments += alignment_results.size();

                // Output alignments
                output_handler(alignment_results, record);

                if constexpr (benchmark_active) {
                    thread_core_benchmark(seeds, anchors, alignment_results, record.id);
                }
            }

#pragma omp critical(statistics)
            {
                anchor_finder_global.m_bm_seeding.Join(anchor_finder.m_bm_seeding);
                anchor_finder_global.m_bm_processing.Join(anchor_finder.m_bm_processing);
                anchor_finder_global.m_bm_pairing.Join(anchor_finder.m_bm_pairing);

                bm_alignment_global.Join(bm_alignment);
                alignment_handler_global.bm_alignment.Join(alignment_handler.bm_alignment);
                alignment_handler_global.bm_seedext.Join(alignment_handler.bm_seedext);
                alignment_handler_global.dummy += alignment_handler.dummy;

                thread_statistics.output_alignments = output_handler.alignments;
                statistics.Join(thread_statistics);

                if constexpr (benchmark_active) {
                    benchmark_global.Join(thread_core_benchmark);
                }



                std::cout << "Tail: " << alignment_handler.total_tail_alignments << std::endl;
                std::cout << "TLen: " << alignment_handler.total_tail_length << std::endl;
                std::cout << "rate: " << (static_cast<double>(alignment_handler.total_tail_length)/static_cast<double>(alignment_handler.total_tail_alignments)) << std::endl;
            }
        }

        if constexpr (benchmark_active) {
            std::cout << "\n------------Alignment benchmarks------------------" << std::endl;
            benchmark_global.WriteRowStats();
            std::cout << "----------------------------------------------------\n" << std::endl;
        }

        std::cout << "---------------Speed benchmarks---------------------" << std::endl;
        anchor_finder_global.m_bm_seeding.PrintResults();
        anchor_finder_global.m_bm_processing.PrintResults();
        anchor_finder_global.m_bm_pairing.PrintResults();
        bm_alignment_global.PrintResults();
        alignment_handler_global.bm_alignment.PrintResults();
        alignment_handler_global.bm_seedext.PrintResults();
        std::cout << alignment_handler_global.dummy << std::endl;
        std::cout << "----------------------------------------------------\n" << std::endl;


        return statistics;
    }


    template<typename KmerHandler, typename AnchorFinder, typename AlignmentHandler, typename OutputHandler, DebugLevel debug, typename AlignmentBenchmark=NoBenchmark>
    requires KmerHandlerConcept<KmerHandler> && AnchorFinderConcept<AnchorFinder>  && AlignmentHandlerConcept<AlignmentHandler>
    static Statistics RunPairedEnd(SeqReaderPE& reader_global, protal::Options const& options, AnchorFinder& anchor_finder_global, AlignmentHandler& alignment_handler_global, OutputHandler& output_handler_global, KmerHandler& kmer_handler_global, AlignmentBenchmark benchmark_global={}) {
        constexpr bool benchmark_active = !std::is_same<AlignmentBenchmark, NoBenchmark>();

        size_t dummy = 0;

        // Set Thread Num
        omp_set_num_threads(options.GetThreads());

        Statistics statistics{};
        Benchmark bm_alignment_global("Alignment handler");

#pragma omp parallel default(none) shared(std::cout, bm_alignment_global, benchmark_global, reader_global, options, dummy, kmer_handler_global, statistics, anchor_finder_global, alignment_handler_global, output_handler_global)
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
            AlignmentBenchmark thread_core_benchmark{ benchmark_global };

            // Intermediate storage objects
            KmerList kmers1;
            SeedList seeds1;
            AlignmentAnchorList anchors1;
            AlignmentResultList alignment_results1;
            KmerList kmers2;
            SeedList seeds2;
            AlignmentAnchorList anchors2;
            AlignmentResultList alignment_results2;

            Benchmark bm_alignment{"Alignment handler"};

            while (reader(record1, record2)) {
                thread_statistics.reads++;

                // Clear intermediate storage objects
                kmers1.clear();
                kmers2.clear();
                seeds1.clear();
                seeds2.clear();
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
                anchor_finder(kmers1, seeds1, anchors1);
                anchor_finder(kmers2, seeds2, anchors2);

                thread_statistics.total_anchors += anchors1.size();
                thread_statistics.total_anchors += anchors2.size();

                // Do Alignment
                bm_alignment.Start();
                alignment_handler(anchors1, alignment_results1, record1.sequence);
                alignment_handler(anchors2, alignment_results2, record2.sequence);
                bm_alignment.Stop();

                thread_statistics.total_alignments += alignment_results1.size();
                thread_statistics.total_alignments += alignment_results2.size();

                // Output alignments
                output_handler(alignment_results1, record1);
                output_handler(alignment_results2, record2);



                if constexpr (benchmark_active) {
                    thread_core_benchmark(seeds1, anchors1, alignment_results1, record1.id);
                    thread_core_benchmark(seeds2, anchors2, alignment_results2, record2.id);
                }

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
                anchor_finder_global.m_bm_seeding.Join(anchor_finder.m_bm_seeding);
                anchor_finder_global.m_bm_processing.Join(anchor_finder.m_bm_processing);
                anchor_finder_global.m_bm_pairing.Join(anchor_finder.m_bm_pairing);

                bm_alignment_global.Join(bm_alignment);
                alignment_handler_global.bm_alignment.Join(alignment_handler.bm_alignment);

                alignment_handler_global.bm_seedext.Join(alignment_handler.bm_seedext);
                alignment_handler_global.dummy += alignment_handler.dummy;

                thread_statistics.output_alignments = output_handler.alignments;
                statistics.Join(thread_statistics);

                if constexpr (benchmark_active) {
                    benchmark_global.Join(thread_core_benchmark);
                }
            }
        }

        if constexpr (benchmark_active) {
            std::cout << "\n------------Alignment benchmarks------------------" << std::endl;
            benchmark_global.WriteRowStats();
            std::cout << "----------------------------------------------------\n" << std::endl;
        }

        std::cout << "---------------Speed benchmarks---------------------" << std::endl;
        anchor_finder_global.m_bm_seeding.PrintResults();
        anchor_finder_global.m_bm_processing.PrintResults();
        anchor_finder_global.m_bm_pairing.PrintResults();
        bm_alignment_global.PrintResults();
        alignment_handler_global.bm_alignment.PrintResults();
        alignment_handler_global.bm_seedext.PrintResults();
        std::cout << alignment_handler_global.dummy << std::endl;
        std::cout << "----------------------------------------------------\n" << std::endl;

        return statistics;
    }

}