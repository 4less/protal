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

    template<typename KmerIterator, typename KmerLookup, typename KmerProcessor, DebugLevel debug, typename ContextFreeMinimizer>
    requires KmerIteratorConcept<KmerIterator> && KmerLookupConcept<KmerLookup>
    static Statistics RunOld(std::string const& read_se_path, protal::Options const& options, KmerLookup& kmer_lookup_global, KmerProcessor& kmer_processor_global, KmerIterator& kmi_global, ContextFreeMinimizer& minimizer_global) {//, GenomeLoader& loader, Seedmap& map) {
        std::cout << "Run classify " << read_se_path << std::endl;

        // Shared
        std::ifstream is(read_se_path, std::ios::in);
        size_t dummy = 0;

        // Set Thread Num
        omp_set_num_threads(options.GetThreads());

        Statistics statistics{};

        std::ofstream varkit_output(options.GetOutputFile());


#pragma omp parallel default(none) shared(std::cout, varkit_output, options, is, dummy, kmer_lookup_global, kmi_global, statistics, minimizer_global, kmer_processor_global)//, loader, map)
    {
        // Private variables
        FastxRecord record;

        // TODO this is a fix to get varkit output as quickly as possible
        VarkitOutputHandler output_handler(varkit_output, 1024*128);

        // Extract variables from kmi_global
        KmerIterator kmi(kmi_global);
        ContextFreeMinimizer minimizer(minimizer_global);

        KmerProcessor processor(kmer_processor_global);

        // IO
        SeqReader reader;
        size_t kmer = 0;
        Statistics thread_statistics;
        thread_statistics.thread_num = omp_get_thread_num();

        KmerList kmers;
        AlignmentResultList alignment_results;

//        std::vector<std::pair<size_t, size_t>> minimizers;

//        LookupList first_anchor;
//        LookupList second_anchor;

        while (reader(is, record)) {
            kmi.SetSequence(const_cast<char *>(record.sequence.c_str()), record.sequence.size());

            if constexpr (debug == DEBUG_VERBOSE) {
                std::cout << "Record: " << record.id << std::endl;
            }

            size_t kmers_total = 0;
            size_t kmers_accepted = 0;

            thread_statistics.reads++;

            kmers.clear();
            alignment_results.clear();

            // Extract all k-mers
            while (kmi(kmer)) {
                kmers_total++;
                // Check if k-mer is minimizer
                // Concept dependant minimizers would not implement ContextFreeMinimizer
                // But would rather implement everthing in the KmerIterator
                if constexpr(protal::IsConceptFreeMinimizer<ContextFreeMinimizer>) {
                    if (!minimizer.IsMinimizer(kmer)) continue;
                }
                kmers.emplace_back( KmerElement(kmer, static_cast<size_t>(kmi.GetPos())) );
                kmers_accepted++;
            }

            thread_statistics.kmers_total += kmers_total;
            thread_statistics.kmers_accepted += kmers_accepted;


            // What to do with that?
            // How to further subdivide?
            // Where to pass output?
            processor(alignment_results, kmers, record);

            dummy += alignment_results.size();
            output_handler(alignment_results, record);

            // Maybe Loop?
            // while (result not good enough)


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
        statistics.Join(thread_statistics);


        std::cout << "total_alignments: " << processor.total_alignments << std::endl;
        std::cout << "alignments.... " << output_handler.alignments << std::endl;

#pragma omp critical(processor_benchmarks)
        processor.PrintSummary();
        processor.PrintBMs();
    }

        std::cout << "Dummy: " << dummy << std::endl;
        is.close();

        return statistics;
    }

    template<typename KmerHandler, typename AnchorFinder, typename AlignmentHandler, DebugLevel debug>
    requires KmerHandlerConcept<KmerHandler> && AnchorFinderConcept<AnchorFinder>  && AlignmentHandlerConcept<AlignmentHandler>
    static Statistics Run(std::string const& read_se_path, protal::Options const& options, AnchorFinder& anchor_finder_global, AlignmentHandler& alignment_handler_global, KmerHandler& kmer_handler_global) {//, GenomeLoader& loader, Seedmap& map) {
        std::cout << "Run classify " << read_se_path << std::endl;

        // Shared
        std::ifstream is(read_se_path, std::ios::in);
        size_t dummy = 0;

        // Set Thread Num
        omp_set_num_threads(options.GetThreads());

        Statistics statistics{};

        std::ofstream varkit_output(options.GetOutputFile());

#pragma omp parallel default(none) shared(std::cout, varkit_output, options, is, dummy, kmer_handler_global, statistics, anchor_finder_global, alignment_handler_global)//, loader, map)
        {
            // Private variables
            FastxRecord record;

            // TODO this is a fix to get varkit output as quickly as possible
            VarkitOutputHandler output_handler(varkit_output, 1024*128);

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

                // Do Alignment
                alignment_handler(anchors, alignment_results, record.sequence);

                statistics.total_alignments += alignment_results.size();

                // Output alignments
                output_handler(alignment_results, record);

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

                // Do Alignment
                alignment_handler(anchors1, alignment_results1, record1.sequence);
                alignment_handler(anchors2, alignment_results2, record2.sequence);

                statistics.total_alignments += alignment_results1.size();
                statistics.total_alignments += alignment_results2.size();

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