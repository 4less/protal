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
#include "ChainAnchorFinder.h"

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

            size_t record_id = omp_get_thread_num();

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
                anchor_finder(kmers, seeds, anchors, record.sequence.length());
                thread_statistics.total_anchors += anchors.size();

//                std::cout << "Anchor length: " << anchors.size() << std::endl;
//                for (auto& anchor : anchors) {
//                    std::cout << anchor.ToString() << std::endl;
//                }

                // Do Alignment
                bm_alignment.Start();
                alignment_handler(anchors, alignment_results, record.sequence);
                bm_alignment.Stop();

//                Utils::Input();
//                continue;

                thread_statistics.total_alignments += alignment_results.size();

                // Output alignments
                output_handler(alignment_results, record);

                if constexpr (benchmark_active) {
                    thread_core_benchmark(seeds, anchors, alignment_results, record.id);
                }

                record_id += options.GetThreads();
            }

#pragma omp critical(statistics)
            {
                anchor_finder_global.m_bm_seeding.Join(anchor_finder.m_bm_seeding);
                anchor_finder_global.m_bm_processing.Join(anchor_finder.m_bm_processing);
                anchor_finder_global.m_bm_pairing.Join(anchor_finder.m_bm_pairing);
                anchor_finder_global.m_bm_sorting_anchors.Join(anchor_finder.m_bm_sorting_anchors);

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
        anchor_finder_global.m_bm_sorting_anchors.PrintResults();
        bm_alignment_global.PrintResults();
        alignment_handler_global.bm_alignment.PrintResults();
        alignment_handler_global.bm_seedext.PrintResults();
        std::cout << alignment_handler_global.dummy << std::endl;
        std::cout << "----------------------------------------------------\n" << std::endl;


        return statistics;
    }

    static void SortAlignmentPairs1(PairedAlignmentResultList &pairs) {
        std::sort(pairs.begin(), pairs.end(), [](PairedAlignment const& a, PairedAlignment const& b) {
            return ScorePairedAlignment(a) > ScorePairedAlignment(b);
        });
    }

    static void SortAlignmentPairs2(PairedAlignmentResultList &pairs) {
        std::sort(pairs.begin(), pairs.end(), [](PairedAlignment const& a, PairedAlignment const& b) {
            return PairedAlignmentComparator(a, b);
        });
    }

    static void SortAlignmentPairs3(PairedAlignmentResultList &pairs) {
        std::sort(pairs.begin(), pairs.end(), [](PairedAlignment const& a, PairedAlignment const& b) {
            return Bitscore(a) > Bitscore(b);
        });
    }

    static void JoinAlignmentPairs(PairedAlignmentResultList &pairs, AlignmentResultList &read1, AlignmentResultList &read2, GenomeLoader& loader) {
//        std::cout << "Join Alignment Pairs -------------- " << std::endl;
//        std::cout << "1: " << read1.size() << ",  2: " << read2.size() << std::endl;
        std::vector<bool> selected2(read2.size(), false);
        for (auto& alignment1 : read1) {
//            std::cout << "Alignment1: " << alignment1.ToString() << std::endl;
            bool paired = false;
            size_t read2_index = 0;
            for (auto& alignment2 : read2) {
                if (alignment1.Taxid() == alignment2.Taxid()) {
                    pairs.emplace_back(PairedAlignment{ alignment1, alignment2 });
                    selected2[read2_index] = true;
                    paired = true;
//                    std::cout << "--> Alignment2: " << alignment2.ToString() << std::endl;
                }
                read2_index++;
            }
//            std::cout << "----" << std::endl;

            if (!paired) {
                pairs.emplace_back(PairedAlignment{ alignment1, AlignmentResult() });
            }
        }
        for (int i = 0; i < selected2.size(); i++) {
            if (!selected2[i]) {
                pairs.emplace_back(PairedAlignment{ AlignmentResult(), read2[i] });
            }
        }
    }

    template<typename KmerHandler, typename AnchorFinder, typename AlignmentHandler, typename OutputHandler, DebugLevel debug, typename AlignmentBenchmark=NoBenchmark>
    requires KmerHandlerConcept<KmerHandler> && AnchorFinderConcept<AnchorFinder>  && AlignmentHandlerConcept<AlignmentHandler>
    static Statistics RunPairedEnd(SeqReaderPE& reader_global, protal::Options const& options, AnchorFinder& anchor_finder_global, AlignmentHandler& alignment_handler_global, OutputHandler& output_handler_global, KmerHandler& kmer_handler_global, GenomeLoader& genome_loader, AlignmentBenchmark benchmark_global={}) {
        constexpr bool benchmark_active = !std::is_same<AlignmentBenchmark, NoBenchmark>();

        constexpr bool output_data = true;

        std::ofstream adoh_os(options.GetOutputPrefix() + ".adata", std::ios::out);
        ProtalAlignmentDataOutputHandler adoh_global(adoh_os, 1024*1024*16, 0.8);

        size_t dummy = 0;

        // Set Thread Num
        omp_set_num_threads(options.GetThreads());

        Statistics statistics{};
        Benchmark bm_alignment_global("Alignment handler");

#pragma omp parallel default(none) shared(std::cout, adoh_global, genome_loader, bm_alignment_global, benchmark_global, reader_global, options, dummy, kmer_handler_global, statistics, anchor_finder_global, alignment_handler_global, output_handler_global)
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

            // Alignment data output handler (set with flag)
            ProtalAlignmentDataOutputHandler adoh(adoh_global);

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

            PairedAlignmentResultList paired_alignment_results;

            size_t record_id = omp_get_thread_num();

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
                paired_alignment_results.clear();

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
                anchor_finder(kmers1, seeds1, anchors1, record1.sequence.length());
                anchor_finder(kmers2, seeds2, anchors2, record2.sequence.length());

                thread_statistics.total_anchors += anchors1.size();
                thread_statistics.total_anchors += anchors2.size();

                // Do Alignment
                bm_alignment.Start();
                alignment_handler(anchors1, alignment_results1, record1.sequence);
                alignment_handler(anchors2, alignment_results2, record2.sequence);
                bm_alignment.Stop();

                thread_statistics.total_alignments += alignment_results1.size();
                thread_statistics.total_alignments += alignment_results2.size();

                // Maybe get rid of this
                // This pairs SE alignments into all possible pairs.
                JoinAlignmentPairs(paired_alignment_results, alignment_results1,
                                   alignment_results2, genome_loader);

//                SortAlignmentPairs1(paired_alignment_results);
//                SortAlignmentPairs2(paired_alignment_results);
                SortAlignmentPairs3(paired_alignment_results);


                {
                    auto [taxonomic_id, gene_id] = KmerUtils::ExtractTaxIdGeneId(record1.id);

                    bool global_label = false;
                    size_t max_al = 5;
                    size_t al_count = 0;

                    auto best_taxid = paired_alignment_results.front().first.IsSet() ? paired_alignment_results.front().first.Taxid() : paired_alignment_results.front().second.Taxid();
                    bool best_true = best_taxid == taxonomic_id;
                    bool ambiguous_best = paired_alignment_results.size() > 1 &&
                                          Score(paired_alignment_results[0]) == Score(paired_alignment_results[1]);


                    if constexpr(false) {
                        if (!best_true && !ambiguous_best) {
                            std::cout << record1.id << std::endl;
                            for (auto &[ar1, ar2]: paired_alignment_results) {
                                if (ar1.IsSet()) std::cout << "1: " << ar1.ToString() << std::endl;
                                if (ar2.IsSet()) std::cout << "2: " << ar2.ToString() << std::endl;
                                std::cout << "Double score: " << Score({ar1, ar2}) << std::endl;
                                std::cout << "Bit score:    " << Bitscore({ar1, ar2}) << std::endl;
                                std::cout << "---" << std::endl;
                            }

                            std::cout << "First read anchors" << std::endl;
                            for (auto &a1: anchors1) {
                                std::cout << a1.ToString() << std::endl;
                            }

                            std::cout << "Second read anchors" << std::endl;
                            for (auto &a2: anchors2) {
                                std::cout << a2.ToString() << std::endl;
                            }

                            std::cout << "First read seeds" << std::endl;
                            for (auto &s1: seeds1) {
                                std::cout << s1.ToString() << std::endl;
                            }
                            std::cout << "Second read seeds" << std::endl;
                            for (auto &s2: seeds2) {
                                std::cout << s2.ToString() << std::endl;
                            }

                            Utils::Input();
                        }
                    }

                }

                // Output alignments
                output_handler(paired_alignment_results, record1, record2, record_id);
//                output_handler(alignment_results1, record1, record_id, true);
//                output_handler(alignment_results2, record2, record_id, false);

                if constexpr(output_data)
                    adoh(paired_alignment_results, record1, record2, record_id);

                if constexpr (benchmark_active) {
                    thread_core_benchmark(seeds1, anchors1, alignment_results1, record1.id);
                    thread_core_benchmark(seeds2, anchors2, alignment_results2, paired_alignment_results, record1.id);
                }

                if constexpr(debug == DEBUG_VERBOSE) {
#pragma omp critical(write)
                    {
                        thread_statistics.WriteStats(std::cout);
                    }
                }
                if constexpr(debug == DEBUG_EXTRAVERBOSE) {

                }
                record_id += options.GetThreads();
            }
#pragma omp critical(statistics)
            {
                anchor_finder_global.m_bm_seeding.Join(anchor_finder.m_bm_seeding);
                anchor_finder_global.m_bm_processing.Join(anchor_finder.m_bm_processing);
                anchor_finder_global.m_bm_pairing.Join(anchor_finder.m_bm_pairing);
                anchor_finder_global.m_bm_sorting_anchors.Join(anchor_finder.m_bm_sorting_anchors);

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

//        adoh_os.close();

        if constexpr (benchmark_active) {
            std::cout << "\n------------Alignment benchmarks------------------" << std::endl;
            benchmark_global.WriteRowStats();
            std::cout << "----------------------------------------------------\n" << std::endl;
        }

        std::cout << "---------------Speed benchmarks---------------------" << std::endl;
        anchor_finder_global.m_bm_seeding.PrintResults();
        anchor_finder_global.m_bm_processing.PrintResults();
        anchor_finder_global.m_bm_pairing.PrintResults();
        anchor_finder_global.m_bm_sorting_anchors.PrintResults();
        bm_alignment_global.PrintResults();
        alignment_handler_global.bm_alignment.PrintResults();
        alignment_handler_global.bm_seedext.PrintResults();
        std::cout << alignment_handler_global.dummy << std::endl;
        std::cout << "----------------------------------------------------\n" << std::endl;

        return statistics;
    }

}