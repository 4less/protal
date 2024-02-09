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
#include "ScoreAlignments.h"
#include "ProgressBar.h"

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
        Benchmark bm_reader_global{"Sequence reader"};
        Benchmark bm_alignment_global{"Alignment handler"};
        Benchmark bm_omp_block{"OMP Loop handler"};
        Benchmark bm_omp_before_loop_global{ "OMP before loop" };

        bm_omp_block.Start();
#pragma omp parallel default(none) shared(std::cout, bm_alignment_global, bm_reader_global, bm_omp_before_loop_global, benchmark_global, options, dummy, kmer_handler_global, statistics, anchor_finder_global, alignment_handler_global, output_handler_global, reader_global)//, loader, map)
        {
            bm_omp_before_loop_global.Start();
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
            Benchmark bm_reader{"Reader"};

            size_t record_id = omp_get_thread_num();

            bm_omp_before_loop_global.Stop();
            bm_reader.Start();
            while (reader(record)) {
                bm_reader.Stop();

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
                anchor_finder(kmers, seeds, anchors, record.sequence);
                thread_statistics.total_anchors += anchors.size();

//                std::cout << "Anchor length: " << anchors.size() << std::endl;
//                for (auto& anchor : anchors) {
//                    std::cout << anchor.ToString() << std::endl;
//                }

                // Do Alignment
                bm_alignment.Start();
                alignment_handler(anchors, alignment_results, record.sequence, options.GetAlignTop(), record.id);
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
                bm_reader.Start();
            }
            bm_reader.Stop();


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
        bm_omp_block.Stop();

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
                if (alignment1.Taxid() == alignment2.Taxid() && alignment1.GeneId() == alignment2.GeneId()) {
                    if (!CorrectOrientation(alignment1, alignment2)) {
                        std::cerr << "Discard wrong orientation: " << std::endl;
                        std::cerr << alignment1.ToString() << std::endl;
                        std::cerr << alignment2.ToString() << std::endl;
                        continue;
                    }

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

//        std::ofstream adoh_os(options.GetOutputPrefix() + ".adata", std::ios::out);
//        ProtalAlignmentDataOutputHandler adoh_global(adoh_os, 1024*1024*16, 0.8);

        size_t dummy = 0;

        // Set Thread Num
        omp_set_num_threads(options.GetThreads());

        Statistics statistics{};
        Benchmark bm_anchor_finder_global{"Seed- and Anchor-finding", 0};
        Benchmark bm_anchor_recovery_global{"Anchor recovery", 0};
        Benchmark bm_alignment_global{"Alignment handler", 0};
        Benchmark bm_alignment_join_sort_global{"Joining alignment pairs and sorting", 0};
        Benchmark bm_output_global{"Output handler", 0};

        Utils::Histogram seed_sizes_global;
        Utils::Histogram anchor_sizes_global;


        Benchmark bm_omp_block{"OMP Loop handler"};
        Benchmark bm_reader_global{"Sequence reader"};
        Benchmark bm_omp_before_loop_global{ "OMP before loop" };

        bm_omp_block.Start();
#pragma omp parallel default(none) shared(std::cout, bm_reader_global, bm_omp_before_loop_global, bm_alignment_join_sort_global, bm_anchor_recovery_global, bm_anchor_finder_global, /*adoh_global,*/ genome_loader, seed_sizes_global, anchor_sizes_global, bm_alignment_global, bm_output_global, benchmark_global, reader_global, options, dummy, kmer_handler_global, statistics, anchor_finder_global, alignment_handler_global, output_handler_global)
        {
            bm_omp_before_loop_global.Start();
            // Private variables
            FastxRecord record1;
            FastxRecord record2;

            // TODO this is a fix to get varkit output as quickly as possible
            OutputHandler output_handler(output_handler_global);

            // Extract variables from kmi_globalF
            KmerHandler kmer_handler(kmer_handler_global);
            AnchorFinder anchor_finder1(anchor_finder_global);
            AnchorFinder anchor_finder2(anchor_finder_global);
            AlignmentHandler alignment_handler(alignment_handler_global);

            // Alignment data output handler (set with flag)
            //            ProtalAlignmentDataOutputHandler adoh(adoh_global);

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

            Benchmark bm_reader{"Sequence reader"};
            Benchmark bm_kmer_extracter{"Retrieve k-mers"};
            Benchmark bm_anchor_finder{"Seed- and Anchor-finding"};
            Benchmark bm_anchor_recovery{"Anchor recovery"};
            Benchmark bm_alignment{"Alignment handler"};
            Benchmark bm_alignment_join_sort{"Joining alignment pairs and sorting"};
            Benchmark bm_output{"Output handler"};

            Utils::Histogram seed_sizes;
            Utils::Histogram anchor_sizes;

            bm_omp_before_loop_global.Stop();
            bm_reader.Start();
            while (reader(record1, record2)) {
                bm_reader.Stop();
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
                bm_kmer_extracter.Start();
                kmer_handler(std::string_view(record1.sequence), kmers1);
                bm_kmer_extracter.Stop();


                if constexpr(KmerStatisticsConcept<KmerHandler>) {
                    thread_statistics.kmers_total += kmer_handler.TotalKmers();
                }
                if constexpr(KmerStatisticsConcept<KmerHandler>) {
                    thread_statistics.kmers_accepted += kmer_handler.TotalMinimizers();
                }

                bm_kmer_extracter.Start();
                kmer_handler(std::string_view(record2.sequence), kmers2);
                bm_kmer_extracter.Stop();

                if constexpr(KmerStatisticsConcept<KmerHandler>) {
                    thread_statistics.kmers_total += kmer_handler.TotalKmers();
                }
                if constexpr(KmerStatisticsConcept<KmerHandler>) {
                    thread_statistics.kmers_accepted += kmer_handler.TotalMinimizers();
                }

                // Calculate Anchors


                bm_anchor_finder.Start();
                anchor_finder1(kmers1, seeds1, anchors1, record1.sequence);
                anchor_finder2(kmers2, seeds2, anchors2, record2.sequence);
                bm_anchor_finder.Stop();

                if (!anchor_finder1.Success() || !anchor_finder2.Success()) {
                    thread_statistics.errors_anchor_finding++;
                }

                //                std::cout << "Anchors1" << std::endl;
                //                for (auto& a : anchors1) {
                //                    std::cout << a.ToString() << std::endl;
                //                }
                //                std::cout << "Anchors2" << std::endl;
                //                for (auto& a : anchors2) {
                //                    std::cout << a.ToString() << std::endl;
                //                }

                bm_anchor_recovery.Start();
                auto recover1 = anchor_finder1.RecoverAnchors(anchors1, anchor_finder2.BestAnchors(), options.GetAlignTop());
                auto recover2 = anchor_finder2.RecoverAnchors(anchors2, anchor_finder1.BestAnchors(), options.GetAlignTop());
                bm_anchor_recovery.Stop();

                //                std::cout << "Anchors1 after" << std::endl;
                //                for (auto& a : anchors1) {
                //                    std::cout << a.ToString() << std::endl;
                //                }
                //                std::cout << "Anchors2 after" << std::endl;
                //                for (auto& a : anchors2) {
                //                    std::cout << a.ToString() << std::endl;
                //                }
                //
                //                if (recover1 || recover2) {
                //                    Utils::Input();
                //                }

                seed_sizes.AddObservation(seeds1.size());
                seed_sizes.AddObservation(seeds2.size());
                anchor_sizes.AddObservation(anchors1.size());
                anchor_sizes.AddObservation(anchors2.size());

                thread_statistics.total_anchors += anchors1.size();
                thread_statistics.total_anchors += anchors2.size();

                // Do Alignment
                bm_alignment.Start();
                alignment_handler(anchors1, alignment_results1, record1.sequence, options.GetAlignTop() + recover1, record1.id);
                alignment_handler(anchors2, alignment_results2, record2.sequence, options.GetAlignTop() + recover2, record2.id);
                bm_alignment.Stop();

                thread_statistics.total_alignments += alignment_results1.size();
                thread_statistics.total_alignments += alignment_results2.size();

                bm_alignment_join_sort.Start();
                // This pairs SE alignments into all possible pairs.
                JoinAlignmentPairs(paired_alignment_results, alignment_results1,
                                   alignment_results2, genome_loader);

                //                SortAlignmentPairs1(paired_alignment_results);
                //                SortAlignmentPairs2(paired_alignment_results);
                SortAlignmentPairs3(paired_alignment_results);
                bm_alignment_join_sort.Stop();

                // Output alignments
                bm_output.Start();
                output_handler(paired_alignment_results, record1, record2, record_id);
                bm_output.Stop();

                if constexpr (benchmark_active) {
                    thread_core_benchmark(seeds1, anchors1, alignment_results1, record1.id);
                    thread_core_benchmark(seeds2, anchors2, alignment_results2, paired_alignment_results, record1.id);
                    //#pragma omp critical(write)
                    //                    thread_core_benchmark.ErrorOutput(seeds1, seeds2, anchors1, anchors2, alignment_results1, alignment_results2, paired_alignment_results, record1, record2);
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

                bm_reader.Start();
            }
            bm_reader.Stop();

#pragma omp critical(statistics)
            {
                anchor_finder_global.m_bm_seeding.Join(anchor_finder1.m_bm_seeding);
                anchor_finder_global.m_bm_seeding.Join(anchor_finder2.m_bm_seeding, false);
                anchor_finder_global.m_bm_reverse_complement.Join(anchor_finder1.m_bm_reverse_complement);
                anchor_finder_global.m_bm_reverse_complement.Join(anchor_finder2.m_bm_reverse_complement, false);
                anchor_finder_global.m_bm_operator.Join(anchor_finder1.m_bm_operator);
                anchor_finder_global.m_bm_operator.Join(anchor_finder2.m_bm_operator, false);
                anchor_finder_global.m_bm_seed_subsetting.Join(anchor_finder1.m_bm_seed_subsetting);
                anchor_finder_global.m_bm_seed_subsetting.Join(anchor_finder2.m_bm_seed_subsetting, false);
                anchor_finder_global.m_bm_processing.Join(anchor_finder1.m_bm_processing);
                anchor_finder_global.m_bm_processing.Join(anchor_finder2.m_bm_processing, false);
                anchor_finder_global.m_bm_pairing.Join(anchor_finder1.m_bm_pairing);
                anchor_finder_global.m_bm_pairing.Join(anchor_finder2.m_bm_pairing, false);
                anchor_finder_global.m_bm_sorting_anchors.Join(anchor_finder1.m_bm_sorting_anchors);
                anchor_finder_global.m_bm_sorting_anchors.Join(anchor_finder2.m_bm_sorting_anchors, false);
                anchor_finder_global.m_bm_recovering_anchors.Join(anchor_finder1.m_bm_recovering_anchors);
                anchor_finder_global.m_bm_recovering_anchors.Join(anchor_finder2.m_bm_recovering_anchors, false);
                anchor_finder_global.m_bm_extend_anchors.Join(anchor_finder1.m_bm_extend_anchors);
                anchor_finder_global.m_bm_extend_anchors.Join(anchor_finder2.m_bm_extend_anchors, false);
                anchor_finder_global.recovered_count += anchor_finder1.recovered_count;
                anchor_finder_global.recovered_count += anchor_finder2.recovered_count;
                anchor_finder_global.total_count += anchor_finder1.total_count;
                anchor_finder_global.total_count += anchor_finder2.total_count;

                reader_global.UpdateSuccess(reader);

                bm_reader_global.Join(bm_reader);

                bm_anchor_finder_global.Join(bm_anchor_finder);
                // bm_anchor_recovery_global.Join(bm_anchor_recovery);
                bm_alignment_global.Join(bm_alignment);
                bm_alignment_join_sort_global.Join(bm_alignment_join_sort);
                bm_output_global.Join(bm_output);
                alignment_handler_global.bm_alignment.Join(alignment_handler.bm_alignment);

                alignment_handler_global.m_bm_alignment.Join(alignment_handler.m_bm_alignment);

                thread_statistics.output_alignments = output_handler.alignments;
                statistics.Join(thread_statistics);

                seed_sizes_global.Join(seed_sizes);
                anchor_sizes_global.Join(anchor_sizes);

                if constexpr (benchmark_active) {
                    benchmark_global.Join(thread_core_benchmark);
                }
            }
        }
        bm_omp_block.Stop();

//        adoh_os.close();

        if constexpr (benchmark_active) {
            if (options.Verbose()) {
                std::cout << "\n------------Alignment benchmarks------------------" << std::endl;
                benchmark_global.WriteRowStats();
                std::cout << "----------------------------------------------------\n" << std::endl;
            }
            benchmark_global.WriteRowStatsToFile(options.GetPrefix(options.GetCurrentIndex()) + "_benchmark.tsv");
        }

        if (options.Verbose()) {
            std::cout << "---------------Speed benchmarks---------------------" << std::endl;
            bm_omp_block.PrintResults();
            bm_omp_before_loop_global.PrintResults();
            std::cout << "-----" << std::endl;
            bm_reader_global.PrintResults();
            bm_anchor_finder_global.PrintResults();
            std::cout << "\t";
            anchor_finder_global.m_bm_operator.PrintResults();
            std::cout << "\t\t";
            anchor_finder_global.m_bm_reverse_complement.PrintResults();
            std::cout << "\t\t";
            anchor_finder_global.m_bm_seeding.PrintResults();
            // anchor_finder_global.m_bm_seed_subsetting.PrintResults();
            std::cout << "\t\t";
            anchor_finder_global.m_bm_processing.PrintResults();
            std::cout << "\t\t";
            anchor_finder_global.m_bm_pairing.PrintResults();
            std::cout << "\t\t";
            anchor_finder_global.m_bm_extend_anchors.PrintResults();
            std::cout << "\t\t";
            anchor_finder_global.m_bm_sorting_anchors.PrintResults();
            std::cout << "\t\t";
            anchor_finder_global.m_bm_recovering_anchors.PrintResults();
            // bm_anchor_recovery_global.PrintResults();
            bm_alignment_global.PrintResults();
            bm_alignment_join_sort_global.PrintResults();
            bm_output_global.PrintResults();
            alignment_handler_global.bm_alignment.PrintResults();
            alignment_handler_global.m_bm_alignment.PrintResults();
            std::cout << "----------------------------------------------------\n" << std::endl;

            std::cout << "Reads that had anchors recovered: "
                      << static_cast<double>(anchor_finder_global.recovered_count) / anchor_finder_global.total_count
                      << std::endl;
        }

        std::ofstream time_os(options.GetPrefix(options.GetCurrentIndex()) + "_runtime.tsv", std::ios::out);
        time_os << anchor_finder_global.m_bm_seeding.GetName() << '\t' << anchor_finder_global.m_bm_seeding.GetDuration(Time::seconds) << '\n';
        time_os << anchor_finder_global.m_bm_processing.GetName() << '\t' << anchor_finder_global.m_bm_processing.GetDuration(Time::seconds) << '\n';
        time_os << anchor_finder_global.m_bm_pairing.GetName() << '\t' << anchor_finder_global.m_bm_pairing.GetDuration(Time::seconds) << '\n';
        time_os << anchor_finder_global.m_bm_sorting_anchors.GetName() << '\t' << anchor_finder_global.m_bm_sorting_anchors.GetDuration(Time::seconds) << '\n';
        time_os << anchor_finder_global.m_bm_extend_anchors.GetName() << '\t' << anchor_finder_global.m_bm_extend_anchors.GetDuration(Time::seconds) << '\n';
        time_os << bm_anchor_finder_global.GetName() << '\t' << bm_anchor_finder_global.GetDuration(Time::seconds) << '\n';
        time_os << bm_anchor_recovery_global.GetName() << '\t' << bm_anchor_recovery_global.GetDuration(Time::seconds) << '\n';
        time_os << bm_alignment_global.GetName() << '\t' << bm_alignment_global.GetDuration(Time::seconds) << '\n';
        time_os << bm_alignment_join_sort_global.GetName() << '\t' << bm_alignment_join_sort_global.GetDuration(Time::seconds) << '\n';
        time_os << bm_output_global.GetName() << '\t' << bm_output_global.GetDuration(Time::seconds) << '\n';
        time_os.close();

        seed_sizes_global.ToTSV(options.GetPrefix(options.GetCurrentIndex()) + "_seedsizes_histogram.tsv");
        anchor_sizes_global.ToTSV(options.GetPrefix(options.GetCurrentIndex()) + "_anchorsizes_histogram.tsv");

        return statistics;
    }

}