//
// Created by fritsche on 14/08/22.
//

#pragma once

#include <iostream>
#include "Utils/Benchmark.h"
#include "Options.h"
#include "Build.h"
#include "Alignment/WFA2Wrapper.h"
#include "Classify.h"
#include "KmerProcessor.h"
#include "AnchorFinder.h"
#include "AlignmentStrategy.h"
#include "gzstream/gzstream.h"

namespace protal {
    template<typename AlignmentBenchmark=NoBenchmark>
    static void RunWrapper(Options& options, AlignmentBenchmark benchmark=NoBenchmark{}) {

        std::cout << "Options:\n" << options.ToString() << std::endl;

        const size_t kmer_size = 15;
        ClosedSyncmer minimizer{kmer_size, 7, 2};
        SimpleKmerHandler iterator{kmer_size, minimizer};

        if (options.BuildMode()) {
            Benchmark bm_build("Run build");
            KmerPutterSM kmer_putter{};
            auto protal_stats = protal::build::Run<SimpleKmerHandler<ClosedSyncmer>, KmerPutterSM, DEBUG_NONE>(
                    options, kmer_putter, iterator);
            bm_build.PrintResults();
            protal_stats.WriteStats(std::cout);

        } else {

            // Benchmark Load Time
            Benchmark bm_load_index("Load Index");
            // Load Index
            Seedmap map;
            std::ifstream idx_in(options.GetIndexFile(), std::ios::binary);
            map.Load(idx_in);
            idx_in.close();
            // Print Load Time
            bm_load_index.PrintResults();


            KmerLookupSM kmer_lookup(map);

            GenomeLoader genomes(options.GetSequenceFile(), options.GetSequenceMapFile());
            WFA2Wrapper aligner(4, 6, 2, options.GetXDrop());

            if (options.PreloadGenomes()) {
                Benchmark bm_preload_genomes("Preload genomes");
                genomes.LoadAllGenomes();
                bm_preload_genomes.PrintResults();
            }

//            using AnchorFinder = SimpleAnchorFinder<KmerLookupSM>;
//            using AnchorFinder = HashMapAnchorFinder<KmerLookupSM>;
            using AnchorFinder = ListAnchorFinder<KmerLookupSM>;
//            using AnchorFinder = NaiveAnchorFinder<KmerLookupSM>;

            using OutputHandler = VarkitOutputHandler;

            using optional_ofstream = std::optional<std::ofstream>;
            std::ofstream varkit_output(options.GetOutputFile(), std::ios::binary);
            std::ofstream sam_output(options.GetOutputFile() + ".sam", std::ios::out);
            optional_ofstream snp_output =
                    options.NoStrains() ? optional_ofstream{} :
                    optional_ofstream{ std::in_place, options.GetOutputFile() + ".snps", std::ios::out };

            OutputHandler output_handler(varkit_output, sam_output, snp_output, 1024*512, 1024*1024*16, 0.8);

            // AnchorFinder
            AnchorFinder anchor_finder(kmer_lookup);

            // AlignmentHandler approach
            SimpleAlignmentHandler alignment_handler(genomes, aligner, kmer_size, options.GetAlignTop(), options.GetMaxScoreAni());

            Benchmark bm_classify("Run classify");
            if (options.PairedMode()) {
//                std::ifstream is1 {options.GetFirstFile(), std::ios::in};
//                std::ifstream is2 {options.GetSecondFile(), std::ios::in};

                igzstream is1 { options.GetFirstFile().c_str() };
                igzstream is2 { options.GetSecondFile().c_str() };
                SeqReaderPE reader{is1, is2};


                auto protal_stats = protal::classify::RunPairedEnd<
                        SimpleKmerHandler<ClosedSyncmer>,
                        AnchorFinder,
                        SimpleAlignmentHandler,
                        OutputHandler,
                        DEBUG_NONE,
                        AlignmentBenchmark>(
                        reader, options, anchor_finder, alignment_handler, output_handler, iterator, benchmark);

                protal_stats.WriteStats();

                is1.close();
                is2.close();

            } else {
//                std::ifstream is {options.GetFirstFile(), std::ios::in};
                igzstream is { options.GetFirstFile().c_str() };
                SeqReader reader{ is };

                auto protal_stats = protal::classify::Run<
                        SimpleKmerHandler<ClosedSyncmer>,
                        AnchorFinder,
                        SimpleAlignmentHandler,
                        OutputHandler,
                        DEBUG_NONE,
                        AlignmentBenchmark>(
                        reader, options, anchor_finder, alignment_handler, output_handler, iterator, benchmark);

                is.close();
                protal_stats.WriteStats();
            }

            // Close output streams;
            varkit_output.close();
            sam_output.close();

            if (snp_output.has_value()) snp_output.value().close();



            bm_classify.PrintResults();

            std::cout << "Loaded genomes: " << genomes.GetLoadedGenomeCount() << std::endl;
        }
    }

    static void Run(int argc, char *argv[]) {
        Benchmark bm_total("Run protal");

        using AlignmentBenchmark = CoreBenchmark;

        auto options = protal::Options::OptionsFromArguments(argc, argv);

        if (options.Help()) {
            options.PrintHelp();
            exit(0);
        }

        // Untangle Template options that ed to be written out specifically.
        if (options.BenchmarkAlignment()) {
            AlignmentBenchmark alignment_benchmark{};
            if (!options.GetBenchmarkAlignmentOutputFile().empty()) {
                alignment_benchmark.SetOutput(options.GetBenchmarkAlignmentOutputFile());
            }
            RunWrapper(options, alignment_benchmark);
        } else {
            RunWrapper(options);
        }

        bm_total.PrintResults();
    }



    static void Test(int argc, char *argv[]) {
        string pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
        string text    = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        WFA2Wrapper aligner(4, 6, 2, 0);
        AlignmentResult result;
        aligner.Alignment(pattern, text, result);
        aligner.PrintAlignment();

    }
}
