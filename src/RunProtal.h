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

namespace protal {
    static void RunOld(int argc, char *argv[]) {
        Benchmark bm_total("Run protal");
        std::cout << "Protal setup" << std::endl;

        auto options = protal::Options::OptionsFromArguments(argc, argv);

        std::cout << "Options:\n" << options.ToString() << std::endl;

        size_t kmer_size = 15;
        CanonicalKmerIterator iterator = CanonicalKmerIterator(kmer_size);
        Syncmer minimizer(15, 7, 2);
//        None dummy_minimizer;



        if (options.BuildMode()) {
//            Benchmark bm_build("Run build");
//            KmerPutterSM kmer_putter{};
//            auto protal_stats = protal::build::Run<protal::CanonicalKmerIterator, KmerPutterSM, DEBUG_NONE>(
//                    options, kmer_putter, iterator);
//            bm_build.PrintResults();
//            protal_stats.WriteStats(std::cout);
        } else {
            // Benchmark Load Time
            Benchmark bm_load_index("Load Index");

            // Load Index
            Seedmap map;

            std::cout << options.GetIndexFile() << std::endl;

            std::ifstream idx_in(options.GetIndexFile(), std::ios::binary);
            map.Load(idx_in);
            idx_in.close();

            // Print Load Time
            bm_load_index.PrintResults();

            KmerLookupSM kmer_lookup(map);

            GenomeLoader genomes(options.GetSequenceFile(), options.GetSequenceMapFile());
            if (options.PreloadGenomes()) {
                std::cout << "Preload genomes" << std::endl;
                genomes.LoadAllGenomes();
            }

//            std::ofstream ofs(options.GetOutputFile());

//            KmerProcessor kmer_processor(kmer_lookup, genomes, kmer_size, ofs);
            KmerProcessor kmer_processor(kmer_lookup, genomes, kmer_size);


            Benchmark bm_classify("Run classify");
            std::string first_file = options.GetFirstFile();

            size_t iterations = 0;

            auto protal_stats = protal::classify::RunOld<protal::CanonicalKmerIterator, KmerLookupSM, KmerProcessor<KmerLookupSM>, DEBUG_NONE, Syncmer>(
                    first_file, options, kmer_lookup, kmer_processor, iterator, minimizer);//, genomes, map);
            protal_stats.WriteStats(std::cout);
            bm_classify.PrintResults();


            std::cout << "Loaded genomes: " << genomes.GetLoadedGenomeCount() << std::endl;
            std::cout << "Dummy: " << kmer_processor.m_dummy << std::endl;
        }

        bm_total.PrintResults();
    }


    static void Run(int argc, char *argv[]) {
        Benchmark bm_total("Run protal");

        auto options = protal::Options::OptionsFromArguments(argc, argv);

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
            std::cout << options.GetIndexFile() << std::endl;
            std::ifstream idx_in(options.GetIndexFile(), std::ios::binary);
            map.Load(idx_in);
            idx_in.close();
            // Print Load Time
            bm_load_index.PrintResults();


            KmerLookupSM kmer_lookup(map);

            GenomeLoader genomes(options.GetSequenceFile(), options.GetSequenceMapFile());
            WFA2Wrapper aligner(4, 6, 2);

            if (options.PreloadGenomes()) {
                std::cout << "Preload genomes" << std::endl;
                genomes.LoadAllGenomes();
            }

//            using AnchorFinder = SimpleAnchorFinder<KmerLookupSM>;
            using AnchorFinder = NaiveAnchorFinder<KmerLookupSM>;

            // new anchor finder approach
//            SimpleAnchorFinder anchor_finder(kmer_lookup);
            AnchorFinder anchor_finder(kmer_lookup);

            // AlignmentHandler approach
            SimpleAlignmentHandler alignment_handler(genomes, aligner, kmer_size);

            Benchmark bm_classify("Run classify");
            if (options.PairedMode()) {
                std::cout << "Paired-end reads" << std::endl;

                std::ifstream is1 {options.GetFirstFile(), std::ios::in};
                std::ifstream is2 {options.GetSecondFile(), std::ios::in};
                SeqReaderPE reader{is1, is2};

                auto protal_stats = protal::classify::RunPairedEnd<
                        SimpleKmerHandler<ClosedSyncmer>,
                        AnchorFinder,
                        SimpleAlignmentHandler,
                        DEBUG_NONE>(
                        reader, options, anchor_finder, alignment_handler, iterator);

                protal_stats.WriteStats();

                is1.close();
                is2.close();

            } else {
                std::cout << "Single-end reads" << std::endl;
                auto protal_stats = protal::classify::Run<
                        SimpleKmerHandler<ClosedSyncmer>,
                        AnchorFinder,
                        SimpleAlignmentHandler,
                        DEBUG_NONE>(
                        options.GetFirstFile(), options, anchor_finder, alignment_handler, iterator);

                protal_stats.WriteStats();
            }

            bm_classify.PrintResults();

            std::cout << "Loaded genomes: " << genomes.GetLoadedGenomeCount() << std::endl;
        }

        bm_total.PrintResults();
    }



    static void Test(int argc, char *argv[]) {
        string pattern = "TCTTTACTCGCGCGTTGGAGAAATACAATAGT";
        string text    = "TCTATACTGCGCGTTTGGAGAAATAAAATAGT";
        WFA2Wrapper aligner(4, 6, 2);
        AlignmentResult result;
        aligner.Alignment(pattern, text, result);
        aligner.PrintAlignment();

    }
}
