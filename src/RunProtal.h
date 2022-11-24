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
#include "ChainAnchorFinder.h"
//#include "AnchorFinder.h"
#include "AlignmentStrategy.h"
#include "gzstream/gzstream.h"
#include "Profiler.h"
#include "Taxonomy.h"

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
//            using AnchorFinder = ListAnchorFinder<KmerLookupSM>;
            using AnchorFinder = ChainAnchorFinder<KmerLookupSM>;
//            using AnchorFinder = NaiveAnchorFinder<KmerLookupSM>;

            using OutputHandler = ProtalOutputHandler;

            std::ofstream sam_output(options.SamFile(), std::ios::out);

            // Maybe change in flag if file should be written.
            genomes.WriteSamHeader(sam_output);

            // AnchorFinder
            AnchorFinder anchor_finder(kmer_lookup);

            // AlignmentHandler approach
            SimpleAlignmentHandler alignment_handler(genomes, aligner, kmer_size, options.GetAlignTop(), options.GetMaxScoreAni(), options.FastAlign());

            Benchmark bm_classify("Run classify");
            if (options.PairedMode()) {
                using OutputHandler = ProtalPairedOutputHandler;
                OutputHandler output_handler(sam_output, 1024*512, 1024*1024*16, 0.8);
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
                        reader, options, anchor_finder, alignment_handler, output_handler, iterator, genomes, benchmark);

                protal_stats.WriteStats();

                is1.close();
                is2.close();

            } else {
                using OutputHandler = ProtalOutputHandler;
                OutputHandler output_handler(sam_output, 1024*512, 1024*1024*16, 0.8);
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
            sam_output.close();


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

        if (!options.SamFile().empty() && options.ProfileOnly()) {
            goto Profile;
        }

        // Untangle Template options that ed to be written out specifically.
        if (options.BenchmarkAlignment()) {
            AlignmentBenchmark alignment_benchmark{};
            if (!options.GetBenchmarkAlignmentOutputFile().empty()) {
                std::cout << "Benchmark to " << options.GetBenchmarkAlignmentOutputFile() << std::endl;
                alignment_benchmark.SetOutput(options.GetBenchmarkAlignmentOutputFile());
            }
            RunWrapper(options, alignment_benchmark);
            if (!options.GetBenchmarkAlignmentOutputFile().empty()) {
                alignment_benchmark.DestroyOutput();
            }
        } else {
            RunWrapper(options);
        }

        /*
         * PROFILER
         */
        std::cout << "Profile? " << options.Profile() << std::endl;
        if (options.Profile()) {
            Profile:


            GenomeLoader genomes(options.GetSequenceFile(), options.GetSequenceMapFile());
            taxonomy::IntTaxonomy taxonomy(options.GetInternalTaxonomyFile());


            std::cout << "Sam file: " << options.SamFile() << std::endl;
            profiler::Profiler profiler(genomes);

            // New Profiler approach
            std::vector<AlignmentPair> unique_pairs;
            std::vector<std::vector<AlignmentPair>> pairs;

            profiler.FromSam(options.SamFile(), pairs, unique_pairs);

            std::cout << "Unique pairs: " << unique_pairs.size() << std::endl;
            std::cout << "Pairs:        " << pairs.size() << std::endl;


//            profiler.FromSam(options.SamFile());
            std::optional<TruthSet> truth = !options.ProfileTruthFile().empty() ?
                                            std::optional<TruthSet>{ protal::GetTruth(options.ProfileTruthFile()) } :
                                            std::optional<TruthSet>{};

            if (options.BenchmarkAlignment()) {
                profiler.TestSNPUtils(pairs);
//                profiler.OutputErrorData(pairs);
            } else if (truth.has_value()) {
//                profiler.TestSNPUtils(pairs);
                profiler.OutputErrorData(truth.value(), pairs);

            }
//            Utils::Input();
            exit(9);

            auto profile = truth.has_value() ? profiler.ProfileWithTruth(truth.value()) : profiler.Profile();

            if (truth.has_value()) {
                profile.AnnotateWithTruth(truth.value());
            }

            profiler::TaxonFilter filter(0.95, 0.7, 70);

            std::cout << "Save profile to " << options.GetOutputPrefix() + ".profile" << std::endl;
            std::ofstream os(options.GetOutputPrefix() + ".profile", std::ios::out);
            profile.WriteSparseProfile(taxonomy, filter, os);
            os.close();

//            std::cout << profile.ToString() << std::endl;
//            for (auto& [key, taxon] : profile.GetTaxa()) {
//                std::cout << taxon.ToString() << std::endl;
//            }
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
