//
// Created by fritsche on 14/08/22.
//

#pragma once

#include "Profiler.h"
#include <iostream>
#include "Utils/Benchmark.h"
#include "Options.h"
#include "Build.h"
#include "Alignment/WFA2Wrapper.h"
#include "Classify.h"
#include "ChainAnchorFinder.h"
#include "Taxonomy.h"
#include "gzstream.h"
#include "AlignmentStrategy.h"
#include "TaxonStatisticsOutput.h"
#include "ProgressBar.h"

// #include "Profiler/AlignmentContainer.h"
// #include "Profiler/Profiler.h"

#include "Profiler/ReadFilter.h"

//#include "AnchorFinder.h"
//#include "AlignmentStrategy.h"
//#include "gzstream/gzstream.h"
//#include "Profiler.h"
//#include "Taxonomy.h"
//#include "TaxonStatisticsOutput.h"
//#include <math.h>
#include <unistd.h>
#include <algorithm>

namespace protal {

    static const std::string PROTAL_LOGO =
            "                                             ,,  \n"
            "`7MM\"\"\"Mq.                   mm            `7MM  \n"
            "  MM   `MM.                  MM              MM  \n"
            "  MM   ,M9 `7Mb,od8 ,pW\"Wq.mmMMmm  ,6\"Yb.    MM  \n"
            "  MMmmdM9    MM' \"'6W'   `Wb MM   8)   MM    MM  \n"
            "  MM         MM    8M     M8 MM    ,pm9MM    MM  \n"
            "  MM         MM    YA.   ,A9 MM   8M   MM    MM  \n"
            ".JMML.     .JMML.   `Ybmd9'  `Mbmo`Moo9^Yo..JMML.\n\n";



    static void PrintLogo(std::ostream& os = std::cout) {
        os << PROTAL_LOGO << std::endl;
    }

    static void PrintProtalInformation(std::ostream& os = std::cout) {
        os << "\nOption --output_names has recently been changed to --prefix."<< std::endl;
    }

    class ProtalDB {
        GenomeLoader m_genomes;
        std::optional<taxonomy::IntTaxonomy> m_taxonomy;

    public:
        ProtalDB(std::string sequence_file, std::string map_file) :
                m_genomes(sequence_file, map_file),
                m_taxonomy() {
        }

        ProtalDB(std::string sequence_file, std::string map_file, std::string hittable_genes_file, std::string unique_kmers_file) :
                m_genomes(sequence_file, map_file),
                m_taxonomy() {
            // m_genomes.LoadHittableGenes(hittable_genes_file);
            m_genomes.LoadUniqueKmers(unique_kmers_file);
        }

        void LoadTaxonomy(std::string file) {
            m_taxonomy = taxonomy::IntTaxonomy(file);
        }

        bool IsTaxonomyLoaded() const {
            return m_taxonomy.has_value();
        }

        taxonomy::IntTaxonomy& GetTaxonomy() {
            if (!m_taxonomy.has_value()) {
                std::cerr << "Taxonomy has not been loaded. Error in code. " << std::endl;
                exit(12);
            }
            return m_taxonomy.value();
        }

        GenomeLoader& GetGenomes() {
            return m_genomes;
        }

        std::optional<taxonomy::IntTaxonomy> GetTaxonomyOptional() {
            return m_taxonomy;
        }
    };

    template<typename AlignmentBenchmark=NoBenchmark>
    static void RunWrapper(Options& options, ProtalDB& db, AlignmentBenchmark benchmark=NoBenchmark{}) {

        const size_t mmer_size = 15;
//        const size_t kmer_size = 27;
        const size_t kmer_size = 31;
        ClosedSyncmer minimizer{mmer_size, 7, 2};
        SimpleKmerHandler iterator{kmer_size, mmer_size, minimizer};


        if (options.BuildMode()) {
            Benchmark bm_build("Run build");
            bm_build.Start();
            KmerPutterSM kmer_putter{};
            auto protal_stats = protal::build::Run<SimpleKmerHandler<ClosedSyncmer>, KmerPutterSM, DEBUG_NONE>(
                    options, kmer_putter, iterator);
            bm_build.PrintResults();
            protal_stats.WriteStats(std::cout);

        } else {

            // Benchmark Load Time
            Benchmark bm_load_index("Load Index");
            bm_load_index.Start();
            // Load Index
            Seedmap map;
            std::ifstream idx_in(options.GetIndexFile(), std::ios::binary);
            map.Load(idx_in);
            idx_in.close();
            bm_load_index.Stop();
            bm_load_index.PrintResults();


            KmerLookupSM kmer_lookup(map, options.GetMaxKeyUbiquity());

            GenomeLoader& genomes = db.GetGenomes();

            WFA2Wrapper2 aligner(4, 6, 2, options.GetXDrop());

            if (options.PreloadGenomes() && !genomes.AllGenomesLoaded()) {
                Benchmark bm_preload_genomes("Preload genomes");
                bm_preload_genomes.Start();
                genomes.LoadAllGenomes();
                bm_preload_genomes.Stop();
                bm_preload_genomes.PrintResults();
            }

//            using AnchorFinder = SimpleAnchorFinder<KmerLookupSM>;
//            using AnchorFinder = HashMapAnchorFinder<KmerLookupSM>;
//            using AnchorFinder = ListAnchorFinder<KmerLookupSM>;
            using AnchorFinder = ChainAnchorFinder<KmerLookupSM>;
//            using AnchorFinder = NaiveAnchorFinder<KmerLookupSM>;

            using OutputHandler = ProtalOutputHandler;




            Benchmark bm_classify("Processing all samples");
            bm_classify.Start();

            if (options.PairedMode()) {
                for (auto index : options.GetRange()) {
                //for (auto index = 0; index < options.GetFileCount(); index++) {

                    Benchmark bm_classify_sample("Aligning reads");
                    bm_classify_sample.Start();

                    // TODO implement logger in protal
                    auto sam = options.SamFile(index);
                    auto dir = std::filesystem::path(sam).parent_path();
                    if (!std::filesystem::create_directories(dir.string()) && !std::filesystem::exists(dir)) {
                        std::cout << "Cannot create directories for this path " << sam << std::endl;
                        exit(2);
                    };


                    std::cout << index << " Process sample " << options.GetSampleId(index) << (std::filesystem::exists(sam) ? " (sam exists)" : " (sam does not exist)") << std::endl;

                    // Avoid aligning files that already exist.
                    if (!options.Force() && std::filesystem::exists(sam)) {
                        std::cout << "Skip " << sam << " continue" << std::endl;
                        continue;
                    }

                    // AnchorFinder
                    AnchorFinder anchor_finder(kmer_lookup, mmer_size, 4, options.GetMaxSeedSize(), genomes);
                    // AlignmentHandler approach
                    SimpleAlignmentHandler alignment_handler(genomes, aligner, kmer_size, options.GetAlignTop(), options.GetMaxScoreAni(), options.FastAlign());



                    options.SetCurrentIndex(index);
                    std::ofstream sam_output(sam, std::ios::out);
//                    ogzstream sam_output(sam.c_str());

                    genomes.WriteSamHeader(sam_output);
                    // {
                    //     Benchmark count_lines_bm{"Count lines for fastq"};
                    //     count_lines_bm.Start();
                    //     igzstream is1 { options.GetFirstFile(index).c_str() };
                    //     igzstream is2 { options.GetSecondFile(index).c_str() };
                    //     auto is1_lines = Utils::CountLines(is1);
                    //     auto is2_lines = Utils::CountLines(is2);
                    //     count_lines_bm.Stop();
                    //     count_lines_bm.PrintResults();
                    //     std::cout << is1_lines << " " << is2_lines << std::endl;
                    //     if (is1_lines != is2_lines) {
                    //         bm_classify_sample.Stop();
                    //         bm_classify.PrintResults();
                    //         std::cerr << "Error: Paired end files have different line counts" << std::endl;
                    //         std::cerr << options.GetFirstFile(index).c_str() << std::endl;
                    //         std::cerr << options.GetSecondFile(index).c_str() << std::endl;
                    //         continue;
                    //     }
                    // }

                    igzstream is1 { options.GetFirstFile(index).c_str() };
                    igzstream is2 { options.GetSecondFile(index).c_str() };
                    SeqReaderPE reader{is1, is2};

                    std::cout << "Process: " << options.GetFirstFile(index) << std::endl;

                    // Main Run Call. This is where the reads are read and alignment happens
                    if (options.GetMAPQDebugOut()) {
                        using OutputHandler = ProtalPairedOutputHandler<true>;

                        OutputHandler output_handler(sam_output, options.GetMaxOut(), 1024*512, 1024*1024*16, genomes, 0.8);
                        auto protal_stats = protal::classify::RunPairedEnd<
                                SimpleKmerHandler<ClosedSyncmer>,
                                AnchorFinder,
                                SimpleAlignmentHandler,
                                OutputHandler,
                                DEBUG_NONE,
                                AlignmentBenchmark>(
                                reader, options, anchor_finder, alignment_handler, output_handler, iterator, genomes, benchmark);
                        if (options.Verbose()) {
                            protal_stats.WriteStats();
                        }
                    } else {
                        using OutputHandler = ProtalPairedOutputHandler<false>;

                        OutputHandler output_handler(sam_output, options.GetMaxOut(), 1024*512, 1024*1024*16, genomes, 0.8);
                        auto protal_stats = protal::classify::RunPairedEnd<
                                SimpleKmerHandler<ClosedSyncmer>,
                                AnchorFinder,
                                SimpleAlignmentHandler,
                                OutputHandler,
                                DEBUG_NONE,
                                AlignmentBenchmark>(
                                reader, options, anchor_finder, alignment_handler, output_handler, iterator, genomes, benchmark);
                        if (options.Verbose()) {
                            protal_stats.WriteStats();
                        }
                    }
                    bm_classify_sample.Stop();
                    bm_classify_sample.PrintResults();

                    is1.close();
                    is2.close();
                    sam_output.close();

                    // CleanUp
                    if (!reader.Success()) {
                        std::cerr << "There was an error reading the fastq files with sample " << options.GetSampleId(index) << " (" << index << ")" << std::endl;
                        std::cerr << options.GetFirstFile(index) << std::endl;
                        std::cerr << options.GetSecondFile(index) << std::endl;
                        std::cerr << "Remove sam file: " << sam << std::endl;
                        std::filesystem::remove(sam);
                    }
                }

            } else {
                // AnchorFinder
                AnchorFinder anchor_finder(kmer_lookup, mmer_size, 4, options.GetMaxSeedSize(), genomes);
                // AlignmentHandler approach
                SimpleAlignmentHandler alignment_handler(genomes, aligner, kmer_size, options.GetAlignTop(), options.GetMaxScoreAni(), options.FastAlign());


                std::ofstream sam_output(options.SamFile(0), std::ios::out);
                genomes.WriteSamHeader(sam_output);
                using OutputHandler = ProtalOutputHandler;
                OutputHandler output_handler(sam_output, 1024*512, 1024*1024*16, genomes, 0.8);
//                std::ifstream is {options.GetFirstFile(), std::ios::in};
//                std::ifstream is {options.GetFirstFile(), std::ios::in};
                igzstream is { options.GetFirstFile(0).c_str() };
                SeqReader reader{ is };

                auto protal_stats = protal::classify::Run<
                        SimpleKmerHandler<ClosedSyncmer>,
                        AnchorFinder,
                        SimpleAlignmentHandler,
                        OutputHandler,
                        DEBUG_NONE,
                        AlignmentBenchmark>(
                        reader, options, anchor_finder, alignment_handler, output_handler, iterator, benchmark);

                protal_stats.WriteStats();
                // Close output streams;
                sam_output.close();
                is.close();
            }
            bm_classify.Stop();
            bm_classify.PrintResults();
        }
    }

    using Profile = profiler::MicrobialProfile;
    using Profiles = std::vector<Profile>;

//     void ProfileWrapper2(Options& options, ProtalDB& db) {
//         std::cout << "ProfileWrapper2" << std::endl;
//         GenomeLoader& genomes = db.GetGenomes();
//
//         if (!db.IsTaxonomyLoaded()) db.LoadTaxonomy(options.GetInternalTaxonomyFile());
//         auto& taxonomy = db.GetTaxonomy();
//
// //        omp_set_num_threads(6);
//
// #pragma omp parallel for default(none) shared(options, cout, taxonomy, genomes, db)
//         for (auto i : options.GetRange()) {
//         //for (auto i = 0; i < options.GetFileCount(); i++) {
//
//             Profiler::AlignmentContainer ac;
//             auto sam = options.SamFile(i);
//
//
// #pragma omp critical(read_sam)
//             ac.LoadSam(sam);
//
//             const auto& alignments = ac.AlignmentPairs();
//             Profiler::ReadFilter filter;
//             Profiler::Profiler<Profiler::ReadFilter> profiler(filter);
//             profiler.Profile(ac);//, db);
//
// //             profiler(genomes);
// //
// //            profiler.FromSam(sam);
//         }
//     }

    Profiles ProfileWrapper(Options& options, ProtalDB& db) {
        GenomeLoader& genomes = db.GetGenomes();

        if (!db.IsTaxonomyLoaded()) db.LoadTaxonomy(options.GetInternalTaxonomyFile());
        auto& taxonomy = db.GetTaxonomy();

        std::vector<profiler::MicrobialProfile> profiles;


        double min_ani = 0.95;
        double min_gene_presence = 0.50; //previously 0.5
        size_t min_total_hits = 60; //previously 70
        size_t min_mean_mapq = 10;




//        using TaxonFilterObj = profiler::TaxonFilterForest;
        using TaxonFilterObj = profiler::TaxonFilterObj;
//        TaxonFilterObj filter(min_ani, min_gene_presence, min_total_hits, min_mean_mapq);
        std::string model = options.GetIndexFolder() + "/random_forest.xml";
        std::cout << "Model: " << model << std::endl;
        TaxonFilterObj filter(model);

        // ProgressBar prog;
        // prog.Reset(options.GetFileCount());

        omp_set_num_threads(options.GetThreads());



        #pragma omp parallel for default(none) firstprivate(filter) shared(options, cout, taxonomy, profiles, genomes, std::cerr)//, bm_read_alignments, bm_profile)
        for (auto i : options.GetRange()) {
#pragma omp critical(print)
            std::cerr << omp_get_thread_num() << " File " << i << " of " << options.GetRange().size() << ":\n\t" << options.SamFile(i) << std::endl;

            auto sam = options.SamFile(i);
            auto sample_name = options.GetSampleId(i);

            if (!Utils::exists(sam)) {
                std::cerr << "Sam file does not exist for sample " << options.GetSampleId(i) << " (" << i << ")" << std::endl;
                profiles.emplace_back(profiler::MicrobialProfile{genomes});
                std::cerr << sam << std::endl;
                continue;
            }

            Benchmark bm_read_alignments{ "Load read alignments" };
            Benchmark bm_profile{ "Profile sample" };

// #pragma omp critical(progress)
//             prog.UpdateAdd(1);

            profiler::Profiler profiler(genomes);
            profiler.SetNoStrain(options.NoStrains());

            // New Profiler approach
            std::vector<AlignmentPair> unique_pairs;
            std::vector<std::vector<AlignmentPair>> pairs;


#pragma omp critical(print)
            std::cout << "Thread " << omp_get_thread_num() << " read sam file " << sam << std::endl;
            if (!Utils::exists(sam)) {
                std::cerr << "File does not exist" << sam << std::endl;
                exit(90);
            }
            std::cout << "Thread " << omp_get_thread_num() << " read now" << std::endl;

#pragma omp critical(load_sam)
            profiler.FromSam(sam);

            std::cout << "Thread " << omp_get_thread_num() << " after Load" << std::endl;

            if (!profiler.HasReads()) {
                std::cerr << "Empty sam file" << std::endl;
                profiles.emplace_back(profiler::MicrobialProfile{genomes});
                continue;
            }


            std::ofstream erro(sam + ".err", std::ios::out);

            profiler.PrintStats();
            bm_profile.Start();
#pragma omp critical(print)
            std::cout << "Thread " << omp_get_thread_num() << " run profile" << std::endl;
            auto profile = profiler.Profile(sample_name, std::optional<std::reference_wrapper<std::ostream>>{erro});
            erro.close();
            bm_profile.Stop();
#pragma omp critical(print)
            std::cout << "Thread " << omp_get_thread_num() << " ";
            bm_profile.PrintResults();

            profile.bm_add_sam.PrintResults();


            std::optional<TruthSet> truth = options.HasProfileTruths() ?
                                            std::optional<TruthSet>{ protal::GetTruth(options.ProfileTruthFile(i), taxonomy) } :
                                            std::optional<TruthSet>{};

            if (options.BenchmarkAlignment()) {
                profiler.TestSNPUtils(pairs);
                profiler.OutputErrorData(pairs);
            } else if (truth.has_value()) {
                std::cout << "Has Truth Available" << std::endl;
                profiler.OutputErrorData(unique_pairs, pairs, &truth.value());
            }

            if (truth.has_value()) {
                std::string truth_output = options.ProfileFile(i) + ".truth_annotated";
                profile.AnnotateWithTruth(truth.value(), filter, truth_output);
            }

            auto dir = std::filesystem::path(options.ProfileFile(i)).parent_path();
            if (!std::filesystem::create_directories(dir.string()) && !std::filesystem::exists(dir)) {
                std::cout << "Cannot create directories for this path " << options.ProfileFile(i) << std::endl;
                exit(2);
            };

            std::ofstream os(options.ProfileFile(i), std::ios::out);
            std::ofstream os_total(options.ProfileFile(i) + ".log", std::ios::out);
            std::ofstream os_dismissed(options.ProfileFile(i) + ".gene.log", std::ios::out);

            profile.WriteSparseProfile(taxonomy, filter, os, &os_total, &os_dismissed);
            os.close();
            os_total.close();
            os_dismissed.close();

            profile.SetName(options.GetSampleId(i));

#pragma omp critical(add_profiles)
            profiles.emplace_back(profile);
        }
        return profiles;
    }

    static std::pair<bool,bool> SharedSNP(VariantBin& a, VariantBin& b) {
        std::sort(a.begin(), a.end(), [](Variant const& va, Variant const& vb) {
            return va.Observations() > vb.Observations();
        });
        std::sort(b.begin(), b.end(), [](Variant const& va, Variant const& vb) {
            return va.Observations() > vb.Observations();
        });

        bool shared_non_ref = false;
        bool shared_ref = false;
        for (auto& var_a : a) {
            if (!var_a.GetValid()) continue;
            for (auto& var_b : b) {
                if (!var_b.GetValid()) continue;

                if (var_a.Match(var_b))  {
                    shared_ref |= var_a.IsReference() && var_b.IsReference();
                    shared_non_ref |= !var_a.IsReference() && !var_b.IsReference();
                }
            }
        }

        return { shared_ref, shared_non_ref };
    }

    static double JukesCantor(size_t mismatch, size_t length) {
        auto strain_identity = static_cast<double>(mismatch) / length;

        auto jc69 = -3.0f/4 * log(1 - (4.0f/3 * strain_identity));
        return jc69 > 1 || jc69 != jc69 ? 1 : jc69;
    }

    static bool IsSNP(VariantBin& bin, size_t min_obs, double min_frequency) {
        for (auto& var : bin) {
            size_t total_obs = std::accumulate(bin.begin(), bin.end(), 0, [](size_t acc, Variant const& var) { return acc + var.Observations(); });
            auto obs = var.Observations();
            double freq = static_cast<double>(obs) / total_obs;
            if (!var.IsReference() && freq >= min_frequency && (obs >= min_obs || var.HasFwdAndRev())) {
                return true;
            }
        }
//
//        for (auto& var : bin) {
//            std::cout << var.ToString() << '\t';
//        }
//        std::cout << endl;
//
//        Utils::Input();
        return false;
    }

    static bool HasReference(VariantBin& bin) {
        return std::any_of(bin.begin(), bin.end(), [](Variant const& var) {
            return var.IsReference();
        });
    }

    static std::pair<size_t, size_t> Distance(size_t shared_length, VariantVec& a, VariantVec& b) {
        const auto pos_getter = [](VariantBin const& bin) { return bin.front().Position(); };
        size_t shared_snps = 0;
        size_t unique_snps_a = 0;
        size_t unique_snps_b = 0;

        size_t min_observations = 2;
        double min_frequency = 0.2;

        int ai = 0, bi = 0;
//        std::cout << "a: " << a.size() << "  b: " << b.size() << std::endl;
        while (ai < a.size() && bi < b.size()) {
            auto pos_a = pos_getter(a[ai]);
            auto pos_b = pos_getter(b[bi]);

            bool has_ref_a = std::any_of(a[ai].begin(), a[ai].end(), [](Variant const& v) {return v.GetValid() && v.IsReference();});
            bool has_ref_b = std::any_of(b[bi].begin(), b[bi].end(), [](Variant const& v) {return v.GetValid() && v.IsReference();});
            bool has_snp_a = std::any_of(a[ai].begin(), a[ai].end(), [](Variant const& v) {return v.GetValid() && !v.IsReference();});
            bool has_snp_b = std::any_of(b[bi].begin(), b[bi].end(), [](Variant const& v) {return v.GetValid() && !v.IsReference();});
//            bool has_snp_b = IsSNP(b[bi], min_observations, min_frequency);

            bool same_pos = (pos_a == pos_b);
            bool advance_a = (pos_a <= pos_b);
            bool advance_b = (pos_b <= pos_a);

//            std::cout << ai << " " << pos_a << "\t->\t" << a[ai].front().ToString() << "\t" << b[bi].front().ToString() << "\t<-\t" << pos_b << " " << bi << std::endl;

            // same_pos:                 required because a shared snp needs variants detected at the same pos.
            // has_snp_a || has_snp_b:   if both have variants at the position but none are considered
            //                           genuine snps, it is assumed that both are reference alleles
            // SharedSNP(a[ai], b[bi]):  Check if they share at least one valid snp, ref or non ref.
            auto [has_shared_ref, has_shared_snp] = same_pos ? SharedSNP(a[ai], b[bi]) : std::pair<bool,bool>(false, false);
            bool is_shared_and_snp = has_shared_snp && !has_shared_ref;
            shared_snps += is_shared_and_snp;
            // has_snp_a: if all variants at this position are dodgy/fp SNPs assume its the reference allele
            // same_po
            // advance var is needed because only the smaller pos can be a snp
            bool is_unique_snp_a = (advance_a && !has_shared_ref && has_snp_a && !has_ref_a && !has_shared_snp);
            bool is_unique_snp_b = (advance_b && !has_shared_ref && has_snp_b && !has_ref_b && !has_shared_snp);
            unique_snps_a += is_unique_snp_a;
            unique_snps_b += is_unique_snp_b;

//            if (is_unique_snp_a || unique_snps_b) {
//                std::cout << std::string(60, '#') << std::endl;
//                std::cout << (is_unique_snp_a ? "SNP IN A" : "SNP IN B") << std::endl;
//
//                std::cout << "____________SNPS A noref:" << !has_ref_a << std::endl;
//                for (auto& var : a[ai]) {
//                    std::cout << var.ToString() << '\t';
//                }
//                std::cout << std::endl;
//                std::cout << "____________SNPS B noref:" << !has_ref_b << std::endl;
//                for (auto& var : b[bi]) {
//                    std::cout << var.ToString() << '\t';
//                }
//                std::cout << std::endl;
//                Utils::Input();
//            }

//            if ((has_snp_a && advance_a && !shared) || (has_snp_b && advance_b && !shared)) {
//                std::cout << "PosA: " << pos_a << "\tPosB:" << pos_b << std::endl;
//                std::cout << "SNPA: " << has_snp_a << "\tSNPB:" << has_snp_b << std::endl;
//                std::cout << "shared: " << shared << std::endl;
//                std::cout << "-------\nA is SNP? " << has_snp_a << std::endl;
//                for (auto& var : a[ai]) {
//                    std::cout << var.ToString() << '\t';
//                }
//                std::cout << endl;
//                std::cout << "B is SNP? " << has_snp_b << std::endl;
//                for (auto& var : b[bi]) {
//                    std::cout << var.ToString() << '\t';
//                }
//                std::cout << endl;
//
//                Utils::Input();
//            }

            ai += advance_a;
            bi += advance_b;
        }

        if (a.size() - ai > 0) {
//            std::cout << "A) SNPS_A: " << unique_snps_a << "\tSNPS_B: " << unique_snps_b << std::endl;
            for (; ai < a.size(); ai++) {
                bool has_var = std::any_of(a[ai].begin(), a[ai].end(), [](Variant const& v) {return v.GetValid() && !v.IsReference();});
                bool has_ref = std::any_of(a[ai].begin(), a[ai].end(), [](Variant const& v) {return v.GetValid() && v.IsReference();});
                unique_snps_a += (has_var && !has_ref);
            }
//            std::cout << "  -->  SNPS_A: " << unique_snps_a << "\tSNPS_B: " << unique_snps_b << std::endl;
        }

        if (b.size() - bi > 0) {
//            std::cout << "B) SNPS_A: " << unique_snps_a << "\tSNPS_B: " << unique_snps_b << std::endl;
            for (; bi < b.size(); bi++) {
                bool has_var = std::any_of(b[bi].begin(), b[bi].end(), [](Variant const& v) {return v.GetValid() && !v.IsReference();});
                bool has_ref = std::any_of(b[bi].begin(), b[bi].end(), [](Variant const& v) {return v.GetValid() && v.IsReference();});
                unique_snps_b += (has_var && !has_ref);
            }
        }

//        if (unique_snps_a + unique_snps_b > 0) {
//            Utils::Input();
//        }

        return std::pair{ unique_snps_a + unique_snps_b, shared_length };

//        return JukesCantor(unique_snps_a + unique_snps_b, shared_length);
    }

    static std::pair<double, size_t> Distance(profiler::Gene& g1, profiler::Gene& g2) {
        auto& strain1 = g1.GetStrainLevel();
        auto& strain2 = g2.GetStrainLevel();
        auto shared_alignment_region = SharedAlignmentRegion::GetSharedAlignmentRegion(strain1, strain2);

        auto& va = shared_alignment_region.variant_handler_a;
        auto& vb = shared_alignment_region.variant_handler_b;
        auto shared_sequence_length = shared_alignment_region.share_range.SequenceLength();


        auto [unique_snps, shared_length] = Distance(shared_sequence_length, va, vb);

        return {unique_snps, shared_sequence_length};
    }

    static std::pair<double, size_t> Distance(profiler::Taxon& t1, profiler::Taxon& t2) {

        double distance_sum = 0;
        size_t total_unique_snps = 0;
        size_t shared_region_sum = 0;

        for (auto& [gid, _] : t1.GetGenes()) {
            if (t2.GetGenes().contains(gid)) {
                //std::cout << "Stuck at " << gid << std::endl;
                auto [unique_snps, shared_region] = Distance(t1.GetGenes().at(gid), t2.GetGenes().at(gid));

                total_unique_snps += unique_snps;
                shared_region_sum += shared_region;
            }
        }

        auto distance = JukesCantor(total_unique_snps, shared_region_sum);

//        if (distance < 0.001) {
//            std::cout << "Distance: " << distance << " with " << total_unique_snps << " snps over " << shared_region_sum << " bases." << std::endl;
//        }

        return { distance, shared_region_sum };
    }


    using StrainResults = tsl::robin_map<size_t, Matrix<double>>;
    using TaxidCounts = tsl::robin_map<uint32_t, uint32_t>;
    using TaxidSet = tsl::robin_set<uint32_t>;
    using TaxidList = std::vector<uint32_t>;
    using OptionalFilter = optional<std::reference_wrapper<const profiler::TaxonFilter>>;
    using SimilarityMatrix = DoubleMatrix;

    TaxidList ExtractTaxa(Profiles const& profiles, std::optional<profiler::TaxonFilterObj> filter= {}) {
        TaxidCounts taxid_counts;
        TaxidSet taxa;
        TaxidList taxid_list;

        for (auto& profile : profiles) {
            for (auto& [id, taxon] : profile.GetTaxa()) {
                if (filter.has_value() && !filter->Pass(taxon)) continue;
                if (!taxid_counts.contains(id)) {
                    taxid_counts.insert({ id, 0 });
                }
                taxid_counts.at(id)++;
            }
        }

        for (auto& [t, c] : taxid_counts) {
            if (c > 1) {
                taxid_list.emplace_back(t);
            }
        }

        std::sort(taxid_list.begin(), taxid_list.end(), [&taxid_counts](uint32_t const& t1, uint32_t const& t2) {
            return taxid_counts.at(t1) > taxid_counts.at(t2);
        });

        return taxid_list;
    }

    static std::pair<double, size_t> GetSimilarity(uint32_t taxid, Profile& profile1, Profile& profile2, Options& options, std::optional<profiler::TaxonFilterObj> filter={}) {
        size_t min_shared_region = 1000;
        if (profile1.GetTaxa().contains(taxid) && profile2.GetTaxa().contains(taxid)) {
            auto& taxon1 = profile1.GetTaxa().at(taxid);
            auto& taxon2 = profile2.GetTaxa().at(taxid);

            if (filter.has_value() && !(filter->Pass(taxon1) && filter->Pass(taxon2))) {
                return { NAN, 0 };
            }

            auto [distance, shared_region] = Distance(taxon1, taxon2);

            auto similarity = 1.0f-distance;

            if (shared_region <= min_shared_region) {
                similarity = NAN;
            }
            return { similarity, shared_region };
        }
        return { NAN, 0 };
    }

    static bool IsRowGood(std::vector<char> row, size_t min_hcov) {
        size_t count_non_n = 0;
        for (auto c : row) {
            count_non_n += (c != 'N' && c != '-');
        }
        return count_non_n >= min_hcov;
    }

    static std::vector<uint32_t> GetInformationVector(MSAVector const& msa_vector) {
        size_t info = 0;
        const auto msa_len = msa_vector.front().size();
        std::vector<uint32_t> information_vector(msa_len, 0);

        for (auto pos = 0; pos < msa_len; pos++) {
            info = 0;
            for (auto const& seq : msa_vector) info += (seq[pos] != 'N' && seq[pos] != '-');
            information_vector[pos] = info;
        }
        return information_vector;
    }

    static MSAVector ProcessMSA(MSAVector const& msa_vector, double position_coverage = 0.5) {
        auto sample_len = msa_vector.size();
        auto info_vector = GetInformationVector(msa_vector);
        MSAVector processed_msa(sample_len, std::vector<char>{});

        auto sample_size_threshold = msa_vector.size() * position_coverage;
        std::cout << "Sample size threshold: " << sample_size_threshold << std::endl;
        for (auto pos = 0; pos < info_vector.size(); pos++) {
            if (info_vector[pos] > sample_size_threshold) {
                for (auto s = 0; s < sample_len; s++) {
                    processed_msa[s].emplace_back(msa_vector[s][pos]);
                }
            }
        }
        return processed_msa;
    }

    static void OutputMSA(MSAVector const& msa, std::vector<std::string> const& names, std::ostream& os=std::cout) {
        for (auto i = 0; i < msa.size(); i++) {
            auto& row = msa[i];

            if (row.empty()) continue;
            os << ">" << names[i] << std::endl;
            os << std::string_view(row.begin(), row.end()) << std::endl;
        }
    }

    static std::vector<size_t> GetProfilesWithTaxon(uint32_t taxid, Profiles& profiles, Options& options, std::optional<profiler::TaxonFilter>& filter) {
        std:vector<size_t> indices;

        for (auto i = 0; i < profiles.size(); i++) {
            auto& profile = profiles[i];
            if (!profile.GetTaxa().contains(taxid)) {
                continue;
            }
            auto& taxon_map = profile.GetTaxa();

            auto& taxon = taxon_map.at(taxid);
            if (filter.has_value() && !filter.value().Pass(taxon)) {
                continue;
            }
            std::cout << taxid << " Passes " << std::endl;

            indices.emplace_back(i);
        }
        return indices;
    }

    static std::vector<uint32_t> SelectGenesForTaxon(uint32_t taxid, std::string name, std::vector<size_t>& selected_profiles, GenomeLoader& loader, Options& options, Profiles& profiles) {
        std::vector<uint32_t> selected_gene_ids;
        std::vector<uint32_t> gene_ids = loader.GetGenome(taxid).GetHittableGenes();

        auto max_gene_id = std::max_element(gene_ids.begin(), gene_ids.end());

        // std::vector<std::vector<double>> per_sample_gene_multiallelic_portion;
        // per_sample_gene_multiallelic_portion.resize(selected_profiles.size(),
        //     std::vector<double>(*max_gene_id + 1, -1) );

        std::vector<std::vector<int>> per_sample_gene_multiallelic_snps;
        per_sample_gene_multiallelic_snps.resize(selected_profiles.size(),
            std::vector<int>(*max_gene_id + 1, -1) );

        std::vector<std::vector<int>> per_sample_gene_snps;
        per_sample_gene_snps.resize(selected_profiles.size(),
            std::vector<int>(*max_gene_id + 1, -1) );

        std::vector<std::vector<int>> per_sample_gene_cov;
        per_sample_gene_cov.resize(selected_profiles.size(),
            std::vector<int>(*max_gene_id + 1, -1) );

        std::vector<std::vector<int>> per_sample_gene_noise;
        per_sample_gene_noise.resize(selected_profiles.size(),
            std::vector<int>(*max_gene_id + 1, -1) );

        for (auto si = 0; si < selected_profiles.size(); si++) {
            std::cout << " _______________Sample: " << si << std::endl;
            auto sample_index = selected_profiles[si];
            auto& taxon = profiles[sample_index].GetTaxa().at(taxid);
            auto& multiallelic_vec = per_sample_gene_multiallelic_snps[si];

            // for (auto& [gid, gene] : taxon.GetGenes()) {
            //     std::cout << "gid: " << gid << " " << gene.ToString() << std::endl;
            // }
            // std::cout << " ------ " << std::endl;
            for (auto& gene_id : gene_ids) {
                if (gene_id >= per_sample_gene_multiallelic_snps[si].size()) {
                    std::cout << gene_id << " >= " << gene_ids.size() << std::endl;
                    exit(123);
                }
                if (!taxon.GetGenes().contains(gene_id)) {
                    per_sample_gene_multiallelic_snps[si][gene_id] = -1;
                    per_sample_gene_snps[si][gene_id] = -1;
                    per_sample_gene_cov[si][gene_id] = 0;
                    continue;
                }

                auto& gene = taxon.GetGene(gene_id);

                auto allele_counts = gene.AlleleSNPCounts(2, 60);
                auto allele_counts_nf = gene.AlleleSNPCounts(0, 0);

                per_sample_gene_multiallelic_snps[si][gene_id] = allele_counts.Multi();
                per_sample_gene_snps[si][gene_id] = allele_counts.AllValid();
                per_sample_gene_cov[si][gene_id] = gene.Coverage();
                per_sample_gene_noise[si][gene_id] = allele_counts.Noisy();
            }
        }

        std::ofstream os(options.GetOutputDir() + '/' + name + ".multiallelic.tsv", std::ios::out);
        std::cout << "OUTPUT: " << (options.GetOutputDir() + '/' + name + ".multiallelic.tsv") << std::endl;
        for (auto si = 0; si < selected_profiles.size(); si++) {
            auto sample_index = selected_profiles[si];
            os << options.GetSampleId(sample_index);

            if (si >= per_sample_gene_multiallelic_snps.size()) exit(244);
            auto& vec = per_sample_gene_multiallelic_snps[si];
            double total = 0;
            for (auto i = 1; i < vec.size(); i++) {
                os << '\t' << vec[i];
                total += vec[i] == -1 ? 0 : vec[i];
            }
            os << '\t' << total << std::endl;
        }
        os.close();

        std::ofstream os2(options.GetOutputDir() + '/' + name + ".total.tsv", std::ios::out);
        std::cout << "OUTPUT: " << (options.GetOutputDir() + '/' + name + ".total.tsv") << std::endl;
        for (auto si = 0; si < selected_profiles.size(); si++) {
            auto sample_index = selected_profiles[si];
            os2 << options.GetSampleId(sample_index);

            if (si >= per_sample_gene_snps.size()) exit(244);
            auto& vec = per_sample_gene_snps[si];
            double total = 0;
            for (auto i = 1; i < vec.size(); i++) {
                os2 << '\t' << vec[i];
                total += vec[i] == -1 ? 0 : vec[i];
            }
            os2 << '\t' << total << std::endl;
        }
        os2.close();

        std::ofstream os3(options.GetOutputDir() + '/' + name + ".cov.tsv", std::ios::out);
        std::cout << "OUTPUT: " << (options.GetOutputDir() + '/' + name + ".cov.tsv") << std::endl;
        for (auto si = 0; si < selected_profiles.size(); si++) {
            auto sample_index = selected_profiles[si];
            os3 << options.GetSampleId(sample_index);

            if (si >= per_sample_gene_snps.size()) exit(244);
            auto& vec = per_sample_gene_cov[si];
            size_t total = 0;
            for (auto i = 1; i < vec.size(); i++) {
                os3 << '\t' << vec[i];
                total += vec[i] == -1 ? 0 : vec[i];
            }
            os3 << '\t' << total << std::endl;
        }
        os3.close();


        for (auto gi = 0; gi < per_sample_gene_noise.front().size(); gi++) {
            std::cout << gi;
            for (auto si = 0; si < selected_profiles.size(); si++) {
                auto sample_index = selected_profiles[si];
                auto sample_name = options.GetSampleId(sample_index);
                std::cout << '\t' << (per_sample_gene_noise[si][gi]/static_cast<double>(per_sample_gene_cov[si][gi])) << "(" << per_sample_gene_noise[si][gi] << "/" << per_sample_gene_cov[si][gi] << ")";
            }
            std::cout << std::endl;
        }


        return gene_ids;
    }

    static void GetMSAForTaxon (uint32_t taxid, std::string taxon_name, GenomeLoader& loader, Options& options, Profiles& profiles, std::ostream* os_meta=nullptr, std::optional<profiler::TaxonFilter> filter={}) {
        auto min_hcov = 5000;
        auto min_cov = 2;
        auto min_qual_sum = 60;
        auto min_samples_with_gene = 3;

        std::vector<size_t> profile_indices = GetProfilesWithTaxon(taxid, profiles, options, filter);

        if (profile_indices.empty()) return;

        MSAVector msa{ profile_indices.size(), std::vector<char>() };

        auto& genome = loader.GetGenome(taxid);
        if (!genome.IsLoaded()) genome.LoadGenomeOMP();
        std::vector<std::string> names;
        std::vector<std::string> partitions;
        size_t partition_start = 0;
        size_t previous_size = 0;

        std::cout << "MULTIALLELIC: " << taxid << " " << taxon_name << std::endl;
        std::vector<uint32_t> selected_genes = SelectGenesForTaxon(taxid, taxon_name, profile_indices, loader, options, profiles);

        ProgressBar prog(120);
        std::cout << "Process " << selected_genes.size() << std::endl;
        for (auto& geneid : selected_genes) {
            std::cout << "GID: " << geneid << std::endl;
            // if (!loader.GetGenome(taxid).IsGeneHittable(geneid)) {
            //     continue;
            // }

            prog.UpdateAdd(1);
            MSASequenceItems items;

            // if (!genome.ValidGene(geneid)) break;
            //
            auto& gene = genome.GetGene(geneid);
            // if (!gene.IsSet()) continue;

            // Check if
            size_t samples_with_gene = 0;

            for (auto i = 0; i < profile_indices.size(); i++) {
                std::cout << "profile: " << i << std::endl;

                auto& profile = profiles[profile_indices[i]];

                auto& taxon_map = profile.GetTaxa();

                if (!taxon_map.contains(taxid)) continue;

                names.emplace_back(profile.GetName());
                auto& genes = profile.GetTaxa().at(taxid).GetGenes();

                if (!genes.contains(geneid)) {
                    items.emplace_back(OptionalMSASequenceItem{});
                } else {
                    samples_with_gene++;
                    auto& gene_obs = genes.at(geneid);
                    auto& strain = gene_obs.GetStrainLevel();
                    auto snps = SharedAlignmentRegion::GetSNPs(strain.GetVariantHandler());
                    auto& region = strain.GetSequenceRangeHandler();
                    items.emplace_back( OptionalMSASequenceItem { { std::move(snps), region } } );

                    if (os_meta) {
#pragma omp critical(metaout)
                        {
                            auto ac = gene_obs.AlleleSNPCounts(min_cov, min_qual_sum);
                            auto length = ranges::count_if(region.CalculateCoverageVector(), [](auto val){ return(val >= 2);});
                            *os_meta << profile.GetName() << '\t';
                            *os_meta << geneid << '\t';
                            *os_meta << "gene" << geneid << '\t';
                            *os_meta << geneid << '\t';
                            *os_meta << (length > 0 ? ac.Multi()/static_cast<double>(length) : 0) << '\t';
                            *os_meta << (length > 0 ? ac.Filtered()/static_cast<double>(length) : 0) << '\t';
                            *os_meta << gene_obs.VerticalCoverage() << std::endl;
                        }
                    }
                }
            }

            if (samples_with_gene > min_samples_with_gene) {
                previous_size = msa.front().size();

                for (auto& opt : items)  {
                    if (!opt.has_value()) continue;
                    auto& [var, shr] = opt.value();
                    if (std::any_of(var.begin(), var.end(), [](std::vector<Variant> const& vv) {
                        return std::any_of(vv.begin(), vv.end(), [](Variant const&  v) {
                            return v.Observations() == 65535;
                        });
                    })) {
                        std::cout << "Wrong variant 2" << std::endl;
                        exit(3);
                    }
                }


                bool result = protal::MSA(items, gene.Sequence(), msa, min_cov, min_qual_sum);

                if (!result) continue;
                if (msa.front().size() > partition_start) {
                    if (partitions.size() > 0) {
                        partitions.back() += std::to_string(previous_size-1);
                    }

                    size_t partition_end = msa.front().size();
                    std::string partition = "DNA, gene";
                    partition += std::to_string(geneid) + " = ";
                    partition += std::to_string(partition_start) + '-';
//                    partition += std::to_string(partition_end-1);
                    partitions.emplace_back(partition);
                    partition_start = msa.front().size();
                }
            }
        }

        bool any_good = std::any_of(msa.begin(), msa.end(), [](MSARow const& row){
            return IsRowGood(row, 5000);
        });
        if (!any_good) {
            std::cout << "No any good" << std::endl;
            return;
        }

        std::ofstream os(options.GetMSAOutput(taxon_name), std::ios::out);
        for (auto i = 0; i < msa.size(); i++) {
            auto& row = msa[i];

            if (!IsRowGood(row, min_hcov)) continue;
            os << ">" << names[i] << std::endl;
            os << std::string_view(row.begin(), row.end()) << std::endl;
        }
        os.close();
        std::cout << " Saved MSA to " << options.GetMSAOutput(taxon_name);

        partitions.back() += std::to_string(msa.front().size());

        std::ofstream os_part(options.GetMSAPartitionOutput(taxon_name), std::ios::out);
        for (auto i = 0; i < partitions.size(); i++) {
            os_part << partitions[i] << std::endl;
        }
        os_part.close();

        std::cout << " Saved Partitions " << std::endl;

        auto processed_msa = protal::ProcessMSA(msa, options.GetMSAMinVCOV());
        std::cout << "Trimmed size: " << processed_msa.front().size() << " with minvcov: " << options.GetMSAMinVCOV() << std::endl;

        os = std::ofstream(options.GetMSAProcessedOutput(taxon_name), std::ios::out);
        std::cout << "Processed output: " << options.GetMSAProcessedOutput(taxon_name) << std::endl;
        OutputMSA(processed_msa, names, os);
        os.close();
    }



    static SimilarityMatrix GetSimilarityMatrixForTaxon(uint32_t taxid, Options& options, Profiles& profiles, std::optional<profiler::TaxonFilterObj> filter={}) {
        SimilarityMatrix matrix;

        size_t min_shared_length = 1000;

        size_t total_combinations = profiles.size() * profiles.size() - profiles.size() / 2;
        size_t count_combinations = 1;

        ProgressBar prog(((profiles.size() * profiles.size()) - profiles.size()) / 2);
        for (auto i = 0; i < profiles.size(); i++) {
            auto& profile1 = profiles[i];
            auto name1 = profile1.GetName();

            for (auto j = i+1; j < profiles.size(); j++) {
                prog.UpdateAdd(1);
                count_combinations++;
//                std::cout << "\rCombination " << count_combinations << " of " << total_combinations;
                auto& profile2 = profiles[j];
                auto name2 = profile2.GetName();

                auto [similarity, shared_length] = GetSimilarity(taxid, profile1, profile2, options, filter);

//                std::cout << name1 << " -- " << name2 <<  " = " << similarity << " over " << shared_length << std::endl;

                if (similarity == NAN || shared_length < min_shared_length) continue;

                if (!matrix.HasName(name1)) {
                    matrix.AddName(name1);
                }
                if (!matrix.HasName(name2)) {
                    matrix.AddName(name2);
                }

//                std::cout << taxid << " " << similarity << " " << shared_length << std::endl;

                matrix.SetValue(name1, name2, similarity, true);
                matrix.SetValue(name1, name2, static_cast<double>(shared_length), false);
            }
        }
        std::cout << std::endl;
        return matrix;
    }



    static void WriteDistanceMatrix(uint32_t id, SimilarityMatrix const& matrix, Options& options, std::string& name) {
        std::ofstream matrix_os(options.GetSimilarityMatrixOutput(name), std::ios::out);
        matrix.PrintMatrix(matrix_os, "\t", 10);
        matrix_os.close();
    }

    static void StrainWrapper2(Options& options, Profiles& profiles, GenomeLoader& loader, taxonomy::IntTaxonomy& taxonomy, std::optional<profiler::TaxonFilterObj> filter={}) {
        Benchmark bm_strain{"Strain-level MSAs"};
        bm_strain.Start();
        std::cout << "Output " << options.GetOutputDir() << std::endl;
        auto dir = std::filesystem::path(options.GetOutputDir());
        if (!std::filesystem::create_directories(dir.string()) && !std::filesystem::exists(dir)) {
            std::cout << "Cannot create directories for this path " << dir << std::endl;
            exit(2);
        };

        auto taxids = ExtractTaxa(profiles, filter);

        auto enable_similarity_matrix = false;


        for (auto& taxid : taxids) {
            std::cout << "Target taxid: " << taxid << " >> " << taxonomy.Get(taxid).scientific_name << std::endl;

            std::string name = taxonomy.Get(taxid).scientific_name;
            std::replace(name.begin(), name.end(), ' ', '_');

            if (enable_similarity_matrix) {
                auto similarities = GetSimilarityMatrixForTaxon(taxid, options, profiles, filter);
                if (!similarities.AnySet()) continue;
                WriteDistanceMatrix(taxid, similarities, options, name);
            }

            std::cout << "GetMSAForTaxon " << taxonomy.Get(taxid).ToString() << std::endl;
            std::ofstream os_meta(options.GetSpeciesMetaOutput(name));
            GetMSAForTaxon(taxid, name, loader, options, profiles, &os_meta);
            os_meta.close();
            std::cout << "Check " << taxid << " " << name << std::endl;
//            Utils::Input();
        }
        bm_strain.Stop();
        bm_strain.PrintResults();
    }


    size_t GetTotalSystemMemory()
    {
        long pages = sysconf(_SC_PHYS_PAGES);
        long page_size = sysconf(_SC_PAGE_SIZE);
        return pages * page_size;
    }

    static void Run(int argc, char *argv[]) {
        PrintLogo();
        PrintProtalInformation();
        std::cout << "Total available memory is " << GetTotalSystemMemory() / (1024 * 1024 * 1024) << "GB" << std::endl;;

        Benchmark bm_total("Run protal");
        bm_total.Start();

        using AlignmentBenchmark = CoreBenchmark;

        auto options = protal::Options::OptionsFromArguments(argc, argv);

        if (options.Help()) {
            options.PrintHelp();
            exit(0);
        }

        std::cout << "Options:\n" << options.ToString() << std::endl;


        // Load protal DB into RAM
        ProtalDB db = options.UniqueKmersFileExists() ?
            ProtalDB(options.GetSequenceFile(), options.GetSequenceMapFile(), options.GetHittableGenesMap(), options.GetUniqueKmersFile()) :
            ProtalDB(options.GetSequenceFile(), options.GetSequenceMapFile());

        // Load fasta sequences of reference into RAM (advised)
        if (options.PreloadGenomes()) {
            std::cout << "Preload genomes" << std::endl;
            Benchmark bm_preload_genomes("Preload genomes");
            bm_preload_genomes.Start();
            db.GetGenomes().LoadAllGenomes();
            bm_preload_genomes.Stop();
            bm_preload_genomes.PrintResults();
        }

        // Skip alignment if files are present. Do not skip if either files are not there or user specified --force
        auto sam_files = options.SamFiles();
        bool all_alignments_exist = std::all_of(sam_files.begin(), sam_files.end(), [](std::string const& file){ return Utils::exists(file); });

        bool skip_alignment = !options.BuildMode() && all_alignments_exist && !options.Force();
        if (!options.BuildMode() && !options.SamFile(0).empty() && (options.ProfileOnly() || skip_alignment)) {
            std::cout << "All alignments are present." << std::endl;
            goto Profile;
        }



        /*
         *  READ ALIGNMENT SECTION
         */
        // Untangle Template options that need to be written out specifically.
        if (options.BenchmarkAlignment()) {
            AlignmentBenchmark alignment_benchmark{};
            if (!options.GetBenchmarkAlignmentOutputFile().empty()) {
                alignment_benchmark.SetOutput(options.GetBenchmarkAlignmentOutputFile());
            }
            RunWrapper(options, db, alignment_benchmark);
            if (!options.GetBenchmarkAlignmentOutputFile().empty()) {
                alignment_benchmark.DestroyOutput();
            }
        } else {
            RunWrapper(options, db);
        }

        /*
         * PROFILER
         */
        std::cout << "Profile? " << options.Profile() << std::endl;
        if (options.Profile()) {
            Profile:

            db.LoadTaxonomy(options.GetInternalTaxonomyFile());

            Benchmark bm_profiling("Profiling");
            bm_profiling.Start();
            // ProfileWrapper2(options, db);
            // exit(9);
            auto profiles = ProfileWrapper(options, db);
            bm_profiling.Stop();
            bm_profiling.PrintResults();


            // Remove later - keep option
            double min_ani = 0.95;
            double min_gene_presence = 0.5; //previously 0.5
            size_t min_total_hits = 60; //previously 70
            size_t min_mean_mapq = 10;

//            profiler::TaxonFilter filter(min_ani, min_gene_presence, min_total_hits, min_mean_mapq);

            std::string model = options.GetIndexFolder() + "/random_forest.xml";
            std::cout << "Model: " << model << std::endl;
            profiler::TaxonFilterObj filter(model);

            // If there is more than one profile, get per taxon output
            if (profiles.size() > 1) {
                auto taxids = ExtractTaxa(profiles);
                auto& taxonomy = db.GetTaxonomy();

                for (auto taxid : taxids) {
                    TaxonStatisticsOutput stats;
                    TaxonStatisticsOutput all_stats;

                    std::string name = taxonomy.Get(taxid).scientific_name;
                    std::replace(name.begin(), name.end(), ' ', '_');

                    std::ofstream os(options.GetOutputDir() + '/' + name + ".statistics.tsv");

                    stats.PrintHeader(os);
                    for (auto& profile : profiles) {
                        auto& taxa = profile.GetTaxa();
                        if (!taxa.contains(taxid)) continue;

                        auto& taxon = taxa.at(taxid);
                        bool accepted = filter.Pass(taxon);
                        stats.PrintLine(os, profile.GetName(), taxon.VerticalCoverage(), taxon.TotalHits(), taxon.TotalLength(), taxon.GetMeanANI(), taxon.GetMeanMAPQ(), accepted);
                    }
                    os.close();
                }
            }

            /*
             * STRAIN PART -  RESOLVE MSAs BETWEEN SAMPLES
             */
            if (!options.NoStrains()) {
                std::cout << "Start StrainWrapper" << std::endl;
                StrainWrapper2(options, profiles, db.GetGenomes(), db.GetTaxonomy(), filter);
            }
        }

        bm_total.Stop();
        bm_total.PrintResults();

        std::cout << "Find the results under:" << std::endl;
        std::cout << options.GetOutputDir() << std::endl;
    }
}
