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
#include "TaxonStatisticsOutput.h"
#include <math.h>

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

            WFA2Wrapper aligner(4, 6, 2, options.GetXDrop());

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




            Benchmark bm_classify("Run classify");
            bm_classify.Start();

            if (options.PairedMode()) {
                for (auto index = 0; index < options.GetFileCount(); index++) {
                    Benchmark bm_classify_sample("Classifying sample");
                    bm_classify_sample.Start();

                    // TODO implement logger in protal
                    auto sam = options.SamFile(index);
                    std::cout << "Process sample " << options.GetSampleId(index) << (std::filesystem::exists(sam) ? " (sam exists)" : " (sam does not exist)") << std::endl;

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
//                    std::ofstream sam_output(sam, std::ios::out);
                    ogzstream sam_output(sam.c_str());
                    genomes.WriteSamHeader(sam_output);
                    igzstream is1 { options.GetFirstFile(index).c_str() };
                    igzstream is2 { options.GetSecondFile(index).c_str() };
                    SeqReaderPE reader{is1, is2};

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

//                auto protal_stats = protal::classify::RunPairedEnd<
//                        SimpleKmerHandler<ClosedSyncmer>,
//                        AnchorFinder,
//                        SimpleAlignmentHandler,
//                        OutputHandler,
//                        DEBUG_NONE,
//                        AlignmentBenchmark>(
//                        reader, options, anchor_finder, alignment_handler, output_handler, iterator, genomes, benchmark);
//                protal_stats.WriteStats();

                    is1.close();
                    is2.close();
                    sam_output.close();

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
//        std::string model = "/usr/users/QIB_fr017/fritsche/Documents/HPC/group/projects/protal/camisim_all_output/output_flex_r214/profiles/mapq4_2/random_forest_cami_r214.pmml";
        TaxonFilterObj filter(model);

        for (auto i = 0; i < options.GetFileCount(); i++) {
            profiler::Profiler profiler(genomes);

            // New Profiler approach
            std::vector<AlignmentPair> unique_pairs;
            std::vector<std::vector<AlignmentPair>> pairs;

            auto sam = options.SamFile(i);
            std::cout << "Read profile from Sam " << sam << " (exists: " << Utils::exists(options.SamFile(i)) << ")" << std::endl;
            profiler.FromSam(sam);

            auto profile = profiler.Profile();

            std::cout << "Print profile" << std::endl;
            std::cout << profile.ToString(db.GetTaxonomy(), filter) << std::endl;


//            profiler.TestSNPUtils(profiler.m_pairs_unique);

//            profiler.FromSam(options.SamFile());

            std::optional<TruthSet> truth = options.HasProfileTruths() ?
                                            std::optional<TruthSet>{ protal::GetTruth(options.ProfileTruthFile(i)) } :
                                            std::optional<TruthSet>{};

            if (options.BenchmarkAlignment()) {
                profiler.TestSNPUtils(pairs);
                profiler.OutputErrorData(pairs);
            } else if (truth.has_value()) {
                profiler.OutputErrorData(unique_pairs, pairs, &truth.value());
            }

            if (truth.has_value()) {
//                profiler::TaxonFilter filter(0.95, 0.7, 70, 10);
                std::string truth_output = options.ProfileFile(i) + ".truth_annotated";
                profile.AnnotateWithTruth(truth.value(), filter, truth_output);
            }


            std::cout << "Save profile to " << options.ProfileFile(i) << std::endl;
            auto dir = std::filesystem::path(options.ProfileFile(i)).parent_path();
            std::cout << "Create directories " << dir.string() << std::endl;
            if (!std::filesystem::create_directories(dir.string()) && !std::filesystem::exists(dir)) {
                std::cout << "Cannot create directories for this path " << options.ProfileFile(i) << std::endl;
                exit(2);
            };

            std::ofstream os(options.ProfileFile(i), std::ios::out);
            std::ofstream os_total(options.ProfileFile(i) + ".log", std::ios::out);
            std::ofstream os_dismissed(options.ProfileFile(i) + ".gene.log", std::ios::out);


            //profile.WriteSparseProfile(taxonomy, filter);
            profile.WriteSparseProfile(taxonomy, filter, os, &os_total, &os_dismissed);
            os.close();
            os_total.close();
            os_dismissed.close();

//            profile.ClearSams();
            profile.SetName(options.GetSampleId(i));
            profiles.emplace_back(profile);

            profiler.m_post_process_bm.PrintResults();
        }

        std::cout << "Profiles: " << profiles.size() << std::endl;
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


    using StrainResults = tsl::robin_map<size_t, DoubleMatrix>;
    using TaxidCounts = tsl::robin_map<uint32_t, uint32_t>;
    using TaxidSet = tsl::robin_set<uint32_t>;
    using TaxidList = std::vector<uint32_t>;
    using OptionalFilter = optional<std::reference_wrapper<const profiler::TaxonFilter>>;
    using SimilarityMatrix = DoubleMatrix;

    TaxidList ExtractTaxa(Profiles const& profiles, std::optional<profiler::TaxonFilter> filter= {}) {
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

    static bool IsRowGood(std::vector<char> row, size_t min_chars) {
        size_t count_non_n = 0;
        for (auto c : row) {
            count_non_n += (c != 'N' && c != '-');
        }
        return count_non_n >= min_chars;
    }


    static void GetMSAForTaxon (uint32_t taxid, std::string taxon_name, GenomeLoader& loader, Options& options, Profiles& profiles, std::optional<profiler::TaxonFilter> filter={}) {


        std::cout << "GET MSA FOR TAXON " << taxid << std::endl;

        auto min_cov = 2;
        auto min_qual_sum = 60;

        std::vector<size_t> profile_indices;
        for (auto i = 0; i < profiles.size(); i++) {
            if (profiles[i].GetTaxa().contains(taxid)) {
                profile_indices.emplace_back(i);
            }
        }
        MSAVector msa{ profile_indices.size(), std::vector<char>() };

        auto& genome = loader.GetGenome(taxid);
        if (!genome.IsLoaded()) genome.LoadGenomeOMP();
        std::vector<std::string> names;
        std::vector<std::string> partitions;
        size_t partition_start = 0;
        size_t previous_size = 0;
        std::cout << "Genes ";
        for (auto geneid = 1; geneid <= 120; geneid++) {
            MSASequenceItems items;

            if (!genome.HasGene(geneid)) break;

            auto& gene = genome.GetGene(geneid);
            if (!gene.IsSet()) continue;



//            std::cout << "found: " << geneid << std::endl;
            bool hasone = false;
            for (auto i = 0; i < profile_indices.size(); i++) {

                auto& profile = profiles[profile_indices[i]];
                auto& taxon = profile.GetTaxa();

                if (!taxon.contains(taxid)) continue;

                names.emplace_back(profile.GetName());
                auto& genes = profile.GetTaxa().at(taxid).GetGenes();

                if (!genes.contains(geneid)) {
                    items.emplace_back(OptionalMSASequenceItem{});
                    continue;
                } else {
                    hasone = true;
                    auto& strain = genes.at(geneid).GetStrainLevel();
                    auto snps = SharedAlignmentRegion::GetSNPs(strain.GetVariantHandler());
                    auto& region = strain.GetSequenceRangeHandler();
                    items.emplace_back( OptionalMSASequenceItem { { std::move(snps), region } } );
                }
            }

            if (hasone) {
                std::cout << " " << geneid;

                previous_size = msa.front().size();
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
        std::cout << " Done ... ";

        bool any_good = std::any_of(msa.begin(), msa.end(), [](MSARow const& row){
            return IsRowGood(row, 5000);
        });
        if (!any_good) return;

        std::ofstream os(options.GetMSAOutput(taxon_name), std::ios::out);
        auto s = 0;
        for (auto i = 0; i < msa.size(); i++) {
            auto& row = msa[i];

            if (!IsRowGood(row, 5000)) continue;
            os << ">" << names[i] << std::endl;
            os << std::string_view(row.begin(), row.end()) << std::endl;
        }
        os.close();

        std::cout << " Saved MSA ";

        partitions.back() += std::to_string(msa.front().size());

        std::ofstream os_part(options.GetMSAPartitionOutput(taxon_name), std::ios::out);
        for (auto i = 0; i < partitions.size(); i++) {
            os_part << partitions[i] << std::endl;
        }
        os_part.close();

        std::cout << " Saved Partitions " << std::endl;
    }


    static SimilarityMatrix GetSimilarityMatrixForTaxon(uint32_t taxid, Options& options, Profiles& profiles, std::optional<profiler::TaxonFilterObj> filter={}) {
        SimilarityMatrix matrix;

        size_t min_shared_length = 1000;

        size_t total_combinations = profiles.size() * profiles.size() - profiles.size() / 2;
        size_t count_combinations = 1;
//        std::cout << "Combination";
        for (auto i = 0; i < profiles.size(); i++) {
            auto& profile1 = profiles[i];
            auto name1 = profile1.GetName();


            for (auto j = i+1; j < profiles.size(); j++) {
                count_combinations++;
//                std::cout << "\rCombination " << count_combinations << " of " << total_combinations;
                auto& profile2 = profiles[j];
                auto name2 = profile2.GetName();

                auto [similarity, shared_length] = GetSimilarity(taxid, profile1, profile2, options, filter);

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

        auto taxids = ExtractTaxa(profiles);

        for (auto& taxid : taxids) {
            std::cout << "Target taxid: " << taxid << " >> " << taxonomy.Get(taxid).scientific_name << std::endl;

            auto similarities = GetSimilarityMatrixForTaxon(taxid, options, profiles, filter);
            std::cout << "Got Similarity matrix" << std::endl;

            if (!similarities.AnySet()) continue;

            if (options.Verbose()) {
                std::cout << "Target taxid: " << taxid << " >> " << taxonomy.Get(taxid).scientific_name << std::endl;
                similarities.PrintMatrix();
                std::cout << "------------" << std::endl;
            }

            std::string name = taxonomy.Get(taxid).scientific_name;
            std::replace(name.begin(), name.end(), ' ', '_');

            WriteDistanceMatrix(taxid, similarities, options, name);

            GetMSAForTaxon(taxid, name, loader, options, profiles);
        }
        std::cout << "Finish StrainWrapper" << std::endl;


    }

    static StrainResults StrainWrapper(Options& options, Profiles& profiles, std::optional<profiler::TaxonFilter> filter={}) {
        std::unordered_set<size_t> taxa1_set;

        auto taxa = ExtractTaxa(profiles, filter);

        StrainResults strain_results;
        auto min_shared_region = 1500;

        auto sample_size = profiles.size();

        double equality_sum = 0;
        size_t equality_count = 0;

        for (auto& profile1 : profiles) {
            auto& taxa_map1 = profile1.GetTaxa();
            taxa1_set = profile1.GetKeySet(filter);

            for (auto& profile2 : profiles) {

                for (auto& [tid, _] : profile2.GetTaxa()) {
                    auto& taxon2 = profile2.GetTaxa().at(tid);

                    bool pass_filter = !filter.has_value() || filter.has_value() && filter->Pass(taxon2);
                    if (!(taxa1_set.contains(tid) && pass_filter)) continue;

                    if (!strain_results.contains(tid)) {
                        strain_results.insert( { tid, DoubleMatrix(sample_size, sample_size, -1.0) });
                        strain_results.at(tid).SetColNames(options.GetPrefixes());
                        strain_results.at(tid).SetRowNames(options.GetPrefixes());
                    }
                    auto& matrix = strain_results.at(tid);



                    auto& taxon1 = profile1.GetTaxa().at(tid);

                    auto [distance, shared_region] = Distance(taxon1, taxon2);

                    auto similarity = 1.0f-distance;
                    if (similarity != similarity) {
                        exit(12);
                    }

                    if (shared_region <= min_shared_region) {
                        similarity = NAN;
                    }


                    if (profile1.GetName() != profile2.GetName() && similarity > 0.999) {
                        equality_sum += similarity;
                        std::cout << "GREP SIMILARITY: " << similarity << std::endl;
                        equality_count++;
                    }


                    matrix.SetValue(profile1.GetName(), profile2.GetName(), similarity, true);
                    if (profile1.GetName() == profile2.GetName()) continue;
                    matrix.SetValue(profile1.GetName(), profile2.GetName(), static_cast<double>(shared_region), false);
                }
            }
        }

        std::cout << std::endl << std::string(60 ,'#') << std::endl;
        std::cout << "StrainLevel Results" << std::endl;
        for (auto& [id, matrix] : strain_results) {
            matrix.PrintMatrix(std::cout, "\t", 6);
            std::ofstream matrix_os(options.GetSimilarityMatrixOutput(std::to_string(id)), std::ios::out);
            std::cout << " Write to : " << options.GetSimilarityMatrixOutput(std::to_string(id)) << std::endl;
            matrix.PrintMatrix(matrix_os, "\t", 10);
            matrix_os.close();
        }
        std::cout << "Equality: " << equality_sum << " (" << equality_count << ")" << std::endl;

        return strain_results;
    }


    static void Run(int argc, char *argv[]) {
        PrintLogo();
        PrintProtalInformation();

        Benchmark bm_total("Run protal");
        bm_total.Start();

        using AlignmentBenchmark = CoreBenchmark;

        auto options = protal::Options::OptionsFromArguments(argc, argv);

        std::cout << "After get options" << std::endl;

        if (options.Help()) {
            options.PrintHelp();
            exit(0);
        }

        std::cout << "Options:\n" << options.ToString() << std::endl;

//        std::unordered_map<std::string, std::string> sample = {
//                {"present_genes", "40"},
//                {"total_hits", "50"},
//                {"unique_hits", "10"},
//                {"expected_gene_presence", "40"},
//                {"mean_ani", "97.9"}
//        };
//        std::string model_path = "/usr/users/QIB_fr017/fritsche/Documents/HPC/group/projects/protal/camisim_all_output/output_flex_r214/profiles/mapq4_2/random_forest_cami_r214.pmml";
//        cpmml::Model model(model_path);
//        std::cout << "Loaded data" << std::endl;
//        auto prediction = stod(model.predict(sample));
//        auto score = model.score(sample);
//        std::cout << "Prediction: " << prediction << std::endl;
//        std::cout << "Score:      " << score.as_double() << std::endl;

        ProtalDB db(options.GetSequenceFile(), options.GetSequenceMapFile());

        if (options.PreloadGenomes()) {
            Benchmark bm_preload_genomes("Preload genomes");
            bm_preload_genomes.Start();
            db.GetGenomes().LoadAllGenomes();
            bm_preload_genomes.Stop();
            bm_preload_genomes.PrintResults();
        }

        bool skip_alignment = Utils::exists(options.SamFile(0)) && !options.Force();
        std::cout << "Skip alignment: " << skip_alignment << std::endl;
        if (!options.BuildMode() && !options.SamFile(0).empty() && (options.ProfileOnly() || skip_alignment)) {
            std::cout << "Goto profile" << std::endl;
            goto Profile;
        }

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
            auto profiles = ProfileWrapper(options, db);
            bm_profiling.Stop();
            bm_profiling.PrintResults();


            // Remove later
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

            if (!options.NoStrains()) {
                std::cout << "Start StrainWrapper" << std::endl;
                StrainWrapper2(options, profiles, db.GetGenomes(), db.GetTaxonomy(), filter);
            }
        }

        bm_total.Stop();
        bm_total.PrintResults();
    }
}
