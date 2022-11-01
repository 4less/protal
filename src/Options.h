//
// Created by fritsche on 14/08/22.
//

#pragma once

#include <cxxopts.hpp>
#include <filesystem>
#include <utility>

namespace protal {
    static const size_t DEFAULT_THREADS = 1;
    static const size_t DEFAULT_ALIGN_TOP = 3;
    static const double DEFAULT_MAX_SCORE_ANI = 0.8;
    static const size_t DEFAULT_X_DROP = 1000;
    static const size_t DEFAULT_OUTPUT_TOP = 3;
    static const bool DEFAULT_NO_STRAIN = false;
    static const std::string PROTAL_DB_ENV_VARIABLE = "PROTAL_DB_PATH";

    static cxxopts::Options CxxOptions() {
        cxxopts::Options options(
                "protal",
                "Protal help text");

        options.positional_help("Help")
                .add_options()
                ("b,build", "Build index from reference file with header format ()")
                ("t,threads", "Specify number of threads to use.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_THREADS)))
                ("d,db", "Path to protal database folder.", cxxopts::value<std::string>())
                ("o,output_file", "Output alignment file.", cxxopts::value<std::string>())
                ("1,first", "Comma separated list of reads. If paired-end, also specify second end.", cxxopts::value<std::string>()->default_value(""))
                ("2,second", "Comma separated list of reads. must have <-1/--first> specified.", cxxopts::value<std::string>()->default_value(""))
                ("3,output", "Comma separated list of output prefixes (optional). If not specified, output file prefixes are generated from the input file names.", cxxopts::value<std::string>()->default_value(""))
                ("h,help", "Print help.")
                ("f,fastalign", "Speed up alignment by allowing approximate alignment between seeds with no indication for indels.")
                ("n,no_strains", "Stray on species level. Do not output SNPs")
                ("0,benchmark_alignment", "Benchmark alignment part of protal based on true taxonomic id and gene id supplied in the read header. Header must fulfill the formatting >taxid_geneid... with the regex: >[0-9]+_[0-9]+([^0-9]+.*)*")
                ("9,benchmark_alignment_output", "Benchmark alignment output. Output is appended to the file.", cxxopts::value<std::string>())
                ("c,align_top", "After seeding, anchor are sorted by quality passed to alignment. <take_top> specifies how many anchors should be aligned starting with the most promising anchor.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_ALIGN_TOP)))
                ("a,max_score_ani", "A max score makes an alignment stop if the alignment diverges too much. This parameter estimates the score for a given ani and is a tradeoff between speed/accuracy. [ Default: " + std::to_string(DEFAULT_MAX_SCORE_ANI) + "]", cxxopts::value<double>()->default_value(std::to_string(DEFAULT_MAX_SCORE_ANI)))
                ("x,x_drop", "Value determines when to cut of branches in the aligment process that are unpromising. [ Default: " + std::to_string(DEFAULT_X_DROP) + "]", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_X_DROP)))
                ("e,output_top", "After alignment, alignments are sorted by score. <output_top> specifies how many alignments should be reported starting with the highest scoring alignment.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_OUTPUT_TOP)))
                ("g,preload_genomes", "Preload complete reference library (can be very memory intensive) instead of dynamic loading. This improves performance. [default off]")
                ("p,profile", "Perform taxonomic profiling.")
                ("q,profile_only", "Provide profile filename (.sam) and only perform profiling based on sam file.", cxxopts::value<std::string>()->default_value(""))
                ("8,profile_truth", "Provide truth file and annotate profile taxa with TP/FP. Format is list of integers (internal ids)", cxxopts::value<std::string>()->default_value(""))
                ("reference", "", cxxopts::value<std::string>()->default_value(""));

        return options;
    }


    class Options {
    private:
        bool m_build = false;
        bool m_profile = false;
        bool m_preload_genomes = false;
        bool m_benchmark_alignment = false;
        bool m_show_help = false;
        bool m_no_strains = false;
        bool m_fastalign = false;
        bool m_profile_only = false;


        std::string m_sequence_file;
        std::string m_database_path;
        std::string m_output_prefix;
        std::string m_first;
        std::string m_second;
        std::string m_sam_file;

        std::string m_profile_truth;

        size_t m_threads = DEFAULT_THREADS;

        std::string m_benchmark_alignment_output;

        size_t m_align_top = DEFAULT_ALIGN_TOP;
        double m_max_score_ani = DEFAULT_MAX_SCORE_ANI;
        size_t m_x_drop = DEFAULT_X_DROP;

        const std::string PROTAL_INDEX_FILE = "index.prx";
        const std::string PROTAL_SEQUENCE_FILE = "reference.fna";
        const std::string PROTAL_SEQUENCE_MAP_FILE = "reference.map";
        const std::string PROTAL_TAXONOMY_FILE = "internal_taxonomy.dmp";
        const std::string PROTAL_TAXONOMY_GTDB_FILE = "internal_taxonomy_gtdb.dmp";
        const std::string PROTAL_TAXONOMY_NCBI_FILE = "internal_taxonomy_ncbi.dmp";

    public:
        Options(bool build, bool profile, bool profile_only, bool no_strains, bool preload_genomes, bool benchmark_alignment,
                std::string benchmark_alignment_output, bool show_help, std::string first,
                std::string second, std::string database_path, std::string output_prefix,
                std::string sequence_file, size_t threads, size_t align_top, double max_score_ani,
                size_t x_drop, bool fastalign, std::string sam_file, std::string profile_truth
                ) :
                m_build(build),
                m_profile(profile),
                m_profile_only(profile_only),
                m_no_strains(no_strains),
                m_preload_genomes(preload_genomes),
                m_show_help(show_help),
                m_first(std::move(first)),
                m_second(std::move(second)),
                m_database_path(std::move(database_path)),
                m_output_prefix(std::move(output_prefix)),
                m_sam_file(sam_file),
                m_sequence_file(std::move(sequence_file)),
                m_threads(threads),
                m_align_top(align_top),
                m_max_score_ani(max_score_ani),
                m_x_drop(x_drop),
                m_fastalign(fastalign),
                m_profile_truth(profile_truth),
                m_benchmark_alignment(benchmark_alignment),
                m_benchmark_alignment_output(benchmark_alignment_output) {

            if (sam_file.empty() && !m_output_prefix.empty()) {
                m_sam_file = m_output_prefix + ".sam";
            }
            if (!output_prefix.empty()) return;

            if (!m_sam_file.empty()) {
                if (m_sam_file.ends_with(".sam")) {
                    m_output_prefix = m_sam_file.substr(0, m_sam_file.length() - 4);
                }
            }
        };

        std::string ToString() const {
            std::string result_str = "";

            result_str += "------ General ------" + std::string(30, '-') + '\n';
            result_str += "build:               " + std::to_string(m_build) + '\n';
            result_str += "no strains:          " + std::to_string(m_no_strains) + '\n';
            result_str += "threads:             " + std::to_string(m_threads) + '\n';
            result_str += "-------- I/O --------" + std::string(30, '-') + '\n';
            result_str += "first:               " + m_first + '\n';
            result_str += "second:              " + m_second + '\n';
            result_str += "db path:             " + m_database_path + '\n';
            result_str += "sequence file:       " + m_sequence_file + '\n';
            result_str += "profile file:        " + m_sam_file + '\n';
            result_str += "output prefix:       " + m_output_prefix + '\n';
            result_str += "preload genomes:     " + std::to_string(m_preload_genomes) + '\n';
            result_str += "----- Alignment -----" + std::string(30, '-') + '\n';
            result_str += "align top:           " + std::to_string(m_align_top) + '\n';
            result_str += "max score ani:       " + std::to_string(m_max_score_ani) + '\n';
            result_str += "x_drop:              " + std::to_string(m_x_drop) + '\n';
            result_str += "fastalign:           " + std::to_string(m_fastalign) + '\n';
            result_str += "---- Dev Options ----" + std::string(30, '-') + '\n';
            result_str += "benchmark alignment: " + std::to_string(m_benchmark_alignment) + '\n';
            result_str += "profile truth:       " + m_profile_truth + '\n';
            result_str += "^ output:            " + m_benchmark_alignment_output + '\n';
            result_str += "---------------------" + std::string(30, '-') + '\n';
            return result_str;
        }

        bool BuildMode() const {
            return m_build;
        }

        bool Profile() const {
            return m_profile;
        }

        bool ProfileOnly() const {
            return m_profile_only;
        }

        std::string SamFile() const {
            return m_sam_file;
        }

        std::string& SamFile() {
            return m_sam_file;
        }

        std::string& ProfileTruthFile() {
            return m_profile_truth;
        }

        bool PreloadGenomes() const {
            return m_preload_genomes;
        }

        bool Help() const {
            return m_show_help;
        }

        bool NoStrains() const {
            return m_no_strains;
        }

        bool BenchmarkAlignment() const {
            return m_benchmark_alignment;
        }

        bool FastAlign() const {
            return m_fastalign;
        }

        std::string GetSequenceFilePath() const {
            return m_sequence_file;
        }

        std::string GetIndexFolder() const {
            return m_database_path;
        }

        std::string GetIndexFile() const {
            return m_database_path + "/" + PROTAL_INDEX_FILE;
        }

        std::string GetInternalTaxonomyFile() const {
            return m_database_path + "/" + PROTAL_TAXONOMY_FILE;
        }

        std::string GetInternalTaxonomyGTDBFile() const {
            return m_database_path + "/" + PROTAL_TAXONOMY_GTDB_FILE;
        }

        std::string GetInternalTaxonomyNCBIFile() const {
            return m_database_path + "/" + PROTAL_TAXONOMY_NCBI_FILE;
        }

        std::string GetOutputPrefix() const {
            return m_output_prefix;
        }

        std::string GetSequenceFile() const {
            return m_database_path + "/" + PROTAL_SEQUENCE_FILE;
        }

        std::string GetSequenceMapFile() const {
            return m_database_path + "/" + PROTAL_SEQUENCE_MAP_FILE;
        }

        std::string GetFirstFile() const {
            return m_first;
        }

        std::string GetSecondFile() const {
            return m_second;
        }

        std::string GetBenchmarkAlignmentOutputFile() const {
            return m_benchmark_alignment_output;
        }

        const std::string& GetBenchmarkAlignmentOutputFile() {
            return m_benchmark_alignment_output;
        }

        bool PairedMode() const {
            return !m_first.empty() && !m_second.empty();
        }

        size_t GetThreads() const {
            return m_threads;
        }

        size_t GetAlignTop() const {
            return m_align_top;
        }

        size_t GetXDrop() const {
            return m_x_drop;
        }

        double GetMaxScoreAni() const {
            return m_max_score_ani;
        }

        void PrintHelp() {
            auto cxx_options = CxxOptions();
            std::cout << cxx_options.help() << std::endl;
        }

        static Options OptionsFromArguments(int argc, char *argv[]) {
            auto cxx_options = CxxOptions();

            cxx_options.parse_positional({ "reference" });
            auto result = cxx_options.parse(argc, argv);

            bool help = result.count("help");
            bool no_strains = result.count("no_strains");
            bool build = result.count("build");
            bool preload_genomes = result.count("preload_genomes");
            bool benchmark_alignment = result.count("benchmark_alignment");
            bool fastalign = result.count("fastalign");
            bool profile = result.count("profile");

            size_t threads = result["threads"].as<size_t>();
            size_t align_top = result["align_top"].as<size_t>();
            size_t x_drop = result["x_drop"].as<size_t>();
            double max_score_ani = result["max_score_ani"].as<double>();


            auto reference = result["reference"].as<std::string>();

            auto first = result.count("first") ? result["first"].as<std::string>() : "";
            auto second = result.count("second") ? result["second"].as<std::string>() : "";

            auto sam_in = result.count("profile_only") ? result["profile_only"].as<std::string>() : "";
            auto profile_truth = result.count("profile_truth") ? result["profile_truth"].as<std::string>() : "";
            bool profile_only = result.count("profile_only");

            auto output_file = result.count("output_file") ? result["output_file"].as<std::string>() : "";
            auto benchmark_alignment_output_file = result.count("benchmark_alignment_output") ? result["benchmark_alignment_output"].as<std::string>() : "";

            std::string db_path;
            if (result.count("db")) {
                db_path = result["db"].as<std::string>();
            } else {
                auto db_path_env = std::getenv(PROTAL_DB_ENV_VARIABLE.c_str());
                if (!db_path_env) {
                    std::cerr << "Error " << PROTAL_DB_ENV_VARIABLE << std::endl;
                } else {
                    db_path = db_path_env;
                }
            }

            return Options(
                    build,
                    profile,
                    profile_only,
                    no_strains,
                    preload_genomes,
                    benchmark_alignment,
                    benchmark_alignment_output_file,
                    help,
                    first,
                    second,
                    db_path,
                    output_file,
                    reference,
                    threads,
                    align_top,
                    max_score_ani,
                    x_drop,
                    fastalign,
                    sam_in,
                    profile_truth);
        }
    };


}