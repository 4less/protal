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
    static const size_t DEFAULT_OUTPUT_TOP = 3;
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
                ("0,benchmark_alignment", "Benchmark alignment part of protal based on true taxonomic id and gene id supplied in the read header. Header must fulfill the formatting >taxid_geneid... with the regex: >[0-9]+_[0-9]+([^0-9]+.*)*")
                ("c,align_top", "After seeding, anchor are sorted by quality passed to alignment. <take_top> specifies how many anchors should be aligned starting with the most promising anchor.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_ALIGN_TOP)))
                ("e,output_top", "After alignment, alignments are sorted by score. <output_top> specifies how many alignments should be reported starting with the highest scoring alignment.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_OUTPUT_TOP)))
                ("p,preload_genomes", "Preload complete reference library (can be very memory intensive) instead of dynamic loading. This improves performance. [default off]")
                ("reference", "", cxxopts::value<std::string>()->default_value(""));

        return options;
    }


    class Options {
    private:
        bool m_build = false;
        bool m_preload_genomes = false;
        bool m_benchmark_alignment = false;
        bool m_show_help = false;

        std::string m_sequence_file;
        std::string m_database_path;
        std::string m_output_file;
        std::string m_first;
        std::string m_second;
        size_t m_threads = DEFAULT_THREADS;

        size_t m_align_top = DEFAULT_ALIGN_TOP;

        const std::string PROTAL_INDEX_FILE = "index.prx";
        const std::string PROTAL_SEQUENCE_FILE = "reference.fna";
        const std::string PROTAL_SEQUENCE_MAP_FILE = "reference.map";

    public:
        Options(bool build, bool preload_genomes, bool benchmark_alignment, bool show_help, std::string first, std::string second, std::string database_path, std::string output_file, std::string sequence_file, size_t threads, size_t align_top) :
            m_build(build),
            m_preload_genomes(preload_genomes),
            m_show_help(show_help),
            m_first(std::move(first)),
            m_second(std::move(second)),
            m_database_path(std::move(database_path)),
            m_output_file(std::move(output_file)),
            m_sequence_file(std::move(sequence_file)),
            m_threads(threads),
            m_align_top(align_top),
            m_benchmark_alignment(benchmark_alignment) {};

        std::string ToString() const {
            std::string result_str = "";

            result_str += "------ General ------" + std::string(30, '-') + '\n';
            result_str += "build:               " + std::to_string(m_build) + '\n';
            result_str += "threads:             " + std::to_string(m_threads) + '\n';
            result_str += "-------- I/O --------" + std::string(30, '-') + '\n';
            result_str += "first:               " + m_first + '\n';
            result_str += "second:              " + m_second + '\n';
            result_str += "db path:             " + m_database_path + '\n';
            result_str += "sequence file:       " + m_sequence_file + '\n';
            result_str += "output file:         " + m_output_file + '\n';
            result_str += "preload genomes:     " + std::to_string(m_preload_genomes) + '\n';
            result_str += "----- Alignment -----" + std::string(30, '-') + '\n';
            result_str += "align top:       " + std::to_string(m_align_top) + '\n';
            result_str += "benchmark alignment: " + std::to_string(m_benchmark_alignment) + '\n';
//            result_str += "seed num:        " + std::to_string(m_seed_num) + '\n';
//            result_str += "advance by:      " + std::to_string(m_advance_by) + '\n';
            result_str += "---------------------" + std::string(30, '-') + '\n';

            return result_str;
        }

        bool BuildMode() const {
            return m_build;
        }

        bool PreloadGenomes() const {
            return m_preload_genomes;
        }

        bool Help() const {
            return m_show_help;
        }

        bool BenchmarkAlignment() const {
            return m_benchmark_alignment;
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

        std::string GetOutputFile() const {
            return m_output_file;
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

        bool PairedMode() const {
            return !m_first.empty() && !m_second.empty();
        }

        size_t GetThreads() const {
            return m_threads;
        }

        size_t GetAlignTop() const {
            return m_align_top;
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
            bool build = result.count("build");
            bool preload_genomes = result.count("preload_genomes");
            bool benchmark_alignment = result.count("benchmark_alignment");

            size_t threads = result["threads"].as<size_t>();
            size_t align_top = result["align_top"].as<size_t>();

            auto reference = result["reference"].as<std::string>();

            auto first = result.count("first") ? result["first"].as<std::string>() : "";
            auto second = result.count("second") ? result["second"].as<std::string>() : "";

            auto output_file = result.count("output_file") ? result["output_file"].as<std::string>() : "";

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
                    preload_genomes,
                    benchmark_alignment,
                    help,
                    first,
                    second,
                    db_path,
                    output_file,
                    reference,
                    threads,
                    align_top);
        }
    };


}