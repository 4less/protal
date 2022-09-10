//
// Created by fritsche on 14/08/22.
//

#pragma once

#include <cxxopts.hpp>
#include <filesystem>
#include <utility>

namespace protal {
    static const size_t DEFAULT_THREADS = 1;
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
                ("h,help", "Print help.")
                ("p,preload_genomes", "Preload complete reference library (can be very memory intensive) instead of dynamic loading. This improves performance. [default off]")
                ("reference", "", cxxopts::value<std::string>()->default_value(""));

        return options;
    }


    class Options {
    private:
        bool m_build = false;
        bool m_preload_genomes = false;

        std::string m_sequence_file;
        std::string m_database_path;
        std::string m_output_file;
        std::string m_first;
        std::string m_second;
        size_t m_threads = DEFAULT_THREADS;

        const std::string PROTAL_INDEX_FILE = "index.prx";
        const std::string PROTAL_SEQUENCE_FILE = "reference.fna";
        const std::string PROTAL_SEQUENCE_MAP_FILE = "reference.map";

    public:
        Options(bool build, bool preload_genomes, std::string first, std::string second, std::string database_path, std::string output_file, std::string sequence_file, size_t threads) :
            m_build(build),
            m_preload_genomes(preload_genomes),
            m_first(std::move(first)),
            m_second(std::move(second)),
            m_database_path(std::move(database_path)),
            m_output_file(std::move(output_file)),
            m_sequence_file(std::move(sequence_file)),
            m_threads(threads) {};

        std::string ToString() const {
            std::string result_str = "";

            result_str += "build:         " + std::to_string(m_build) + '\n';
            result_str += "threads:       " + std::to_string(m_threads) + '\n';
            result_str += "db path:       " + m_database_path + '\n';
            result_str += "sequence file: " + m_sequence_file;
            return result_str;
        }

        bool BuildMode() const {
            return m_build;
        }

        bool PreloadGenomes() const {
            return m_preload_genomes;
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

        static Options OptionsFromArguments(int argc, char *argv[]) {
            auto cxx_options = CxxOptions();

            cxx_options.parse_positional({ "reference" });
            auto result = cxx_options.parse(argc, argv);

            bool help = result.count("help");
            bool build = result.count("build");
            bool preload_genomes = result.count("preload_genomes");
            size_t threads = result["threads"].as<size_t>();

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
                    std::cerr << "Error " << std::endl;
                } else {
                    db_path = db_path_env;
                }
            }

            return Options(
                    build,
                    preload_genomes,
                    first,
                    second,
                    db_path,
                    output_file,
                    reference,
                    threads);
        }
    };


}