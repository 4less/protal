//
// Created by fritsche on 14/08/22.
//

#pragma once

#include <cxxopts.hpp>
#include <filesystem>
#include <utility>
#include "LineSplitter.h"
#include <fstream>
#include <regex>


namespace protal {
    static const size_t DEFAULT_THREADS = 1;
    static const size_t DEFAULT_ALIGN_TOP = 3;
    static const double DEFAULT_MAX_SCORE_ANI = 0.9;
    static const double DEFAULT_MSA_MIN_VCOV = 0.5;
    static const size_t DEFAULT_X_DROP = 1000;
    static const size_t DEFAULT_MAX_KEY_UBIQUITY = 256;
    static const size_t DEFAULT_MAX_SEED_SIZE = 128;
    static const size_t DEFAULT_MAX_OUT = 1;
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
                ("v,verbose", "Have verbose program output")
                ("t,threads", "Specify number of threads to use.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_THREADS)))
                ("d,db", "Path to protal database folder.", cxxopts::value<std::string>())
                ("o,outdir", "Comma separated list of output prefixes (optional). If not specified, output file prefixes are generated from the input file names.", cxxopts::value<std::string>())
                ("1,first", "Comma separated list of reads. If paired-end, also specify second end.", cxxopts::value<std::string>()->default_value(""))
                ("2,second", "Comma separated list of reads. must have <-1/--first> specified.", cxxopts::value<std::string>()->default_value(""))
                ("3,prefix", "Comma separated list of output prefixes (optional). If not specified, output file prefixes are generated from the input file names.", cxxopts::value<std::string>()->default_value(""))
                ("4,map", "For larger datasets you can define parameters -1, -2, -3 and -o in a file.", cxxopts::value<std::string>()->default_value(""))
                ("5,map_help", "Get help how to format the map file.")
                ("6,map_range", "If you specified a map file with -4 or --map you can also pass a range to protal to run protal only on a subset. The first entry is 1, the end is inclusive. e.g.: 1-10. If the end open or larger than the number of entries in the map file, the last entry in the map file is selected as end.", cxxopts::value<std::string>()->default_value("1-"))
                ("7,mapq_debug_output", "Output mapq debug info to stderr")
                ("8,profile_truth", "Provide truth file and annotate profile taxa with TP/FP. Format is list of integers (internal ids)", cxxopts::value<std::string>()->default_value(""))
                ("9,benchmark_alignment_output", "Benchmark alignment output. Output is appended to the file.", cxxopts::value<std::string>())
                ("h,help", "Print help.")
                ("z,force", "Force redo alignment even if sam files exists.")
                ("f,fastalign", "Speed up alignment by allowing approximate alignment between seeds with no indication for indels.")
                ("n,no_strains", "Stray on species level. Do not output SNPs")
                ("0,benchmark_alignment", "Benchmark alignment part of protal based on true taxonomic id and gene id supplied in the read header. Header must fulfill the formatting >taxid_geneid... with the regex: >[0-9]+_[0-9]+([^0-9]+.*)*")

                ("c,align_top", "After seeding, anchor are sorted by quality passed to alignment. <take_top> specifies how many anchors should be aligned starting with the most promising anchor.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_ALIGN_TOP)))
                ("m,max_out", "Maximum alignments that should be outputted", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_MAX_OUT)))
                ("u,max_key_ubiquity", "Max key ubiquity. Best matching Flexkey count for seed must be lower or equal", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_MAX_KEY_UBIQUITY)))
                ("s,max_seed_size", "Max seed size after which seeding is stopped.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_MAX_SEED_SIZE)))
                ("a,max_score_ani", "A max score makes an alignment stop if the alignment diverges too much. This parameter estimates the score for a given ani and is a tradeoff between speed/accuracy. [ Default: " + std::to_string(DEFAULT_MAX_SCORE_ANI) + "]", cxxopts::value<double>()->default_value(std::to_string(DEFAULT_MAX_SCORE_ANI)))
                ("x,x_drop", "Value determines when to cut of branches in the aligment process that are unpromising. [ Default: " + std::to_string(DEFAULT_X_DROP) + "]", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_X_DROP)))
                ("e,output_top", "After alignment, alignments are sorted by score. <output_top> specifies how many alignments should be reported starting with the highest scoring alignment.", cxxopts::value<size_t>()->default_value(std::to_string(DEFAULT_OUTPUT_TOP)))
                ("g,preload_genomes_off", "Do not preload complete reference library (reference.fna and reference.map in protal index folder) and instead do dynamic loading. This usually decreases performance but saves memory.")
                ("k,msa_min_vcov", "Protal outputs two MSAs. The processed MSA is condensed horizontally such that each position in the MSA is covered by at least msa_min_cov percent of the sequences with bases that are neither '-' nor 'N'", cxxopts::value<double>()->default_value(std::to_string(DEFAULT_MSA_MIN_VCOV)))
                ("p,profile", "Perform taxonomic profiling.")
                ("i,full_reference", "", cxxopts::value<std::string>()->default_value(""))
                ("q,profile_only", "Provide profile filename (.sam) and only perform profiling based on sam file.", cxxopts::value<std::string>()->default_value(""))
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
        bool m_force = false;
        bool m_verbose = false;
        bool m_gzip_sam = true;

        bool m_mapq_debug_out = false;

        size_t m_current_index = 0;


        std::string m_sequence_file;
        std::string m_full_sequence_file;
        std::string m_database_path;
        std::string m_output_dir;
        std::string m_map_file;

        std::vector<std::string> m_first_list;
        std::vector<std::string> m_second_list;
        std::vector<std::string> m_prefix_list;
        std::vector<std::string> m_sam_list;
        std::vector<std::string> m_profile_list;
        std::vector<std::string> m_sampleid_list;
        std::vector<std::string> m_profile_truth_list;

        std::vector<size_t> m_range;

        std::string m_profile_truth;

        size_t m_threads = DEFAULT_THREADS;

        std::string m_benchmark_alignment_output;

        size_t m_align_top = DEFAULT_ALIGN_TOP;
        double m_max_score_ani = DEFAULT_MAX_SCORE_ANI;
        size_t m_x_drop = DEFAULT_X_DROP;
        size_t m_max_key_ubiquity = DEFAULT_MAX_KEY_UBIQUITY;
        size_t m_max_seed_size = DEFAULT_MAX_SEED_SIZE;
        size_t m_max_out = DEFAULT_MAX_OUT;
        double m_msa_min_vcov = DEFAULT_MSA_MIN_VCOV;

    public:
        static inline const std::string PROTAL_INDEX_FILE = "index.prx";
        static inline const std::string PROTAL_SEQUENCE_FILE = "reference.fna";
        static inline const std::string PROTAL_SEQUENCE_MAP_FILE = "reference.map";
        static inline const std::string PROTAL_HITTABLE_GENES_FILE = "species_gene_mask.tsv";
        static inline const std::string PROTAL_UNIQUE_KMER_FILE = "unique_kmers.tsv";
        static inline const std::string PROTAL_TAXONOMY_FILE = "internal_taxonomy.dmp";

        // DEPRECATED
        static inline const std::string PROTAL_TAXONOMY_GTDB_FILE = "internal_taxonomy_gtdb.dmp";
        static inline const std::string PROTAL_TAXONOMY_NCBI_FILE = "internal_taxonomy_ncbi.dmp";

        static inline const std::string MAP_SAMPLEID = "#SAMPLEID";
        static inline const std::string MAP_FIRST_READ = "FIRST";
        static inline const std::string MAP_SECOND_READ = "SECOND";
        static inline const std::string MAP_SAM = "SAM";
        static inline const std::string MAP_PROFILE = "PROFILE";
        static inline const std::string MAP_PROFILE_TRUTH = "PROFILE_TRUTH";
        static inline const std::string MAP_HEADER_TO_SPECIES = "HEADER_TO_SPECIES";
        static inline const std::string MAP_PREFIX = "PREFIX";

        static inline const std::string MAP_VAR_INPUT_DIR = "#INPUT_DIR";
        static inline const std::string MAP_VAR_OUTPUT_DIR = "#OUTPUT_DIR";
        static inline const std::string MAP_VAR_SAM_OUTPUT_DIR = "#SAM_OUTPUT_DIR";
        static inline const std::string MAP_VAR_PROFILE_OUTPUT_DIR = "#PROFILE_OUTPUT_DIR";
        static inline const std::string MAP_VAR_STRAIN_OUTPUT_DIR = "#STRAIN_OUTPUT_DIR";

        const size_t MAP_SAMPLE_ID_COL = 0;

        Options() : m_show_help(true) {}

        Options(bool build, bool profile, bool profile_only, bool no_strains, bool preload_genomes, bool benchmark_alignment,
                std::string benchmark_alignment_output, bool show_help, bool mapq_debug_output,
                std::vector<std::string>& first_list, std::vector<std::string>& second_list, std::vector<std::string>& samplename_list,
                std::string database_path, std::vector<std::string>& output_prefix_list, std::string full_sequence_file,
                std::string sequence_file, std::string map_file, std::string& output_dir, size_t threads, size_t align_top, size_t max_out, double max_score_ani,
                double msa_min_vcov, size_t x_drop, size_t max_key_ubiquity, size_t max_seed_size, bool fastalign, std::vector<std::string>& sam_file_list,
                std::vector<std::string>& profile_file_list, std::vector<std::string>& profile_truth_list, std::string profile_truth,
                bool force, bool verbose, std::vector<size_t> range) :
                m_build(build),
                m_profile(profile),
                m_profile_only(profile_only),
                m_no_strains(no_strains),
                m_preload_genomes(preload_genomes),
                m_show_help(show_help),
                m_first_list(std::move(first_list)),
                m_second_list(std::move(second_list)),
                m_database_path(std::move(database_path)),
                m_prefix_list(std::move(output_prefix_list)),
                m_output_dir(output_dir),
                m_sam_list(sam_file_list),
                m_profile_list(profile_file_list),
                m_profile_truth_list(profile_truth_list),
                m_sequence_file(std::move(full_sequence_file)),
                m_full_sequence_file(std::move(sequence_file)),
                m_map_file(std::move(map_file)),
                m_threads(threads),
                m_align_top(align_top),
                m_max_out(max_out),
                m_max_score_ani(max_score_ani),
                m_msa_min_vcov(msa_min_vcov),
                m_x_drop(x_drop),
                m_max_key_ubiquity(max_key_ubiquity),
                m_max_seed_size(max_seed_size),
                m_fastalign(fastalign),
                m_profile_truth(std::move(profile_truth)),
                m_benchmark_alignment(benchmark_alignment),
                m_benchmark_alignment_output(benchmark_alignment_output),
                m_mapq_debug_out(mapq_debug_output),
                m_force(force),
                m_verbose(verbose),
                m_range(range) {
            if (samplename_list.empty()) {
                m_sampleid_list = m_prefix_list;
            } else {
                m_sampleid_list = samplename_list;
            }

            if (m_profile_list.empty()) {
                for (auto i = 0; i < m_prefix_list.size(); i++) {
                    auto& prefix = m_prefix_list[i];
                    m_profile_list.emplace_back(prefix + ".profile");
                    std::cout << m_profile_list.size() << " " << prefix + ".profile" << std::endl;
                }
            }

//            if (m_sam_list.empty() && !m_prefix_list.empty()) {
//                for (auto i = 0; i < m_sam_list.size(); i++) {
//                    m_sam_list[i].
//                }
//                m_sam_file = m_output_prefix + ".sam";
//            }
//            if (!output_prefix.empty()) return;
//
//            if (!m_sam_file.empty()) {
//                if (m_sam_file.ends_with(".sam")) {
//                    m_output_prefix = m_sam_file.substr(0, m_sam_file.length() - 4);
//                }
//            }
        };

        std::string ToString() const {
            std::string result_str = "";

            std::string first_list_str = m_first_list.empty() ? "" : m_first_list.front();
            for (auto i = 1; i < m_first_list.size(); i++) first_list_str += ", " + m_first_list[i];
            std::string second_list_str = m_second_list.empty() ? "" : m_second_list.front();
            for (auto i = 1; i < m_second_list.size(); i++) second_list_str += ", " + m_second_list[i];
            std::string sam_list_str = m_sam_list.empty() ? "" : m_sam_list.front();
            for (auto i = 1; i < m_sam_list.size(); i++) sam_list_str += ", " + m_sam_list[i];
            std::string prefix_list_str = m_prefix_list.empty() ? "" : m_prefix_list.front();
            for (auto i = 1; i < m_prefix_list.size(); i++) prefix_list_str += ", " + m_prefix_list[i];
            std::string profile_list_str = m_profile_list.empty() ? "" : m_profile_list.front();
            for (auto i = 1; i < m_profile_list.size(); i++) profile_list_str += ", " + m_profile_list[i];

            for (auto& file : m_profile_list) std::cout << file << std::endl;

            result_str += "------ General ------" + std::string(30, '-') + '\n';
            result_str += "build:               " + std::to_string(m_build) + '\n';
            result_str += "no strains:          " + std::to_string(m_no_strains) + '\n';
            result_str += "threads:             " + std::to_string(m_threads) + '\n';
            result_str += "-------- I/O --------" + std::string(30, '-') + '\n';
            result_str += "first:               " + (first_list_str.length() > 50 ? std::to_string(m_first_list.size()) + " files" : first_list_str) + '\n';
            result_str += "second:              " + (second_list_str.length() > 50 ? std::to_string(m_second_list.size()) + " files" : second_list_str) + '\n';
            result_str += "db path:             " + m_database_path + '\n';
            result_str += "sequence file:       " + m_sequence_file + '\n';
            result_str += "sam file:            " + (sam_list_str.length() > 50 ? std::to_string(m_sam_list.size()) + " files" : sam_list_str) + '\n';
            result_str += "profile file:        " + (profile_list_str.length() > 50 ? std::to_string(m_profile_list.size()) + " files" : profile_list_str) + '\n';
            result_str += "output prefix:       " + (sam_list_str.length() > 50 ? std::to_string(m_sam_list.size()) + " files" : sam_list_str) + '\n';
            result_str += "Output dir:          " + m_output_dir + '\n';
            result_str += "preload genomes:     " + std::to_string(m_preload_genomes) + '\n';
            result_str += "----- Alignment -----" + std::string(30, '-') + '\n';
            result_str += "align top:           " + std::to_string(m_align_top) + '\n';
            result_str += "max key ubiquity:    " + std::to_string(m_max_key_ubiquity) + '\n';
            result_str += "max seed size:       " + std::to_string(m_max_seed_size) + '\n';
            result_str += "max score ani:       " + std::to_string(m_max_score_ani) + '\n';
            result_str += "x-drop:              " + std::to_string(m_x_drop) + '\n';
            result_str += "fastalign:           " + std::to_string(m_fastalign) + '\n';
            result_str += "max out:             " + std::to_string(m_max_out) + '\n';
            result_str += "---- Dev Options ----" + std::string(30, '-') + '\n';
            result_str += "benchmark alignment: " + std::to_string(m_benchmark_alignment) + '\n';
            result_str += "profile truth:       " + std::to_string(HasProfileTruths()) + '\n';
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

//        std::string SamFile() const {
//            return m_sam_file;
//        }
//
//        std::string& SamFile() {
//            return m_sam_file;
//        }

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

        std::string GetFullSequenceFilePath() const {
            return m_full_sequence_file;
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

        std::string GetSequenceFile() const {
            return m_database_path + "/" + PROTAL_SEQUENCE_FILE;
        }

        std::string GetSequenceMapFile() const {
            return m_database_path + "/" + PROTAL_SEQUENCE_MAP_FILE;
        }

        std::string GetHittableGenesMap() const {
            return m_database_path + "/" + PROTAL_HITTABLE_GENES_FILE;
        }

        std::string GetUniqueKmersFile() const {
            return m_database_path + "/" + PROTAL_UNIQUE_KMER_FILE;
        }

        bool HittableGenesMapExists() const {
            return Utils::exists(GetHittableGenesMap());
        }

        bool UniqueKmersFileExists() const {
            return Utils::exists(GetUniqueKmersFile());
        }

        void SetCurrentIndex(size_t i) {
            if (i >= m_prefix_list.size()) {
                std::cerr << "SetCurrentIndex to " << i << " not possible." << std::endl;
                std::cerr << "Max is " << m_prefix_list.size() << std::endl;
                exit(9);
            }
            m_current_index = i;
        }

        size_t GetCurrentIndex() const {
            return m_current_index;
        }

        std::vector<size_t> GetRange() const {
            return m_range;
        }

//        std::string GetOutputPrefix() const {
//            return m_output_prefix;
//        }
//
//        std::string GetFirstFile() const {
//            return m_first;
//        }
//
//        std::string GetSecondFile() const {
//            return m_second;
//        }

        std::vector<std::string> GetFirstFiles() const {
            return m_first_list;
        }

        std::vector<std::string> GetSecondFiles() const {
            return m_second_list;
        }

        std::vector<std::string> GetPrefixes() const {
            return m_prefix_list;
        }

        size_t GetFileCount() const {
            return m_prefix_list.size();
        }

        std::string GetFirstFile(int index) const {
            if (index > m_first_list.size()) {
                std::cerr << "Cannot access index " << index << " of first files (Length: " << m_first_list.size() << ")" << std::endl;
                exit(33);
            }
            return m_first_list[index];
        }

        std::string GetSecondFile(int index) const {
            if (index > m_second_list.size()) {
                std::cerr << "Cannot access index " << index << " of second files (Length: " << m_second_list.size() << ")" << std::endl;
                exit(33);
            }
            return m_second_list[index];
        }

        std::string SamFile(int index) const {
            if (index >= m_sam_list.size()) {
                std::cerr << "Cannot access index " << index << " of sam files (Length: " << m_sam_list.size() << ")" << std::endl;
                exit(33);
            }
            return m_sam_list[index];
        }

        std::vector<std::string> SamFiles() {
            return m_sam_list;
        }

        std::string ProfileFile(int index) const {
            if (index >= m_profile_list.size()) {
                std::cerr << "Cannot access index " << index << " of profile files (Length: " << m_profile_list.size() << ")" << std::endl;
                exit(33);
            }
            return m_profile_list[index];
        }

        std::string ProfileTruthFile(int index) const {
            if (index >= m_profile_truth_list.size()) {
                std::cerr << "Cannot access index " << index << " of profile files (Length: " << m_profile_truth_list.size() << ")" << std::endl;
                exit(33);
            }
            return m_profile_truth_list[index];
        }

        std::string GetPrefix(int index) const {
            if (index > m_prefix_list.size()) {
                std::cerr << "Cannot access index " << index << " of prefix files (Length: " << m_prefix_list.size() << ")" << std::endl;
                exit(33);
            }
            return m_prefix_list[index];
        }

        std::string GetSampleId(int index) const {
            if (index > m_sampleid_list.size()) {
                std::cerr << "Cannot access index " << index << " of sample names (Length: " << m_sampleid_list.size() << ")" << std::endl;
                exit(33);
            }
            return m_sampleid_list[index];
        }

        std::string GetOutputDir() const {
            return m_output_dir;
        }

        std::string GetSimilarityMatrixOutput(std::string species_name) const {
            return m_output_dir + '/' + species_name + ".tsv";
        }

        std::string GetMSAOutput(std::string species_name) const {
            return m_output_dir + '/' + species_name + ".msa.fna";
        }

        std::string GetSpeciesMetaOutput(std::string species_name) const {
            return m_output_dir + '/' + species_name + ".meta.tsv";
        }

        std::string GetMSAProcessedOutput(std::string species_name) const {
            return m_output_dir + '/' + species_name + ".processed.msa.fna";
        }

        std::string GetMSAPartitionOutput(std::string species_name) const {
            return m_output_dir + '/' + species_name + ".partition.txt";
        }

        std::string GetBenchmarkAlignmentOutputFile() const {
            return m_benchmark_alignment_output;
        }

        const std::string& GetBenchmarkAlignmentOutputFile() {
            return m_benchmark_alignment_output;
        }

        bool PairedMode() const {
            return !m_first_list.empty() && !m_second_list.empty();
        }

        bool Force() const {
            return m_force;
        }

        bool Verbose() const {
            return m_verbose;
        }

        size_t GetThreads() const {
            return m_threads;
        }

        bool GetMAPQDebugOut() const {
            return m_mapq_debug_out;
        }

        bool HasProfileTruths() const {
            return m_profile_truth_list.size() == m_profile_list.size();
        }

        auto GetMSAMinVCOV() {
            return m_msa_min_vcov;
        }

        size_t GetAlignTop() const {
            return m_align_top;
        }

        size_t GetMaxOut() const {
            return m_max_out;
        }

        size_t GetMaxSeedSize() const {
            return m_max_seed_size;
        }

        size_t GetXDrop() const {
            return m_x_drop;
        }

        size_t GetMaxKeyUbiquity() const {
            return m_max_key_ubiquity;
        }

        double GetMaxScoreAni() const {
            return m_max_score_ani;
        }

        void PrintHelp() {
            auto cxx_options = CxxOptions();
            std::cout << cxx_options.help() << std::endl;
        }


        static bool LoadFromMap(std::string map_path, std::string& strain_output_dir, std::vector<std::string>& prefix_list,
                                std::vector<std::string>& first_list, std::vector<std::string>& second_list,
                                std::vector<std::string>& sam_list, std::vector<std::string>& profile_list,
                                std::vector<std::string>& samplenames_list, std::vector<std::string>& profile_truth_list) {
            using namespace std::filesystem;

            if (!std::filesystem::exists(map_path)) {
                std::cerr << "Map file " << map_path << " does not exist." << std::endl;
                return false;
            }
            std::ifstream is(map_path, std::ios::in);

            std::string prefix_output_dir = "";
            std::string sam_output_dir = "";
            std::string profile_output_dir = "";
            std::string strain_output_dir2 = "";
            std::string input_dir = "";

            LineSplitter splitter;

            int sample_id_column = 0;
            int first_column = -1;
            int second_column = -1;
            int prefix_column = -1;
            int sam_column = -1;
            int profile_column = -1;
            int profile_truth_column = -1;
            int header_to_species_column = -1;

            bool header = true;
            size_t line_num = 0;
            for (std::string line; std::getline(is, line);) {
                if (line.empty()) continue;
                splitter.Split(line);
                auto& tokens = splitter.Tokens();
                // Header
                if (line.starts_with('#')) {
                    // Check if header is expected
                    if (!header) {
                        std::cerr << "Line " << line_num << ": Did not expect header line but line starts with #" << std::endl;
                        return false;
                    }

                    // Check if variable definition or header
                    if (!tokens.empty() && tokens[0] == MAP_SAMPLEID) {
                        splitter.Split(line);

                        for (auto i = 1; i < splitter.Tokens().size(); i++) {
                            auto& token = splitter.Tokens()[i];
                            if (token == MAP_SAMPLEID) continue;
                            if (token == MAP_FIRST_READ) {
                                if (first_column != -1) {
                                    std::cerr << "Column '" << MAP_FIRST_READ << "' is defined twice" << std::endl;
                                    return false;
                                }
                                first_column = i;
                            }
                            if (token == MAP_SECOND_READ) {
                                if (second_column != -1) {
                                    std::cerr << "Column '" << MAP_SECOND_READ << "' is defined twice" << std::endl;
                                    return false;
                                }
                                second_column = i;
                            }
                            if (token == MAP_PREFIX) {
                                if (prefix_column != -1) {
                                    std::cerr << "Column '" << MAP_PREFIX << "' is defined twice" << std::endl;
                                    return false;
                                }
                                prefix_column = i;
                            }
                            if (token == MAP_SAM) {
                                if (sam_column != -1) {
                                    std::cerr << "Column '" << MAP_SAM << "' is defined twice" << std::endl;
                                    return false;
                                }
                                sam_column = i;
                            }
                            if (token == MAP_PROFILE) {
                                if (profile_column != -1) {
                                    std::cerr << "Column '" << MAP_PROFILE << "' is defined twice" << std::endl;
                                    return false;
                                }
                                profile_column = i;
                            }
                            if (token == MAP_PROFILE_TRUTH) {
                                if (profile_truth_column != -1) {
                                    std::cerr << "Column '" << MAP_PROFILE_TRUTH << "' is defined twice" << std::endl;
                                    return false;
                                }
                                profile_truth_column = i;
                            }
                            if (token == MAP_HEADER_TO_SPECIES) {
                                if (header_to_species_column != -1) {
                                    std::cerr << "Column '" << MAP_HEADER_TO_SPECIES << "' is defined twice" << std::endl;
                                    return false;
                                }
                                header_to_species_column = i;
                            }
                        }
                        header = false;
                    }

                    if (!tokens.empty() && tokens[0] == MAP_VAR_INPUT_DIR) {
                        if (tokens.size() < 2) {
                            std::cerr << "Line " << line_num << ": Expected value for key " << MAP_VAR_INPUT_DIR << std::endl;
                            return false;
                        }
                        input_dir = tokens[1];
                    }
                    if (!tokens.empty() && tokens[0] == MAP_VAR_STRAIN_OUTPUT_DIR) {
                        if (tokens.size() < 2) {
                            std::cerr << "Line " << line_num << ": Expected value for key " << MAP_VAR_STRAIN_OUTPUT_DIR << std::endl;
                            return false;
                        }
                        strain_output_dir = tokens[1];
                        strain_output_dir2 = tokens[1];
                    }
                    if (!tokens.empty() && tokens[0] == MAP_VAR_SAM_OUTPUT_DIR) {
                        if (tokens.size() < 2) {
                            std::cerr << "Line " << line_num << ": Expected value for key " << MAP_VAR_SAM_OUTPUT_DIR << std::endl;
                            return false;
                        }
                        sam_output_dir = tokens[1];
                    }
                    if (!tokens.empty() && tokens[0] == MAP_VAR_PROFILE_OUTPUT_DIR) {
                        if (tokens.size() < 2) {
                            std::cerr << "Line " << line_num << ": Expected value for key " << MAP_VAR_PROFILE_OUTPUT_DIR << std::endl;
                            return false;
                        }
                        profile_output_dir = tokens[1];
                    }
                    if (!tokens.empty() && tokens[0] == MAP_VAR_OUTPUT_DIR) {
                        if (tokens.size() < 2) {
                            std::cerr << "Line " << line_num << ": Expected value for key " << MAP_VAR_OUTPUT_DIR << std::endl;
                            return false;
                        }
                        prefix_output_dir = tokens[1];
                    }

                } else {
                    if (prefix_column == -1 || first_column == -1 || second_column == 1) {
                        std::cerr << "The columns must be specified: " << MAP_PREFIX << ", " << MAP_FIRST_READ << ", " << MAP_SECOND_READ << std::endl;
                        return false;
                    }

                    if (header) {
                        std::cerr << "Line " << line_num << ": Expected header line starting with #" << std::endl;
                        return false;
                    }

                    splitter.Split(line);
                    auto sample_id = tokens[prefix_column];
                    auto prefix_path = path(prefix_output_dir).append(tokens[prefix_column]);
                    auto first_path = path(input_dir).append(tokens[first_column]);
                    auto second_path = path(input_dir).append(tokens[second_column]);

                    // mandatory
                    samplenames_list.emplace_back(sample_id);
                    prefix_list.emplace_back(prefix_path);
                    first_list.emplace_back(first_path);
                    second_list.emplace_back(second_path);

                    // optional
                    if (sam_column != -1) {
                        auto sam_path = path(sam_output_dir).append(tokens[sam_column]);
                        sam_list.emplace_back(sam_path);
                    }
                    if (profile_column != -1) {
                        auto profile_path = path(profile_output_dir).append(tokens[profile_column]);
                        profile_list.emplace_back(profile_path);
                    }
                    if (profile_truth_column != -1) {
                        auto profile_truth_path = tokens[profile_truth_column];
                        profile_truth_list.emplace_back(profile_truth_path);
                    }
                }

                line_num++;
            }


            // Fill in gaps
            if (profile_list.empty()) {
                if (prefix_list.empty()) {
                    std::cout << "column PREFIX must be specified if column PROFILE is not specified." << std::endl;
                }
                for (auto prefix : prefix_list) {
                    profile_list.emplace_back(prefix + ".profile");
                }
            }
            if (sam_list.empty()) {
                if (prefix_list.empty()) {
                    std::cout << "column PREFIX must be specified if column SAM is not specified." << std::endl;
                }
                for (auto prefix : prefix_list) {
                    sam_list.emplace_back(prefix + ".sam");
                }
            }


//            std::cout << "strain_output_dir " << strain_output_dir << std::endl;
//            std::cout << "sam_output_dir " << sam_output_dir << std::endl;
//            std::cout << "profile_output_dir " << sam_output_dir << std::endl;
//            std::cout << "input_dir " << input_dir << std::endl;
//
//            std::cout << "\nPrefixes" << std::endl;
//            for (auto& prefix : prefix_list) {
//                std::cout << prefix << std::endl;
//            }
//
//            std::cout << "\nFirsts" << std::endl;
//            for (auto& first : first_list) {
//                std::cout << first << std::endl;
//            }
//
//            std::cout << "\nSeconds" << std::endl;
//            for (auto& second : second_list) {
//                std::cout << second << std::endl;
//            }
//
//            std::cout << "\nSams" << std::endl;
//            for (auto& sam : sam_list) {
//                std::cout << sam << std::endl;
//            }
//
//            std::cout << "\nProfiles" << std::endl;
//            for (auto& profile : profile_list) {
//                std::cout << profile << std::endl;
//            }


            is.close();
            return true;
        }

        bool PrepareAndCheckValidity(bool force_read_check=false) {
            std::vector<std::string> error_log;
            std::vector<std::string> warning_log;

            if (!std::filesystem::exists(GetSequenceFile())) {
                error_log.emplace_back("Sequence file does not exist: " + GetSequenceFile());
            }
            if (!std::filesystem::exists(GetSequenceMapFile())) {
                error_log.emplace_back("Sequence map file does not exist: " + GetSequenceMapFile());
            }
            if (!std::filesystem::exists(GetInternalTaxonomyFile())) {
                error_log.emplace_back("Taxonomy file does not exist: " + GetInternalTaxonomyFile());
            }
            if (!m_build && !std::filesystem::exists(GetIndexFile())) {
                error_log.emplace_back("Index file does not exist: " + GetIndexFile());
            }

            // Length of files
            bool valid_lengths1 =
                    m_first_list.size() == m_second_list.size() &&
                    m_first_list.size() == m_prefix_list.size();

            bool valid_lengths2 =
                    m_prefix_list.size() == m_sam_list.size();


            if (!valid_lengths1 && !m_first_list.empty()) {
                std::string error =
                        "You must provide equal amounts of items in options -1, -2 and -3. "
                        "Provided are -1 (" +
                        std::to_string(m_first_list.size()) +
                        "), -2 (" +
                        std::to_string(m_second_list.size()) +
                        "), and -3 (" +
                        std::to_string(m_prefix_list.size()) +
                        ")";
                error_log.emplace_back(error);
            }

            if (!valid_lengths2 && m_profile_only) {
                std::string error =
                        "You must provide equal amounts of items in options -3 and --profile_only";
                error_log.emplace_back(error);
            }

            if (m_first_list.size() != m_second_list.size() || m_first_list.size() != m_sam_list.size()) {
                std::cerr << "First:  " << m_first_list.size() << std::endl;
                std::cerr << "Second: " << m_second_list.size() << std::endl;
                std::cerr << "Sam:    " << m_sam_list.size() << std::endl;
                std::cerr << "lists different sizes" << std::endl;
                exit(9);
            }

            // Check files
            for (auto i = 0; i < m_first_list.size(); i++) {

                auto first = m_first_list[i];
                auto second = m_second_list[i];
                auto sam = m_sam_list[i];

                auto first_exists = std::filesystem::exists(first);
                auto second_exists = std::filesystem::exists(second);
                auto sam_exists = std::filesystem::exists(sam);

                if (!first_exists || !second_exists) {
                    if (!force_read_check && sam_exists) {
                        if (!first_exists) {
                            std::string warning = "-1 file does not exist: " + first + " ( But .sam file does )";
                            warning_log.emplace_back(warning);
                        }
                        if (!second_exists) {
                            std::string warning = "-2 file does not exist: " + second + " ( But .sam file does )";
                            warning_log.emplace_back(warning);
                        }
                    } else {
                        if (!first_exists) {
                            std::string error = "-1 file does not exist: " + first;
                            error_log.emplace_back(error);
                        }
                        if (!second_exists) {
                            std::string error = "-2 file does not exist: " + second;
                            error_log.emplace_back(error);
                        }
                    }
                }
            }

            if (!warning_log.empty()) {
                for (auto& line : error_log) {
                    std::cerr << "Warning: " << line << std::endl;
                }
            }

            if (!error_log.empty()) {
                for (auto& line : error_log) {
                    std::cerr << "Error: " << line << std::endl;
                }
                return false;
            }
            return true;
        }

        static std::vector<size_t> ProcessRange(std::string const& s, size_t const& total) {
            std::regex pattern("^\\d+(-\\d+)?(,\\d+(-(\\d+)?)?)*$");
            if (!std::regex_match(s, pattern)) {
                std::cout << "String " << s << " does not match the regex." << std::endl;
                exit(9);
            }

            std::vector<size_t> samples;
            std::vector<std::string> tokens;
            std::vector<std::string> subtokens;

            Utils::split(tokens, s, ",");

            for (auto& token : tokens) {
                if (token.empty()) continue;
                if (token.find('-') != std::string::npos) {
                    Utils::split(subtokens, token, "-");

                    auto begin = stoll(subtokens[0]);
                    auto end = !subtokens[1].empty() ? stoll(subtokens[1]) + 1 : total + 1;

                    if (begin >= total || end - 1 > total) {
                        std::cout << "There are only " << total << " samples." << std::endl;
                        exit(9);
                    }

                    for (auto i = begin; i < end; i++) {
                        samples.emplace_back(i-1);
                    }
                } else {
                    samples.emplace_back(stoll(token)-1);
                }
            }

            std::unordered_set<size_t> as_set(samples.begin(), samples.end());
            if (as_set.size() != samples.size()) {
                std::cout << "Ranges and numbers may not overlap." << std::endl;
                exit(9);
            }

            return samples;
        }

        static Options OptionsFromArguments(int argc, char *argv[]) {
            auto cxx_options = CxxOptions();


            if (argc <= 1) {
                return {};
            }

            cxx_options.parse_positional({ "reference" });
            auto result = cxx_options.parse(argc, argv);

            bool help = result.count("help");

            if (help) {
                return {};
            }

            bool no_strains = result.count("no_strains");
            bool build = result.count("build");
            bool preload_genomes_off = result.count("preload_genomes_off");
            bool benchmark_alignment = result.count("benchmark_alignment");
            bool fastalign = result.count("fastalign");
            bool profile = result.count("profile");
            bool verbose = result.count("verbose");

            size_t threads = result["threads"].as<size_t>();
            size_t align_top = result["align_top"].as<size_t>();
            size_t x_drop = result["x_drop"].as<size_t>();
            size_t max_key_ubiquity = result["max_key_ubiquity"].as<size_t>();
            size_t max_seed_size = result["max_seed_size"].as<size_t>();
            double max_score_ani = result["max_score_ani"].as<double>();
            double msa_min_vcov = result["msa_min_vcov"].as<double>();
            size_t max_out = result["max_out"].as<size_t>();


            auto reference = result["reference"].as<std::string>();
            auto full_reference = result["full_reference"].as<std::string>();

            auto map_file = result.count("map") ? result["map"].as<std::string>() : "";
            auto first = result.count("first") ? result["first"].as<std::string>() : "";
            auto second = result.count("second") ? result["second"].as<std::string>() : "";
            auto sam_in = result.count("profile_only") ? result["profile_only"].as<std::string>() : "";
            auto universal_prefix = result.count("prefix") ? result["prefix"].as<std::string>() : "";
            auto output_dir = result.count("outdir") ? result["outdir"].as<std::string>() : "";

            auto range_arg = result.count("map_range") ? result["map_range"].as<std::string>() : "";
            std::vector<size_t> range;


            std::vector<std::string> first_list;
            std::vector<std::string> second_list;
            std::vector<std::string> prefix_list;
            std::vector<std::string> sam_list;
            std::vector<std::string> profile_list;
            std::vector<std::string> samplenames_list;
            std::vector<std::string> profile_truth_list;

            if (!map_file.empty()) {
                output_dir = "";
                LoadFromMap(map_file, output_dir, prefix_list, first_list, second_list, sam_list, profile_list, samplenames_list, profile_truth_list);

                if (range_arg != "") {
                    range = ProcessRange(range_arg, first_list.size());
                }
            } else {
                LineSplitter::Split(first, ",", first_list);
                LineSplitter::Split(second, ",", second_list);
                LineSplitter::Split(universal_prefix, ",", prefix_list);
                LineSplitter::Split(sam_in, ",", sam_list);
            }


            auto profile_truth = result.count("profile_truth") ? result["profile_truth"].as<std::string>() : "";
            bool profile_only = result.count("profile_only");
            bool mapq_debug_output = result.count("mapq_debug_output");
            bool force = result.count("force");
            bool gzip_sam = !result.count("gzip_sam");

            if (profile_only) {
                if (!first_list.empty()) {
                    std::cerr << "Warning: profile only selected but first list is not empty. List is cleared" << std::endl;
                }
                if (!second_list.empty()) {
                    std::cerr << "Warning: profile only selected but second list is not empty. List is cleared" << std::endl;
                }
                first_list.clear();
                second_list.clear();

                if (!prefix_list.empty() && prefix_list.size() != sam_list.size()) {
                    std::cerr << "Warning: If profile only is selected, the user can either specify output prefixes for all sam files or for non. Sizes are not equal." << std::endl;
                    exit(34);
                }

                if (prefix_list.empty()) {
                    prefix_list.resize(sam_list.size());
                    for (auto i = 0; i < sam_list.size(); i++) {
                        auto& sam_file = sam_list[i];
                        if (sam_file.ends_with(".sam")) {
                            prefix_list[i] = sam_file.substr(0, sam_file.length() - 4);
                        } else if (sam_file.ends_with(".sam.gz")) {
                            prefix_list[i] = sam_file.substr(0, sam_file.length() - 7);
                        } else {
                            std::cerr << sam_file << " does not end on .sam" << std::endl;
                            exit(35);
                        }
                    }
                }


            } else {
                if (sam_list.empty()) {
                    for (auto i = 0; i < prefix_list.size(); i++) {
                        sam_list.emplace_back(prefix_list[i] + (gzip_sam ? ".sam.gz" : ".sam"));
                    }
                }
            }

            if (profile && profile_list.empty()) {
                for (auto i = 0; i < prefix_list.size(); i++) {

                    profile_list.emplace_back(prefix_list[i] + ".profile");
                    // std::cout << profile_list[i] << std::endl;
                }
            }

            auto benchmark_alignment_output_file = result.count("benchmark_alignment_output") ? result["benchmark_alignment_output"].as<std::string>() : "";

            std::string db_path;
            if (result.count("db")) {
                db_path = result["db"].as<std::string>();
            } else {
                auto db_path_env = std::getenv(PROTAL_DB_ENV_VARIABLE.c_str());
                std::cout << "Get DB from environment variable $" << PROTAL_DB_ENV_VARIABLE << std::endl;
                if (!db_path_env) {
                    std::cerr << "Error " << PROTAL_DB_ENV_VARIABLE << std::endl;
                } else {
                    db_path = db_path_env;
                }
            }

            if (range.empty()) {
                range.resize(prefix_list.size());
                std::iota(range.begin(), range.end(), 0);
            }

            auto options = Options(
                    build,
                    profile,
                    profile_only,
                    no_strains,
                    !preload_genomes_off,
                    benchmark_alignment,
                    benchmark_alignment_output_file,
                    help,
                    mapq_debug_output,
                    first_list,
                    second_list,
                    samplenames_list,
                    db_path,
                    prefix_list,
                    reference,
                    full_reference,
                    map_file,
                    output_dir,
                    threads,
                    align_top,
                    max_out,
                    max_score_ani,
                    msa_min_vcov,
                    x_drop,
                    max_key_ubiquity,
                    max_seed_size,
                    fastalign,
                    sam_list,
                    profile_list,
                    profile_truth_list,
                    profile_truth,
                    force,
                    verbose,
                    range);


            if (!options.PrepareAndCheckValidity()) {
                std::cerr << "Exit Program" << std::endl;
                exit(30);
            }

            return options;
        }
    };


}