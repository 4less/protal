// SPDX-License-Identifier: MIT
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include <cxxopts.hpp>
#include "RandomForest/MetagenomeSimulator.h"

namespace fs = std::filesystem;
using protal::sim::AbundanceDistribution;
using protal::sim::ArtIlluminaOptions;
using protal::sim::MetagenomeSimulator;
using protal::sim::ProfileDesignOptions;

static std::optional<fs::path> find_executable_path(const std::string& arg0) {
    std::error_code ec;
    fs::path direct = arg0;
    if (!direct.empty() && fs::exists(direct, ec)) {
        fs::path canonical = fs::canonical(direct, ec);
        if (ec) {
            canonical = fs::absolute(direct);
        }
        return canonical;
    }
    const char* path_env = std::getenv("PATH");
    if (!path_env) {
        return std::nullopt;
    }
    std::string path_list(path_env);
    std::size_t start = 0;
    while (start <= path_list.size()) {
        auto end = path_list.find(':', start);
        std::string dir = path_list.substr(start, end == std::string::npos ? std::string::npos : end - start);
        if (!dir.empty()) {
            fs::path candidate = fs::path(dir) / arg0;
            if (fs::exists(candidate, ec)) {
                fs::path canonical = fs::canonical(candidate, ec);
                if (ec) {
                    canonical = fs::absolute(candidate);
                }
                return canonical;
            }
        }
        if (end == std::string::npos) {
            break;
        }
        start = end + 1;
    }
    return std::nullopt;
}

static std::string derive_prefix_from_r1(const fs::path& r1) {
    std::string base = r1.filename().string();
    auto strip_suffix = [&](const std::string& suf) {
        if (base.size() >= suf.size() && base.compare(base.size() - suf.size(), suf.size(), suf) == 0) {
            base.erase(base.size() - suf.size());
        }
    };
    strip_suffix(".gz");
    strip_suffix(".fastq");
    strip_suffix(".fq");
    strip_suffix(".fasta");
    strip_suffix(".fa");
    strip_suffix(".FASTQ");
    strip_suffix(".FQ");
    if (base.size() > 3 && base.compare(base.size() - 3, 3, "_R1") == 0) {
        base.erase(base.size() - 3);
    } else if (base.size() > 2 && base.compare(base.size() - 2, 2, "_1") == 0) {
        base.erase(base.size() - 2);
    }
    return base;
}

static void write_protal_metafile(
    const std::vector<protal::sim::SampleOutput>& samples,
    const fs::path& output_dir,
    const fs::path& reads_dir,
    const fs::path& metafile_path,
    const std::vector<fs::path>& profile_truth_paths,
    const std::optional<fs::path>& output_dir_override) {
    if (samples.empty()) {
        return;
    }
    if (!metafile_path.parent_path().empty()) {
        fs::create_directories(metafile_path.parent_path());
    }
    auto to_abs = [](const fs::path& p) -> fs::path {
        std::error_code ec;
        auto c = fs::canonical(p, ec);
        return ec ? fs::absolute(p) : c;
    };
    fs::path out_abs = to_abs(output_dir_override ? *output_dir_override : output_dir);
    fs::path reads_abs = to_abs(reads_dir);

    std::ofstream out(metafile_path);
    if (!out) {
        throw std::runtime_error("Unable to write protal metafile: " + metafile_path.string());
    }
    out << "#OUTPUT_DIR\t" << out_abs.string() << "\n";
    out << "#INPUT_DIR\t" << reads_abs.string() << "\n";
    out << "#SAMPLEID\tFIRST\tSECOND\tSAM\tPREFIX\tPROFILE\tPROFILE_TRUTH\n";
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const auto& sample = samples[i];
        auto first = sample.read1_path.filename().string();
        auto second = sample.read2_path.filename().string();
        auto prefix = derive_prefix_from_r1(sample.read1_path);
        fs::path truth_abs = to_abs(profile_truth_paths[i]);
        out << sample.sample_name << '\t' << first << '\t' << second << '\t' << prefix << ".sam.gz"
            << '\t' << prefix << '\t' << prefix << ".profile" << '\t' << truth_abs.string() << '\n';
    }
}

struct CliOptions {
    fs::path genome_table;
    fs::path output_dir;
    std::size_t samples{1};
    std::string sample_prefix{"sample"};
    std::uint64_t total_read_pairs{100'000};
    std::size_t species_per_sample{10};
    std::size_t species_per_sample_min{0};
    AbundanceDistribution distribution{AbundanceDistribution::PoissonLognormal};
    double alpha{2.0};
    int nb_r{5};
    double nb_p{0.5};
    double pln_mu{0.0};
    double pln_sigma{1.3};
    std::string strain_probabilities;
    std::string include_species;
    std::string genus_counts;
    std::string taxon_counts;
    ArtIlluminaOptions art;
    std::optional<std::uint64_t> seed;
    bool plot_png{false};
    bool test_mode{false};
    bool keep_tmp{false};
    bool pick_random_demand_if_fail{false};
    int threads{1};
    std::string pigz_path{"pigz"};
    std::optional<fs::path> protal_metafile_output_dir;
    fs::path strain_sharing_file;
};

static std::vector<protal::sim::StrainSharingSpec> parse_strain_sharing_file(const fs::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open strain sharing file: " + path.string());

    std::vector<protal::sim::StrainSharingSpec> specs;
    std::string line;
    std::size_t line_no = 0;
    while (std::getline(in, line)) {
        ++line_no;
        if (line.empty() || line[0] == '#') continue;

        // tab-split
        std::vector<std::string> fields;
        std::size_t start = 0;
        while (true) {
            auto pos = line.find('\t', start);
            std::string tok = (pos == std::string::npos) ? line.substr(start) : line.substr(start, pos - start);
            // trim whitespace
            auto nb = [](unsigned char c) { return !std::isspace(c); };
            tok.erase(tok.begin(), std::find_if(tok.begin(), tok.end(), nb));
            tok.erase(std::find_if(tok.rbegin(), tok.rend(), nb).base(), tok.end());
            fields.push_back(tok);
            if (pos == std::string::npos) break;
            start = pos + 1;
        }

        if (fields.size() < 4) {
            throw std::runtime_error("strain_sharing_file line " + std::to_string(line_no) +
                ": expected at least 4 tab-separated columns (SPECIES, SAMPLE_FRACTION, N_STRAINS, MIN_OCCURRENCE)");
        }

        // Catch common mistake: decimal point in an integer column usually means a tab was
        // replaced by a space, shifting all subsequent columns by one.
        auto require_integer_field = [&](const std::string& val, const char* col_name) -> std::size_t {
            if (val.find('.') != std::string::npos) {
                throw std::runtime_error(
                    "strain_sharing_file line " + std::to_string(line_no) + ": " + col_name +
                    " must be an integer, got \"" + val + "\" (decimal point found — "
                    "check for spaces instead of tabs or a missing column)");
            }
            return static_cast<std::size_t>(std::stoull(val));
        };

        protal::sim::StrainSharingSpec spec;
        spec.species         = fields[0];
        spec.sample_fraction = std::stod(fields[1]);
        spec.n_strains       = require_integer_field(fields[2], "N_STRAINS");
        spec.min_occurrence  = require_integer_field(fields[3], "MIN_OCCURRENCE");
        spec.min_vcov        = (fields.size() >= 5 && !fields[4].empty()) ? std::stod(fields[4]) : 0.0;
        if (fields.size() >= 6 && !fields[5].empty()) {
            spec.conspecific_strain_probabilities = protal::sim::parse_strain_probabilities(fields[5]);
        }

        if (spec.sample_fraction < 0.0 || spec.sample_fraction > 1.0) {
            throw std::runtime_error("strain_sharing_file line " + std::to_string(line_no) +
                ": SAMPLE_FRACTION must be in [0,1]");
        }
        specs.push_back(std::move(spec));
    }
    return specs;
}

static cxxopts::Options build_cxxopts() {
    cxxopts::Options options("simulate_metagenomes", "Simulate shotgun metagenomes from a genome table");

    options.add_options("I/O")
        ("genome_table", "TSV with genome name, GTDB taxonomy, FASTA path (.gz ok)", cxxopts::value<std::string>())
        ("o,output_dir",  "Output directory for FASTQs and manifest", cxxopts::value<std::string>());

    options.add_options("Sampling")
        ("n,samples",           "Number of metagenome samples", cxxopts::value<std::size_t>()->default_value("1"))
        ("sample_prefix",       "Prefix for sample names", cxxopts::value<std::string>()->default_value("sample"))
        ("total_read_pairs",    "Read pairs per sample", cxxopts::value<std::uint64_t>()->default_value("100000"))
        ("species_per_sample",  "Number of species per sample, or an inclusive range e.g. 20-80", cxxopts::value<std::string>()->default_value("10"))
        ("distribution",        "Abundance model: power_law | negative_binomial | poisson_lognormal", cxxopts::value<std::string>()->default_value("poisson_lognormal"))
        ("alpha",               "Power law alpha", cxxopts::value<double>()->default_value("2.0"))
        ("nb_r",                "Negative binomial r", cxxopts::value<int>()->default_value("5"))
        ("nb_p",                "Negative binomial p", cxxopts::value<double>()->default_value("0.5"))
        ("pln_mu",              "Poisson-lognormal mean (log-scale)", cxxopts::value<double>()->default_value("0.0"))
        ("pln_sigma",           "Poisson-lognormal sigma (log-scale)", cxxopts::value<double>()->default_value("1.3"))
        ("strains_per_species", "Probabilities for adding 2nd, 3rd, ... strains per species, e.g. \"0.4,0.2,0.1\"", cxxopts::value<std::string>()->default_value(""))
        ("include_species",     "Comma-separated species to force-include in each sample", cxxopts::value<std::string>()->default_value(""))
        ("genus",               "Comma-separated genus:count pairs, e.g. \"g__A:10,g__B:2\"", cxxopts::value<std::string>()->default_value(""))
        ("taxon",               "Comma-separated taxon:count pairs, e.g. \"d__Archaea:10\"", cxxopts::value<std::string>()->default_value(""))
        ("pick_random_demand_if_fail", "If --genus/--taxon demand more species than available, cap and fill randomly instead of failing")
        ("strain_sharing_file", "TSV file for cross-sample strain sharing. Columns (tab-separated): "
                                "SPECIES  SAMPLE_FRACTION  N_STRAINS  MIN_OCCURRENCE  MIN_VCOV  CONSPECIFIC_STRAINS. "
                                "Lines starting with # are ignored. MIN_VCOV and CONSPECIFIC_STRAINS are optional. "
                                "CONSPECIFIC_STRAINS uses the same probability format as --strains_per_species.",
                                cxxopts::value<std::string>()->default_value(""))
        ("seed",                "RNG seed (default: random)", cxxopts::value<std::uint64_t>());

    options.add_options("ART")
        ("art_path",        "art_illumina executable", cxxopts::value<std::string>()->default_value("art_illumina"))
        ("read_length",     "Read length", cxxopts::value<int>()->default_value("150"))
        ("fragment_mean",   "Fragment mean", cxxopts::value<int>()->default_value("350"))
        ("fragment_stdev",  "Fragment stdev", cxxopts::value<int>()->default_value("50"))
        ("sequencer",       "ART sequencer profile", cxxopts::value<std::string>()->default_value("HS25"))
        ("extra_art_args",  "Extra ART arguments, e.g. \"--qprof1 q1 --qprof2 q2\"", cxxopts::value<std::string>()->default_value(""));

    options.add_options("General")
        ("t,threads",       "Threads for ART/pigz", cxxopts::value<int>()->default_value("1"))
        ("pigz_path",       "Path to pigz", cxxopts::value<std::string>()->default_value("pigz"))
        ("protal_metafile", "Write a Protal meta file (output_dir/protal.meta) but set OUTPUT_DIR to <path>", cxxopts::value<std::string>())
        ("test",            "Generate profiles/manifests but skip read simulation (fast dry run)")
        ("keep_tmp",        "Keep the individual per-genome reads")
        ("plot_png",        "Generate barplot PNG of species abundances")
        ("h,help",          "Print help.");

    return options;
}

static CliOptions parse_cli(int argc, char** argv) {
    auto cxx = build_cxxopts();

    if (argc <= 1) {
        std::cout << cxx.help({"I/O", "Sampling", "ART", "General"}) << std::endl;
        std::exit(0);
    }

    auto result = cxx.parse(argc, argv);

    if (result.count("help")) {
        std::cout << cxx.help({"I/O", "Sampling", "ART", "General"}) << std::endl;
        std::exit(0);
    }

    if (!result.count("genome_table") || !result.count("output_dir")) {
        std::cerr << "Error: --genome_table and --output_dir are required.\n\n";
        std::cout << cxx.help({"I/O", "Sampling", "ART", "General"}) << std::endl;
        std::exit(1);
    }

    CliOptions opts;
    opts.genome_table      = result["genome_table"].as<std::string>();
    opts.output_dir        = result["output_dir"].as<std::string>();
    opts.samples           = result["samples"].as<std::size_t>();
    opts.sample_prefix     = result["sample_prefix"].as<std::string>();
    opts.total_read_pairs  = result["total_read_pairs"].as<std::uint64_t>();
    {
        const std::string sps_arg = result["species_per_sample"].as<std::string>();
        const auto dash = sps_arg.find('-');
        if (dash != std::string::npos) {
            opts.species_per_sample_min = std::stoull(sps_arg.substr(0, dash));
            opts.species_per_sample     = std::stoull(sps_arg.substr(dash + 1));
            if (opts.species_per_sample_min > opts.species_per_sample) {
                throw std::runtime_error("--species_per_sample range min (" +
                    std::to_string(opts.species_per_sample_min) + ") exceeds max (" +
                    std::to_string(opts.species_per_sample) + ")");
            }
        } else {
            opts.species_per_sample = std::stoull(sps_arg);
        }
    }
    opts.alpha             = result["alpha"].as<double>();
    opts.nb_r              = result["nb_r"].as<int>();
    opts.nb_p              = result["nb_p"].as<double>();
    opts.pln_mu            = result["pln_mu"].as<double>();
    opts.pln_sigma         = result["pln_sigma"].as<double>();
    opts.strain_probabilities = result["strains_per_species"].as<std::string>();
    opts.include_species   = result["include_species"].as<std::string>();
    opts.genus_counts      = result["genus"].as<std::string>();
    opts.taxon_counts      = result["taxon"].as<std::string>();
    opts.threads           = result["threads"].as<int>();
    opts.pigz_path         = result["pigz_path"].as<std::string>();
    opts.test_mode         = result.count("test") > 0;
    opts.keep_tmp          = result.count("keep_tmp") > 0;
    opts.plot_png          = result.count("plot_png") > 0;
    opts.pick_random_demand_if_fail = result.count("pick_random_demand_if_fail") > 0;

    if (result.count("seed")) {
        opts.seed = result["seed"].as<std::uint64_t>();
    }
    if (result.count("protal_metafile")) {
        opts.protal_metafile_output_dir = result["protal_metafile"].as<std::string>();
    }
    {
        const std::string ssf = result["strain_sharing_file"].as<std::string>();
        if (!ssf.empty()) opts.strain_sharing_file = ssf;
    }

    const std::string dist_str = result["distribution"].as<std::string>();
    if (dist_str == "power_law") {
        opts.distribution = AbundanceDistribution::PowerLaw;
    } else if (dist_str == "negative_binomial") {
        opts.distribution = AbundanceDistribution::NegativeBinomial;
    } else if (dist_str == "poisson_lognormal") {
        opts.distribution = AbundanceDistribution::PoissonLognormal;
    } else {
        throw std::runtime_error("Unknown distribution: " + dist_str);
    }

    opts.art.art_path      = result["art_path"].as<std::string>();
    opts.art.read_length   = result["read_length"].as<int>();
    opts.art.fragment_mean = result["fragment_mean"].as<int>();
    opts.art.fragment_stdev = result["fragment_stdev"].as<int>();
    opts.art.sequencer     = result["sequencer"].as<std::string>();

    const std::string extra = result["extra_art_args"].as<std::string>();
    if (!extra.empty()) {
        std::istringstream iss(extra);
        std::string token;
        while (iss >> token) opts.art.extra_args.push_back(token);
    }

    return opts;
}

int main(int argc, char** argv) {
    CliOptions cli;
    try {
        cli = parse_cli(argc, argv);
    } catch (const std::exception& ex) {
        std::cerr << "Error parsing arguments: " << ex.what() << '\n';
        return 1;
    }

    try {
        auto genomes = protal::sim::read_genome_table(cli.genome_table);
        if (cli.seed) {
            std::cerr << "Using fixed seed: " << *cli.seed << '\n';
        }

        ProfileDesignOptions profile{};
        profile.species_per_sample     = cli.species_per_sample;
        profile.species_per_sample_min = cli.species_per_sample_min;
        profile.distribution = cli.distribution;
        profile.powerlaw_alpha = cli.alpha;
        profile.negative_binomial_r = cli.nb_r;
        profile.negative_binomial_p = cli.nb_p;
        profile.pln_mu = cli.pln_mu;
        profile.pln_sigma = cli.pln_sigma;
        profile.strain_probabilities = protal::sim::parse_strain_probabilities(cli.strain_probabilities);
        profile.include_species = protal::sim::parse_species_list(cli.include_species);
        profile.genus_species_counts = protal::sim::parse_genus_selection(cli.genus_counts);
        profile.taxon_species_counts = protal::sim::parse_taxon_selection(cli.taxon_counts);
        profile.total_read_pairs = cli.total_read_pairs;
        profile.pick_random_demand_if_fail = cli.pick_random_demand_if_fail;
        if (!cli.strain_sharing_file.empty()) {
            profile.strain_sharing = parse_strain_sharing_file(cli.strain_sharing_file);
            std::cerr << "Loaded " << profile.strain_sharing.size()
                      << " strain sharing spec(s) from " << cli.strain_sharing_file << '\n';
        }

        std::uint64_t seed = cli.seed ? *cli.seed : std::random_device{}();
        cli.art.threads = std::max(1, cli.threads);
        MetagenomeSimulator simulator(std::move(genomes), cli.art, seed, cli.pigz_path);

        auto samples =
            simulator.simulate_samples(profile, cli.samples, cli.sample_prefix, cli.output_dir, cli.test_mode, cli.keep_tmp);

        auto combined_manifest_path = cli.output_dir / "manifest.tsv";
        protal::sim::write_combined_manifest(samples, combined_manifest_path);

        fs::path manifests_dir = cli.output_dir / "manifests";
        fs::create_directories(manifests_dir);
        for (const auto& sample : samples) {
            protal::sim::write_sample_manifest(sample, manifests_dir / (sample.sample_name + ".tsv"));
        }

        protal::sim::write_abundance_matrix(samples, cli.output_dir / "abundance_matrix.tsv");

        std::cout << "Wrote " << samples.size() << " samples to " << cli.output_dir << '\n';

        if (cli.plot_png) {
            fs::path plots_dir = cli.output_dir / "plots";
            fs::create_directories(plots_dir);
            fs::path script_path = fs::path("scripts/plot_abundances.R");
            if (!fs::exists(script_path)) {
                // Try locating next to the executable.
                if (auto exe_path = find_executable_path(argv[0])) {
                    script_path = exe_path->parent_path() / "plot_abundances.R";
                }
            }
            if (!fs::exists(script_path)) {
                throw std::runtime_error("plot_abundances.R not found; please run from repo root");
            }
            for (const auto& sample : samples) {
                fs::path manifest_path = manifests_dir / (sample.sample_name + ".tsv");
                fs::path plot_out = plots_dir / (sample.sample_name + ".png");
                std::stringstream cmd;
                cmd << "Rscript \"" << script_path.string() << "\" \"" << manifest_path.string() << "\" \""
                    << plot_out.string() << "\"";
                int rc = std::system(cmd.str().c_str());
                if (rc != 0) {
                    throw std::runtime_error("Plot generation failed for " + sample.sample_name + " with exit code " +
                                             std::to_string(rc));
                }
                std::cout << "Generated plot: " << plot_out << '\n';
            }
        }

        if (cli.protal_metafile_output_dir) {
            fs::path reads_dir = cli.output_dir / "reads";
            fs::path goldstd_dir = cli.output_dir / "protal_goldstd";
            fs::create_directories(goldstd_dir);

            std::vector<fs::path> truth_paths;
            truth_paths.reserve(samples.size());
            for (const auto& sample : samples) {
                fs::path truth_path = goldstd_dir / (sample.sample_name + ".profile_truth");
                truth_paths.push_back(truth_path);
                std::unordered_set<std::string> seen;
                std::ofstream truth_out(truth_path);
                if (!truth_out) {
                    throw std::runtime_error("Unable to write profile truth: " + truth_path.string());
                }
                for (const auto& assignment : sample.assignments) {
                    if (seen.insert(assignment.genome.taxonomy).second) {
                        truth_out << assignment.genome.taxonomy << '\n';
                    }
                }
            }

            fs::path metafile_path = cli.output_dir / "protal.meta";
            write_protal_metafile(
                samples, cli.output_dir, reads_dir, metafile_path, truth_paths, cli.protal_metafile_output_dir);
            std::cout << "Wrote Protal metafile: " << metafile_path << '\n';
        }
    } catch (const std::exception& ex) {
        std::cerr << "Simulation failed: " << ex.what() << '\n';
        return 1;
    }
    return 0;
}
