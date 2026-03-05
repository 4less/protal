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
    std::size_t genomes_per_sample{0};  // deprecated fallback
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
    int threads{1};
    std::string pigz_path{"pigz"};
    std::optional<fs::path> protal_metafile_output_dir;
};

void print_usage() {
    std::cout << "simulate_metagenomes --genome-table <file.tsv> --output-dir <dir> [options]\n"
              << "Required:\n"
              << "  --genome-table <file.tsv>       TSV with genome name, GTDB taxonomy, FASTA path (.gz ok)\n"
              << "  --output-dir <dir>              Output directory for FASTQs and manifest\n"
              << "Options:\n"
              << "  --samples <int>                 Number of metagenome samples (default: 1)\n"
              << "  --sample-prefix <str>           Prefix for sample names (default: sample)\n"
              << "  --total-read-pairs <int>        Read pairs per sample (default: 100000)\n"
              << "  --species-per-sample <int>      Number of species per sample (default: 10)\n"
              << "  --distribution <power_law|negative_binomial|poisson_lognormal>  Abundance model (default: poisson_lognormal)\n"
              << "  --alpha <float>                 Power law alpha (default: 2.0)\n"
              << "  --nb-r <int>                    Negative binomial r (default: 5)\n"
              << "  --nb-p <float>                  Negative binomial p (default: 0.5)\n"
              << "  --pln-mu <float>                Poisson-lognormal mean (log-scale) (default: 0.0)\n"
              << "  --pln-sigma <float>             Poisson-lognormal sigma (log-scale) (default: 1.3)\n"
              << "  --strains-per-species \"0.4,0.2,0.1\"  Probabilities for adding 2nd, 3rd, ... strains per species\n"
              << "  --include-species \"SpeciesA,SpeciesB\" Comma-separated list of species to force-include in each sample\n"
              << "  --genus \"g__A:10,g__B:2\"       Comma-separated genus:count pairs; randomly pick <count> species per genus\n"
              << "  --taxon \"d__Archaea:10\"        Comma-separated taxon:count pairs; randomly pick <count> species per taxon\n"
              << "  --test                          Generate profiles/manifests but skip read simulation (fast dry run)\n"
              << "  --keep-tmp                      Keep the individual reads\n"
              << "  --art-path <path>               art_illumina executable (default: art_illumina)\n"
              << "  --read-length <int>             Read length (default: 150)\n"
              << "  --fragment-mean <int>           Fragment mean (default: 350)\n"
              << "  --fragment-stdev <int>          Fragment stdev (default: 50)\n"
              << "  --sequencer <id>                ART sequencer profile (default: HS25)\n"
              << "  --extra-art-args \"--qprof1 q1 --qprof2 q2\"    Extra ART arguments\n"
              << "  --seed <int>                    RNG seed (default: random)\n"
              << "  --threads <int>                 Threads for ART/pigz (default: 1)\n"
              << "  --pigz-path <path>              Path to pigz (default: pigz)\n"
              << "  --protal_metafile <path>        Write a Protal meta file (output_dir/protal.meta) but set OUTPUT_DIR to <path>\n"
              << "  --plot-png                      Generate barplot PNG of species abundances\n"
              << "  --help                          Show this message\n";
}

bool parse_cli(int argc, char** argv, CliOptions& opts, std::string& err) {
    auto require_value = [&](int& i) -> std::string {
        if (i + 1 >= argc) {
            throw std::runtime_error("Missing value for argument " + std::string(argv[i]));
        }
        return std::string(argv[++i]);
    };

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage();
            return false;
        } else if (arg == "--genome-table") {
            opts.genome_table = require_value(i);
        } else if (arg == "--output-dir") {
            opts.output_dir = require_value(i);
        } else if (arg == "--samples") {
            opts.samples = static_cast<std::size_t>(std::stoull(require_value(i)));
        } else if (arg == "--sample-prefix") {
            opts.sample_prefix = require_value(i);
        } else if (arg == "--total-read-pairs") {
            opts.total_read_pairs = std::stoull(require_value(i));
        } else if (arg == "--species-per-sample") {
            opts.species_per_sample = static_cast<std::size_t>(std::stoull(require_value(i)));
        } else if (arg == "--distribution") {
            std::string val = require_value(i);
            if (val == "power_law") {
                opts.distribution = AbundanceDistribution::PowerLaw;
            } else if (val == "negative_binomial") {
                opts.distribution = AbundanceDistribution::NegativeBinomial;
            } else if (val == "poisson_lognormal") {
                opts.distribution = AbundanceDistribution::PoissonLognormal;
            } else {
                throw std::runtime_error("Unknown distribution: " + val);
            }
        } else if (arg == "--alpha") {
            opts.alpha = std::stod(require_value(i));
        } else if (arg == "--nb-r") {
            opts.nb_r = std::stoi(require_value(i));
        } else if (arg == "--nb-p") {
            opts.nb_p = std::stod(require_value(i));
        } else if (arg == "--pln-mu") {
            opts.pln_mu = std::stod(require_value(i));
        } else if (arg == "--pln-sigma") {
            opts.pln_sigma = std::stod(require_value(i));
        } else if (arg == "--strains-per-species") {
            opts.strain_probabilities = require_value(i);
        } else if (arg == "--include-species") {
            opts.include_species = require_value(i);
        } else if (arg == "--genus") {
            opts.genus_counts = require_value(i);
        } else if (arg == "--taxon") {
            opts.taxon_counts = require_value(i);
        } else if (arg == "--art-path") {
            opts.art.art_path = require_value(i);
        } else if (arg == "--read-length") {
            opts.art.read_length = std::stoi(require_value(i));
        } else if (arg == "--fragment-mean") {
            opts.art.fragment_mean = std::stoi(require_value(i));
        } else if (arg == "--fragment-stdev") {
            opts.art.fragment_stdev = std::stoi(require_value(i));
        } else if (arg == "--sequencer") {
            opts.art.sequencer = require_value(i);
        } else if (arg == "--extra-art-args") {
            std::istringstream iss(require_value(i));
            std::string token;
            while (iss >> token) {
                opts.art.extra_args.push_back(token);
            }
        } else if (arg == "--seed") {
            opts.seed = std::stoull(require_value(i));
        } else if (arg == "--threads") {
            opts.threads = std::stoi(require_value(i));
        } else if (arg == "--pigz-path") {
            opts.pigz_path = require_value(i);
        } else if (arg == "--protal_metafile") {
            opts.protal_metafile_output_dir = require_value(i);
        } else if (arg == "--plot-png") {
            opts.plot_png = true;
        } else if (arg == "--test") {
            opts.test_mode = true;
        } else if (arg == "--keep-tmp") {
            opts.keep_tmp = true;
        } else {
            err = "Unknown argument: " + arg;
            return false;
        }
    }

    if (opts.genome_table.empty() || opts.output_dir.empty()) {
        err = "Required: --genome-table and --output-dir";
        return false;
    }
    return true;
}

int main(int argc, char** argv) {
    CliOptions cli;
    std::string error;
    try {
        bool parsed = parse_cli(argc, argv, cli, error);
        if (!parsed) {
            if (!error.empty()) {
                std::cerr << error << "\n\n";
            }
            return error.empty() ? 0 : 1;
        }
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
        // Species count: if user didn't provide species-per-sample, fall back to previous flag value.
        if (cli.species_per_sample == 0 && cli.genomes_per_sample != 0) {
            profile.species_per_sample = cli.genomes_per_sample;
        } else {
            profile.species_per_sample = cli.species_per_sample;
        }
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
