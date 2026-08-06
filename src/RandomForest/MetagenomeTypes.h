// SPDX-License-Identifier: MIT
#pragma once

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace protal::sim {

struct GenomeRecord {
    std::string name;
    std::string taxonomy;
    std::filesystem::path fasta_path;
    std::optional<std::uint64_t> genome_length;
};

struct GenomeAssignment {
    GenomeRecord genome;
    std::string species;
    std::uint64_t read_pairs{};
    double relative_abundance{0.0};
    double vertical_coverage{0.0};
    std::uint64_t genome_length{0};
    // ART's --rndSeed for this genome. Recorded in the manifest so that a replay
    // (--from_manifest) reproduces the reads themselves, not just the composition.
    std::optional<std::uint64_t> art_seed;
};

enum class AbundanceDistribution {
    PowerLaw,
    NegativeBinomial,
    PoissonLognormal
};

// Cross-sample strain sharing config for a single species.
struct StrainSharingSpec {
    std::string species;            // species name (resolved against genome table)
    double sample_fraction{1.0};    // fraction of all samples that include this species [0,1]
    std::size_t n_strains{1};       // total distinct strains drawn for this species across all samples
    std::size_t min_occurrence{1};  // each drawn strain must appear in >= this many samples
    double min_vcov{0.0};           // minimum vertical coverage per strain per sample (0 = no floor)
    // Probabilities of adding a 2nd, 3rd, ... conspecific strain to a sample from the forced pool.
    // Same semantics as --strains_per_species but scoped to the n_strains drawn for this species.
    std::vector<double> conspecific_strain_probabilities;
};

struct ProfileDesignOptions {
    std::uint64_t total_read_pairs{100'000};
    std::size_t species_per_sample{10};
    std::size_t species_per_sample_min{0};  // if > 0, count is drawn uniformly from [min, species_per_sample] per sample
    AbundanceDistribution distribution{AbundanceDistribution::PoissonLognormal};
    double powerlaw_alpha{2.0};
    int negative_binomial_r{5};
    double negative_binomial_p{0.5};
    double pln_mu{0.0};
    double pln_sigma{1.3};
    std::vector<std::string> include_species;  // species that must be present in each sample
    std::unordered_map<std::string, std::size_t> genus_species_counts;  // requested species counts per genus
    std::vector<double> strain_probabilities;  // probabilities for adding 2nd, 3rd, ... strain of a species
    std::unordered_map<std::string, std::size_t> taxon_species_counts;  // requested species counts per taxon token
    bool pick_random_demand_if_fail{false};  // if true, cap genus/taxon demands to species_per_sample instead of failing

    // Cross-sample strain sharing: consumed by MetagenomeSimulator before the per-sample loop.
    std::vector<StrainSharingSpec> strain_sharing;

    // Per-sample fields populated by MetagenomeSimulator from strain_sharing pre-assignment:
    // Maps species key -> specific GenomeRecord list to use (bypasses normal strain picking).
    std::unordered_map<std::string, std::vector<GenomeRecord>> forced_strains;
    // Maps species key -> minimum relative abundance floor derived from min_vcov.
    std::unordered_map<std::string, double> species_min_abundance;
};

struct ArtIlluminaOptions {
    std::string art_path{"art_illumina"};
    int read_length{150};
    int fragment_mean{350};
    int fragment_stdev{50};
    std::string sequencer{"HS25"};
    std::vector<std::string> extra_args;
    int threads{1};
};

struct SampleOutput {
    std::string sample_name;
    std::filesystem::path read1_path;
    std::filesystem::path read2_path;
    std::vector<GenomeAssignment> assignments;
};

}  // namespace protal::sim
