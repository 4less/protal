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
    std::uint64_t genome_length{0};
};

enum class AbundanceDistribution {
    PowerLaw,
    NegativeBinomial
};

struct ProfileDesignOptions {
    std::uint64_t total_read_pairs{100'000};
    std::size_t species_per_sample{10};
    AbundanceDistribution distribution{AbundanceDistribution::PowerLaw};
    double powerlaw_alpha{2.0};
    int negative_binomial_r{5};
    double negative_binomial_p{0.5};
    std::vector<std::string> include_species;  // species that must be present in each sample
    std::unordered_map<std::string, std::size_t> genus_species_counts;  // requested species counts per genus
    std::vector<double> strain_probabilities;  // probabilities for adding 2nd, 3rd, ... strain of a species
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
