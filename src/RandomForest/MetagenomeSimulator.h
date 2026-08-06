// SPDX-License-Identifier: MIT
#pragma once

#include <filesystem>
#include <optional>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "../Utilities/Benchmark.h"
#include "ArtIlluminaWrapper.h"
#include "CommunityProfileDesigner.h"
#include "MetagenomeTypes.h"

namespace protal::sim {

std::vector<GenomeRecord> read_genome_table(const std::filesystem::path& tsv_path);
std::vector<double> parse_strain_probabilities(const std::string& text);
std::vector<std::string> parse_species_list(const std::string& text);
std::unordered_map<std::string, std::size_t> parse_genus_selection(const std::string& text);
std::unordered_map<std::string, std::size_t> parse_taxon_selection(const std::string& text);
void write_combined_manifest(const std::vector<SampleOutput>& samples, const std::filesystem::path& manifest_path);
void write_sample_manifest(const SampleOutput& sample, const std::filesystem::path& manifest_path);
void write_abundance_matrix(const std::vector<SampleOutput>& samples, const std::filesystem::path& matrix_path);

// Reads a manifest written by write_combined_manifest or write_sample_manifest back
// into a design: which genome at which read depth in which sample, plus the ART seed
// and fasta path when the manifest carries them. Sample and row order are preserved.
std::vector<SampleOutput> read_manifest(const std::filesystem::path& manifest_path);

// Fills in fasta paths for manifests written before the fasta_path column existed,
// looking each genome up in a genome table by name. Throws naming the genomes it
// cannot resolve.
void resolve_manifest_fasta_paths(std::vector<SampleOutput>& samples, const std::vector<GenomeRecord>& genomes);
class MetagenomeSimulator {
public:
    MetagenomeSimulator(
        std::vector<GenomeRecord> genomes,
        ArtIlluminaOptions art_options = {},
        std::uint64_t seed = std::random_device{}(),
        std::string pigz_path = "pigz");

    std::vector<SampleOutput> simulate_samples(
        const ProfileDesignOptions& profile_options,
        std::size_t sample_count,
        const std::string& sample_prefix,
        const std::filesystem::path& output_dir,
        bool skip_reads = false,
        bool keep_tmp = false);

    // Re-simulates a design read back from a manifest, bypassing the profile designer
    // entirely. Rows carrying an art_seed reproduce their reads exactly; rows without
    // one keep the recorded composition and get fresh read realizations.
    std::vector<SampleOutput> replay_samples(
        std::vector<SampleOutput> design,
        const std::filesystem::path& output_dir,
        bool skip_reads = false,
        bool keep_tmp = false);

private:
    std::vector<GenomeRecord> genomes_;
    ArtIlluminaWrapper art_;
    CommunityProfileDesigner designer_;
    std::mt19937_64 rng_;
    std::string pigz_path_;

    SampleOutput simulate_single(
        const ProfileDesignOptions& profile_options,
        const std::string& sample_name,
        const std::filesystem::path& output_dir,
        const std::unordered_map<std::string, std::uint64_t>& genome_lengths,
        std::uint64_t paired_read_length,
        bool skip_reads,
        bool keep_tmp);

    void render_sample(
        SampleOutput& sample,
        const std::filesystem::path& output_dir,
        std::uint64_t paired_read_length,
        bool skip_reads,
        bool keep_tmp);
};

}  // namespace protal::sim
