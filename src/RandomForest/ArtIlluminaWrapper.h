// SPDX-License-Identifier: MIT
#pragma once

#include <filesystem>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "MetagenomeTypes.h"

namespace protal::sim {

class ArtIlluminaWrapper {
public:
    explicit ArtIlluminaWrapper(ArtIlluminaOptions options = {});

    std::pair<std::filesystem::path, std::filesystem::path> simulate_read_pairs(
        const GenomeRecord& genome,
        std::uint64_t read_pairs,
        std::uint64_t genome_length,
        const std::filesystem::path& output_prefix,
        unsigned int art_seed,
        const std::filesystem::path& temp_dir) const;

    const ArtIlluminaOptions& options() const { return options_; }

    // The seed from an -rs/--rndSeed in extra_args, which overrides the one this
    // wrapper is handed. Callers record it so the manifest names the seed ART really
    // used, not the one that was ignored.
    std::optional<std::uint64_t> seed_override() const;

private:
    ArtIlluminaOptions options_;

    std::filesystem::path ensure_fasta(const std::filesystem::path& fasta, const std::filesystem::path& temp_dir) const;
    void decompress_gzip(const std::filesystem::path& gz_path, const std::filesystem::path& output_path) const;
    void run_command(const std::vector<std::string>& args, const std::vector<std::pair<std::string, std::string>>& env) const;
};

}  // namespace protal::sim
