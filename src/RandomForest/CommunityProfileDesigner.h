// SPDX-License-Identifier: MIT
#pragma once

#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include "MetagenomeTypes.h"

namespace protal::sim {

std::string extract_species(const std::string& taxonomy);

class CommunityProfileDesigner {
public:
    explicit CommunityProfileDesigner(std::vector<GenomeRecord> genomes);

    std::vector<GenomeAssignment> design_profile(
        const ProfileDesignOptions& options,
        std::mt19937_64& rng) const;

private:
    std::vector<GenomeRecord> genomes_;

    std::unordered_map<std::string, std::vector<GenomeRecord>> group_by_species() const;
    std::vector<double> draw_weights(std::size_t count, const ProfileDesignOptions& options, std::mt19937_64& rng) const;
};

}  // namespace protal::sim
