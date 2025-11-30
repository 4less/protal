// SPDX-License-Identifier: MIT
#include "CommunityProfileDesigner.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_map>

namespace protal::sim {

std::string extract_species(const std::string& taxonomy) {
    std::string current;
    std::string best = "unknown_species";
    for (char c : taxonomy) {
        if (c == ';') {
            if (current.rfind("s__", 0) == 0 && current.size() > 3) {
                return current.substr(3);
            }
            if (!current.empty()) {
                best = current.size() > 3 ? current.substr(3) : current;
            }
            current.clear();
        } else {
            current.push_back(c);
        }
    }
    if (!current.empty()) {
        if (current.rfind("s__", 0) == 0 && current.size() > 3) {
            return current.substr(3);
        }
        best = current.size() > 3 ? current.substr(3) : current;
    }
    return best.empty() ? "unknown_species" : best;
}

CommunityProfileDesigner::CommunityProfileDesigner(std::vector<GenomeRecord> genomes)
    : genomes_(std::move(genomes)) {
    if (genomes_.empty()) {
        throw std::invalid_argument("CommunityProfileDesigner requires at least one genome");
    }
}

std::unordered_map<std::string, std::vector<GenomeRecord>> CommunityProfileDesigner::group_by_species() const {
    std::unordered_map<std::string, std::vector<GenomeRecord>> grouped;
    for (const auto& genome : genomes_) {
        grouped[extract_species(genome.taxonomy)].push_back(genome);
    }
    return grouped;
}

std::vector<double> CommunityProfileDesigner::draw_weights(
    std::size_t count, const ProfileDesignOptions& options, std::mt19937_64& rng) const {
    std::vector<double> weights;
    weights.reserve(count);
    if (options.distribution == AbundanceDistribution::PowerLaw) {
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        const double alpha = options.powerlaw_alpha > 0.0 ? options.powerlaw_alpha : 1.0;
        for (std::size_t i = 0; i < count; ++i) {
            double u = std::clamp(uniform(rng), 1e-12, 0.999999999999);
            weights.push_back(std::pow(1.0 - u, -1.0 / alpha));
        }
    } else if (options.distribution == AbundanceDistribution::NegativeBinomial) {
        const int r = std::max(1, options.negative_binomial_r);
        const double p = std::clamp(options.negative_binomial_p, 1e-6, 0.999999);
        std::negative_binomial_distribution<int> nb(r, p);
        for (std::size_t i = 0; i < count; ++i) {
            weights.push_back(static_cast<double>(nb(rng) + 1));
        }
    }
    return weights;
}

std::vector<GenomeAssignment> CommunityProfileDesigner::design_profile(
    const ProfileDesignOptions& options, std::mt19937_64& rng) const {
    auto grouped = group_by_species();
    std::vector<GenomeRecord> picked;
    picked.reserve(options.genomes_per_sample);

    // Honor user-specified strains per species first.
    for (const auto& [species, requested] : options.strains_per_species) {
        auto it = grouped.find(species);
        if (it == grouped.end() || it->second.empty()) {
            continue;
        }
        const auto& genomes = it->second;
        std::vector<std::size_t> indices(genomes.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), rng);
        const std::size_t to_take = std::min<std::size_t>(requested, genomes.size());
        for (std::size_t i = 0; i < to_take; ++i) {
            picked.push_back(genomes[indices[i]]);
        }
    }

    // Fill the remaining slots with random genomes across all species.
    std::vector<std::size_t> all_indices(genomes_.size());
    std::iota(all_indices.begin(), all_indices.end(), 0);
    std::shuffle(all_indices.begin(), all_indices.end(), rng);
    for (std::size_t idx : all_indices) {
        if (picked.size() >= options.genomes_per_sample) {
            break;
        }
        const auto& genome = genomes_[idx];
        // Avoid picking duplicates.
        if (std::find_if(
                picked.begin(), picked.end(), [&](const GenomeRecord& g) { return g.name == genome.name; }) ==
            picked.end()) {
            picked.push_back(genome);
        }
    }

    if (picked.empty()) {
        throw std::runtime_error("No genomes selected for profile design");
    }

    auto weights = draw_weights(picked.size(), options, rng);
    const double sum_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
    std::vector<double> rel_abundances;
    rel_abundances.reserve(picked.size());
    for (double w : weights) {
        rel_abundances.push_back(sum_weights > 0.0 ? w / sum_weights : 1.0 / picked.size());
    }

    // Convert relative abundances to read counts so the total matches total_read_pairs.
    std::vector<std::uint64_t> counts;
    counts.reserve(picked.size());
    for (double rel : rel_abundances) {
        auto c = static_cast<std::uint64_t>(std::llround(rel * static_cast<double>(options.total_read_pairs)));
        counts.push_back(c == 0 ? 1 : c);
    }
    // Adjust to match total exactly.
    std::int64_t diff = static_cast<std::int64_t>(options.total_read_pairs) -
                        static_cast<std::int64_t>(std::accumulate(counts.begin(), counts.end(), std::uint64_t{0}));
    if (diff != 0 && !counts.empty()) {
        std::uniform_int_distribution<std::size_t> pick(0, counts.size() - 1);
        while (diff != 0) {
            std::size_t idx = pick(rng);
            if (diff > 0) {
                ++counts[idx];
                --diff;
            } else if (counts[idx] > 1) {
                --counts[idx];
                ++diff;
            }
        }
    }

    std::vector<GenomeAssignment> assignments;
    assignments.reserve(picked.size());
    for (std::size_t i = 0; i < picked.size(); ++i) {
        assignments.push_back(
            GenomeAssignment{picked[i], extract_species(picked[i].taxonomy), counts[i], rel_abundances[i]});
    }
    return assignments;
}

}  // namespace protal::sim
