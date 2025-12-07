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
    std::vector<std::pair<std::string, std::vector<GenomeRecord>>> selected_species;
    selected_species.reserve(options.species_per_sample);

    // Shuffle species order to pick initial strain per species.
    std::vector<std::string> species_order;
    species_order.reserve(grouped.size());
    for (const auto& [spec, _] : grouped) species_order.push_back(spec);
    std::shuffle(species_order.begin(), species_order.end(), rng);

    for (const auto& species : species_order) {
        if (selected_species.size() >= options.species_per_sample) break;
        const auto& genomes = grouped.at(species);
        if (genomes.empty()) continue;
        std::vector<std::size_t> idx(genomes.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::shuffle(idx.begin(), idx.end(), rng);

        // Always pick one strain for the species.
        std::vector<GenomeRecord> strains;
        strains.push_back(genomes[idx[0]]);

        // Additional strains decided by probabilities sequence.
        std::size_t next_idx = 1;
        for (double p : options.strain_probabilities) {
            if (next_idx >= idx.size()) break;
            if (std::uniform_real_distribution<double>(0.0, 1.0)(rng) <= p) {
                strains.push_back(genomes[idx[next_idx]]);
                ++next_idx;
            } else {
                break;  // stop adding further strains for this species
            }
        }
        selected_species.emplace_back(species, std::move(strains));
    }

    if (selected_species.empty()) {
        throw std::runtime_error("No genomes selected for profile design");
    }

    // Species-level abundance weights
    auto weights = draw_weights(selected_species.size(), options, rng);
    const double sum_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
    std::vector<double> species_rel;
    species_rel.reserve(selected_species.size());
    for (double w : weights) {
        species_rel.push_back(sum_weights > 0.0 ? w / sum_weights : 1.0 / selected_species.size());
    }

    // Convert species abundances to read counts matching total_read_pairs.
    std::vector<std::uint64_t> species_counts;
    species_counts.reserve(selected_species.size());
    for (double rel : species_rel) {
        auto c = static_cast<std::uint64_t>(std::llround(rel * static_cast<double>(options.total_read_pairs)));
        species_counts.push_back(c == 0 ? 1 : c);
    }
    std::int64_t diff = static_cast<std::int64_t>(options.total_read_pairs) -
                        static_cast<std::int64_t>(std::accumulate(species_counts.begin(), species_counts.end(), std::uint64_t{0}));
    if (diff != 0 && !species_counts.empty()) {
        std::uniform_int_distribution<std::size_t> pick(0, species_counts.size() - 1);
        while (diff != 0) {
            std::size_t idx = pick(rng);
            if (diff > 0) {
                ++species_counts[idx];
                --diff;
            } else if (species_counts[idx] > 1) {
                --species_counts[idx];
                ++diff;
            }
        }
    }

    // Distribute species counts across strains uniformly at random (Dirichlet via random weights).
    std::vector<GenomeAssignment> assignments;
    for (std::size_t i = 0; i < selected_species.size(); ++i) {
        const auto& species = selected_species[i].first;
        const auto& strains = selected_species[i].second;
        const std::size_t n_strains = strains.size();
        std::vector<double> w(n_strains);
        std::uniform_real_distribution<double> uni(0.0, 1.0);
        for (double& v : w) v = uni(rng);
        double sumw = std::accumulate(w.begin(), w.end(), 0.0);
        if (sumw == 0.0) sumw = static_cast<double>(n_strains);
        std::vector<std::uint64_t> strain_counts;
        strain_counts.reserve(n_strains);
        for (double v : w) {
            auto c = static_cast<std::uint64_t>(std::llround(species_counts[i] * (v / sumw)));
            strain_counts.push_back(c == 0 ? 1 : c);
        }
        // Adjust per-species counts to match species_counts[i]
        std::int64_t diff_strain = static_cast<std::int64_t>(species_counts[i]) -
                                   static_cast<std::int64_t>(std::accumulate(strain_counts.begin(), strain_counts.end(), std::uint64_t{0}));
        if (diff_strain != 0) {
            std::uniform_int_distribution<std::size_t> pick(0, n_strains - 1);
            while (diff_strain != 0) {
                std::size_t idx = pick(rng);
                if (diff_strain > 0) {
                    ++strain_counts[idx];
                    --diff_strain;
                } else if (strain_counts[idx] > 1) {
                    --strain_counts[idx];
                    ++diff_strain;
                }
            }
        }
        for (std::size_t j = 0; j < n_strains; ++j) {
            assignments.push_back(GenomeAssignment{strains[j], species, strain_counts[j], 0.0, 0});
        }
    }
    return assignments;
}

}  // namespace protal::sim
