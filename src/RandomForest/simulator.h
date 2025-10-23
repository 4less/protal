#ifndef PROTAL_RANDOMFOREST_METAGENOME_SIMULATOR_H
#define PROTAL_RANDOMFOREST_METAGENOME_SIMULATOR_H

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <random>
#include <algorithm>
#include <cstdint>
#include <utility>
#include <fstream>

struct SimulationOptions {
    // Total number of reads to simulate (sequencing depth)
    uint64_t sequencing_depth = 0;

    // Per-base sequencing error rate (0.0 - 1.0)
    double error_rate = 0.0;

    // Mean insert size for paired-end reads
    uint32_t insert_size = 300;

    // Read length for single-end / each mate in paired-end
    uint32_t read_length = 150;
};

class MetagenomeSimulator {
public:
    MetagenomeSimulator() = default;

    // Load accession -> path from a file where each line is:
    // accession<TAB>path
    // Returns true on success (file opened and at least one mapping loaded).
    bool load_accession_to_path(const std::string& filepath) {
        m_accession_to_path.clear();
        std::ifstream input(filepath);
        if (!input.is_open()) return false;
        std::string line;
        size_t loaded = 0;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#') continue;
            auto tab = line.find(   '\t');
            if (tab == std::string::npos) continue;
            std::string accession = line.substr(0, tab);
            std::string path = line.substr(tab + 1);
            if (accession.empty() || path.empty()) continue;
            m_accession_to_path[accession] = path;
            ++loaded;
        }
        return loaded > 0;
    }

    // Load accession -> species from a file where each line is:
    // accession<TAB>species
    // This will populate both m_accession_to_species_ and the reverse
    // m_species_to_accessions_.
    bool load_accession_to_species_and_reverse(const std::string& filepath) {
        m_accession_to_species.clear();
        m_species_to_accessions.clear();
        std::ifstream input(filepath);
        if (!input.is_open()) return false;
        std::string line;
        size_t loaded = 0;
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#') continue;
            auto tab = line.find('\t');
            if (tab == std::string::npos) continue;
            std::string accession = line.substr(0, tab);
            std::string species = line.substr(tab + 1);
            if (accession.empty() || species.empty()) continue;
            m_accession_to_species[accession] = species;
            m_species_to_accessions[species].insert(accession);
            ++loaded;
        }
        return loaded > 0;
    }

    // Add a genome accession mapped to a path and a species name
    void add_genome(const std::string& accession, const std::string& path, const std::string& species) {
        m_accession_to_path[accession] = path;
        m_accession_to_species[accession] = species;
        m_species_to_accessions[species].insert(accession);
    }

    // Remove a genome accession if present
    void remove_genome(const std::string& accession) {
        auto it = m_accession_to_species.find(accession);
        if (it != m_accession_to_species.end()) {
            const std::string& species = it->second;
            auto sit = m_species_to_accessions.find(species);
            if (sit != m_species_to_accessions.end()) {
                sit->second.erase(accession);
                if (sit->second.empty()) m_species_to_accessions.erase(sit);
            }
            m_accession_to_species.erase(it);
        }
        m_accession_to_path.erase(accession);
    }

    // Accessors to maps
    const std::unordered_map<std::string, std::string>& accession_to_path() const { return m_accession_to_path; }
    const std::unordered_map<std::string, std::string>& accession_to_species() const { return m_accession_to_species; }
    const std::unordered_map<std::string, std::unordered_set<std::string>>& species_to_accessions() const { return m_species_to_accessions; }

    // Generate abundances for `num_genomes` genomes.
    // Returns a pair of vectors:
    //  - vector<uint64_t> : total read count per selected genome (same order)
    //  - vector<double>   : relative abundance per selected genome (sums to 1.0)
    //
    // Selection of genomes is random without replacement from available accessions.
    std::tuple<std::vector<uint64_t>, std::vector<double>, std::vector<double>>
    generate_abundances(size_t num_genomes, uint64_t total_reads, const std::vector<uint64_t>& genome_lengths) const {
        std::vector<std::string> accessions;
        const double nb_r = 5.0; // shape parameter for negative binomial
        const double nb_p = 0.5; // probability parameter
        accessions.reserve(m_accession_to_path.size());
        for (const auto& kv : m_accession_to_path) accessions.push_back(kv.first);

        if (accessions.empty() || num_genomes == 0 || total_reads == 0 || genome_lengths.size() < num_genomes) {
            return {{}, {}, {}};
        }

        if (num_genomes > accessions.size()) num_genomes = accessions.size();

        // random selection without replacement
        std::random_device rd;
        std::mt19937_64 rng(rd());
        std::shuffle(accessions.begin(), accessions.end(), rng);
        accessions.resize(num_genomes);

        // Generate counts using negative binomial distribution
        std::vector<uint64_t> counts(num_genomes);
        std::vector<double> rel_abundances(num_genomes);
        std::vector<double> sequencing_depths(num_genomes);
        std::negative_binomial_distribution<uint64_t> nb(nb_r, nb_p);
        uint64_t total_count = 0;

        // Generate initial counts
        for (size_t i = 0; i < num_genomes; ++i) {
            counts[i] = nb(rng);
            total_count += counts[i];
        }

        // Scale counts to match desired total_reads
        double scale = static_cast<double>(total_reads) / static_cast<double>(total_count);
        
        uint64_t scaled_total = 0;
        for (size_t i = 0; i < num_genomes; ++i) {
            counts[i] = static_cast<uint64_t>(counts[i] * scale);
            scaled_total += counts[i];
            rel_abundances[i] = static_cast<double>(counts[i]) / static_cast<double>(total_reads);
            // Calculate sequencing depth: (number of reads * read length) / genome length
            sequencing_depths[i] = static_cast<double>(counts[i]) / static_cast<double>(genome_lengths[i]);
        }

        // Distribute any remaining reads due to rounding
        uint64_t remaining = total_reads - scaled_total;
        while (remaining > 0) {
            for (size_t i = 0; i < num_genomes && remaining > 0; ++i) {
                counts[i]++;
                rel_abundances[i] = static_cast<double>(counts[i]) / static_cast<double>(total_reads);
                sequencing_depths[i] = static_cast<double>(counts[i]) / static_cast<double>(genome_lengths[i]);
                remaining--;
            }
        }

        return {counts, rel_abundances, sequencing_depths};
    }

    // High-level simulate function that uses SimulationOptions.
    // For now it delegates to generate_abundances using opts.sequencing_depth.
    // Returns the same pair as generate_abundances.
    std::tuple<std::vector<uint64_t>, std::vector<double>, std::vector<double>>
    simulate(const SimulationOptions& opts, size_t num_genomes, const std::vector<uint64_t>& genome_lengths) {
        return generate_abundances(num_genomes, opts.sequencing_depth, genome_lengths);
    }

private:
    std::unordered_map<std::string, std::string> m_accession_to_path;
    std::unordered_map<std::string, std::string> m_accession_to_species;
    std::unordered_map<std::string, std::unordered_set<std::string>> m_species_to_accessions;
};

#endif // PROTAL_RANDOMFOREST_METAGENOME_SIMULATOR_H