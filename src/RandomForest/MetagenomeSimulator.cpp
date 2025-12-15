// SPDX-License-Identifier: MIT
#include "MetagenomeSimulator.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <zlib.h>

#include "../Utilities/Benchmark.h"

namespace fs = std::filesystem;

namespace protal::sim {

static std::string lowercase(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

std::vector<GenomeRecord> read_genome_table(const fs::path& tsv_path) {
    std::cout << "Load genome table" << std::endl;
    std::ifstream in(tsv_path);
    if (!in) {
        throw std::runtime_error("Unable to open genome table: " + tsv_path.string());
    }
    std::vector<GenomeRecord> genomes;
    std::string line;
    std::size_t line_no = 0;
    while (std::getline(in, line)) {
        ++line_no;
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::istringstream iss(line);
        std::string name, taxonomy, fasta_path;
        if (!std::getline(iss, name, '\t') || !std::getline(iss, taxonomy, '\t') ||
            !std::getline(iss, fasta_path, '\t')) {
            throw std::runtime_error("Malformed line " + std::to_string(line_no) + " in " + tsv_path.string());
        }
        if (line_no == 1) {
            auto lname = lowercase(name);
            auto ltax = lowercase(taxonomy);
            auto lpath = lowercase(fasta_path);
            const bool looks_like_header =
                (lname.find("name") != std::string::npos || lname.find("genome") != std::string::npos ||
                 lname.find("accession") != std::string::npos) &&
                (ltax.find("tax") != std::string::npos || ltax.find("taxonomy") != std::string::npos) &&
                (lpath.find("path") != std::string::npos || lpath.find("fasta") != std::string::npos);
            if (looks_like_header) {
                continue;  // header row
            }
        }
        genomes.push_back(GenomeRecord{std::move(name), std::move(taxonomy), fs::path(fasta_path)});
    }
    if (genomes.empty()) {
        throw std::runtime_error("Genome table is empty: " + tsv_path.string());
    }
    return genomes;
}

std::vector<double> parse_strain_probabilities(const std::string& text) {
    std::vector<double> probs;
    std::size_t start = 0;
    while (start < text.size()) {
        auto end = text.find(',', start);
        if (end == std::string::npos) {
            end = text.size();
        }
        std::string token = text.substr(start, end - start);
        token.erase(token.begin(), std::find_if(token.begin(), token.end(), [](unsigned char c) { return !std::isspace(c); }));
        token.erase(std::find_if(token.rbegin(), token.rend(), [](unsigned char c) { return !std::isspace(c); }).base(), token.end());
        if (!token.empty()) {
            probs.push_back(std::stod(token));
        }
        start = end + 1;
    }
    return probs;
}

void write_sample_manifest(const SampleOutput& sample, const fs::path& manifest_path) {
    if (!manifest_path.parent_path().empty()) {
        fs::create_directories(manifest_path.parent_path());
    }
    std::ofstream out(manifest_path);
    if (!out) {
        throw std::runtime_error("Unable to write manifest: " + manifest_path.string());
    }
    out << "sample\tgenome\tspecies\ttaxonomy\tgenome_length\tread_pairs\trelative_abundance\tfastq_r1\tfastq_r2\n";
    for (const auto& assignment : sample.assignments) {
        out << sample.sample_name << '\t' << assignment.genome.name << '\t' << assignment.species << '\t'
            << assignment.genome.taxonomy << '\t' << assignment.genome_length << '\t' << assignment.read_pairs << '\t'
            << assignment.relative_abundance << '\t' << sample.read1_path.string() << '\t' << sample.read2_path.string()
            << '\n';
    }
}

void write_combined_manifest(const std::vector<SampleOutput>& samples, const fs::path& manifest_path) {
    if (!manifest_path.parent_path().empty()) {
        fs::create_directories(manifest_path.parent_path());
    }
    std::ofstream out(manifest_path);
    if (!out) {
        throw std::runtime_error("Unable to write manifest: " + manifest_path.string());
    }
    out << "sample\tgenome\tspecies\ttaxonomy\tgenome_length\tread_pairs\trelative_abundance\tfastq_r1\tfastq_r2\n";
    for (const auto& sample : samples) {
        for (const auto& assignment : sample.assignments) {
            out << sample.sample_name << '\t' << assignment.genome.name << '\t' << assignment.species << '\t'
                << assignment.genome.taxonomy << '\t' << assignment.genome_length << '\t' << assignment.read_pairs
                << '\t' << assignment.relative_abundance << '\t' << sample.read1_path.string()
                << '\t' << sample.read2_path.string() << '\n';
        }
    }
}

MetagenomeSimulator::MetagenomeSimulator(
    std::vector<GenomeRecord> genomes, ArtIlluminaOptions art_options, std::uint64_t seed, std::string pigz_path)
    : genomes_(std::move(genomes)),
      art_(std::move(art_options)),
      designer_(genomes_),
      rng_(seed),
      pigz_path_(std::move(pigz_path)) {}

static std::uint64_t read_genome_length(const fs::path& fasta_path) {
    auto accumulate_length = [](auto&& getter, auto&& handle) -> std::uint64_t {
        std::string line;
        std::uint64_t total = 0;
        while (getter(handle, line)) {
            if (!line.empty() && line[0] == '>') {
                continue;
            }
            for (char c : line) {
                if (std::isalpha(static_cast<unsigned char>(c))) {
                    ++total;
                }
            }
        }
        return total;
    };

    if (fasta_path.extension() == ".gz") {
        gzFile input = gzopen(fasta_path.string().c_str(), "rb");
        if (!input) {
            throw std::runtime_error("Unable to open compressed fasta: " + fasta_path.string());
        }
        std::string buffer;
        buffer.resize(8192);
        auto getter = [&](gzFile file, std::string& out) -> bool {
            char* res = gzgets(file, buffer.data(), static_cast<int>(buffer.size()));
            if (!res) {
                return false;
            }
            out.assign(res);
            while (!out.empty() && (out.back() == '\n' || out.back() == '\r')) {
                out.pop_back();
            }
            return true;
        };
        std::uint64_t len = accumulate_length(getter, input);
        gzclose(input);
        return len;
    }

    std::ifstream in(fasta_path);
    if (!in) {
        throw std::runtime_error("Unable to open fasta: " + fasta_path.string());
    }
    auto getter = [](std::ifstream& file, std::string& out) -> bool {
        return static_cast<bool>(std::getline(file, out));
    };
    return accumulate_length(getter, in);
}

static std::unordered_map<std::string, std::uint64_t> build_length_cache(const std::vector<GenomeRecord>& genomes) {
    std::unordered_map<std::string, std::uint64_t> lengths;
    lengths.reserve(genomes.size());
    for (const auto& genome : genomes) {
        lengths.emplace(genome.name, read_genome_length(genome.fasta_path));
    }
    return lengths;
}

static void append_fastq(const fs::path& src, std::ofstream& dst) {
    std::ifstream in(src, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Unable to open FASTQ chunk: " + src.string());
    }
    dst << in.rdbuf();
}

SampleOutput MetagenomeSimulator::simulate_single(
    const ProfileDesignOptions& profile_options,
    const std::string& sample_name,
    const fs::path& output_dir,
    const std::unordered_map<std::string, std::uint64_t>& genome_lengths,
    std::uint64_t paired_read_length) {
    auto assignments = designer_.design_profile(profile_options, rng_);
    fs::path reads_dir = output_dir / "reads";
    fs::create_directories(reads_dir);
    const fs::path sample_prefix = reads_dir / sample_name;
    const fs::path r1_path = sample_prefix.string() + "_R1.fq";
    const fs::path r2_path = sample_prefix.string() + "_R2.fq";
    fs::path temp_dir = output_dir / (sample_name + "_tmp");
    fs::create_directories(temp_dir);

    std::ofstream r1_out(r1_path, std::ios::binary);
    std::ofstream r2_out(r2_path, std::ios::binary);
    if (!r1_out || !r2_out) {
        throw std::runtime_error("Unable to create output FASTQ files for " + sample_name);
    }

    for (auto& assignment : assignments) {
            auto it_len = genome_lengths.find(assignment.genome.name);
            if (it_len == genome_lengths.end()) {
                throw std::runtime_error("Missing genome length for " + assignment.genome.name);
            }
            const auto genome_len = it_len->second;
            fs::path genome_prefix = temp_dir / assignment.genome.name;
            auto [fq1, fq2] =
                art_.simulate_read_pairs(assignment.genome, assignment.read_pairs, genome_len, genome_prefix, rng_, temp_dir);
            append_fastq(fq1, r1_out);
            append_fastq(fq2, r2_out);
            const double bases = static_cast<double>(assignment.read_pairs * paired_read_length);
            const double coverage = bases / static_cast<double>(genome_len);
        assignment.relative_abundance = coverage;  // normalized later
        assignment.genome_length = genome_len;
    }

    // Compress reads with pigz
    auto compress_with_pigz = [&](const fs::path& fq_path) {
        std::stringstream cmd;
        cmd << pigz_path_ << " -p " << std::max(1, art_.options().threads) << " -f \"" << fq_path.string() << "\"";
        int rc = std::system(cmd.str().c_str());
        if (rc != 0) {
            throw std::runtime_error("pigz failed on " + fq_path.string() + " with code " + std::to_string(rc));
        }
    };

    compress_with_pigz(r1_path);
    compress_with_pigz(r2_path);

    fs::remove_all(temp_dir);
    fs::path r1_gz = r1_path;
    r1_gz += ".gz";
    fs::path r2_gz = r2_path;
    r2_gz += ".gz";

    return SampleOutput{sample_name, r1_gz, r2_gz, std::move(assignments)};
}

std::vector<SampleOutput> MetagenomeSimulator::simulate_samples(
    const ProfileDesignOptions& profile_options,
    std::size_t sample_count,
    const std::string& sample_prefix,
    const fs::path& output_dir) {
    if (sample_count == 0) {
        return {};
    }
    const auto genome_lengths = build_length_cache(genomes_);
    const auto paired_read_length = static_cast<std::uint64_t>(art_.options().read_length) * 2ULL;
    std::vector<SampleOutput> outputs;
    outputs.reserve(sample_count);
    for (std::size_t i = 0; i < sample_count; ++i) {
        std::ostringstream name;
        name << sample_prefix << "_" << (i + 1);
        protal::Benchmark timer("simulate_metagenome " + name.str());
        timer.Start();
        auto sample = simulate_single(profile_options, name.str(), output_dir, genome_lengths, paired_read_length);
        timer.Stop();

        double coverage_sum = 0.0;
        for (const auto& assignment : sample.assignments) {
            coverage_sum += assignment.relative_abundance;
        }
        for (auto& assignment : sample.assignments) {
            assignment.relative_abundance =
                coverage_sum > 0.0 ? assignment.relative_abundance / coverage_sum : 0.0;
        }

        timer.PrintResults();
        outputs.push_back(std::move(sample));
    }
    return outputs;
}

void write_abundance_matrix(const std::vector<SampleOutput>& samples, const fs::path& matrix_path) {
    if (samples.empty()) {
        return;
    }
    if (!matrix_path.parent_path().empty()) {
        fs::create_directories(matrix_path.parent_path());
    }

    // Collect all species
    std::unordered_map<std::string, std::vector<double>> species_to_samples;
    const std::size_t sample_count = samples.size();
    for (const auto& sample : samples) {
        for (const auto& assignment : sample.assignments) {
            species_to_samples[assignment.species].resize(sample_count, 0.0);
        }
    }
    // Sum relative abundances per species per sample
    for (std::size_t idx = 0; idx < samples.size(); ++idx) {
        for (const auto& assignment : samples[idx].assignments) {
            species_to_samples[assignment.species][idx] += assignment.relative_abundance;
        }
    }

    std::ofstream out(matrix_path);
    if (!out) {
        throw std::runtime_error("Unable to write abundance matrix: " + matrix_path.string());
    }
    out << "species";
    for (const auto& sample : samples) {
        out << '\t' << sample.sample_name;
    }
    out << '\n';
    for (const auto& [species, values] : species_to_samples) {
        out << species;
        for (double val : values) {
            out << '\t' << val;
        }
        out << '\n';
    }
}

}  // namespace protal::sim
