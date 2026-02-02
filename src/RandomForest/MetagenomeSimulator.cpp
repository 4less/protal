// SPDX-License-Identifier: MIT
#include "MetagenomeSimulator.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <zlib.h>

#include "../Utilities/Benchmark.h"

namespace fs = std::filesystem;

namespace protal::sim {

static std::string lowercase(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

static std::string trim_ws(std::string text) {
    auto not_space = [](unsigned char c) { return !std::isspace(c); };
    text.erase(text.begin(), std::find_if(text.begin(), text.end(), not_space));
    text.erase(std::find_if(text.rbegin(), text.rend(), not_space).base(), text.end());
    return text;
}

static std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (start <= line.size()) {
        auto pos = line.find('\t', start);
        if (pos == std::string::npos) {
            fields.emplace_back(trim_ws(line.substr(start)));
            break;
        }
        fields.emplace_back(trim_ws(line.substr(start, pos - start)));
        start = pos + 1;
    }
    return fields;
}

static std::vector<std::string> split_semicolon(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (start <= line.size()) {
        auto pos = line.find(';', start);
        if (pos == std::string::npos) {
            fields.emplace_back(line.substr(start));
            break;
        }
        fields.emplace_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return fields;
}

static bool has_rank_prefix(const std::string& token) {
    if (token.size() < 3 || token[1] != '_' || token[2] != '_') {
        return false;
    }
    switch (token[0]) {
        case 'd':
        case 'k':
        case 'p':
        case 'c':
        case 'o':
        case 'f':
        case 'g':
        case 's':
            return true;
        default:
            return false;
    }
}

static std::string strip_rank_prefix(const std::string& token) {
    if (has_rank_prefix(token) && token.size() > 3) {
        return token.substr(3);
    }
    return token;
}

static bool taxonomy_has_unclassified(const std::string& taxonomy) {
    for (const auto& token : split_semicolon(taxonomy)) {
        if (token.empty()) {
            continue;
        }
        if (strip_rank_prefix(token) == "Unclassified") {
            return true;
        }
    }
    return false;
}

static std::string make_unclassified_lineage(const std::string& genome_name) {
    return "k__" + genome_name + ";p__" + genome_name + ";c__" + genome_name + ";o__" + genome_name +
           ";f__" + genome_name + ";g__" + genome_name + ";s__" + genome_name;
}

std::vector<GenomeRecord> read_genome_table(const fs::path& tsv_path) {
    std::cout << "Load genome table" << std::endl;
    std::ifstream in(tsv_path);
    if (!in) {
        throw std::runtime_error("Unable to open genome table: " + tsv_path.string());
    }
    std::vector<GenomeRecord> genomes;
    int name_idx = 0;
    int tax_idx = 1;
    int path_idx = 2;
    int len_idx = -1;
    bool require_length = false;
    bool header_checked = false;
    std::string line;
    std::size_t line_no = 0;
    while (std::getline(in, line)) {
        ++line_no;
        if (line.empty() || line[0] == '#') {
            continue;
        }

        auto fields = split_tab(line);
        if (!header_checked) {
            int header_name = -1;
            int header_tax = -1;
            int header_path = -1;
            int header_len = -1;
            for (int i = 0; i < static_cast<int>(fields.size()); ++i) {
                auto lf = lowercase(fields[i]);
                if (header_name == -1 &&
                    (lf.find("name") != std::string::npos || lf.find("genome") != std::string::npos ||
                     lf.find("accession") != std::string::npos)) {
                    header_name = i;
                }
                if (header_tax == -1 &&
                    (lf.find("tax") != std::string::npos || lf.find("taxonomy") != std::string::npos)) {
                    header_tax = i;
                }
                if (header_path == -1 &&
                    (lf.find("path") != std::string::npos || lf.find("fasta") != std::string::npos ||
                     lf.find("file") != std::string::npos)) {
                    header_path = i;
                }
                if (header_len == -1 && lf.find("length") != std::string::npos) {
                    header_len = i;
                }
            }
            const bool looks_like_header = header_name != -1 && header_tax != -1 && header_path != -1;
            if (looks_like_header) {
                name_idx = header_name;
                tax_idx = header_tax;
                path_idx = header_path;
                len_idx = header_len;
                require_length = header_len != -1;
                header_checked = true;
                continue;  // header row
            }
            len_idx = static_cast<int>(fields.size()) > 3 ? 3 : -1;  // fallback: fourth column if present
            header_checked = true;
        }

        if (fields.size() <= static_cast<std::size_t>(std::max({name_idx, tax_idx, path_idx}))) {
            throw std::runtime_error("Malformed line " + std::to_string(line_no) + " in " + tsv_path.string());
        }

        std::string name = fields[name_idx];
        std::string taxonomy = fields[tax_idx];
        if (taxonomy_has_unclassified(taxonomy)) {
            taxonomy = make_unclassified_lineage(name);
        }
        std::string fasta_path = fields[path_idx];
        std::optional<std::uint64_t> provided_length;
        if (len_idx >= 0 && static_cast<std::size_t>(len_idx) < fields.size()) {
            if (!fields[len_idx].empty()) {
                try {
                    provided_length = std::stoull(fields[len_idx]);
                } catch (const std::exception&) {
                    throw std::runtime_error("Invalid genome_length on line " + std::to_string(line_no) + " in " +
                                             tsv_path.string());
                }
            } else if (require_length) {
                throw std::runtime_error("Missing genome_length on line " + std::to_string(line_no) + " in " +
                                         tsv_path.string());
            }
        } else if (require_length) {
            throw std::runtime_error("Missing genome_length column on line " + std::to_string(line_no) + " in " +
                                     tsv_path.string());
        }
        genomes.push_back(GenomeRecord{
            std::move(name),
            std::move(taxonomy),
            fs::path(std::move(fasta_path)),
            provided_length});
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

std::vector<std::string> parse_species_list(const std::string& text) {
    std::vector<std::string> species;
    std::size_t start = 0;
    while (start < text.size()) {
        auto end = text.find(',', start);
        if (end == std::string::npos) {
            end = text.size();
        }
        std::string token = text.substr(start, end - start);
        token.erase(token.begin(),
                    std::find_if(token.begin(), token.end(), [](unsigned char c) { return !std::isspace(c); }));
        token.erase(std::find_if(token.rbegin(), token.rend(), [](unsigned char c) { return !std::isspace(c); }).base(),
                    token.end());
        if (!token.empty()) {
            species.push_back(token);
        }
        start = end + 1;
    }
    return species;
}

std::unordered_map<std::string, std::size_t> parse_genus_selection(const std::string& text) {
    std::unordered_map<std::string, std::size_t> genus_counts;
    auto trim = [](std::string s) {
        s.erase(s.begin(),
                std::find_if(s.begin(), s.end(), [](unsigned char c) { return !std::isspace(c); }));
        s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char c) { return !std::isspace(c); }).base(), s.end());
        return s;
    };

    std::size_t start = 0;
    while (start < text.size()) {
        auto end = text.find(',', start);
        if (end == std::string::npos) {
            end = text.size();
        }
        std::string token = trim(text.substr(start, end - start));
        if (!token.empty()) {
            auto sep = token.find(':');
            if (sep == std::string::npos) {
                throw std::runtime_error("Invalid --genus entry (expected genus:count): " + token);
            }
            std::string genus = trim(token.substr(0, sep));
            if (genus.rfind("g__", 0) == 0 && genus.size() > 3) {
                genus = genus.substr(3);  // allow g__ prefix from GTDB-style names
            }
            std::string count_str = trim(token.substr(sep + 1));
            if (genus.empty() || count_str.empty()) {
                throw std::runtime_error("Invalid --genus entry (missing genus or count): " + token);
            }
            std::size_t count = 0;
            try {
                count = static_cast<std::size_t>(std::stoull(count_str));
            } catch (const std::exception&) {
                throw std::runtime_error("Invalid --genus count for entry: " + token);
            }
            if (count > 0) {
                genus_counts[genus] += count;
            }
        }
        start = end + 1;
    }
    return genus_counts;
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
    size_t read_genomes = 0;
    for (const auto& genome : genomes) {
        if ((read_genomes % 100) == 0) {
            std::cout << "genomes processed: " << read_genomes << std::endl;
        }
        const auto len = genome.genome_length ? *genome.genome_length : read_genome_length(genome.fasta_path);
        lengths.emplace(genome.name, len);
        ++read_genomes;
    }
    return lengths;
}

//static void append_fastq(const fs::path& src, std::ofstream& dst) {
//    std::ifstream in(src, std::ios::binary);
//    if (!in) {
//        throw std::runtime_error("Unable to open FASTQ chunk: " + src.string());
//    }
//    dst << in.rdbuf();
//}

//static void append_fastq(const fs::path& src, std::ofstream& dst) {
//    std::ifstream in(src, std::ios::binary);
//    if (!in) {
//        throw std::runtime_error("Unable to open FASTQ chunk: " + src.string());
//    }
//
//    dst << in.rdbuf();
//
//    if (!dst) {
//        throw std::runtime_error("Write failed while appending " + src.string());
//    }
//
//    dst.clear();  // ← THIS IS THE CRITICAL FIX
//}

static void append_fastq(const fs::path& src, std::ofstream& dst)
{
    std::ifstream in(src, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Unable to open FASTQ chunk: " + src.string());
    }

    constexpr std::size_t bufsize = 1 << 20; // 1 MB
    std::vector<char> buffer(bufsize);

    while (in) {
        in.read(buffer.data(), buffer.size());
        std::streamsize n = in.gcount();
        if (n > 0) {
            dst.write(buffer.data(), n);
            if (!dst) {
                throw std::runtime_error("Write failed while appending " + src.string());
            }
        }
    }

    dst.flush();
    if (!dst) {
        throw std::runtime_error("Flush failed after appending " + src.string());
    }
}


SampleOutput MetagenomeSimulator::simulate_single(
        const ProfileDesignOptions& profile_options,
        const std::string& sample_name,
        const fs::path& output_dir,
        const std::unordered_map<std::string, std::uint64_t>& genome_lengths,
        std::uint64_t paired_read_length,
        bool skip_reads,
        bool keep_tmp) 
{
    auto assignments = designer_.design_profile(profile_options, rng_);
    fs::path reads_dir = output_dir / "reads";
    fs::create_directories(reads_dir);
    const fs::path sample_prefix = reads_dir / sample_name;
    fs::path r1_path = sample_prefix.string() + (skip_reads ? "_R1.fq.gz" : "_R1.fq");
    fs::path r2_path = sample_prefix.string() + (skip_reads ? "_R2.fq.gz" : "_R2.fq");
    fs::path temp_dir = output_dir / (sample_name + "_tmp");

    if (!skip_reads) {
        fs::create_directories(temp_dir);
    }

    std::ofstream r1_out;
    std::ofstream r2_out;
    if (!skip_reads) {
        r1_out.open(r1_path, std::ios::binary | std::ios::app);
        r2_out.open(r2_path, std::ios::binary | std::ios::app);
        if (!r1_out || !r2_out) {
            throw std::runtime_error("Unable to create output FASTQ files for " + sample_name);
        }
    }

    for (auto& assignment : assignments) {
            auto it_len = genome_lengths.find(assignment.genome.name);
            if (it_len == genome_lengths.end()) {
                throw std::runtime_error("Missing genome length for " + assignment.genome.name);
            }
            const auto genome_len = it_len->second;
            if (!skip_reads) {
                fs::path genome_prefix = temp_dir / assignment.genome.name;
                auto [fq1, fq2] = art_.simulate_read_pairs(
                    assignment.genome, assignment.read_pairs, genome_len, genome_prefix, rng_, temp_dir);
                append_fastq(fq1, r1_out);
                append_fastq(fq2, r2_out);
            }
            const double bases = static_cast<double>(assignment.read_pairs * paired_read_length);
            const double coverage = bases / static_cast<double>(genome_len);
        assignment.relative_abundance = coverage;  // normalized later
        assignment.genome_length = genome_len;
    }

    fs::path r1_gz = r1_path;
    fs::path r2_gz = r2_path;
    if (!skip_reads) {
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

        if (!keep_tmp) {
            fs::remove_all(temp_dir);
        }
        r1_gz = r1_path;
        r1_gz += ".gz";
        r2_gz = r2_path;
        r2_gz += ".gz";
    } else {
        // Create empty placeholder files for manifests to point to.
        std::ofstream placeholder1(r1_gz, std::ios::binary);
        std::ofstream placeholder2(r2_gz, std::ios::binary);
        (void)placeholder1;
        (void)placeholder2;
    }

    return SampleOutput{sample_name, r1_gz, r2_gz, std::move(assignments)};
}

std::vector<SampleOutput> MetagenomeSimulator::simulate_samples(
    const ProfileDesignOptions& profile_options,
    std::size_t sample_count,
    const std::string& sample_prefix,
    const fs::path& output_dir,
    bool skip_reads,
    bool keep_tmp)
{
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
        auto sample =
            simulate_single(profile_options, name.str(), output_dir, genome_lengths, paired_read_length, skip_reads, keep_tmp);
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
