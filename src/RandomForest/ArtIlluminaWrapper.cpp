// SPDX-License-Identifier: MIT
#include "ArtIlluminaWrapper.h"

#include <array>
#include <cerrno>
#include <csignal>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <cstdio>
#include <optional>
#include <unordered_map>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>
#include <zlib.h>

namespace fs = std::filesystem;

namespace protal::sim {

struct ArtOverrideInfo {
    std::unordered_map<std::string, std::optional<std::string>> values;

    bool has(const std::string& flag) const {
        return values.find(flag) != values.end();
    }

    std::optional<std::string> get(const std::string& flag) const {
        auto it = values.find(flag);
        if (it == values.end()) {
            return std::nullopt;
        }
        return it->second;
    }
};

static ArtOverrideInfo parse_art_overrides(const std::vector<std::string>& extra_args) {
    ArtOverrideInfo info;
    struct FlagAlias {
        const char* alias;
        const char* canonical;
        bool takes_value;
    };
    static const FlagAlias kAliases[] = {
        {"-ss", "-ss", true},
        {"--seqSys", "-ss", true},
        {"-i", "-i", true},
        {"--in", "-i", true},
        {"-l", "-l", true},
        {"--len", "-l", true},
        {"-f", "-f", true},
        {"--fcov", "-f", true},
        {"-m", "-m", true},
        {"--mflen", "-m", true},
        {"-s", "-s", true},
        {"--sdev", "-s", true},
        {"-rs", "-rs", true},
        {"--rndSeed", "-rs", true},
        {"-o", "-o", true},
        {"--out", "-o", true},
        {"-p", "-p", false},
        {"--paired", "-p", false},
        {"-na", "-na", false},
        {"--noALN", "-na", false},
        {"-sam", "-sam", false},
        {"--samout", "-sam", false},
        {"-1", "-1", true},
        {"--qprof1", "-1", true},
        {"-2", "-2", true},
        {"--qprof2", "-2", true}
    };

    for (std::size_t i = 0; i < extra_args.size(); ++i) {
        const std::string& token = extra_args[i];
        for (const auto& spec : kAliases) {
            const std::string flag(spec.alias);
            if (token == flag) {
                std::optional<std::string> value;
                if (spec.takes_value && i + 1 < extra_args.size()) {
                    value = extra_args[i + 1];
                    ++i;
                }
                info.values[spec.canonical] = value;
                break;
            }
            if (spec.takes_value) {
                const std::string prefix = flag + "=";
                if (token.rfind(prefix, 0) == 0) {
                    info.values[spec.canonical] = token.substr(prefix.size());
                    break;
                }
            }
        }
    }
    return info;
}

ArtIlluminaWrapper::ArtIlluminaWrapper(ArtIlluminaOptions options)
    : options_(std::move(options)) {}

std::optional<std::uint64_t> ArtIlluminaWrapper::seed_override() const {
    const auto overrides = parse_art_overrides(options_.extra_args);
    if (!overrides.has("-rs")) {
        return std::nullopt;
    }
    const auto value = overrides.get("-rs");
    if (!value || value->empty()) {
        return std::nullopt;
    }
    try {
        return std::stoull(*value);
    } catch (const std::exception&) {
        return std::nullopt;  // non-numeric: let ART complain about its own argument
    }
}

void ArtIlluminaWrapper::decompress_gzip(const fs::path& gz_path, const fs::path& output_path) const {
    gzFile input = gzopen(gz_path.string().c_str(), "rb");
    if (!input) {
        throw std::runtime_error("Failed to open compressed genome: " + gz_path.string());
    }
    std::ofstream output(output_path, std::ios::binary);
    if (!output) {
        gzclose(input);
        throw std::runtime_error("Failed to create decompressed genome: " + output_path.string());
    }

    std::array<char, 8192> buffer{};
    int bytes_read = 0;
    while ((bytes_read = gzread(input, buffer.data(), static_cast<unsigned int>(buffer.size()))) > 0) {
        output.write(buffer.data(), bytes_read);
    }
    gzclose(input);

    if (!output) {
        throw std::runtime_error("Failed to write decompressed genome: " + output_path.string());
    }
}

fs::path ArtIlluminaWrapper::ensure_fasta(const fs::path& fasta, const fs::path& temp_dir) const {
    if (!fs::exists(fasta)) {
        throw std::runtime_error("Genome FASTA not found: " + fasta.string());
    }
    if (fasta.extension() == ".gz") {
        fs::create_directories(temp_dir);
        fs::path decompressed = temp_dir / fasta.stem();
        decompress_gzip(fasta, decompressed);
        return decompressed;
    }
    return fasta;
}

void ArtIlluminaWrapper::run_command(
    const std::vector<std::string>& args, const std::vector<std::pair<std::string, std::string>>& env) const {
    if (args.empty()) {
        throw std::invalid_argument("Command list is empty");
    }

    std::cerr << "[simulate_metagenomes] ART command:";
    for (const auto& arg : args) {
        std::cerr << ' ' << arg;
    }
    std::cerr << '\n';

    std::vector<char*> argv;
    argv.reserve(args.size() + 1);
    for (const auto& arg : args) {
        argv.push_back(const_cast<char*>(arg.c_str()));
    }
    argv.push_back(nullptr);

    pid_t pid = fork();
    if (pid == -1) {
        throw std::runtime_error("fork() failed: " + std::string(std::strerror(errno)));
    }
    if (pid == 0) {
        for (const auto& [key, val] : env) {
            ::setenv(key.c_str(), val.c_str(), 1);
        }
        execvp(argv[0], argv.data());
        std::fprintf(stderr, "Failed to exec %s: %s\n", argv[0], std::strerror(errno));
        _exit(127);
    }

    int status = 0;
    if (waitpid(pid, &status, 0) == -1) {
        throw std::runtime_error("waitpid() failed: " + std::string(std::strerror(errno)));
    }
    if (WIFEXITED(status) && WEXITSTATUS(status) != 0) {
        std::ostringstream oss;
        oss << args[0] << " exited with status " << WEXITSTATUS(status);
        throw std::runtime_error(oss.str());
    }
    if (WIFSIGNALED(status)) {
        std::ostringstream oss;
        oss << args[0] << " terminated by signal " << WTERMSIG(status);
        throw std::runtime_error(oss.str());
    }
}

std::pair<fs::path, fs::path> ArtIlluminaWrapper::simulate_read_pairs(
    const GenomeRecord& genome,
    std::uint64_t read_pairs,
    std::uint64_t genome_length,
    const fs::path& output_prefix,
    unsigned int art_seed,
    const fs::path& temp_dir) const {
    if (read_pairs == 0 || genome_length == 0) {
        throw std::invalid_argument("read_pairs and genome_length must be greater than zero");
    }

    const auto overrides = parse_art_overrides(options_.extra_args);
    auto log_override = [&](const std::string& flag, const std::optional<std::string>& value) {
        std::cerr << "[simulate_metagenomes] ART override: " << flag;
        if (value && !value->empty()) {
            std::cerr << " -> " << *value;
        }
        std::cerr << '\n';
    };

    fs::path fasta = ensure_fasta(genome.fasta_path, temp_dir);
    fs::path output_prefix_used = output_prefix;
    if (overrides.has("-o")) {
        auto override_value = overrides.get("-o");
        if (!override_value || override_value->empty()) {
            throw std::runtime_error("ART override -o requires a value");
        }
        output_prefix_used = fs::path(*override_value);
        log_override("-o", override_value);
    }
    if (!output_prefix_used.parent_path().empty()) {
        fs::create_directories(output_prefix_used.parent_path());
    }

    const double coverage =
        (static_cast<double>(read_pairs) * 2.0 * static_cast<double>(options_.read_length)) /
        static_cast<double>(genome_length);

    std::vector<std::string> cmd;
    cmd.reserve(24 + options_.extra_args.size());
    cmd.push_back(options_.art_path);
    if (overrides.has("-ss")) {
        log_override("-ss", overrides.get("-ss"));
    } else {
        cmd.push_back("-ss");
        cmd.push_back(options_.sequencer);
    }
    if (overrides.has("-i")) {
        log_override("-i", overrides.get("-i"));
    } else {
        cmd.push_back("-i");
        cmd.push_back(fasta.string());
    }
    if (overrides.has("-p")) {
        log_override("-p", std::nullopt);
    } else {
        cmd.push_back("-p");
    }
    if (overrides.has("-l")) {
        log_override("-l", overrides.get("-l"));
    } else {
        cmd.push_back("-l");
        cmd.push_back(std::to_string(options_.read_length));
    }
    if (overrides.has("-f")) {
        log_override("-f", overrides.get("-f"));
    } else {
        cmd.push_back("-f");
        cmd.push_back(std::to_string(coverage));
    }
    if (overrides.has("-m")) {
        log_override("-m", overrides.get("-m"));
    } else {
        cmd.push_back("-m");
        cmd.push_back(std::to_string(options_.fragment_mean));
    }
    if (overrides.has("-s")) {
        log_override("-s", overrides.get("-s"));
    } else {
        cmd.push_back("-s");
        cmd.push_back(std::to_string(options_.fragment_stdev));
    }
    if (overrides.has("-sam")) {
        std::cerr << "[simulate_metagenomes] ART override: -na suppressed by -sam\n";
    } else if (overrides.has("-na")) {
        log_override("-na", std::nullopt);
    } else {
        cmd.push_back("-na");
    }

    cmd.insert(cmd.end(), options_.extra_args.begin(), options_.extra_args.end());

    if (overrides.has("-rs")) {
        log_override("-rs", overrides.get("-rs"));
    } else {
        cmd.emplace_back("-rs");
        cmd.emplace_back(std::to_string(art_seed));
    }

    if (!overrides.has("-o")) {
        cmd.emplace_back("-o");
        cmd.emplace_back(output_prefix_used.string());
    }

    std::vector<std::pair<std::string, std::string>> env;
    if (options_.threads > 0) {
        env.emplace_back("OMP_NUM_THREADS", std::to_string(options_.threads));
    }

    run_command(cmd, env);

    fs::path fq1 = output_prefix_used.string() + "1.fq";
    fs::path fq2 = output_prefix_used.string() + "2.fq";
    return {fq1, fq2};
}

}  // namespace protal::sim
