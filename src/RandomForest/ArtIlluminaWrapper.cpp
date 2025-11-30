// SPDX-License-Identifier: MIT
#include "ArtIlluminaWrapper.h"

#include <array>
#include <cerrno>
#include <csignal>
#include <cstring>
#include <fstream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <cstdio>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>
#include <zlib.h>

namespace fs = std::filesystem;

namespace protal::sim {

ArtIlluminaWrapper::ArtIlluminaWrapper(ArtIlluminaOptions options)
    : options_(std::move(options)) {}

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
    std::mt19937_64& rng,
    const fs::path& temp_dir) const {
    if (read_pairs == 0 || genome_length == 0) {
        throw std::invalid_argument("read_pairs and genome_length must be greater than zero");
    }

    fs::path fasta = ensure_fasta(genome.fasta_path, temp_dir);
    if (!output_prefix.parent_path().empty()) {
        fs::create_directories(output_prefix.parent_path());
    }

    const double coverage =
        (static_cast<double>(read_pairs) * 2.0 * static_cast<double>(options_.read_length)) /
        static_cast<double>(genome_length);

    std::vector<std::string> cmd{
        options_.art_path,
        "-ss",
        options_.sequencer,
        "-i",
        fasta.string(),
        "-p",
        "-l",
        std::to_string(options_.read_length),
        "-f",
        std::to_string(coverage),
        "-m",
        std::to_string(options_.fragment_mean),
        "-s",
        std::to_string(options_.fragment_stdev),
        "-na"};

    cmd.insert(cmd.end(), options_.extra_args.begin(), options_.extra_args.end());

    const auto seed = static_cast<unsigned int>(rng());
    cmd.emplace_back("-rs");
    cmd.emplace_back(std::to_string(seed));

    cmd.emplace_back("-o");
    cmd.emplace_back(output_prefix.string());

    std::vector<std::pair<std::string, std::string>> env;
    if (options_.threads > 0) {
        env.emplace_back("OMP_NUM_THREADS", std::to_string(options_.threads));
    }

    run_command(cmd, env);

    fs::path fq1 = output_prefix.string() + "1.fq";
    fs::path fq2 = output_prefix.string() + "2.fq";
    return {fq1, fq2};
}

}  // namespace protal::sim
