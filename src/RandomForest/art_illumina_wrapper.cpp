#include "art_illumina_wrapper.h"
#include <sstream>
#include <cstdlib>

bool ArtIlluminaOptions::is_valid() const {
    // Check required fields
    if (ref_seq_file.empty() || output_prefix.empty() || quality_profile.empty()) {
        return false;
    }

    // Check that either fold_coverage or total_reads is set
    if (fold_coverage <= 0 && total_reads <= 0) {
        return false;
    }

    // Check read length
    if (read_length <= 0) {
        return false;
    }

    // For paired-end reads, check fragment length settings
    if (paired_end && (fragment_length <= 0 || fragment_std <= 0)) {
        return false;
    }

    return true;
}

std::string ArtIlluminaWrapper::build_command(const ArtIlluminaOptions& options) {
    std::stringstream cmd;
    
    // Base command
    cmd << "art_illumina";

    // Required parameters
    cmd << " -i " << options.ref_seq_file;
    cmd << " -o " << options.output_prefix;
    cmd << " -l " << options.read_length;
    cmd << " -q " << options.quality_profile;

    // Coverage or total reads
    if (options.fold_coverage > 0) {
        cmd << " -f " << options.fold_coverage;
    } else if (options.total_reads > 0) {
        cmd << " -c " << options.total_reads;
    }

    // Paired-end specific parameters
    if (options.paired_end) {
        cmd << " -p";
        cmd << " -m " << options.fragment_length;
        cmd << " -s " << options.fragment_std;
    }

    return cmd.str();
}

bool ArtIlluminaWrapper::run(const ArtIlluminaOptions& options) {
    if (!options.is_valid()) {
        return false;
    }

    std::string command = build_command(options);
    int result = std::system(command.c_str());

    return result == 0;
}