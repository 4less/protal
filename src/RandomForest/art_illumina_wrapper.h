// art_illumina_wrapper.h
#ifndef ART_ILLUMINA_WRAPPER_H
#define ART_ILLUMINA_WRAPPER_H

#include <string>

class ArtIlluminaOptions {
public:
    std::string ref_seq_file;
    std::string output_prefix;
    std::string quality_profile;
    int read_length = 0;
    bool paired_end = false;
    int fragment_length = 0;
    int fragment_std = 0;
    double fold_coverage = 0;
    int total_reads = 0;

    bool is_valid() const;
};

class ArtIlluminaWrapper {
public:
    static bool run(const ArtIlluminaOptions& options);

private:
    static std::string build_command(const ArtIlluminaOptions& options);
};

#endif // ART_ILLUMINA_WRAPPER_H
