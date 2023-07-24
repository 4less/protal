//
// Created by fritsche on 21/02/23.
//

#pragma once

#include "Variant.h"
#include "SNP.h"
#include "robin_map.h"
#include "SNPUtils.h"
#include "GenomeLoader.h"

namespace protal {
    class VariantHandler {
        uint16_t m_quality_version_offset = 36;
        Variants m_variants;
        const std::string& m_reference;

    public:
        VariantHandler(const std::string& reference) : m_reference(reference) {};

        bool HasVariantBin(VariantPos position) {
            return m_variants.contains(position);
        }

        VariantBin& GetVariantBin(VariantPos position) {
            return m_variants.at(position);
        }

        bool HasReference(VariantBin& variant_bin, VariantPos pos) {
            return std::any_of(variant_bin.begin(), variant_bin.end(), [](Variant const& var) {
                return var.IsReference();
            });
        }

        Variant& GetVariant(VariantBin& variant_bin, VariantPos pos, Base snp, Base ref) {
            for (auto& variant :  variant_bin) {
                if (variant.Match(pos, snp, ref)) {
                    return variant;
                }
            }
            variant_bin.emplace_back( Variant(pos, snp, ref) );
            return variant_bin.back();
        }


        const std::string& GetReference() const {
            return m_reference;
        }

        std::string GetReference() {
            return m_reference;
        }

        Variant& GetVariant(VariantBin& variant_bin, VariantType type, VariantPos pos, Base ref, std::string& structural) {
            for (auto& variant :  variant_bin) {
                if (variant.Match(type, pos, structural)) {
                    return variant;
                }
            }
            variant_bin.emplace_back( Variant(type, pos, ref, structural) );
            return variant_bin.back();
        }

        const Variants& GetVariants() const {
            return m_variants;
        }

        void AddSNP(VariantPos position, Base snp, Base ref, bool on_forward, Qual quality) {
            auto& variant_bin = HasVariantBin(position) ? GetVariantBin(position) : m_variants[position];
            auto& variant = GetVariant(variant_bin, position, snp, ref);
            variant.AddObservation(quality, on_forward);
        }

        void AddINDEL(VariantType type, VariantPos position, Base ref, std::string&& structural, bool on_forward, Qual quality) {
            auto& variant_bin = HasVariantBin(position) ? GetVariantBin(position) : m_variants[position];
            auto& variant = GetVariant(variant_bin, type, position, ref, structural);
            variant.AddObservation(quality, on_forward);
        }

        bool ExtractVariants(SamEntry const& sam, size_t read_id = 0) {
            int qpos = 0;
            int rpos = sam.m_pos - 1;

            int count = 0;
            char op = ' ';

            std::string query = "";
            std::string ref = "";

            std::string align = "";
            std::string cigar = "";
            int cpos = 0;
            bool faulty = false;



            SNP snp;
            snp.readid = read_id;
            snp.orientation = sam.IsReversed();
            bool is_fwd = !sam.IsReversed();

            bool output = false;
//            bool output = sam.m_rname == "54758_110";

            if (output) {
                std::cout << sam.ToString() << std::endl;
                PrintAlignment(sam, m_reference, std::cout);
            }
            while (NextCompressedCigar(cpos, sam.m_cigar, count, op)) {
                if (op == 'M') {
                    for (auto i = 0; i < count; i++) {
                        if (sam.m_seq[qpos + i] != m_reference[rpos+i] && sam.m_seq[qpos + i] != 'N' && m_reference[rpos+i] != 'N') {
                            std::cerr << sam.m_seq[qpos + i] << " " << m_reference[rpos+i] << std::endl;
                            faulty = true;
                        }
                    }
                }

                if (op == 'I') {
                    // Mean over qualities in insertion
                    auto qual_sum = 0;
                    for (auto i = 0; i < count; i++) {
                        qual_sum += sam.m_qual[qpos + i] - m_quality_version_offset;
                    }
                    AddINDEL(VariantType::INS, rpos, m_reference[rpos], sam.m_seq.substr(qpos, count), is_fwd, qual_sum/count);
                    if (output) {
                        std::cout << rpos << " Insertion of " << count << " (" << sam.m_seq.substr(qpos, count) << ") at " << qpos << ", " << rpos << std::endl;
                    }
                } else if (op == 'D') {
                    // Deletion has  no quality because it is not present in read. (Maybe incorporate following bases)
                    AddINDEL(VariantType::DEL, rpos, m_reference[rpos], m_reference.substr(rpos, count), is_fwd, 0);
                    if (output) {
                        std::cout << rpos << " Deletion of  " << count << " (" << m_reference.substr(rpos, count) << ") at "
                                  << qpos << ", " << rpos << std::endl;
                    }
                } else if (op == 'X') {
                    for (auto i = 0; i < count; i++) {
//                        std::cout << "SNP: " << sam.m_seq[qpos + i] << " " << static_cast<uint16_t>(sam.m_qual[qpos + i] - m_quality_version_offset) << std::endl;

                        AddSNP(rpos + i, sam.m_seq[qpos + i], m_reference[rpos + i], is_fwd, sam.m_qual[qpos + i] - m_quality_version_offset);
                        if (output) {
                            std::cout << rpos + i << " " << sam.m_seq[qpos + i] << " -> " << m_reference[rpos + i] << " (" << qpos+i << ", " << rpos+i
                                      << ") Qual: " << static_cast<int>(sam.m_qual[qpos + i]-33) << std::endl;
                        }
                    }
                }

                qpos += (op != 'D') * count;
                rpos += (!(op == 'I' || op == 'S')) * count;
            }

            if (output) {
                std::cout << std::flush << std::endl;
                Utils::Input();
            }

            if (faulty) {
                std::cerr << "FAULTY ---------------------------------" << std::endl;
                PrintAlignment(sam, m_reference, std::cerr);
                std::cerr << sam.m_seq << std::endl;
                std::cerr << "Faulty sam: \n" << sam.ToString() << std::endl;
                return false;
            }
            return true;
        }

        static std::string VariantBinToString(const VariantBin& variant_bin) {
            std::string str;
            for (auto& variant : variant_bin) {
                str += variant.ToString() + '\t';
            }
            return str;
        }

        bool FilterSNPs(Variant& var, size_t coverage, size_t min_observations, size_t min_observations_fwdrev, double min_frequency, size_t min_avg_quality) {
            auto observations = var.Observations();
            auto frequency = static_cast<double>(observations) / coverage;
            auto mean_qual = var.MeanQuality();

            bool valid =
                    (observations >= min_observations ||
                    var.HasFwdAndRev() && observations >= min_observations_fwdrev) &&
                    frequency >= min_frequency &&
                    mean_qual >= min_avg_quality;

            return valid;
        }

        void FilterSNPs(VariantBin& bin, size_t coverage, size_t min_observations, size_t min_observations_fwdrev, double min_frequency, size_t min_avg_quality) {
            for (auto& var : bin) {
                auto is_valid = FilterSNPs(var, coverage, min_observations, min_observations_fwdrev, min_frequency, min_avg_quality);
                var.SetValid(is_valid);
            }
        }

        void PostProcessSNPBin(VariantBin& bin, size_t coverage, size_t min_observations=5, size_t min_observations_fwdrev=3, double min_frequency=0.2, size_t min_avg_quality=15) {
            auto var_pos = bin.front().Position();

            // Total variant var_observations
            auto var_observations = std::accumulate(bin.begin(), bin.end(), 0, [](size_t acc, Variant const & a) {
                return acc + a.Observations();
            });

            // Sort SNPs based on
            std::sort(bin.begin(), bin.end(), [](Variant const& a, Variant const& b) {
                return a.Observations() > b.Observations();
            });

            // Average quality for bin.
            auto qual = std::accumulate(bin.begin(), bin.end(), 0, [](size_t acc, Variant const & a) {
                return acc + a.MeanQuality();
            });
            qual /= bin.size();

            // No reads carrying reference allele (all variants)
            if (coverage == var_observations) {
                return;
            }

            // Reference allele (also stored in insertions, deleteions)
            char ref = bin.front().Reference();

            auto& variant = GetVariant(bin, var_pos, ref, ref);
            variant.SetObservations(coverage - var_observations);

            FilterSNPs(bin, coverage, min_observations, min_observations_fwdrev, min_frequency, min_avg_quality);
        }

        void PostProcessSNPs(std::vector<uint16_t>& coverage, size_t min_observations=2, size_t min_observations_fwdrev=2, double min_frequency=0.2, size_t min_avg_quality=15) {
            // Iterate all variant positions.
            for (auto& [variant_pos, variant_bin] : m_variants) {
                auto& bin = m_variants.at(variant_pos);
                auto cov = coverage[variant_pos];

                PostProcessSNPBin(bin, cov, min_observations, min_observations_fwdrev, min_frequency, min_avg_quality);
            }
        }

        void AddVariantsFromSam(SamEntry const& sam, size_t read_id = 0) {
            ExtractVariants(sam, read_id);
        }

    };
}