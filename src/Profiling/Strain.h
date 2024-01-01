//
// Created by fritsche on 21/02/23.
//

#pragma once

#include "Matrix.h"
#include <optional>
#include <unordered_set>
#include "VariantHandler.h"
#include "SequenceRangeHandler.h"

namespace protal {
    static Variant& GetConsensusCall(VariantBin& bin) {
        std::sort(bin.begin(), bin.end(), [](Variant const& a, Variant const& b) {
            return a.QualitySum() > b.QualitySum();
        });
        return bin.front();
    }


    struct SharedAlignmentRegion;
    using VariantVec = std::vector<VariantBin>;

    // This should be in alignment utils but there are weird circular dependencies
    static size_t AlignmentLengthRef(const std::string& cigar) {
        if (cigar.empty()) return 0;

        int score = 0;
        int count = 0;
        char c;
        int digit_start = -1;

        size_t ins = 0;
        size_t del = 0;
        size_t soft = 0;
        size_t mismatch = 0;
        size_t match = 0;

        for (auto i = 0; i < cigar.length(); i++) {
            if (!std::isdigit(cigar[i])) {
                count = stoi(cigar.substr(digit_start, i-digit_start));
                c = cigar[i];
                digit_start = -1;
                if (c == 'X') mismatch += count;
                if (c == 'M') match += count;
                if (c == 'I') ins += count;
                if (c == 'D') del += count;
            } else if (digit_start == -1) {
                digit_start = i;
            }
        }
        return del + mismatch + match;
    }

    class StrainLevelContainer {
        friend SharedAlignmentRegion;
        // Needs two things, variant handler which stores and handles SNPs and INDELS
        // And sequence range handler which stores which regions are covered by the alignment
        VariantHandler m_variant_handler;
        SequenceRangeHandler m_sequence_range_handler;

        std::vector<SNP> m_snp_tmp;
        const Gene& m_reference;

    public:
        StrainLevelContainer(Gene const& reference) :
                m_reference(reference), m_variant_handler(reference.Sequence()) {
        };

        void AddToVariants(SamEntry const& sam, size_t read_id) {
            m_variant_handler.AddVariantsFromSam(sam, read_id);

        }

        const VariantHandler& GetVariantHandler() const {
            return m_variant_handler;
        }

        const SequenceRangeHandler& GetSequenceRangeHandler() const {
            return m_sequence_range_handler;
        }

        SequenceRangeHandler& GetSequenceRangeHandler() {
            return m_sequence_range_handler;
        }

        void PostProcess(size_t min_observations=2, size_t min_observations_fwdrev=2, double min_frequency=0.2, size_t min_avg_quality=15) {
            auto cov = m_sequence_range_handler.CalculateCoverageVector2();
            m_variant_handler.PostProcessSNPs(cov, min_observations, min_observations_fwdrev, min_frequency, min_avg_quality);
        }

        void AddToSequenceRange(SamEntry const& sam, size_t read_id) {
            size_t start = sam.m_pos-1;
            size_t length = AlignmentLengthRef(sam.m_cigar);
            size_t end = start + length;

            SequenceRange query_range(start, end);
            ReadInfo rinfo;
            rinfo.read_id = read_id;
            rinfo.length = length;
            rinfo.start = start;
            rinfo.forward = !Flag::IsRead1ReverseComplement(sam.m_flag);
            query_range.AddReadInfo(rinfo);
            auto range_it = m_sequence_range_handler.FindSequenceRange(start, end);

            auto& ranges = m_sequence_range_handler.GetRanges();
            if (range_it == ranges.end() || *range_it != query_range) {
                // Range has no overlap with existing range and thus add as new
                ranges.insert(range_it, query_range);

            } else if (*range_it == query_range) {
                // Range has overlap with existing range, include in existing
                range_it->Union(query_range);
                // If updated range now overlaps with next range merge.
                if (range_it+1 != ranges.end() && *(range_it+1) == *range_it) {
                    auto del_it = range_it + 1;
                    range_it->Union(*del_it);
                    ranges.erase(del_it);
                }
            }
        }

        void AddSam(SamEntry const& sam, size_t read_id) {
//            m_sequence_range_handler.Add(sam.m_pos, sam.m_cigar.length());
            AddToSequenceRange(sam, read_id);
            AddToVariants(sam, read_id);
        }
    };

    class SharedAlignmentRegion {
    public:
        const SequenceRangeHandler share_range;
        VariantVec variant_handler_a;
        VariantVec variant_handler_b;
//        const VariantVec variant_handler_a;
//        const VariantVec variant_handler_b;

        SharedAlignmentRegion(SequenceRangeHandler const&& shared_range,
                              VariantVec const&& variant_handler_a,
                              VariantVec const&& variant_handler_b) :
                share_range(shared_range),
                variant_handler_a(variant_handler_a),
                variant_handler_b(variant_handler_b) {};

        SharedAlignmentRegion(SequenceRangeHandler const& shared_range,
                              VariantVec const& variant_handler_a,
                              VariantVec const& variant_handler_b) :
                share_range(shared_range),
                variant_handler_a(variant_handler_a),
                variant_handler_b(variant_handler_b) {};

        static VariantVec GetSNPs(VariantHandler const& a) {
            VariantVec var;

            for (auto& [pos, varbin] : a.GetVariants()) {
                var.emplace_back(varbin);

                for (auto& var : varbin) {
                    if (var.Observations() == 65535) {
                        std::cout << var.ToString() << std::endl;
                        std::cout << "Faulty Variant" << std::endl;
                        exit(3);
                    }
                }
            }

            std::sort(var.begin(), var.end(), [](VariantBin const& var1, VariantBin const& var2) {
                return var1.front().Position() < var2.front().Position();
            });

            return var;
        }

        static VariantVec GetSNPsInRegion(CoverageVec const& coverage, VariantHandler const& a) {
            VariantVec var;

//            std::cout << "GREP GetSNPsInRegion: " << a.GetVariants().size() << std::endl;
            for (auto& [pos, varbin] : a.GetVariants()) {
//                std::cout << "Pos: " << pos << " cov: " << coverage[pos] << std::endl;
                if (coverage[pos] != 0) {
                    var.emplace_back(varbin);
                    varbin.front().Position();
                }
            }
            std::sort(var.begin(), var.end(), [](VariantBin const& var1, VariantBin const& var2) {
                return var1.front().Position() < var2.front().Position();
            });

            return var;
        }

        static SharedAlignmentRegion GetSharedAlignmentRegion(StrainLevelContainer& a, StrainLevelContainer& b) {
            size_t min_cov = 3;
            size_t min_len = 30;
            SequenceRangeHandler ranges_a, ranges_b;
//            a.GetSequenceRangeHandler().CalculateCoverageVector();
//            b.GetSequenceRangeHandler().CalculateCoverageVector();
            VariantHandler variants_a(a.GetVariantHandler().GetReference()), variants_b(b.GetVariantHandler().GetReference());

//            std::cout << "GetSharedAlignmentRegion: " << a.GetVariantHandler().GetVariants().size() << " " << b.GetVariantHandler().GetVariants().size() << std::endl;
            auto intersection = a.GetSequenceRangeHandler().Intersect(b.GetSequenceRangeHandler(), min_cov, 10);

            auto variant_vec_a = GetSNPsInRegion(intersection.GetCoverageVector(), a.GetVariantHandler());
            auto variant_vec_b = GetSNPsInRegion(intersection.GetCoverageVector(), b.GetVariantHandler());

            return SharedAlignmentRegion(intersection, variant_vec_a, variant_vec_b);
        };
    };



    using DoubleMatrix = Matrix<double>;


    static bool VariantPass(Variant const& call, VariantBin const& bin, size_t min_var_qual_sum=50, size_t min_var_cov=3) {
        if (call.QualitySum() < min_var_qual_sum || call.Observations() < min_var_cov) return false;
//        std::cout << VariantHandler::VariantBinToMinimalString(bin);
//        std::cout << " -- PASS --> " << bin.front().ToString() << std::endl;
//        size_t total = std::accumulate(bin.begin(), bin.end(), 0, [](size_t acc, Variant const& v) { return acc + v.Observations(); });
//        double ratio = static_cast<double>(bin.front().Observations())/total;
//        if (!bin.front().IsReference() && bin.size() > 1 && ratio < 0.8) {
//            std::cout << VariantHandler::VariantBinToString(bin);
//            std::cout << "            (" << bin.front().Observations() << "/" << total << ")" << std::endl;
//            Utils::Input();
//        }
        return true;
    }

    using MSAVector = std::vector<std::vector<char>>;
    using MSARow = std::vector<char>;
    static void AddInsertionGap(MSARow& msa_row, size_t ins_count) {
        for (auto j = ins_count; j > 0; j--) msa_row.emplace_back('-');
    }

    using VariantVecRef = std::reference_wrapper<VariantVec>;
    using SequenceRangeHandlerRef = std::reference_wrapper<SequenceRangeHandler>;
    using OptionalMSASequenceItem = std::optional<std::pair<VariantVec, SequenceRangeHandlerRef>>;
    using MSASequenceItems = std::vector<OptionalMSASequenceItem>;
    using CoverageVecs = std::vector<CoverageVec>;
    using OptionalVariant = std::optional<Variant>;
    using OptionalVariantBin = std::optional<VariantBin>;

    static auto FindUnequalIndex(std::vector<OptionalVariant> const& column, std::vector<bool> const& column_pass) {
        auto count = std::count_if(column.begin(), column.end(), [](OptionalVariant const& v) { return v.has_value(); });
        if (count <= 1) return -1;

        char baseline = 'X';
        for (auto i = 0; i < column.size(); i++) {
            if (column_pass[i] && column[i].has_value() && column[i].value().GetVariant() != 'N') {
                if (baseline == 'X') {
                    baseline = column[i].value().GetVariant();
                } else {
                    if (column[i].value().GetVariant() != baseline)
                        return i;
                }
            }
        }

//        auto baseline = std::find_if(column.begin(), column.end(), [](OptionalVariant const& v) { return v.has_value() && v.value().GetVariant() != 'N'; });
//        for (auto& v : column) {
//            if (v.has_value() && v.value().GetVariant() != 'N' && !v.value().Match(baseline->value())) {
//                return true;
//            }
//        }
        return -1;
    }

    static auto WrongPos(std::vector<OptionalVariant> const& column) {
        auto count = std::count_if(column.begin(), column.end(), [](OptionalVariant const& v) { return v.has_value(); });
        if (count <= 1) return false;
        auto baseline = std::find_if(column.begin(), column.end(), [](OptionalVariant const& v) { return v.has_value(); });
        auto pos = baseline->value().Position();
        auto correct = std::all_of(column.begin(), column.end(), [pos](OptionalVariant const& var) {
            return !var.has_value() || var.value().Position() == pos;
        });
        std::cout << pos << ", " << count << " -> " << correct << std::endl;
        return !correct;
    }

    static auto IsVariantBinBad(VariantBin const& bin) {
        auto baseline = bin.front().Position();
        for (auto& v : bin) {
            if (baseline != v.Position()) return true;
        }
        return false;
    }

    static auto PrintColumn(std::vector<OptionalVariant> const& column) {
        std::cout << "Bad column" << std::endl;
        for (auto& var : column) {
            if (var.has_value())
                std::cout << var->ToString() << std::endl;
            else std::cout << " no " << std::endl;
        }
    }



    static bool MSA2(MSASequenceItems const& items, std::string const& reference, MSAVector& msa, uint32_t min_cov, uint32_t min_qual_sum, bool ignore_insertions) {
        // Get Coverages
        CoverageVecs covs(items.size(), std::vector<uint16_t>());
        for (auto i = 0; i < items.size(); i++) {
            auto& item = items[i];
            if (!item.has_value()) continue;
            covs[i] = item->second.get().CalculateCoverageVector2();
        }

        // Store all variants for column
        std::vector<OptionalVariant> column( items.size(), OptionalVariant{} );
        std::vector<int> index(0, items.size());

        // Iterate all positions
        for (auto rpos = 0; rpos < reference.size(); rpos++) {
            std::fill(column.begin(), column.end(), OptionalVariant{});

            for (auto i = 0; i < items.size(); i++) {
                auto& item = items[i];
                if (!item.has_value()) continue;

                auto& cov = covs[i];
                column[i] = Variant();
            }
        }
        return true;
    }


    static bool MSA(MSASequenceItems const& items, std::string const& reference, MSAVector& msa, uint32_t min_cov, uint32_t min_qual_sum) {
        if (msa.size() != items.size()) {
            std::cerr << msa.size() << " != " << items.size() << " <- items" << std::endl;
            std::cerr << "Msa object must be of the same length as items" << std::endl;
            return false;
        }
        if (items.empty()) return false;

        CoverageVecs covs;
        for (auto i = 0; i < items.size(); i++) {
            auto& optional_item = items[i];
            auto& msa_item = msa[i];
            covs.emplace_back(optional_item.has_value() ?
                              optional_item.value().second.get().CalculateCoverageVector2() :
                              CoverageVec());

            // if (optional_item.has_value()) {
            //     auto& [var, srh] = optional_item.value();
            //     std::cout << (!var.empty() ? var.front().front().ToString() + " " + var.back().front().ToString() : "empty") << std::endl;
            //     std::cout << srh.get().ToString() << std::endl;
            // }

            msa_item.reserve(msa_item.size() + covs.back().size());
        }

        int valid_bases = 0;

        for (auto& cv : covs) {
            for (auto& c : cv) {
                valid_bases += (c >= min_cov);
            }
        }
        if (valid_bases == 0) return false;

        std::vector<bool> has_variant(reference.length(), false);

//        std::cout << "Loop items and check if pos has variant " << std::endl;
        // Extract information if position has variant or not.
        for (auto& item : items) {
            if (!item.has_value()) {
                continue;
            }
            for (const auto& var : item.value().first) {
                has_variant[var.front().Position()] = true;
            }
        }

        std::vector<size_t> indices(items.size(), 0);
        std::vector<uint16_t> pause_timer(items.size(), 0);

        std::vector<OptionalVariant> column( items.size(), OptionalVariant{} );
        std::vector<bool> column_pass( items.size(), false );

        bool had_indel = false;


//        std::cout << "Loop over: " << reference.length() << std::endl;
        for (size_t rpos = 0; rpos < reference.length(); rpos++) {
            std::fill(column.begin(), column.end(), OptionalVariant{});
//            std::cout << "rpos: " << rpos << std::endl;
            char ref = reference[rpos];
            auto max_ins = 0;

            std::vector<std::string> outs(items.size(),"");

            // Iterate All ITems
            for (auto i = 0; i < items.size(); i++) {
                auto& cov = covs[i];
                if (cov.size() > reference.size()) {
                    std::cerr << "Coverage values are faulty: " << cov.size() << " > " << reference.size() << std::endl;
                    std::cerr << "abort." << std::endl;
                    exit(4);
                }
                outs[i] += std::to_string(i) + '\t' + std::to_string(rpos < cov.size() ? cov[rpos] : -1) + '\t';
                outs[i] += std::to_string(rpos) + '\t' + std::to_string(indices[i]) + '\t';
                if (!items[i].has_value() || (rpos < cov.size() && cov[rpos] == 0)) {
                    outs[i] += "A\t";
                    column[i] = OptionalVariant();
                    continue;
                }

                const VariantVec& const_variants = items[i].value().first;
                if (std::any_of(const_variants.begin(), const_variants.end(), [](std::vector<Variant> const& vv) {
                    return std::any_of(vv.begin(), vv.end(), [](Variant const&  v) {
                        return v.Observations() == 65535;
                    });
                })) {
                    std::cout << "Wrong variant" << std::endl;
                    exit(3);
                }
                VariantVec& variants = const_cast<VariantVec&>(const_variants);



                // update index;
                while (indices[i] < variants.size() && variants[indices[i]].front().Position() < rpos) indices[i]++;
//                std::cout << (indices[i] == variants.size() || variants[indices[i]].front().Position() != rpos) << std::endl;
                if (indices[i] == variants.size() || variants[indices[i]].front().Position() != rpos) {
                    outs[i] += "B(" + std::to_string(variants.size()) + ", " + std::to_string(indices[i] < variants.size() ? variants[indices[i]].front().Position() : -1) + ")\t";
                    column[i] = OptionalVariant();
                } else {
                    auto& variant = variants[indices[i]];
                    auto& call = GetConsensusCall(variant);
                    bool pass = VariantPass(call, variant, min_qual_sum);

                    outs[i] += std::to_string(pass);
                    outs[i] += "\t";

                    column_pass[i] = pass;
                    column[i] = call;

                    if (pass && column[i]->IsINS()) max_ins = column[i]->GetStructuralSize() > max_ins ? column[i]->GetStructuralSize() : max_ins;
                }
            }



            bool current_had_indel = false;
            bool bad = false;
            size_t before = 0;

            bool stop = false;


            // Iterate column
            for (auto i = 0; i < column.size(); i++) {
                auto& msa_row = msa[i];
                auto& cov = covs[i];
                auto& var = column[i];
                auto var_pass = column_pass[i];

                before = msa_row.size();

//                if (var.has_value()) std::cout << var->ToString() << std::endl;

                if (pause_timer[i] > 0) {
                    msa_row.emplace_back('-');
                    pause_timer[i]--;
                    for (auto j = max_ins; j > 0; j--) msa_row.emplace_back('-');

                    if (i > 0 && msa[i].size() != msa[i-1].size()) {
                        bad = true;
                    }
                    if (msa_row.size() == before) {
                        std::cout << "Hey this is wrong!! Pause" << std::endl;
                        if (var.has_value()) std::cout << var->ToString() << std::endl;
                        std::cout << "Pause timer " << pause_timer[i] << std::endl;
                    }
                    continue;
                }

                bool no_cov = (rpos < cov.size() && cov[rpos] == 0) || rpos >= cov.size();
                bool lower_min_cov = rpos < cov.size() && cov[rpos] < min_cov;
                bool sufficient_cov = rpos < cov.size() && cov[rpos] >= min_cov;

                const char LACKING_COVERAGE = '-'; // Changed from 'N'
                const char VARIANT_NO_PASS = 'N';
                const char REFERENCE_NO_PASS = 'N';

                if (!items[i].has_value() || (rpos < cov.size() && cov[rpos] < min_cov) || rpos >= cov.size()) {
                    for (auto j = max_ins; j > 0; j--) msa_row.emplace_back(LACKING_COVERAGE); // changed from '-'
                    msa_row.emplace_back(LACKING_COVERAGE);
                    outs[i] += "A";
                } else if (cov[rpos] > 0) {
                    outs[i] += "B";
                    if (!var.has_value()) {
                        // NO VARIANT: ---------------------------------------------------------------------------------
                        outs[i] += "C";
                        for (auto j = max_ins; j > 0; j--) msa_row.emplace_back('-');
                        msa_row.emplace_back(column_pass[i] ? ref : REFERENCE_NO_PASS);
                    } else {
                        // VARIANT: ------------------------------------------------------------------------------------
                        outs[i] += "D";
//                        std::cout << var->ToString() << " pass: " << var_pass << std::endl;
                        if (!var_pass) {
                            outs[i] += "E";
                            // NO PASS: IGNORE COLUMN ------------------------------------------------------------------
                            // Variant does not pass - add 'N' for ambiguous base.
                            AddInsertionGap(msa_row, max_ins);
                            msa_row.emplace_back(VARIANT_NO_PASS);
                        } else if (var->IsSNP()) {
                            // PASS: SNP -------------------------------------------------------------------------------
                            outs[i] += "F";
                            // Variant passes and is SNP
                            AddInsertionGap(msa_row, max_ins);
                            msa_row.emplace_back(var->GetVariant());
                        } else if (var->IsINS()) {
                            // PASS: INSERTION -------------------------------------------------------------------------
                            outs[i] += "G";
                            // Variant passes and is Insertion
                            had_indel = true;
                            current_had_indel=true;
                            for (auto c: var->GetStructural()) msa_row.emplace_back(c);
                            AddInsertionGap(msa_row, max_ins - var->GetStructuralSize());
                            msa_row.emplace_back(ref);
                        } else {
                            // PASS: DELETION --------------------------------------------------------------------------
                            outs[i] += "H";
                            // Variant passes and is Deletion
                            had_indel = true;
                            current_had_indel=true;
                            AddInsertionGap(msa_row, max_ins);
                            pause_timer[i] = column[i]->GetStructuralSize() - 1;
                            msa_row.emplace_back('-');
                        }

                    }
                }


                if (msa_row.size() == before) {
                    std::cout << "Hey this is wrong!! " << std::endl;
                    std::cout << "items[i].has_value(): " << items[i].has_value() << std::endl;
                    std::cout << "rpos: " << rpos << std::endl;
                    std::cout << "cov.size(): " << cov.size() << std::endl;
                    if (rpos < cov.size()) {
                        std::cout << "cov[rpos]: " << cov[rpos] << std::endl;
                    }
                    if (var.has_value()) std::cout << var->ToString() << std::endl;
                }

                if (i > 0 && msa[i].size() != msa[i-1].size()) {
                    std::cout << "BAD IN " << i << std::endl;
                    bad = true;
                }

                auto row_str = "";
                for (auto& c : msa_row) row_str += c;
            }

            constexpr bool debug = false;
            stop = true;
            if constexpr (debug) {
                auto unequal_index = FindUnequalIndex(column, column_pass);
                if (unequal_index != -1) {

                    auto& item = items[unequal_index].value();
                    auto& [var_handler, range_handler] = item;

                    std::cout << range_handler.get().GetRange(rpos).ToVerboseString() << std::endl;
                    for (auto& var : var_handler) {
                        std::cout << "Bin- " << var.size() << std::endl;
                        std::cout << VariantHandler::VariantBinToString(var) << std::endl;
                    }

                    std::cout << "Rpos: " << rpos << " " << column[unequal_index]->Position() << std::endl;
                    for (auto i = 0; i < items.size(); i++) {
                        outs[i] += '\t' + std::to_string((rpos < covs[i].size() ? covs[i][rpos] : -1)) + '\t';
                        outs[i] += (column[i].has_value() ? column[i]->ToString() : "NULL") + '\t';
                        outs[i] += std::to_string(column_pass[i]) + '\t';
                    }
                    stop = true;

                    char first = 'X';
                    for (auto i = 0; i < items.size(); i++) {
                        if (msa[i].back() == 'N' || msa[i].back() == '-') continue;
                        if (first == 'X') first = msa[i].back();
                        else if (first != msa[i].back()) stop = true;
                    }
                }

                if (stop) {
                    std::unordered_set<char> obs;
                    for (auto i = 0; i < msa.size(); i++) {
                        char last = msa[i].back();
                        if (last != 'N' && last != '-') obs.insert(last);
                        outs[i] += '\t';
                        outs[i] += last;
//                        std::cout << outs[i] << std::endl;
                    }
                    if (obs.size() > 1) {
                        std::cout << "---" << std::endl;
                        for (auto i = 0; i < msa.size(); i++) {
                            std::cout << column_pass[i] << "\t" << outs[i] << std::endl;
                        }
                        for (auto& e : obs) std::cout << e << std::endl;
                        Utils::Input();

                    }
                }
            }


            if (bad) {
                std::cout << "MSA" << std::endl;
                size_t first_len = msa.front().size();
                size_t show = 200;
                size_t start = first_len < show ? 0 : first_len - show;
                auto index = 0;
                for (auto& row : msa) {
                    std::cout << index++ << " " << row.size() << std::string_view(row.begin() + start, row.end()) << std::endl;
                }
                std::cout << "CURRENT INDEL" << std::endl;
                Utils::Input();
            }
        }

        return true;

//        if (had_indel) {
//            std::cout << "MSA" << std::endl;
//            for (auto& row : msa) {
//                for (auto& col : row) std::cout << col;
//                std::cout << std::endl;
//            }
//            Utils::Input();
//        }
    }
}