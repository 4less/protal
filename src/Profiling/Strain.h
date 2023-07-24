//
// Created by fritsche on 21/02/23.
//

#pragma once

#include "SequenceRangeHandler.h"
#include "VariantHandler.h"
#include "SamHandler.h"
#include "SNP.h"
#include "SNPUtils.h"
#include "GenomeLoader.h"
#include "AlignmentUtils.h"
#include "robin_map.h"

namespace protal {
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
            auto cov = m_sequence_range_handler.CalculateCoverageVector();
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

        void GetSharedRegions(StrainLevelContainer const& other) {

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
            a.GetSequenceRangeHandler().CalculateCoverageVector();
            b.GetSequenceRangeHandler().CalculateCoverageVector();
            VariantHandler variants_a(a.GetVariantHandler().GetReference()), variants_b(b.GetVariantHandler().GetReference());

//            std::cout << "GetSharedAlignmentRegion: " << a.GetVariantHandler().GetVariants().size() << " " << b.GetVariantHandler().GetVariants().size() << std::endl;
            auto intersection = a.GetSequenceRangeHandler().Intersect(b.GetSequenceRangeHandler(), min_cov, 10);

            auto variant_vec_a = GetSNPsInRegion(intersection.GetCoverageVector(), a.GetVariantHandler());
            auto variant_vec_b = GetSNPsInRegion(intersection.GetCoverageVector(), b.GetVariantHandler());

            return SharedAlignmentRegion(intersection, variant_vec_a, variant_vec_b);
        };
    };


    template<typename T>
    class Matrix {
        using NamesDict = tsl::robin_map<std::string, size_t>;
        using NamesList = std::vector<std::string>;
        using RowType = std::vector<T>;
        using MatrixType = std::vector<RowType>;

        NamesDict m_row_names_to_index;
        NamesDict m_col_names_to_index;
        NamesList m_row_names;
        NamesList m_col_names;
        MatrixType m_matrix;

        bool any_set = false;

    public:
        Matrix<T>(size_t row_size, size_t col_size, T default_value) :
                m_matrix(row_size, RowType(col_size, default_value) ) {};

        Matrix<T>(){};

        void SetNames(NamesList& from, NamesList& to_list, NamesDict& to_dict) {
            to_list = from;
            to_dict.clear();
            for (auto i = 0; i < to_list.size(); i++) {
                to_dict.insert({to_list[i], i});
            }
        }

        void AddName(std::string new_name) {
            m_row_names_to_index.insert( { new_name, m_row_names.size() });
            m_col_names_to_index.insert( { new_name, m_col_names.size() });
            m_row_names.emplace_back(new_name);
            m_col_names.emplace_back(new_name);

            m_matrix.resize(m_matrix.size() + 1, RowType(m_matrix.size(), NAN));
            std::for_each(m_matrix.begin(), m_matrix.end(), [](RowType& row) {
                row.resize(row.size() + 1, NAN);
            });
        }

        void SetColNames(NamesList& col_names) {
            SetNames(col_names, m_col_names, m_col_names_to_index);
        }
        void SetRowNames(NamesList& row_names) {
            SetNames(row_names, m_row_names, m_row_names_to_index);
        }
        void SetColNames(NamesList&& col_names) {
            SetNames(col_names, m_col_names, m_col_names_to_index);
        }
        void SetRowNames(NamesList&& row_names) {
            SetNames(row_names, m_row_names, m_row_names_to_index);
        }

        std::string GetName(NamesList const& names, size_t index) const {
            if (names.empty()) {
                return std::to_string(index);
            } else {
                return names[index];
            }
        }

        bool HasName(std::string& name) {
            return m_row_names_to_index.contains(name);
        }

        bool HasIdx(uint32_t index) {
            return index < m_row_names.size();
        }

        bool AnySet() const {
            return any_set;
        }

        void SetValue(size_t row_idx, size_t col_idx, T value) {
            if (row_idx >= m_matrix.size()) {
                std::cerr << "row_idx out of bounds." << row_idx << std::endl;
                exit(9);
            }
            if (m_matrix.empty() || col_idx >= m_matrix.front().size()) {
                std::cerr << "col_idx out of bounds." << col_idx << std::endl;
                exit(9);
            }
            m_matrix[row_idx][col_idx] = value;
            any_set = true;
        }

        void SetValue(std::string& row_name, std::string& col_name, T value) {
            auto row_idx_find = m_row_names_to_index.find(row_name);
            auto col_idx_find = m_col_names_to_index.find(col_name);

            if (row_idx_find == m_row_names_to_index.end()) {
                std::cerr << "row has no name " << row_name << std::endl;
                exit(9);
            }
            if (col_idx_find == m_col_names_to_index.end()) {
                std::cerr << "col has no name " << col_name << std::endl;
                exit(9);
            }
            SetValue(row_idx_find->second, col_idx_find->second, value);
        }

        void SetValue(size_t row_idx, size_t col_idx, T value, bool lower_triangle) {
            if (row_idx >= m_matrix.size()) {
                std::cerr << "row_idx out of bounds." << row_idx << std::endl;
                exit(9);
            }
            if (m_matrix.empty() || col_idx >= m_matrix.front().size()) {
                std::cerr << "col_idx out of bounds." << col_idx << std::endl;
                exit(9);
            }

            if (lower_triangle) {
                if (row_idx < col_idx) std::swap(row_idx, col_idx);
                m_matrix[row_idx][col_idx] = value;
            } else {
                if (row_idx > col_idx) std::swap(row_idx, col_idx);
                m_matrix[row_idx][col_idx] = value;
            }
        }

        void SetValue(std::string& row_name, std::string& col_name, T value, bool lower_triangle) {
            auto row_idx_find = m_row_names_to_index.find(row_name);
            auto col_idx_find = m_col_names_to_index.find(col_name);

            if (row_idx_find == m_row_names_to_index.end()) {
                std::cerr << "row has no name " << row_name << std::endl;
                exit(9);
            }
            if (col_idx_find == m_col_names_to_index.end()) {
                std::cerr << "col has no name " << col_name << std::endl;
                exit(9);
            }
            SetValue(row_idx_find->second, col_idx_find->second, value, lower_triangle);
            any_set = true;
        }


        void PrintMatrix(std::ostream& os=std::cout, std::string const& sep="\t", size_t precision=4) const {
            os << fixed;
            for (auto col_idx = 0; col_idx < m_col_names.size(); col_idx++) {
                os << sep << GetName(m_col_names, col_idx);
            }
            os << std::endl;
            for (auto row_idx = 0; row_idx < m_row_names.size(); row_idx++) {
                os << GetName(m_row_names, row_idx);

                for (auto col_idx = 0; col_idx < m_col_names.size(); col_idx++) {
                    os << sep << std::setprecision(precision) << m_matrix[row_idx][col_idx];
                }
                os << std::endl;
            }
        };


    };
    using DoubleMatrix = Matrix<double>;


    static Variant& GetConsensusCall(VariantBin& bin) {
        std::sort(bin.begin(), bin.end(), [](Variant const& a, Variant const& b) {
            return a.QualitySum() > b.QualitySum();
        });
        return bin.front();
    }

    static bool VariantPass(Variant const& call, VariantBin const& bin, size_t min_var_qual_sum=50) {
        if (call.QualitySum() < min_var_qual_sum) return false;
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
    static bool MSA(MSASequenceItems const& items, std::string const& reference, MSAVector& msa, uint32_t min_cov, uint32_t min_qual_sum) {
        if (msa.size() != items.size()) {
            std::cerr << msa.size() << " != " << items.size() << " <- items" << std::endl;
            std::cerr << "Msa object must be of the same length as items" << std::endl;
            return false;
        }
        if (items.empty()) return false;

//        std::cout << "Extract coverages " << std::endl;
        CoverageVecs covs;
        for (auto i = 0; i < items.size(); i++) {
            auto& optional_item = items[i];
            auto& msa_item = msa[i];
            covs.emplace_back(optional_item.has_value() ?
                              optional_item.value().second.get().CalculateCoverageVector() :
                              CoverageVec());

//            if (optional_item.has_value()) {
//                std::cout << i << " cov_length: " <<  covs.back().size() << std::endl;
//            } else {
//                std::cout << i << " no cov" << std::endl;
//            }

            msa_item.reserve(msa_item.size() + covs.back().size());
        }

        int valid_bases = 0;
        for (auto cv : covs) {
            for (auto c : cv) {
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
//        std::vector<size_t> msa_indices(0, items.size());
//        std::vector<OptionalVariantBin> column;
        std::vector<OptionalVariant> column( items.size(), OptionalVariant{} );
        std::vector<bool> column_pass( items.size(), false );

        bool had_indel = false;


//        std::cout << "Loop over: " << reference.length() << std::endl;
        for (size_t rpos = 0; rpos < reference.length(); rpos++) {
//            std::cout << "rpos: " << rpos << std::endl;
            char ref = reference[rpos];
            auto max_ins = 0;

            for (auto i = 0; i < items.size(); i++) {
                auto& cov = covs[i];
                if (!items[i].has_value() || (rpos < cov.size() && cov[rpos] == 0)) {
                    column[i] = OptionalVariant();
                    continue;
                }
                const VariantVec& const_variants = items[i].value().first;
                VariantVec& variants = const_cast<VariantVec&>(const_variants);



                // update index;
                while (indices[i] < variants.size() && variants[indices[i]].front().Position() < rpos) indices[i]++;
//                std::cout << (indices[i] == variants.size() || variants[indices[i]].front().Position() != rpos) << std::endl;
                if (indices[i] == variants.size() || variants[indices[i]].front().Position() != rpos) {
                    column[i] = OptionalVariant();
                } else {
                    auto& variant = variants[indices[i]];
                    auto& call = GetConsensusCall(variant);
                    bool pass = VariantPass(call, variant, min_qual_sum);

                    column_pass[i] = pass;
                    column[i] = call;

                    if (pass && column[i]->IsINS()) max_ins = column[i]->GetStructuralSize() > max_ins ? column[i]->GetStructuralSize() : max_ins;
                }
            }

            bool current_had_indel = false;
            bool bad = false;
            size_t before = 0;
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

                if (!items[i].has_value() || (rpos < cov.size() && cov[rpos] < min_cov) || rpos >= cov.size()) {
                    for (auto j = max_ins; j > 0; j--) msa_row.emplace_back('N'); // changed from '-'
                    msa_row.emplace_back('N');
                } else if (cov[rpos] > 0) {
                    if (!var.has_value()) {
                        for (auto j = max_ins; j > 0; j--) msa_row.emplace_back('-');
                        msa_row.emplace_back(ref);
                    } else {
//                        std::cout << var->ToString() << " pass: " << var_pass << std::endl;
                        if (!var_pass) {
                            // Variant does not pass - add 'N' for ambiguous base.
                            AddInsertionGap(msa_row, max_ins);
                            msa_row.emplace_back('N');
                        } else if (var->IsSNP()) {
                            // Variant passes and is SNP
                            AddInsertionGap(msa_row, max_ins);
                            msa_row.emplace_back(var->GetVariant());
                        } else if (var->IsINS()) {
                            // Variant passes and is Insertion
                            had_indel = true;
                            current_had_indel=true;
                            for (auto c: var->GetStructural()) msa_row.emplace_back(c);
                            AddInsertionGap(msa_row, max_ins - var->GetStructuralSize());
                            msa_row.emplace_back(ref);
                        } else {
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