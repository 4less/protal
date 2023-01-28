//
// Created by fritsche on 01/09/22.
//

#pragma once

#include "WFA2Wrapper.h"
#include <cstdint>
#include <string>
#include <vector>
#include "Constants.h"
#include "BufferedOutput.h"
#include "VarkitInterface.h"
#include "KmerLookup.h"
#include "FastxReader.h"
#include "SamHandler.h"
#include "SNP.h"
#include <htslib/sam.h>
#include "AlignmentUtils.h"

namespace protal {
    class AlignmentResult {
        static constexpr int DEFAULT_ALIGNMENT_SCORE = INT32_MIN;
        static constexpr int DEFAULT_VALUE = UINT32_MAX;

        /* Only stores crucial alignment information, no information about header sequence etc
         * This is to save space and it is not necessary as this information can be passed down by another function.
         */
        int m_alignment_score = DEFAULT_ALIGNMENT_SCORE;

        uint32_t m_taxid = DEFAULT_VALUE;
        uint32_t m_geneid = DEFAULT_VALUE;
        int32_t m_genepos = INT32_MAX;

        bool m_forward = false;
        std::string m_cigar;
        std::string m_compressed_cigar;



    public:
        AlignmentResult(int alignment_score, std::string& cigar, uint32_t taxid, uint32_t geneid, int32_t genepos, bool forward) :
                m_alignment_score(alignment_score),
                m_cigar(cigar),
                m_taxid(taxid),
                m_geneid(geneid),
                m_genepos(genepos),
                m_forward(forward) {};

        AlignmentResult(int alignment_score, std::string&& cigar, uint32_t taxid, uint32_t geneid, int32_t genepos, bool forward) :
                m_alignment_score(alignment_score),
                m_cigar(std::move(cigar)),
                m_taxid(taxid),
                m_geneid(geneid),
                m_genepos(genepos),
                m_forward(forward) {};

        AlignmentResult() {};

        std::string& Cigar() {
            return m_cigar;
        }

        std::string Cigar() const {
            return m_cigar;
        }

        uint32_t Taxid() const {
            return m_taxid;
        }

        uint32_t GeneId() const {
            return m_geneid;
        }

        int32_t GenePos() const {
            return m_genepos;
        }

        bool Forward() const {
            return m_forward;
        }

        int AlignmentScore() const {
            return m_alignment_score;
        }


        void Set(int score, std::string& cigar) {
            m_cigar = cigar;
            m_alignment_score = score;
        }

        void Set(int score, std::string&& cigar) {
            m_cigar = cigar;
            m_alignment_score = score;
        }

        void Set(size_t taxid, size_t geneid, int32_t abs_pos, bool forward) {
            m_taxid = taxid;
            m_geneid = geneid;
            m_genepos = abs_pos;
            m_forward = forward;
        }


        void Reset() {
            m_alignment_score = DEFAULT_ALIGNMENT_SCORE;
            m_cigar.clear();
        }

        bool IsSet() const {
            return m_alignment_score != DEFAULT_ALIGNMENT_SCORE;
        }

        std::string ToString() {
            std::string str = "";
            str += std::to_string(AlignmentScore()) + '\t';
            str += Cigar() + '\t';
            str += std::to_string(Taxid()) + '\t';
            str += std::to_string(GeneId()) + '\t';
            str += std::to_string(GenePos()) + '\t';
            str += std::to_string(Forward());
            return str;
        }
    };

    static int ScorePairedAlignment(AlignmentResult const& a, AlignmentResult const& b) {
        int score = 0;
        if (a.IsSet() && b.IsSet()) {
            return (a.AlignmentScore() + b.AlignmentScore()) / 2.2;
        } else if (a.IsSet()) {
            return a.AlignmentScore() - 2;
        } else if (b.IsSet()) {
            return b.AlignmentScore() - 2;
        }
        return score;
    }

    static int ScorePairedAlignment(PairedAlignment const& pair) {
        return ScorePairedAlignment(pair.first, pair.second);
    }

    static double Score(AlignmentInfo const& info1, AlignmentInfo const& info2) {
        double score_sum = std::abs(info1.alignment_score + info2.alignment_score) * 0.25;
        double length_sum = static_cast<double>(info1.alignment_length + info2.alignment_length);

        return score_sum/length_sum;
    }

    static double Score(AlignmentInfo const& info) {
        double score_sum = std::abs(info.alignment_score) * 0.25;
        double length_sum = static_cast<double>(info.alignment_length);

        return score_sum/length_sum;
    }

    static double Score(PairedAlignment const& a) {
        AlignmentInfo info1;
        AlignmentInfo info2;
        GetAlignmentInfo(info1, a.first.Cigar());
        GetAlignmentInfo(info2, a.second.Cigar());

        return Score(info1, info2);
    }


    static int Bitscore(PairedAlignment const& a) {
        AlignmentInfo info1;
        AlignmentInfo info2;
        GetAlignmentInfo(info1, a.first.Cigar());
        GetAlignmentInfo(info2, a.second.Cigar());

        auto score1 = CigarScore(a.first.Cigar(), 1, 3, 1, 2);
        auto score2 = CigarScore(a.second.Cigar(), 1, 3, 1, 2);

        return score1 + score2;
    }

    static bool PairedAlignmentComparator(PairedAlignment const& a, PairedAlignment const& b) {
        // Return true if a > b
        AlignmentInfo info_a_1;
        AlignmentInfo info_a_2;
        AlignmentInfo info_b_1;
        AlignmentInfo info_b_2;

        GetAlignmentInfo(info_a_1, a.first.Cigar());
        GetAlignmentInfo(info_a_2, a.second.Cigar());
        GetAlignmentInfo(info_b_1, b.first.Cigar());
        GetAlignmentInfo(info_b_2, b.second.Cigar());

        double score_a = 0.0;
        double score_b = 0.0;

        bool both = a.first.IsSet() && a.second.IsSet() &&
                b.first.IsSet() && b.second.IsSet();

        if (both) {
//            std::cout << "Compare both " << std::endl;
            score_a = Score(info_a_1, info_a_2);
            score_b = Score(info_b_1, info_b_2);
        } else if (a.first.IsSet() && b.first.IsSet()) {
//            std::cout << "Compare firsts" << std::endl;
            score_a = Score(info_a_1);
            score_b = Score(info_b_1);
        } else if (a.second.IsSet() && b.second.IsSet()) {
//            std::cout << "Compare seconds" << std::endl;
            score_a = Score(info_a_2);
            score_b = Score(info_b_2);
        } else {
//            std::cout << "Compare first/second" << std::endl;
            score_a = a.first.IsSet() ? Score(info_a_2) : Score(info_a_1);
            score_b = b.first.IsSet() ? Score(info_b_2) : Score(info_b_1);
        }



//        std::cout << "Comparator: " << score_a << " > " << score_b << " = " << (score_a > score_b) << std::endl;

        return score_a < score_b;
    }

    static void ArtoSAM(SamEntry &sam, AlignmentResult const& ar, AlignmentInfo &info, FastxRecord &record) {
        sam.m_qname = record.id;
        sam.m_flag = 0;
        sam.m_rname = std::to_string(ar.Taxid()) + "_" + std::to_string(ar.GeneId());
        sam.m_pos = ar.GenePos() + info.gene_start_offset;
        sam.m_mapq = ar.AlignmentScore();
        sam.m_cigar = info.compressed_cigar;
        sam.m_rnext = "*";
        sam.m_pnext = 0;
        sam.m_tlen = ar.Cigar().length();
        sam.m_seq = ar.Forward() ? record.sequence : KmerUtils::ReverseComplement(record.sequence);
        sam.m_qual = record.quality;
    }

    using SNPList = std::vector<SNP>;

    class ProtalAlignmentDataOutputHandler {
    private:
        std::ostream& m_os;
        BufferedStringOutput m_output;
        AlignmentInfo m_info;
        double m_min_cigar_ani;

        std::string m_non_existent_alignment = ',' + std::to_string(INT32_MIN) + ',' + std::to_string(INT32_MIN) + ',' + std::to_string(INT32_MIN) + ',' + std::to_string(INT32_MIN) + ',' + std::to_string(INT32_MIN) + ",0";
    public:
        ProtalAlignmentDataOutputHandler(std::ostream& os, size_t buffer_capacity, double min_cigar_ani) :
            m_os(os),
            m_output(buffer_capacity),
            m_min_cigar_ani(min_cigar_ani) {}

        ProtalAlignmentDataOutputHandler(ProtalAlignmentDataOutputHandler const& other) :
                m_os(other.m_os),
                m_output(other.m_output.Capacity()),
                m_min_cigar_ani(other.m_min_cigar_ani) {}

        ~ProtalAlignmentDataOutputHandler() {
#pragma omp critical(adata_output)
            m_output.Write(m_os);
        }

        void operator () (PairedAlignmentResultList& alignment_results, FastxRecord& record1, FastxRecord& record2, size_t read_id=0, bool first_pair=true) {
            if (alignment_results.empty()) return;

            auto &best = alignment_results.front();


            if (CigarANI(best.first.Cigar()) < m_min_cigar_ani) {
                return;
            }



            /* ##############################################################
             * Output alignments.
             */

            bool first = true;

            std::string alignment_data_str = "";

//            std::cout << record1.id << std::endl;

            auto [taxonomic_id, gene_id] = KmerUtils::ExtractTaxIdGeneId(record1.id);

            bool global_label = false;
            size_t max_al = 5;
            size_t al_count = 0;

            auto best_taxid = alignment_results.front().first.IsSet() ? alignment_results.front().first.Taxid() : alignment_results.front().second.Taxid();
            bool best_true = best_taxid == taxonomic_id;
            bool ambiguous_best = alignment_results.size() > 1 &&
                    Score(alignment_results[0]) == Score(alignment_results[1]);


//            if (!best_true && !ambiguous_best) {
//                std::cout << record1.id << std::endl;
//                for (auto &[ar1, ar2]: alignment_results) {
//                    if (ar1.IsSet()) std::cout << "1: " << ar1.ToString() << std::endl;
//                    if (ar2.IsSet()) std::cout << "2: " << ar2.ToString() << std::endl;
//                    std::cout << "Double score: " << Score({ ar1, ar2 }) << std::endl;
//                    std::cout << "Bit score:    " << Bitscore({ ar1, ar2 }) << std::endl;
//                    std::cout << "---" << std::endl;
//
//                }
//                Utils::Input();
//            }

            for (auto &[ar1, ar2]: alignment_results) {
//                std::cout << al_count+1 << std::endl;
//                std::cerr << ar1.ToString() << "\n" << ar2.ToString() << std::endl;

                bool both = ar1.IsSet() && ar2.IsSet();

                int alignment_length = 0;
                int alignment_score = 0;
                int alignment_mismatch_quality = 0;
                int adjusted_score = 0;

                int pred = ar1.IsSet() ? ar1.Taxid() : ar2.Taxid();
                bool label = pred == taxonomic_id;
                global_label |= label;

//                if (ar1.IsSet()) {
//                    std::cout << "    " << ar1.Taxid() << " " << ar1.GeneId() << std::endl;
//                } else {
//                    std::cout << "    " << ar2.Taxid() << " " << ar2.GeneId() << std::endl;
//                }

                if (ar1.IsSet()) {
                    GetExtendedAlignmentInfo(m_info, ar1.Cigar(), record1.sequence, record2.quality);

                    alignment_length += m_info.alignment_length;
                    alignment_score += m_info.alignment_score;
                    alignment_mismatch_quality += m_info.mismatch_quality_sum;

                }
                if (ar2.IsSet()) {
                    GetExtendedAlignmentInfo(m_info, ar2.Cigar(), record2.sequence, record2.quality);

                    alignment_length += m_info.alignment_length;
                    alignment_score += m_info.alignment_score;
                    alignment_mismatch_quality += m_info.mismatch_quality_sum;
                }
                if (both) {

                }



                if (al_count < max_al) {
                    if (al_count > 0) alignment_data_str += ',';
                    alignment_data_str += std::to_string(pred);
                    alignment_data_str += ',' + std::to_string(alignment_score);
                    alignment_data_str += ',' + std::to_string(alignment_length);
                    alignment_data_str += ',' + std::to_string(alignment_mismatch_quality);
                    alignment_data_str += ',' + std::to_string(both);
                    alignment_data_str += ',' + std::to_string(label);
                }

                al_count++;
            }



            for (auto i = al_count; i < max_al; i++) {
                alignment_data_str += m_non_existent_alignment;
            }

            alignment_data_str += ',' + std::to_string(global_label);
            alignment_data_str += ',' + record1.id + '\n';
//            std::cout << alignment_data_str << std::endl;


            if (!m_output.Write(alignment_data_str)) {
#pragma omp critical(adata_output)
                m_output.Write(m_os);
            }

//            Utils::Input();

            // output stats

        }
     };

    /*
     * Protal Output Handler
     */
    class ProtalOutputHandler {
    private:
        std::ostream& m_sam_os;
        BufferedStringOutput m_sam_output;
        SamEntry m_sam;

        AlignmentInfo m_info;

        double m_min_cigar_ani = 0;
    public:
        size_t alignments = 0;

        ProtalOutputHandler(std::ostream& sam_os, size_t varkit_buffer_capacity, size_t sam_buffer_capacity, double min_cigar_ani=0.0f) :
                m_sam_os(sam_os),
                m_sam_output(sam_buffer_capacity),
                m_min_cigar_ani(min_cigar_ani) {}

        ProtalOutputHandler(ProtalOutputHandler const& other) :
                m_sam_os(other.m_sam_os),
                m_sam_output(other.m_sam_output.Capacity()),
                m_min_cigar_ani(other.m_min_cigar_ani) {}

        ~ProtalOutputHandler() {
#pragma omp critical(sam_output)
            m_sam_output.Write(m_sam_os);
        }


        void operator () (AlignmentResultList& alignment_results, FastxRecord& record, size_t read_id=0, bool first_pair=true) {
            if (alignment_results.empty()) return;
//            bool multiple_best = alignment_results.size() > 1 && alignment_results[0].AlignmentScore() == alignment_results[1].AlignmentScore();
//
//            if (multiple_best) {
//                return;
//            }

            auto& best = alignment_results.front();
            if (protal::CigarANI(best.Cigar()) < m_min_cigar_ani) {
                return;
            }

            /* ##############################################################
             * Output alignments.
             */
            bool first = true;
            for (auto& ar :  alignment_results) {
                GetExtendedAlignmentInfo(m_info, ar.Cigar(), record.sequence, record.quality);
                ArtoSAM(m_sam, ar, m_info, record);
                Flag::SetPairedEnd(m_sam.m_flag, false, false, first_pair, !first_pair);
                Flag::SetRead2Unmapped(m_sam.m_flag, true);
                Flag::SetRead1ReverseComplement(m_sam.m_flag, ar.Forward());
                Flag::SetNotPrimaryAlignment(m_sam.m_flag, !first);

                alignments++;
                if (!m_sam_output.Write(m_sam.ToString())) {
#pragma omp critical(sam_output)
                    m_sam_output.Write(m_sam_os);
                }
                first = false;
            }
        }
    };


    /*
     * Protal Output Handler
     */
    class ProtalPairedOutputHandler {
    private:
        std::ostream& m_sam_os;
        BufferedStringOutput m_sam_output;
        SamEntry m_sam1;
        SamEntry m_sam2;

        AlignmentInfo m_info;

        size_t m_max_out = 1;
        double m_min_cigar_ani = 0;
    public:
        size_t alignments = 0;

        ProtalPairedOutputHandler(std::ostream& sam_os, size_t max_out, size_t varkit_buffer_capacity, size_t sam_buffer_capacity, double min_cigar_ani=0.0f) :
                m_sam_os(sam_os),
                m_sam_output(sam_buffer_capacity),
                m_min_cigar_ani(min_cigar_ani),
                m_max_out(max_out) {}

        ProtalPairedOutputHandler(ProtalPairedOutputHandler const& other) :
                m_sam_os(other.m_sam_os),
                m_sam_output(other.m_sam_output.Capacity()),
                m_min_cigar_ani(other.m_min_cigar_ani),
                m_max_out(other.m_max_out) {}

        ~ProtalPairedOutputHandler() {
#pragma omp critical(sam_output)
            m_sam_output.Write(m_sam_os);
        }


        void operator () (PairedAlignmentResultList& alignment_results, FastxRecord& record1, FastxRecord& record2, size_t read_id=0, bool first_pair=true) {
            if (alignment_results.empty()) return;

            auto& best = alignment_results.front();

            if (CigarANI(best.first.Cigar()) < m_min_cigar_ani) {
                return;
            }

            /* ##############################################################
             * Output alignments.
             */

            bool first = true;

            std::string alignment_data_str = "";

            size_t output_counter = 0;
            for (auto& [ar1, ar2] :  alignment_results) {
                bool both = ar1.IsSet() && ar2.IsSet();

                int alignment_length = 0;
                int alignment_score = 0;
                int adjusted_score = 0;

                if (ar1.IsSet()) {

                    GetExtendedAlignmentInfo(m_info, ar1.Cigar(), record1.sequence, record2.quality);
                    ArtoSAM(m_sam1, ar1, m_info, record1);
                    Flag::SetPairedEnd(m_sam1.m_flag, true, both, true);
                    Flag::SetRead2Unmapped(m_sam1.m_flag, !ar2.IsSet());
                    Flag::SetRead1ReverseComplement(m_sam1.m_flag, ar1.Forward());
                    Flag::SetNotPrimaryAlignment(m_sam1.m_flag, !first);

//                    std::cout << "1: " << m_info.ToString() << std::endl;

                    alignment_length += m_info.alignment_length;
                    alignment_score += m_info.alignment_score;

                }
                if (ar2.IsSet()) {
                    auto len = std::count_if(ar2.Cigar().begin(), ar2.Cigar().end(), [](char c) {
                        return c != 'I';
                    });

                    if (len != record2.sequence.length()) {
                        std::cout << len << std::endl;
                        std::cout << record2.sequence << std::endl;
                        std::cout << ar2.Cigar() << std::endl;
                        exit(11);
                    }

                    GetExtendedAlignmentInfo(m_info, ar2.Cigar(), record2.sequence, record2.quality);
                    ArtoSAM(m_sam2, ar2, m_info, record2);
                    Flag::SetPairedEnd(m_sam2.m_flag, true, both, false, true);
                    Flag::SetRead1Unmapped(m_sam2.m_flag, !ar1.IsSet());
                    Flag::SetRead2ReverseComplement(m_sam2.m_flag, ar2.Forward());
                    Flag::SetNotPrimaryAlignment(m_sam2.m_flag, !first);

//                    std::cout << "2: " << m_info.ToString() << std::endl;

                    alignment_length += m_info.alignment_length;
                    alignment_score += m_info.alignment_score;
                }
                if (both) {
                    m_sam1.m_rnext = m_sam2.m_rname;
                    m_sam2.m_rnext = m_sam1.m_rname;
                    m_sam1.m_pnext = m_sam2.m_pos;
                    m_sam2.m_pnext = m_sam1.m_pos;
                }

                first = false;
                alignments++;
                if (ar1.IsSet() && !m_sam_output.Write(m_sam1.ToString())) {
#pragma omp critical(sam_output)
                    m_sam_output.Write(m_sam_os);
                }
                if (ar2.IsSet() && !m_sam_output.Write(m_sam2.ToString())) {
#pragma omp critical(sam_output)
                    m_sam_output.Write(m_sam_os);
                }
                if (++output_counter == m_max_out) {
                    break;
                }
            }
        }
    };
}