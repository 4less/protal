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
//#include <htslib/sam.h>
#include "AlignmentUtils.h"
#include "SNPUtils.h"

namespace protal {
    static bool CorrectOrientation(AlignmentResult const& a1, AlignmentResult const& a2) {
        if (!a1.IsSet() || !a2.IsSet()) {
            std::cout << "CorrectOrientation both alignments needs to be set" << std::endl;
        }
        return a1.Forward() != a2.Forward();
    }

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
        return Score(a.first.GetAlignmentInfo(), a.second.GetAlignmentInfo());
    }


    static bool PairedAlignmentComparator(PairedAlignment const& a, PairedAlignment const& b) {
        // Return true if a > b
        AlignmentInfo info_a_1 = a.first.GetAlignmentInfo();
        AlignmentInfo info_a_2 = a.second.GetAlignmentInfo();
        AlignmentInfo info_b_1 = b.first.GetAlignmentInfo();
        AlignmentInfo info_b_2 = b.second.GetAlignmentInfo();

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
        sam.m_qname = record.id.substr(0, record.id.length()-2);
        sam.m_flag = 0;
        sam.m_rname = std::to_string(ar.Taxid()) + "_" + std::to_string(ar.GeneId());
        sam.m_pos = info.gene_alignment_start + 1;
        sam.m_mapq = ar.AlignmentScore();
        sam.m_cigar = info.compressed_cigar;
        sam.m_rnext = "*";
        sam.m_pnext = 0;
        sam.m_tlen = ar.Cigar().length();
        sam.m_seq = ar.Forward() ? record.sequence : KmerUtils::ReverseComplement(record.sequence);
        sam.m_qual = record.quality;
        if (!ar.Forward()) reverse(sam.m_qual.begin(), sam.m_qual.end());
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
                    auto& info = ar1.GetAlignmentInfo();
                    alignment_length += info.alignment_length;
                    alignment_score += info.Score();
                    alignment_mismatch_quality += 0;

                }
                if (ar2.IsSet()) {
                    auto& info = ar2.GetAlignmentInfo();
                    alignment_length += info.alignment_length;
                    alignment_score += info.Score();
                    alignment_mismatch_quality += 0;
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
        GenomeLoader& m_genomes;

        double m_min_cigar_ani = 0;
    public:
        size_t alignments = 0;

        ProtalOutputHandler(std::ostream& sam_os, size_t varkit_buffer_capacity, size_t sam_buffer_capacity, GenomeLoader& genomes, double min_cigar_ani=0.0f) :
                m_sam_os(sam_os),
                m_sam_output(sam_buffer_capacity),
                m_min_cigar_ani(min_cigar_ani),
                m_genomes(genomes) {}

        ProtalOutputHandler(ProtalOutputHandler const& other) :
                m_sam_os(other.m_sam_os),
                m_sam_output(other.m_sam_output.Capacity()),
                m_min_cigar_ani(other.m_min_cigar_ani),
                m_genomes(other.m_genomes) {}

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

                auto& info = ar.GetAlignmentInfo();
                ArtoSAM(m_sam, ar, info, record);
                Flag::SetPairedEnd(m_sam.m_flag, false, false, first_pair, !first_pair);
                Flag::SetReadUnmapped(m_sam.m_flag, true);
                Flag::SetReadReverseComplement(m_sam.m_flag, ar.Forward());
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
    template<bool DEBUG=false>
    class ProtalPairedOutputHandler {
    private:
        std::ostream& m_sam_os;
//        ogzstream& ogz;
        BufferedStringOutput m_sam_output;
        SamEntry m_sam1;
        SamEntry m_sam2;

        AlignmentInfo m_info;

        GenomeLoader& m_genomes;

        size_t m_max_out = 1;
        double m_min_cigar_ani = 0;
    public:
        size_t alignments = 0;

        ProtalPairedOutputHandler(std::ostream& sam_os, size_t max_out, size_t varkit_buffer_capacity, size_t sam_buffer_capacity, GenomeLoader& genomes, double min_cigar_ani=0.0f) :
                m_sam_os(sam_os),
                m_sam_output(sam_buffer_capacity),
                m_min_cigar_ani(min_cigar_ani),
                m_max_out(max_out),
                m_genomes(genomes) {}

        ProtalPairedOutputHandler(ProtalPairedOutputHandler const& other) :
                m_sam_os(other.m_sam_os),
                m_sam_output(other.m_sam_output.Capacity()),
                m_min_cigar_ani(other.m_min_cigar_ani),
                m_max_out(other.m_max_out),
                m_genomes(other.m_genomes) {}

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


            int mapq = MAPQv1(alignment_results);

            auto& any = best.first.IsSet() ? best.first : best.second;

            constexpr bool debug = true;

            if constexpr (DEBUG) {
                auto [mapqs, s1, s2] = MAPQv1Debug(alignment_results);
                std::vector<std::string> tokens;
                LineSplitter::Split(record1.id, "-", tokens);
                auto true_ref = tokens[0];

                auto ref = std::to_string(any.Taxid()) + "_" + std::to_string(any.GeneId());
                bool is_correct = ref == true_ref;
#pragma omp critical(errout)
                std::cerr << int(is_correct) << '\t' << mapqs << '\t' << s1 << '\t' << s2 << '\t' << record1.id << std::endl;
            }


            SNPList snps;

            bool valid1 = false;
            bool valid2 = false;

            size_t output_counter = 0;
            for (auto& [ar1, ar2] :  alignment_results) {
                bool both = ar1.IsSet() && ar2.IsSet();

                int alignment_length = 0;
                int alignment_score = 0;
                int adjusted_score = 0;

                if (ar1.IsSet()) {
                    auto& info = ar1.GetAlignmentInfo();
                    ArtoSAM(m_sam1, ar1, info, record1);
                    Flag::SetPairedEnd(m_sam1.m_flag, true, both, true);
                    Flag::SetReadUnmapped(m_sam1.m_flag, false);
                    Flag::SetReadReverseComplement(m_sam1.m_flag, !ar1.Forward());
                    Flag::SetNotPrimaryAlignment(m_sam1.m_flag, !first);


                    alignment_length += info.alignment_length;
                    alignment_score += info.Score();
                    m_sam1.m_mapq = first ? mapq : 0;

                    valid1 = ExtractSNPs(m_sam1, m_genomes.GetGenome(ar1.Taxid()).GetGene(ar1.GeneId()).Sequence(), snps, ar1.Taxid(), ar1.GeneId(), 0);
                }
                if (ar2.IsSet()) {
                    auto len = std::count_if(ar2.Cigar().begin(), ar2.Cigar().end(), [](char c) {
                        return c != 'D';
                    });

                    if (len != record2.sequence.length()) {
#pragma omp critical(debug_out)
                        {
                            std::cout << " --------- Problemo --------- " << std::endl;
                            std::cout << "Len: " << len << std::endl;
                            std::cout << record2.sequence << std::endl;
                            std::cout << ar2.Cigar() << std::endl;

                        }
//                        exit(11);
                    }

                    auto& info = ar2.GetAlignmentInfo();
                    ArtoSAM(m_sam2, ar2, info, record2);
                    Flag::SetPairedEnd(m_sam2.m_flag, true, both, false, true);
                    Flag::SetReadUnmapped(m_sam2.m_flag, false);
                    Flag::SetReadReverseComplement(m_sam2.m_flag, !ar2.Forward());
                    Flag::SetNotPrimaryAlignment(m_sam2.m_flag, !first);

                    alignment_length += info.alignment_length;
                    alignment_score += info.alignment_score;
                    m_sam2.m_mapq = first ? mapq : 0;

                    valid2 = ExtractSNPs(m_sam2, m_genomes.GetGenome(ar2.Taxid()).GetGene(ar2.GeneId()).Sequence(), snps, ar2.Taxid(), ar2.GeneId(), 0);
                }
                if (both) {
                    m_sam1.m_rnext = "=";
                    m_sam2.m_rnext = "=";
                    m_sam1.m_pnext = m_sam2.m_pos;
                    m_sam2.m_pnext = m_sam1.m_pos;
                    Flag::SetMateReverseComplement(m_sam1.m_flag, (FLAG_t)!ar2.Forward());
                    Flag::SetMateReverseComplement(m_sam2.m_flag, (FLAG_t)!ar1.Forward());
                    Flag::SetMateUnmapped(m_sam1.m_flag, false);
                    Flag::SetMateUnmapped(m_sam2.m_flag, false);
                    Flag::SetPairBothAlign(m_sam1.m_flag, true);
                    Flag::SetPairBothAlign(m_sam2.m_flag, true);
                }

                if ((ar1.IsSet() && !valid1) || (ar1.IsSet() && !valid2)) {
                    if (ar1.IsSet() && !valid1) {
#pragma omp critical(err_out)
                        {
                            std::cerr << record1.to_string() << std::endl;
                            std::cerr << m_sam1.ToString() << std::endl;
                            std::string const& reference = m_genomes.GetGenome(ar1.Taxid()).GetGene(ar1.GeneId()).Sequence();
                            PrintAlignment(m_sam1, reference, std::cerr);
                        }
                    }
                    if (ar2.IsSet() && !valid2) {
#pragma omp critical(err_out)
                        {
                            std::cerr << record2.to_string() << std::endl;
                            std::cerr << m_sam2.ToString() << std::endl;
                            auto& reference = m_genomes.GetGenome(ar2.Taxid()).GetGene(ar2.GeneId()).Sequence();
                            PrintAlignment(m_sam2, reference, std::cerr);
                        }
                    }
                    return;
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